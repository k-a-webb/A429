__author__ = 'kawebb'

"""
sex -c phot_t7.sex ../release/CB68_J_sub.fits -CATALOG_NAME CB68_J_sex_t7.txt
sex -c phot_t7.sex ../release/CB68_H_sub.fits -CATALOG_NAME CB68_H_sex_t7.txt
sex -c phot_t7.sex ../release/CB68_Ks_sub.fits -CATALOG_NAME CB68_Ks_sex_t7.txt

python phot_curves.py -img CB68/CB68_J_sub.fits -sf phot_t7/CB68_J_sex_t7.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_J_mags_t7.txt --outimage CB68/CB68_J_mags_t7.png
python phot_curves.py -img CB68/CB68_H_sub.fits -sf phot_t7/CB68_H_sex_t7.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_H_mags_t7.txt --outimage CB68/CB68_H_mags_t7.png
python phot_curves.py -img CB68/CB68_Ks_sub.fits -sf phot_t7/CB68_Ks_sex_t7.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_Ks_mags_t7.txt --outimage CB68/CB68_Ks_mags_t7.png

python phot_curves.py -img CB68/CB68_J_sub.fits -sf phot_t20/CB68_J_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_J_mags_t20.txt --outimage CB68/CB68_J_mags_t20.png
python phot_curves.py -img CB68/CB68_H_sub.fits -sf phot_t20/CB68_H_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_H_mags_t20.txt --outimage CB68/CB68_H_mags_t20.png
python phot_curves.py -img CB68/CB68_Ks_sub.fits -sf phot_t20/CB68_Ks_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_Ks_mags_t20.txt --outimage CB68/CB68_Ks_mags_t20.png
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
from astropy.io import fits
from matplotlib.colors import LogNorm
from scipy import optimize


def main():
    parser = argparse.ArgumentParser(
        description='For a given image, and list of star positions from source extractor, return list of stars ' \
                    'which are not saturated, and plot curve of growths')
    parser.add_argument("--image", '-img',
                        action="store",
                        default=None,
                        help="Fits image of star field")
    parser.add_argument("--apertures", '-aps',
                        nargs='*',
                        action='store',
                        default=None,
                        help='Apertures of photometry corresponding to sexfils')
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default=None,
                        help='Source extractor file for J,H, or Ks bands.')
    parser.add_argument("--outfile", '-of',
                        action='store',
                        default=None,
                        help='Outfile of matched magnitudes.')
    parser.add_argument("--outimage",
                        action='store',
                        default=None,
                        help='Output plot of curve of growths.')
    parser.add_argument("--toplot",
                        action='store',
                        default=False,
                        help='Wheter or not to plot figures.')
    parser.add_argument("--r_fed",
                        action='store',
                        default=30,
                        help='Feducial radius.')
    args = parser.parse_args()

    if args.image is None:
        raise Exception, 'ERROR: No image file given'
    if args.apertures is None:
        raise Exception, 'ERROR: No apertures specified'
    if args.sexfile is None:
        raise Exception, 'ERROR: No source extractor files given'

    imgdata, imghdr = read_fits(args.image)  # Access fits data to read locations of saturated stars

    if not os.path.exists(args.outfile):

        phottable, names = read_sex(args.sexfile, aperture=args.apertures)  # Import source extractor file
        unsat_idxs = detect_satur(imgdata, phottable['x_pix'].values, phottable['y_pix'].values,
                                  phottable['kron_radius'].values)
        phottable = phottable[unsat_idxs]  # splice phot. table according to index array
        phottable.reset_index(inplace=True,
                              drop=True)  # drop index column (inherent to pandas dataframe) as out of order

        if args.outfile is None:
            outfile = 'test.txt'
        else:
            outfile = args.outfile

        phottable[names].to_csv(outfile, index=False)
    else:
        phottable = pd.read_csv(args.outfile)

    curve_of_growths(phottable, np.array(args.apertures, dtype=np.int), int(args.r_fed), args.outimage)

    if args.toplot:
        plt.imshow(imgdata, norm=LogNorm())
        plt.scatter(phottable['x_pix'].values, phottable['y_pix'].values, marker='x', color='k')
        # plt.colorbar()
        plt.show()


def curve_of_growths(phottable, apertures, r_feducial, outimage=None):
    """
    Plot curve of growth, relative to a given feducial radius
    """

    # make numpy ndarrays from the dataframes, easier to plot
    magtable = np.ones((len(apertures), len(phottable)))
    magerrtable = np.ones((len(apertures), len(phottable)))
    for i, ap in enumerate(apertures):
        magtable[i, :] = phottable['mag_aper_{}'.format(ap)].values
        magerrtable[i, :] = phottable['magerr_aper_{}'.format(ap)].values

    # determine the feducial aperture by closest value to given feducial radius
    iap = np.argmin(abs(r_feducial - np.array(apertures)))
    ap_fed = apertures[iap]
    print 'Given feducial radius {}, closest aperture in data set {}'.format(r_feducial, ap_fed)

    # calculate the weighted average and error, same dimension as apertures
    average = apertures * 0.
    yerr = apertures * 0.
    for i in range(len(apertures)):
        weight = (magerrtable[i, :]) ** (-2)
        average[i] = np.average(magtable[i, :], weights=weight)
        yerr[i] = np.sqrt(1/np.sum(weight))  # sqrt(1/sum(w))

    average_rel = average * 0.
    for i in range(len(apertures)):
        average_rel[i] = average[i] - average[iap]

    expo = lambda x, amp, alpha: amp*np.exp(-1*x*alpha)
    popt, pcov = optimize.curve_fit(expo, apertures, average_rel, p0=[1.5, .2], sigma=yerr, absolute_sigma=True)


    fit = expo(apertures, popt[0], popt[1])

    ax = plt.gca()
    ax.plot(apertures, average_rel, '-b', label=r'$\Delta mag$')
    ax.plot(apertures, fit, '-r', label=r'$\Delta mag$ fit')
    ax.errorbar(apertures, average_rel, yerr=yerr)
    ax.plot(apertures, apertures * 0., '-k')
    ax.invert_yaxis()
    ax.set_xlabel('apertures')
    ax.set_ylabel(r'$\Delta$ mag relative to aperture {}'.format(ap_fed))
    if outimage is not None:
        plt.savefig(outimage)
    plt.show()


def detect_satur(fdata, x, y, r, nz=1.):
    """
    Identify saturation by looking for at least 'nz' zeros in a radius 'r' around the measured coordinate
    """

    idx = []
    n = np.arange(len(x))

    for i in n:
        zeros = []
        rr = int(r[i])
        # r**2 = (xx - x0)**2 + (yy - y0)**2
        for xx in range(int(x[i]) - rr, int(x[i]) + rr):
            yy_p = np.sqrt(rr ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            yy_m = -np.sqrt(rr ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            for yy in range(int(yy_m), int(yy_p)):
                if fdata[yy, xx] == 0.:
                    zeros.append(1.)
        if len(zeros) >= nz:
            idx.append(i)

    ix = np.invert(np.in1d(n.ravel(), idx).reshape(n.shape))
    np.any(ix) is False, '  WARNING: No unsaturated stars detected'
    return ix


def read_sex(sfile, skiplines=12, aperture=None):
    assert os.path.exists(sfile), 'ERROR: Source extractor file {} does not exist'.format(sfile)

    names = []
    for ap in aperture:
        names.append('flux_aper_{}'.format(ap))
    for ap in aperture:
        names.append('fluxerr_aper_{}'.format(ap))
    for ap in aperture:
        names.append('mag_aper_{}'.format(ap))
    for ap in aperture:
        names.append('magerr_aper_{}'.format(ap))
    names = np.concatenate(
        (names, ['kron_radius', 'x_pix', 'y_pix', 'x_wcs', 'y_wcs', 'fwhm_image', 'fwhm_world', 'flags']))

    #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
    #  24 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
    #  47 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #  70 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
    #  93 KRON_RADIUS            Kron apertures in units of A or B
    #  94 X_IMAGE                Object position along x                                    [pixel]
    #  95 Y_IMAGE                Object position along y                                    [pixel]
    #  96 X_WORLD                Barycenter position along world x axis                     [deg]
    #  97 Y_WORLD                Barycenter position along world y axis                     [deg]
    #  98 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
    #  99 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
    # 100 FLAGS                  Extraction flags

    # read in file as pandas dataframe
    table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python', names=names)

    # drop all rows with saturated stars
    table = table[table[names[2].format(aperture)] != 99.]

    return table, names


def read_fits(ffile):
    assert os.path.exists(ffile), 'ERROR: Fits file {} does not exist'.format(ffile)

    with fits.open(ffile) as hdu:
        # hdu.info()
        # Filename: L1544_J_sub.fits
        # No.    Name         Type      Cards   Dimensions   Format
        # 0    PRIMARY     PrimaryHDU      54   (5475, 5493)   float32
        fdata = hdu[0].data
        fhdr = hdu[0].header
    return fdata, fhdr


if __name__ == '__main__':
    main()
