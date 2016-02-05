__author__ = 'kawebb'

"""
sex -c phot_t7.sex ../release/CB68_J_sub.fits -CATALOG_NAME CB68_J_sex_t7.txt
sex -c phot_t7.sex ../release/CB68_H_sub.fits -CATALOG_NAME CB68_H_sex_t7.txt
sex -c phot_t7.sex ../release/CB68_Ks_sub.fits -CATALOG_NAME CB68_Ks_sex_t7.txt

python phot_curves.py -img CB68/CB68_J_sub.fits -sf phot_t20/CB68_J_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_J_mags_t20.txt --outcog CB68/CB68_J_mags_t20.png --outmags CB68/CB68_J_diff.png
python phot_curves.py -img CB68/CB68_H_sub.fits -sf phot_t20/CB68_H_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_H_mags_t20.txt --outcog CB68/CB68_H_mags_t20.png --outmags CB68/CB68_H_diff.png
python phot_curves.py -img CB68/CB68_Ks_sub.fits -sf phot_t20/CB68_Ks_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_Ks_mags_t20.txt --outcog CB68/CB68_Ks_mags_t20.png --outmags CB68/CB68_Ks_diff.png

python phot_curves.py -img CB68/CB68_J_sub.fits -sf phot_t20/CB68_J_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outfile CB68/CB68_J_mags_t20.txt --outcog CB68/CB68_J_mags_t20.png --r_fed 30 --outhlr CB68/CB68_J_halflight.png
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
from astropy.io import fits
from matplotlib.colors import LogNorm
from scipy import optimize, stats


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
    parser.add_argument("--outcog",
                        action='store',
                        default=None,
                        help='Output plot of curve of growths.')
    parser.add_argument("--outhlr",
                        action='store',
                        default=None,
                        help='Output plot of half right radius.')
    parser.add_argument("--outmags",
                        action='store',
                        default=None,
                        help='Output plot of mag aper and auto difference.')
    parser.add_argument("--toplot",
                        action='store',
                        default=False,
                        help='Wheter or not to plot figures.')
    parser.add_argument("--r_fed",
                        action='store',
                        default=None,
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
        phottable.reset_index(inplace=True, drop=True)  # drop index column as now out of order

        if args.outfile is None:
            outfile = 'test.txt'
        else:
            outfile = args.outfile

        phottable[names].to_csv(outfile, index=False)
    else:
        phottable = pd.read_csv(args.outfile)

    # curve_of_growths(phottable, np.array(args.apertures, dtype=np.int), args.r_fed, args.outcog)
    # aper_auto_check(phottable, args.apertures, args.outmags)

    if args.r_fed is None:
        args.r_fed = np.mean(phottable['fwhm_image'].values) * 10.
    iap = np.argmin(abs(int(args.r_fed) - np.array(args.apertures, dtype=np.int)))
    aper = args.apertures[iap]
    half_light(phottable['mag_aper_{}'.format(aper)].values, phottable['flux_radius'].values, args.outhlr)

    if args.toplot:
        plt.imshow(imgdata, norm=LogNorm())
        plt.scatter(phottable['x_pix'].values, phottable['y_pix'].values, marker='x', color='k')
        # plt.colorbar()
        plt.show()





def aper_auto_check(phottable, apertures, outfig=None):
    """
    Plot difference between mag_auto and mag_aper for the apertures
    """

    mag_auto = phottable['mag_auto'].values
    magerr_auto = phottable['magerr_auto'].values
    mag_auto_avg = np.average(mag_auto, weights=magerr_auto ** (-2))
    magerr_auto_avg = np.sqrt(1 / np.sum(magerr_auto ** (-2)))

    mag_aper_avg = np.zeros(len(apertures))
    magerr_aper_avg = np.zeros(len(apertures))
    for i, ap in enumerate(apertures):
        mag_aper = phottable['mag_aper_{}'.format(ap)].values
        magerr_aper = phottable['magerr_aper_{}'.format(ap)].values
        mag_aper_avg[i] = np.average(mag_aper, weights=magerr_aper ** (-2))
        magerr_aper_avg[i] = np.sqrt(1 / np.sum(magerr_aper ** (-2)))

    diff = np.subtract(mag_auto_avg, mag_aper_avg)
    differr = np.sqrt(magerr_auto_avg ** 2 + magerr_aper_avg ** 2)

    plt.plot(apertures, diff)
    plt.errorbar(apertures, diff, yerr=differr)
    plt.xlabel('apertures')
    plt.ylabel(r'$mag_{auto} - mag_{aper}$')
    if outfig is not None:
        plt.savefig(outfig)
    plt.show()


def half_light(x, y, outfig=None):
    """
    Plot flux radius as a function of magnitude
    """

    linear = lambda x, m, b: m * x + b

    weight = x ** (-1)
    p, cov = optimize.curve_fit(linear, x, y, p0=[1., np.min(y)], sigma=weight, absolute_sigma=True)
    perr = np.sqrt(np.diag(cov))
    fit = p[0] * x + p[1]

    idx = (y <= fit + np.std(y))*(fit - np.std(y) <= y)
    yy = y[np.where(idx)]
    xx = x[np.where(idx)]

    weight = xx ** (-1)
    p, cov = optimize.curve_fit(linear, xx, yy, p0=[1., np.min(yy)], sigma=weight, absolute_sigma=True)
    perr = np.sqrt(np.diag(cov))
    fit = p[0] * xx + p[1]
    fiterr = 3*np.std(yy)

    plt.scatter(x, y, marker='.', alpha=0.3, color='g')
    plt.plot(xx, fit, '-k', label="linear fit")
    plt.plot(xx, fit - fiterr, '-b', label="linear fit")
    plt.plot(xx, fit + fiterr, '-b', label="linear fit")
    plt.xlabel('mag')
    plt.ylabel('half light radius')
    plt.ylim(0, 5)
    if outfig is not None:
        plt.savefig(outfig)
    plt.show()

    idx_bool =  np.zeros_like(x, dtype=bool)
    idx_bool[np.where(idx)] = True

    return idx_bool


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
    if r_feducial is None:
        r_feducial = np.average(phottable['fwhm_image'].values, weights=(magerrtable[i, :]) ** (-2)) * 5.
    iap = np.argmin(abs(int(r_feducial) - np.array(apertures)))
    ap_fed = apertures[iap]

    print 'Given feducial radius {}, closest aperture in data set {}'.format(r_feducial, ap_fed)

    # calculate the weighted average and error, same dimension as apertures
    average = apertures * 0.
    yerr = apertures * 0.
    for i in range(len(apertures)):
        weight = (magerrtable[i, :]) ** (-2)
        average[i] = np.average(magtable[i, :], weights=weight)
        yerr[i] = np.sqrt(1 / np.sum(weight))  # sqrt(1/sum(w))

    average_rel = average * 0.
    for i in range(len(apertures)):
        average_rel[i] = average[i] - average[iap]

    expo = lambda x, amp, alpha, b: amp * np.exp(-1 * x * alpha) + b
    popt, pcov = optimize.curve_fit(expo, apertures, average_rel, p0=[1.5, .2, 0.], sigma=yerr, absolute_sigma=True)

    fit = expo(apertures, popt[0], popt[1], popt[2])

    ax = plt.gca()
    ax.plot(apertures, average_rel, '-b', label=r'$\Delta mag$')
    ax.plot(apertures, fit, '-r', label=r'$\Delta mag$ fit')
    ax.errorbar(apertures, average_rel, yerr=yerr)
    ax.plot(apertures, apertures * 0., '-k')
    ax.invert_yaxis()
    ax.set_xlabel('apertures')
    ax.set_ylabel(r'$\Delta$ mag relative to aperture {}'.format(ap_fed))
    plt.legend(loc=4)
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


def read_sex(sfile, skiplines=15, aperture=None):
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
    names = np.concatenate((names, ['mag_auto', 'magerr_auto', 'kron_radius', 'flux_radius', 'x_pix', 'y_pix', 'x_wcs',
                                    'y_wcs', 'fwhm_image', 'fwhm_world', 'flags']))

    #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
    #  24 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
    #  47 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #  70 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
    #  93 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
    #  94 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
    #  95 KRON_RADIUS            Kron apertures in units of A or B
    #  96 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
    #  97 X_IMAGE                Object position along x                                    [pixel]
    #  98 Y_IMAGE                Object position along y                                    [pixel]
    #  99 X_WORLD                Barycenter position along world x axis                     [deg]
    # 100 Y_WORLD                Barycenter position along world y axis                     [deg]
    # 101 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
    # 102 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
    # 103 FLAGS                  Extraction flags

    # read in file as pandas dataframe
    table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python', names=names)

    # drop all rows with saturated stars
    for ap in aperture:
        table = table[table['mag_aper_{}'.format(ap)] != 99.]

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
