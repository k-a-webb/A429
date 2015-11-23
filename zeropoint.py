__author__ = 'kawebb'

"""
Procedure:
Determine the magnitude zeropoint of our CFHT observations by comparing bright
unsaturated stars to the same stars in the 2mass catalog.
Apply this zeropoint to all of the stars measured.
"""

import os
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits


def main():
    parser = argparse.ArgumentParser(
        description='For a given image, and list of star positions, return list of stars which are not saturated')
    parser.add_argument("--image", '-img',
                        action="store",
                        default='/Users/kawebb/Dropbox/A429/L1544/L1544_J_sub.fits',
                        help="Fits image of star field")
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default='/Users/kawebb/Dropbox/A429/L1544/L1544_J_sex_2m.txt',
                        help='List of star positions as output by sextractor for 2mass comparison.')
    parser.add_argument("--usout", '-us',
                        action='store',
                        default=None,
                        help='List of unsaturated stars in image, will be created if not specified')
    parser.add_argument("--tmass", '-tm',
                        action='store',
                        default='/Users/kawebb/Dropbox/A429/2mass_L1544.tbl',
                        help='File of 2Mass stars in region of interest')
    parser.add_argument("--tmout",
                        action='store',
                        default=None,
                        help='File of 2Mass stars in the image')
    parser.add_argument("--toplot", '-plt',
                        action='store',
                        default=False,
                        help='Whether or not to plot saturated stars over top of the image')
    parser.add_argument("--band", '-b',
                        action='store',
                        default='J',
                        help='What colour band the image is in')
    args = parser.parse_args()

    uslist = remove_satur(args.image, args.sexfile, args.usout, args.toplot)
    tmstars = parse_2mass(args.tmass, uslist, args.tmout, args.band)
    offset = compare_magnitudes(uslist, tmstars, args.toplot)

    plot_plots(args.image, uslist, tmstars, args.toplot)


def plot_plots(ffile, uslist, tmstars, toplot):

    if not toplot:
        return

    fdata, fhdr = read_fits(ffile)

    # Unsaturated stars
    plt.imshow(fdata, origin='lower', norm=matplotlib.colors.LogNorm())
    plt.plot(uslist.x.values, uslist.y.values, 'x')
    plt.title('Unsaturated stars')
    plt.show()

    idx = np.zeros(len(tmstars))
    for i in tmstars.index:
        sig = 0.005
        stars = []
        while len(stars) == 0:
            ras = uslist.query('({} < x_wcs) & (x_wcs < {})'.format(tmstars.ra[i] - sig, tmstars.ra[i] + sig))
            stars = ras.query('({} < y_wcs) & (y_wcs < {})'.format(tmstars.dec[i] - sig, tmstars.dec[i] + sig))
            sig += 0.001
        idx[i] = stars.index[0]
    us2mstars = uslist.loc[idx, :]

    # Unsaturated stars and 2Mass stars
    plt.imshow(fdata, origin='lower', norm=matplotlib.colors.LogNorm())
    plt.plot(uslist.x.values, uslist.y.values, 'x', label='Unsaturated stars')
    plt.plot(us2mstars.x.values, us2mstars.y.values, 'o', label='2Mass stars')
    plt.legend()
    plt.title('Saturated stars')
    plt.show()


def compare_magnitudes(uslist, tmstars, toplot=False):
    """
    Comparing magnitudes from source extractor and 2mass.
    Calculate the magnitude offset between the two to find correction to zeropoint in source extractor
    """

    print '  Comparing colours of 2Mass stars and stars in the image'
    tmmag = tmstars.j_m.values
    sexmag = uslist.mag_aper.values
    sexmagerr = uslist.magerr_aper.values

    wgt = (sexmagerr)**(-2)

    # Weighting as function of magnitude
    # df = pd.DataFrame(data={'mag':sexmag, 'weight':wgt})
    # dfs = df.sort('mag')
    # plt.plot(dfs.mag.values, dfs.weight.values)
    # plt.show()

    offset = tmmag-sexmag
    wavg = np.average(offset, weights=wgt)
    sig = np.var(offset)

    print 'weighted avg offset, unweighted avg offset, variance: {} {} {}'.format(wavg, np.mean(offset), sig)
    print 'weighted average, unweighted average: {} {}'.format(np.average(sexmag, weights=wgt), np.mean(sexmag))

    return wavg


def parse_2mass(tmfile, uslist, tmout=None, band='J'):
    """
    For a given list of coordinates returned by source extractor, find the same star in the
    2Mass catalog, and save the cropped catalog list to a file
    """

    if tmout is None:
        tmout = tmfile.split('.')[0] + '_{}.tbl'.format(band)
    if os.path.exists(tmout):
        print '  2Mass list has already been parsed for this image: {}'.format(tmout)
        return read_2mass(tmout, 1)

    assert os.path.exists(tmfile), '2Mass file does not exist: {}'.format(tmfile)
    print '  Identifying 2Mass catalog stars with Unsaturated stars in CFHT image'

    tmlist = read_2mass(tmfile, skiplines=55)

    idx = np.zeros(len(uslist))
    for ii,i in enumerate(uslist.index):
        sig = 0.005
        stars = []
        while len(stars) == 0:
            ras = tmlist.query('({} < ra) & (ra < {})'.format(uslist.x_wcs[i] - sig, uslist.x_wcs[i] + sig))
            stars = ras.query('({} < dec) & (dec < {})'.format(uslist.y_wcs[i] - sig, uslist.y_wcs[i] + sig))
            sig += 0.001
        if len(stars) > 0:
            idx[ii] = stars.index[0]
        else:
            print 'No 2Mass star found at {},{}'.format(uslist.x_wcs[i],uslist.y_wcs[i])
            idx[ii] = np.nan()

    tmstars = tmlist.loc[idx, :]

    tmstars.to_csv(tmout, index=False, sep=' ')
    return tmstars


def remove_satur(ffile, sfile, ofile=None, toplot=False):
    """
    Detect saturated stars in the image, and remove them from the list returned by source extractor
    In our images, saturated stars have low (0) values at the center rather than a very high value
    and are therefore not identified by source extractor automatticaly.
    """

    if ofile is None:
        ofile = sfile.strip('.txt') + '_unsat.txt'

    if os.path.exists(ofile):
        print '  List of unsaturated stars already exists: {}'.format(ofile)
        return read_sex(ofile, skiplines=1)

    print "  Creating list of unsaturated stars for image {}, writing to {}".format(ffile, ofile)

    stable = read_sex(sfile, skiplines=13)
    fdata, fhdr = read_fits(ffile)

    # Iterate through stars
    ix = detect_satur(fdata, stable['x'].values, stable['y'].values, stable['b'].values)
    us_stable = stable[ix]
    # us_stable.to_csv(ofile, index=False, sep=' ')

    if toplot:
        plt.imshow(fdata, origin='lower', norm=matplotlib.colors.LogNorm())
        s_stable = stable[np.invert(ix)]
        plt.plot(s_stable.x.values, s_stable.y.values, 'x')
        plt.title('Saturated stars')
        plt.show()

    return us_stable


def detect_satur(fdata, x, y, r, nz=1.):
    """
    Identify saturation by looking for at least 'nz' zeros in a radius 'r' around the measured coordinate
    """

    idx = []
    n = np.arange(len(x))

    for i in n:
        zeros = []
        rr = int(r[i]) + 1
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
    return ix


def read_sex(sfile, skiplines=13):
    #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
    #   2 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
    #   3 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   4 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
    #   5 THRESHOLD              Detection threshold above background                       [count]
    #   6 X_IMAGE                Object position along x                                    [pixel]
    #   7 Y_IMAGE                Object position along y                                    [pixel]
    #   8 X_WORLD                Barycenter position along world x axis                     [deg]
    #   9 Y_WORLD                Barycenter position along world y axis                     [deg]
    #  10 B_IMAGE                Profile RMS along minor axis                               [pixel]
    #  11 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
    #  12 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
    #  13 FLAGS                  Extraction flags
    table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python',
                        names=['flux_aper', 'fluxerr_aper', 'mag_aper', 'magerr_aper', 'thresh',
                               'x', 'y', 'x_wcs', 'y_wcs', 'b', 'fhm_image', 'fwhm_world', 'flag'])
    return table


def read_fits(ffile):
    with fits.open(ffile) as hdu:
        # hdu.info()
        # Filename: L1544_J_sub.fits
        # No.    Name         Type      Cards   Dimensions   Format
        # 0    PRIMARY     PrimaryHDU      54   (5475, 5493)   float32
        fdata = hdu[0].data
        fhdr = hdu[0].header
    return fdata, fhdr


def read_2mass(tmfile, skiplines=55):
    # ra, dec, j_m, j_cmsig, j_msigcom, j_snr, h_m, h_cmsig, h_msigcom, h_snr, k_m, k_cmsig, k_msigcom, k_snr
    #   rd_flg, dist, angle, j_h, h_k, j_k
    table = pd.read_csv(tmfile, skiprows=skiplines, sep=r"\s*", engine='python',
                        names=['ra', 'dec', 'j_m', 'j_cmsig', 'j_msigcom', 'j_snr', 'h_m', 'h_cmsig', 'h_msigcom',
                               'h_snr', 'k_m', 'k_cmsig', 'k_msigcom', 'k_snr', 'rd_flg', 'dist', 'angle', 'j_h',
                               'h_k', 'j_k'])
    return table


if __name__ == '__main__':
    main()
