
"""
Preform an analysis on the initial photometry to determine the optimal parameters.
Input consists of a fits image with corresponding source extractor photometric catalogue

  This script requires, at minimum, the source extractor values:
  'mag_aper', 'magerr_aper', 'kron_radius', 'x_image', 'y_image', 'fwhm_image' (if r_fed not specified)
  to use aper_auto_check: 'mag_auto', 'mag_err_auto'
  to use half_light: 'flux_radius'

  The apertures used with source extractor must be specified, or the default values are used:
  [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]

To calculate the curve of growth, the photometry catalogue is parsed for stars with accurate
  measurements (i.e. mag < 99.), and the location of each star is compared to the image, where
  stars that are saturated (have cores with 0's) are removed. A curve of growth is then produced
  for each remaining star, and the final curve is the weighted average (weight =  1/sigma^2).
  This curve is shown relative to a chosen feducial radius, default 30, and fitted to an exponential
  curve. The plot is saved if the parameter 'outcurvegrowth' is specified.

Additional operations:
  aper_auto_check: calculation of the difference between the measurement of the magnitude from source
  extractors 'mag_aper' for the fiducial radius specified and the 'mag_auto'

  half_light: fits a linear function to the distribution of 'flux_radius' vs. 'mag_aper', removes objects
  outside of 1 standard deviation, fits a linear function again, and removes objects outside 3 standard
  deviations. The fitting is used twice as the large scatter of 'flux_radius' due to highly extended
  objects effects the linear fitting process.


Example code:

Run image through source extractor to obtain photometric catalogues:
sex -c phot_t7.sex ../release/CB68_J_sub.fits -CATALOG_NAME CB68_J_sex_t7.txt
sex -c phot_t7.sex ../release/CB68_H_sub.fits -CATALOG_NAME CB68_H_sex_t7.txt
sex -c phot_t7.sex ../release/CB68_Ks_sub.fits -CATALOG_NAME CB68_Ks_sex_t7.txt

For every aperture selected in source extractor, calculate the average magnitude, and display as a curve of growth:
python phot_curves.py -img CB68/CB68_J_sub.fits -sf phot_t20/CB68_J_sex_t20.txt --outcog CB68/CB68_J_mags_t20.png --outmags CB68/CB68_J_diff.png
python phot_curves.py -img CB68/CB68_H_sub.fits -sf phot_t20/CB68_H_sex_t20.txt --outcog CB68/CB68_H_mags_t20.png --outmags CB68/CB68_H_diff.png
python phot_curves.py -img CB68/CB68_Ks_sub.fits -sf phot_t20/CB68_Ks_sex_t20.txt --outcog CB68/CB68_Ks_mags_t20.png --outmags CB68/CB68_Ks_diff.png

python phot_curves.py -img CB68/CB68_J_sub.fits -sf phot_t20/CB68_J_sex_t20.txt -aps 5 10 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 --outcog CB68/CB68_J_mags_t20.png --r_fed 30 --outhlr CB68/CB68_J_halflight.png
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from astropy.io import fits, ascii
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
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default=None,
                        help='Source extractor file for J,H, or Ks bands.')

    # optional parameters
    parser.add_argument("--apertures", '-aps',
                        nargs='*',
                        action='store',
                        default=None,
                        help='Apertures of photometry corresponding to sexfils')
    parser.add_argument("--r_fed",
                        action='store',
                        default=None,
                        help='Feducial radius.')

    # optional methods
    parser.add_argument("--outcurvegrowth",
                        action='store',
                        default=None,
                        help='Plot curve of growth.')
    parser.add_argument("--outhalflight",
                        action='store',
                        default=None,
                        help='Plot half right radius.')
    parser.add_argument("--outcomparemag",
                        action='store',
                        default=None,
                        help='Plot difference between MAG_AUTO and MAG_APER.')
    parser.add_argument("--toplot",
                        action='store',
                        default=False,
                        help='Plot identified star coordinates over image.')

    args = parser.parse_args()

    # raise errors if missing input files
    if args.image is None:
        raise Exception, 'ERROR: No image file given'
    if args.sexfile is None:
        raise Exception, 'ERROR: No source extractor files given'
    if args.outfile is None:
        raise Exception, 'ERROR: No output files given'

    # raise warnings if missing input values
    if args.apertures is None:
        apertures = [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
        print 'WARNING: No apertures specified, using default {}'.format(apertures)
    else:
        apertures = args.apertures
    if args.r_fed is None:
        r_fed = 30
        print 'WARNING: No apertures specified, using default {}'.format(r_fed)
    else:
        r_fed = args.r_fed

    # remove saturated stars from catalogue, output as new catalogue
    phottable = remove_saturated(args.image, args.sexfile)

    # calculate curve of growths
    if args.outcurvegrowth is not None:
        curve_of_growths(phottable, np.array(apertures, dtype=np.int), r_fed, args.outcurvegrowth)

    # compare MAG_APER with MAG_AUTO
    if args.outcomparemag is not None:
        aper_auto_check(phottable, apertures, args.outcomparemag)

    # calculate half light radius
    if args.outhalflight is not None:
        iap = np.argmin(abs(int(r_fed) - np.array(apertures, dtype=np.int)))
        half_light(phottable['MAG_APER_{}'.format(iap + 1)], phottable['FLUX_RADIUS'], args.outhalflight)

    # plot coordinates of identified stars over image
    if args.toplot:
        imgdata, imghdr = read_fits(args.image)
        plt.imshow(imgdata, norm=LogNorm())
        plt.scatter(phottable['X_IMAGE'], phottable['Y_IMAGE'], marker='x', color='k')
        plt.show()


def remove_saturated(image, sexfile):
    imgdata, imghdr = read_fits(image)  # Access fits data to read locations of saturated stars

    phottable = ascii.read(sexfile, format='sextractor')
    unsat_idxs = detect_satur(imgdata, phottable['X_IMAGE'], phottable['Y_IMAGE'], phottable['KRON_RADIUS'])
    phottable = phottable[unsat_idxs]

    phottable = phottable[np.where(phottable['MAG_APER'] < 99.)]

    return phottable


def detect_satur(img_data, x, y, r, max_zeros=1.):
    """
    Identify saturation by looking for at least 'nz' zeros in a radius 'r' around the measured coordinate
    """

    idx = []
    for i in range(len(x)):
        zeros = []
        rr = int(r[i])
        # r**2 = (xx - x0)**2 + (yy - y0)**2
        for xx in range(int(x[i]) - rr, int(x[i]) + rr):
            yy_p = np.sqrt(rr ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            yy_m = -np.sqrt(rr ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            for yy in range(int(yy_m), int(yy_p)):
                if img_data[yy, xx] == 0.:
                    zeros.append(1.)
        if len(zeros) < max_zeros:
            idx.append(i)
    return idx


def curve_of_growths(phottable, apertures, r_feducial=None, outimage=None):
    """
    Plot curve of growth, relative to a given feducial radius
    """

    # make numpy ndarrays from the dataframes, easier to plot
    magtable = np.ones((len(apertures), len(phottable)))
    magerrtable = np.ones((len(apertures), len(phottable)))
    magtable[0, :] = np.array(phottable['MAG_APER'])
    magerrtable[0, :] = np.array(phottable['MAGERR_APER'])

    for i in range(1,len(apertures)):
        magtable[i, :] = np.array(phottable['MAG_APER_{}'.format(i)])
        magerrtable[i, :] = np.array(phottable['MAGERR_APER_{}'.format(i)])

    # determine the feducial aperture by closest value to given feducial radius
    if r_feducial is None:
        r_feducial = np.average(phottable['FWHM_IMAGE'], weights=(magerrtable[i, :]) ** (-2)) * 5.
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

    aps_fine = np.linspace(np.min(apertures), np.max(apertures), 100.)
    fit = expo(aps_fine, popt[0], popt[1], popt[2])

    data_fitted = expo(apertures, *popt)
    # determine "goodness of fit" by chi square statistic
    chi2 = np.sum((average_rel- data_fitted)**2 / data_fitted)

    ax = plt.gca()
    ax.plot(aps_fine, fit, '-r', label=r'exponential fit, $\chi^2=${:<2.3f}'.format(chi2))
    ax.errorbar(apertures, average_rel, yerr=yerr, fmt='.', label=r'$\Delta$ magnitude')
    ax.plot(aps_fine, aps_fine*0., ':k')
    ax.invert_yaxis()
    ax.set_xlabel('aperture')
    ax.set_ylabel(r'$\Delta$ magnitude relative to aperture {}'.format(ap_fed))
    plt.legend(loc=4)
    if outimage is not None:
        plt.savefig(outimage)
    plt.show()


def aper_auto_check(phottable, apertures, outfig=None):
    """
    Plot difference between mag_auto and mag_aper for the apertures
    """

    mag_auto = np.array(phottable['MAG_AUTO'])
    magerr_auto = np.array(phottable['MAGERR_AUTO'])
    mag_auto_avg = np.average(mag_auto, weights=magerr_auto ** (-2))
    magerr_auto_avg = np.sqrt(1 / np.sum(magerr_auto ** (-2)))

    mag_aper_avg = np.zeros(len(apertures))
    magerr_aper_avg = np.zeros(len(apertures))
    for i in range(len(apertures)):
        mag_aper = np.array(phottable['MAG_APER_{}'.format(i+1)])
        magerr_aper = np.array(phottable['MAGERR_APER_{}'.format(i+1)])
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

    # plt.plot(apertures, diff)
    plt.scatter(phottable['magerr_auto'].values, phottable['magerr_aper_{}'.format(ap)].values)
    plt.xlabel('magerr_auto_avg')
    plt.ylabel(r'magerr_aper_avg')
    if outfig is not None:
        plt.savefig(outfig)
    plt.show()


def half_light(x, y, outfig=None, toplot=True):
    """
    Plot flux radius as a function of magnitude
    """

    linear = lambda x, m, b: m * x + b

    weight = x ** (-1)
    p, cov = optimize.curve_fit(linear, x, y, p0=[1., np.min(y)], sigma=weight, absolute_sigma=True)
    perr = np.sqrt(np.diag(cov))
    fit = p[0] * x + p[1]

    idx = (y <= fit + np.std(y)) * (fit - np.std(y) <= y)
    yy = y[np.where(idx)]
    xx = x[np.where(idx)]

    weight = xx ** (-1)
    p, cov = optimize.curve_fit(linear, xx, yy, p0=[1., np.min(yy)], sigma=weight, absolute_sigma=True)
    perr = np.sqrt(np.diag(cov))
    fit = p[0] * xx + p[1]
    fiterr = 3 * np.std(yy)

    if toplot:
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

    idx_bool = np.zeros_like(x, dtype=bool)
    idx_bool[np.where(idx)] = True

    return idx


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
