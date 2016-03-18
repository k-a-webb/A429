__author__ = 'kawebb'

"""
Procedure:
Determine the magnitude zeropoint of our CFHT observations by comparing bright unsaturated stars to the same stars
in the 2mass catalog.

Input includes either: image and corresponding source extractor file
 as well as: a 2mass file listing photometric values for bright stars in the region of the image
                - this can be retreived at: http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd
             the aperture to calculate the zeropoint at. Any float may be specified, and the aperture
               closest to an aperture used with source extractor will be used.
               Default source extractor apertures are:
               [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
             the band of the photometry, used to specify the magnitudes from the 2mass file for comparison to the
               measured magnitudes for the calculation of the zeropoint
 optional:  an outfile to print a few of the calculated values to

To calculate the zeropoint correction to be entered into source extractor for the correct calculation of magnitudes,
 the process is as follows:
 if source extractor and file given:
   parse photometric cataloque to remove saturated stars
 else:
   read in photmetry file output from phot_curves, which has already been parsed
 use the cloud names, and shape paremeters specified in the header to remove stars which are in the region of the
   cloud core
 remove stars which do not have a star light PSF using the phot_curves half_light function.
 Identify the stars in our image remaining with those from the 2mass catalogue, uses a nearest neighbour approach
   with a maximum separation of 0.005 degrees (hardcoded into query_seperation function)
 subtract the magnitudes measured from those in the 2mass catalogue, calculate a weighted average


Example commands:
python zeropoint.py -img L429/L429_J_sub.fits -sf phot_t20/L429_J_sex_t20.txt -band j -ap 30 -tm L429/2mass_L429.tbl -sf phot_t20/L429_J_sex_t20.txt --cloud L429
python zeropoint.py -img L429/L429_KS_sub.fits -sf phot_t20/L429_Ks_sex_t20.txt -band k -ap 28 -tm L429/2mass_L429.tbl -sf phot_t20/L429_Ks_sex_t20.txt --cloud L429
python zeropoint.py -img L429/L429_H_sub.fits -sf phot_t20/L429_H_sex_t20.txt -band h -ap 26 -tm L429/2mass_L429.tbl -sf phot_t20/L429_H_sex_t20.txt --cloud L429

Rerun the photomerty with the zeropoint correction.
sex -c phot_t3.sex ../release/CB68_J_sub.fits -CATALOG_NAME CB68_J_sex_t3_ap30.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 30-0.192943311392
sex -c phot_t3.sex ../release/CB68_Ks_sub.fits -CATALOG_NAME CB68_Ks_sex_t3_ap24.txt -PHOT_APERTURES 24 -MAG_ZEROPOINT 30+0.605746118244
sex -c phot_t3.sex ../release/CB68_H_sub.fits -CATALOG_NAME CB68_H_sex_t3_ap40.txt -PHOT_APERTURES 40 -MAG_ZEROPOINT 30+0.471278827929

"""

import numpy as np
import argparse
from astropy.io import ascii
import phot_curves

_CLOUDS = ['CB68', 'L429', 'L1521E', 'L1544', 'L1552']
_CENTERS = np.array([[2610.6899, 2868.1778], [2496.3232, 2158.1909], [2532.6025, 2753.7333], [2345.069, 2855.932],
                        [2710.337, 2593.2019]])  # in pixels
_SIZES = np.array([[1382.2207, 1227.3661], [1869.7259, 2345.7605], [1880.5105, 1777.3177], [1788.7782, 1570.9142],
                      [1639.7134, 177.3117]])  # in pixels

def main():
    parser = argparse.ArgumentParser(
        description='For a given image, and list of star positions, return list of stars which are not saturated')
    parser.add_argument("--image", '-img',
                        action="store",
                        default=None,
                        help="Fits image of star field")
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default=None,
                        help='List of star positions as output by sextractor for 2mass comparison.')
    parser.add_argument("--tmass", '-tm',
                        action='store',
                        default=None,
                        help='File of 2Mass stars in region of interest')
    parser.add_argument("--outfile", '-out',
                        action='store',
                        default=None,
                        help='Whether or not to plot saturated stars over top of the image')
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=None,
                        help='Aperture of photometry')
    parser.add_argument("--band", '-band',
                        action='store',
                        default=None,
                        help='Band of photometry')
    parser.add_argument("--cloud",
                        action='store',
                        default=None,
                        help='Band of photometry')
    args = parser.parse_args()

    if args.image is None:
        raise Exception, 'ERROR: No image file given'
    if args.sexfile is None:
        print 'Warning: No souce extractor file given'
    if args.outfile is None:
        print 'Warning: No output file given'
    if args.tmass is None:
        raise Exception, 'ERROR: No 2Mass file given'
    if args.band is None:
        raise Exception, 'ERROR: No band given'
    assert args.band in ['j', 'k', 'h'], 'ERROR: band must be "j", "k", or "h"'

    apertures = [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
    iap = np.argmin(np.array(apertures) - int(args.aperture)) + 1

    for i, cloud in enumerate(_CLOUDS):
        if cloud in args.cloud:
            center = _CENTERS[i]
            size = _SIZES[i]

    zeropoint(args.tmass, iap, args.band, args.image, args.sexfile, args.outfile, center, size)


def zeropoint(twomassfile, iap, band, image, sexfile, outfile=None, center=(None,None),
              size=(None,None)):

    # read in catalogue of unsaturated stars
    photfile = phot_curves.remove_saturated(image, sexfile)

    # Removed reddened stars from zeropoint calculation
    if center[0] is not None:
        idx_reddened = id_region(photfile['X_IMAGE'], photfile['Y_IMAGE'], center, size)  # Find reddened
        unreddened = photfile[~idx_reddened]
    else:
        unreddened = photfile

    # Remove objects which are not stars by PSF
    idx_stars = phot_curves.half_light(unreddened['MAG_APER_{}'.format(iap)], unreddened['FLUX_RADIUS'], toplot=False)
    stars = unreddened[idx_stars]

    # Compare magnitudes with same in 2Mass catalogue
    twomass, stars_tm = parse_2mass(twomassfile, stars)
    offset = compare_magnitudes(stars_tm['MAG_APER_{}'.format(iap)], twomass['{}_m'.format(band)],
                                stars_tm['MAGERR_APER_{}'.format(iap)], outfile)


def id_region(x, y, center, size):
    """
    For a given centerpoint and size of a rectangular region (in pixels), remove stars in the region from the list
    """

    y_min = center[1] - size[1] / 2.
    y_max = center[1] + size[1] / 2.
    x_min = center[0] - size[0] / 2.
    x_max = center[0] + size[0] / 2.

    # idxx = np.array(x[np.where(x > x_min)] < x_max)
    # idxy = np.array(y[np.where(y > y_min)] < y_max)
    # idx = np.intersect1d(np.arange(len(x))[idxx], np.arange(len(y))[idxy])

    idxx = np.arange(len(x))[(x_min < x) * (x < x_max)]
    idxy = np.arange(len(y))[(y_min < y) * (y < y_max)]
    idx = np.intersect1d(idxx, idxy)
    return idx


def parse_2mass(twomass_file, stars):
    """
    For a given list of coordinates returned by source extractor, find the same star in the
    2Mass catalog, and save the cropped catalog list to a file
    """

    print '  Identifying 2Mass catalog stars with stars in CFHT image'

    twomass = ascii.read(twomass_file)

    idx_stars = []
    idx_2mass = []
    for i in range(len(stars)):
        max_separation = 0.005
        imin, sep_min = query_separation((twomass['ra'], twomass['dec']), (stars['X_WORLD'][i], stars['Y_WORLD'][i]))
        if sep_min < max_separation:
            idx_stars.append(i)
            idx_2mass.append(imin)

    return twomass[idx_2mass], stars[idx_stars]


def query_separation((x1, y1), (x2, y2)):
    sep = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
    return np.argmin(sep), np.min(sep)


def compare_magnitudes(mag_sex, mag_tm, magerr, outfile=None):
    """
    Comparing magnitudes from source extractor and 2mass.
    Calculate the magnitude offset between the two to find correction to zeropoint in source extractor
    """

    wgt = (magerr) ** (-2)

    offset = mag_tm - mag_sex
    wavg_off = np.average(offset, weights=wgt)
    var_off = np.var(offset)
    sig_off = np.sqrt(np.sum(wgt))

    wavg_mag = np.average(mag_sex, weights=wgt)
    var_mag = np.var(mag_sex)

    print '{} & {} & {} & {} & {}'.format(wavg_off, var_off, sig_off, wavg_mag, var_mag)
    if outfile is not None:
        with open(outfile, 'a') as of:
            of.write(
                '& {} & {} & {} & {} & {} \\\\ \n'.format(wavg_off, var_off, sig_off, wavg_mag, var_mag))

    return wavg_off


if __name__ == '__main__':
    main()
