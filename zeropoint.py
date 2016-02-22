__author__ = 'kawebb'

"""
Procedure:
Determine the magnitude zeropoint of our CFHT observations by comparing bright
unsaturated stars to the same stars in the 2mass catalog.
Apply this zeropoint to all of the stars measured.

python zeropoint.py -img CB68/CB68_J_sub.fits -ph phot_t20/CB68_J_mags_t20.txt -band j -ap 30 -tm CB68/2mass_CB68.tbl -out CB68/CB68_offsets.txt  -sf phot_t20/CB68_J_sex_t20.txt
python zeropoint.py -img CB68/CB68_H_sub.fits -ph phot_t20/CB68_H_mags_t20.txt -band h -ap 30 -tm CB68/2mass_CB68.tbl -out CB68/CB68_offsets.txt -sf phot_t20/CB68_H_sex_t20.txt
python zeropoint.py -img CB68/CB68_Ks_sub.fits -ph phot_t20/CB68_Ks_mags_t20.txt -band k -ap 24 -tm CB68/2mass_CB68.tbl -out CB68/CB68_offsets.txt -sf phot_t20/CB68_Ks_sex_t20.txt

sex -c phot_t3.sex ../release/CB68_J_sub.fits -CATALOG_NAME CB68_J_sex_t3_ap30.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 30-0.192943311392
sex -c phot_t3.sex ../release/CB68_Ks_sub.fits -CATALOG_NAME CB68_Ks_sex_t3_ap24.txt -PHOT_APERTURES 24 -MAG_ZEROPOINT 30+0.605746118244
sex -c phot_t3.sex ../release/CB68_H_sub.fits -CATALOG_NAME CB68_H_sex_t3_ap40.txt -PHOT_APERTURES 40 -MAG_ZEROPOINT 30+0.471278827929

"""

import os
import numpy as np
import argparse
from astropy.io import ascii
import phot_curves


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
    parser.add_argument("--photfile", '-ph',
                        action='store',
                        default=None,
                        help='File of photometry')
    parser.add_argument("--band", '-band',
                        action='store',
                        default=None,
                        help='Band of photometry')
    args = parser.parse_args()

    if args.image is None:
        raise Exception, 'ERROR: No image file given'
    if args.sexfile is None:
        print 'Warning: No souce extractor file given'
    if args.outfile is None:
        print 'Warning: No out file file given'
    if args.photfile is None and args.sexfile is None:
        raise Exception, 'ERROR: No souce extractor file or photometry file given'
    if args.tmass is None:
        raise Exception, 'ERROR: No 2Mass file given'
    if args.band is None:
        raise Exception, 'ERROR: No band given'

    apertures = [5, 10, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
    iap = np.argmin(np.array(apertures) - int(args.aperture)) + 1

    zeropoint(args.image, args.sexfile, args.outfile, args.photfile, args.tmass, iap, args.band)


def zeropoint(image, sexfile, outfile, photfile, twomassfile, iap, band):
    clouds = ['CB68', 'L429', 'L1521E', 'L1544', 'L1552']
    centers = np.array([[2610.6899, 2868.1778], [2496.3232, 2158.1909], [2532.6025, 2753.7333], [2345.069, 2855.932],
                        [2710.337, 2593.2019]])  # in pixels
    sizes = np.array([[1382.2207, 1227.3661], [1869.7259, 2345.7605], [1880.5105, 1777.3177], [1788.7782, 1570.9142],
                      [1639.7134, 177.3117]])  # in pixels
    for i, cloud in enumerate(clouds):
        if cloud in image:
            center = centers[i]
            size = sizes[i]

    # read in catalogue of unsaturated stars
    if not os.path.exists(photfile):
        photfile = phot_curves.remove_saturated(image, sexfile, photfile)
    else:
        photfile = ascii.read(photfile)

    # Removed reddened stars from zeropoint calculation
    idx_reddened = id_region(photfile['X_IMAGE'], photfile['Y_IMAGE'], center, size)  # Find reddened
    unreddened = photfile[~idx_reddened]

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

    print '{} {} {} {} {}'.format(wavg_off, var_off, sig_off, wavg_mag, var_mag)
    if outfile is not None:
        with open(outfile, 'a') as of:
            of.write(
                '& {} & {} & {} & {} & {} \\\\ \n'.format(wavg_off, var_off, sig_off, wavg_mag, var_mag))

    return wavg_off


if __name__ == '__main__':
    main()
