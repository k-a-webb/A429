__author__ = 'kawebb'

"""
Procedure:
Determine the magnitude zeropoint of our CFHT observations by comparing bright
unsaturated stars to the same stars in the 2mass catalog.
Apply this zeropoint to all of the stars measured.

python zeropoint.py -img CB68/CB68_J_sub.fits -ph CB68/CB68_J_mags_t20.txt -ap 30 -tm CB68/2mass_CB68.tbl -out CB68/CB68_offsets.txt
python zeropoint.py -img CB68/CB68_H_sub.fits -ph CB68/CB68_H_mags_t20.txt -ap 30 -tm CB68/2mass_CB68.tbl -out CB68/CB68_offsets.txt
python zeropoint.py -img CB68/CB68_Ks_sub.fits -ph CB68/CB68_Ks_mags_t20.txt -ap 24 -tm CB68/2mass_CB68.tbl -out CB68/CB68_offsets.txt
"""

import os
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import match_phot

TEMP_FILE = '/Users/kawebb/a429/temp.txt'


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

    clouds = ['CB68', 'L429', 'L1521E', 'L1544', 'L1552']
    centers = np.array([[2610.6899, 2868.1778], [2496.3232, 2158.1909], [2532.6025, 2753.7333], [2345.069, 2855.932],
                        [2710.337, 2593.2019]])  # in pixels
    sizes = np.array([[1382.2207, 1227.3661], [1869.7259, 2345.7605], [1880.5105, 1777.3177], [1788.7782, 1570.9142],
                      [1639.7134, 177.3117]])  # in pixels
    for i, cloud in enumerate(clouds):
        if cloud in args.image:
            center = centers[i]
            size = sizes[i]

    if args.photfile is not None:
        unsat_list = pd.read_csv(args.photfile)
    else:
        unsat_list = remove_satur(args.image, args.sexfile,
                                  args.toplot)  # create list of unsatured stars identified by sex

    bkg_list = remove_region(unsat_list, center, size)  # Remove stars reddened by cloud from zeropoint calculation
    twomass_list, star_list = parse_2mass(args.tmass, bkg_list)
    offset = compare_magnitudes(star_list, twomass_list, args.aperture, args.outfile)


def compare_magnitudes(uslist, tmstars, aperture=None, outfile=None):
    """
    Comparing magnitudes from source extractor and 2mass.
    Calculate the magnitude offset between the two to find correction to zeropoint in source extractor
    """

    tmmag = tmstars.j_m.values
    sexmag = uslist['mag_aper_{}'.format(aperture)].values
    sexmagerr = uslist['magerr_aper_{}'.format(aperture)].values

    wgt = (sexmagerr) ** (-2)

    offset = tmmag - sexmag
    wavg_off = np.average(offset, weights=wgt)
    var_off = np.var(offset)
    sig_off = np.sqrt(np.sum(wgt))

    wavg_mag = np.average(sexmag, weights=wgt)
    var_mag = np.var(sexmag)

    print '{} {} {} {} {} {}'.format(aperture, wavg_off, var_off, sig_off, wavg_mag, var_mag)
    if outfile is not None:
        with open(outfile, 'a') as of:
            of.write('{} {} {} {} {} {}\n'.format(aperture, wavg_off, var_off, sig_off, wavg_mag, var_mag))

    return wavg_off


def parse_2mass(twomass_file, star_list):
    """
    For a given list of coordinates returned by source extractor, find the same star in the
    2Mass catalog, and save the cropped catalog list to a file
    """

    print '  Identifying 2Mass catalog stars with Unsaturated stars in CFHT image'

    twomass_list = read_2mass(twomass_file)

    idx_stars = []
    idx_2mass = []
    for ii, i in enumerate(star_list.index):
        max_separation = 0.005
        imin, sep_min = query_separation(twomass_list, star_list.x_wcs[i], star_list.y_wcs[i])
        if sep_min < max_separation:
            idx_stars.append(i)
            idx_2mass.append(imin)
            # else:
            #     print '  {},{}  {},{} {}'.format(star_list.x_wcs[i], star_list.y_wcs[i], twomass_list['ra'][i_min],
            #                                      twomass_list['dec'][i_min], sep_min)

    return twomass_list.loc[idx_2mass, :], star_list.loc[idx_stars, :]


def query_separation(table, ra, dec):
    sep = np.sqrt((ra - table.ra.values) ** 2 + (dec - table.dec.values) ** 2)
    return np.argmin(sep), np.min(sep)


def remove_region(starlist, center, size):
    """
    For a given centerpoint and size of a rectangular region (in pixels), remove stars in the region from the list
    """

    in_x = starlist.query('{} < x_pix < {}'.format(center[0] - size[0] / 2, center[0] + size[0] / 2))
    in_xy = in_x.query('{} < y_pix < {}'.format(center[1] - size[1] / 2, center[1] + size[1] / 2))

    starlist = starlist.drop(starlist.index[in_xy.index.values])
    starlist.reset_index(inplace=True, drop=True)  # drop index column (inherent to pandas dataframe) as out of order

    return starlist


def remove_satur(ffile, sfile, toplot=False):
    """
    Detect saturated stars in the image, and remove them from the list returned by source extractor
    In our images, saturated stars have low (0) values at the center rather than a very high value
    and are therefore not identified by source extractor automatticaly.
    """

    stable = match_phot.read_sex(sfile)
    fdata, fhdr = match_phot.read_fits(ffile)

    # Iterate through stars
    ix = match_phot.detect_satur(fdata, stable['x_pix'].values, stable['y_pix'].values, stable['kron_radius'].values)
    unsat_stable = stable[ix]

    if toplot:
        plt.imshow(fdata, origin='lower', norm=matplotlib.colors.LogNorm())
        sat_stable = unsat_stable[np.invert(ix)]
        plt.plot(sat_stable.x.values, sat_stable.y.values, 'x')
        plt.title('Saturated stars')
        plt.show()

    return unsat_stable


def read_2mass(tmfile, skiplines=55):
    assert os.path.exists(tmfile), 'ERROR: 2Mass file {} does not exist'.format(tmfile)
    # ra, dec, j_m, j_cmsig, j_msigcom, j_snr, h_m, h_cmsig, h_msigcom, h_snr, k_m, k_cmsig, k_msigcom, k_snr
    #   rd_flg, dist, angle, j_h, h_k, j_k
    table = pd.read_csv(tmfile, skiprows=skiplines, sep=r"\s*", engine='python',
                        names=['ra', 'dec', 'j_m', 'j_cmsig', 'j_msigcom', 'j_snr', 'h_m', 'h_cmsig', 'h_msigcom',
                               'h_snr', 'k_m', 'k_cmsig', 'k_msigcom', 'k_snr', 'rd_flg', 'dist', 'angle', 'j_h',
                               'h_k', 'j_k'])
    return table


if __name__ == '__main__':
    main()
