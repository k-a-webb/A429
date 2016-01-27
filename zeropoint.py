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
import match_phot

TEMP_FILE = '/Users/kawebb/a429/temp.txt'


def main():
    parser = argparse.ArgumentParser(
        description='For a given image, and list of star positions, return list of stars which are not saturated')
    parser.add_argument("--image", '-img',
                        action="store",
                        default='/Users/kawebb/a429/CB68/CB68_H_sub.fits',
                        help="Fits image of star field")
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default='/Users/kawebb/a429/CB68/CB68_H_sex_2m.txt',
                        help='List of star positions as output by sextractor for 2mass comparison.')
    parser.add_argument("--tmass", '-tm',
                        action='store',
                        default='/Users/kawebb/a429/CB68/2mass_CB68.tbl',
                        help='File of 2Mass stars in region of interest')
    parser.add_argument("--toplot", '-plt',
                        action='store',
                        default=False,
                        help='Whether or not to plot saturated stars over top of the image')
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=None,
                        help='Aperture of photometry')
    args = parser.parse_args()

    unsat_list = remove_satur(args.image, args.sexfile, args.toplot)  # create list of unsatured stars identified by sex
    twomass_list, star_list = parse_2mass(args.tmass, unsat_list)
    offset = compare_magnitudes(star_list, twomass_list, args.toplot, args.aperture, TEMP_FILE)

    # plot_plots(args.image, uslist, tmstars, args.toplot)


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


def compare_magnitudes(uslist, tmstars, toplot=False, aperture=None, outfile=None):
    """
    Comparing magnitudes from source extractor and 2mass.
    Calculate the magnitude offset between the two to find correction to zeropoint in source extractor
    """

    print '  Comparing colours of 2Mass stars and stars in the image'
    tmmag = tmstars.j_m.values
    sexmag = uslist.mag_aper.values
    sexmagerr = uslist.magerr_aper.values

    wgt = (sexmagerr) ** (-2)

    # Weighting as function of magnitude
    # df = pd.DataFrame(data={'mag':sexmag, 'weight':wgt})
    # dfs = df.sort('mag')
    # plt.plot(dfs.mag.values, dfs.weight.values)
    # plt.show()

    try:
        offset = tmmag - sexmag
        wavg = np.average(offset, weights=wgt)
        sig = np.var(offset)
    except:
        offset = 0. * tmmag
        wavg = 0.
        sig = 0.

    print '{} {} {} {} {} {}'.format(aperture, wavg, np.mean(offset), sig, np.average(sexmag, weights=wgt),
                                     np.mean(sexmag))

    if outfile is not None:
        with open(outfile, 'a') as of:
            of.write('{} {} {} {} {} {}\n'.format(aperture, wavg, np.mean(offset), sig, np.average(sexmag, weights=wgt),
                                                  np.mean(sexmag)))

    return wavg


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
