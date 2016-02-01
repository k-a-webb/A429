__author__ = 'kawebb'

"""
Match stars from J,H,Ks bands, plot colour magnitude diagrams

python plot_cmd.py -is CB68/CB68_J_sub.fits CB68/CB68_H_sub.fits CB68/CB68_Ks_sub.fits -sfs phot_t3/CB68_J_sex_t3_ap30.txt phot_t3/CB68_H_sex_t3_ap40.txt phot_t3/CB68_Ks_sex_t3_ap24.txt -ofile phot_t3/CB68_allmags.txt -ofig CB68/CB68_cmd.png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phot_curves
import os
from scipy.spatial import KDTree


def main():
    parser = argparse.ArgumentParser(
        description='For photometry output from sextractor for 3 bands, compile magnitude arrays into singe file')
    parser.add_argument("--sexfiles", '-sfs',
                        nargs='*',
                        action='store',
                        default=None,
                        help='Source extractor files for J,H,Ks bands.')
    parser.add_argument("--maxsep", '-ms',
                        action='store',
                        default=0.0001,
                        help='Maximum angular separation for nearest neighbour search.')
    parser.add_argument("--outfile", '-ofile',
                        action='store',
                        default=None,
                        help='Output file.')
    parser.add_argument("--outfig", '-ofig',
                        action='store',
                        default=None,
                        help='Output figure.')
    parser.add_argument("--images", '-is',
                        nargs='*',
                        action='store',
                        default=None,
                        help='Fits images for J,H,Ks bands.')
    parser.add_argument("--toplot",
                        action='store',
                        default=False,
                        help='Plot figures.')

    args = parser.parse_args()

    if args.sexfiles is None:
        raise Exception, 'No source extractor files given'
    if args.images is None:
        raise Exception, 'No image files given'
    if args.outfile is None:
        print 'Warning: No output file specified'
    if args.outfig is None:
        print 'Warning: No output figure specified'

    sexfile_j, sexfile_h, sexfile_ks = args.sexfiles
    fits_j, fits_h, fits_ks = args.images

    # if not os.path.exists(str(args.outfile)):

    # remove stars saturated in the fits images, marked with 0 center
    phot_j = loop_remove_satur(fits_j, sexfile_j, band='J')
    phot_ks = loop_remove_satur(fits_ks, sexfile_ks, band='Ks')
    phot_h = loop_remove_satur(fits_h, sexfile_h, band='H')

    # match the stars in each source extractor file by nearest neighbour
    mags = sort_by_coords_tree(phot_j, phot_ks, phot_h, float(args.maxsep))
    # mags.to_csv(args.outfile, index=False)
    # else:
    #     mags = pd.read_csv(args.outfile)

    if args.toplot:
        # plot to make sure matching the same points
        plt.scatter(mags['x_pix_J'].values, mags['y_pix_J'].values, color='k', label='J', marker='x')
        plt.scatter(mags['x_pix_Ks'].values, mags['y_pix_Ks'].values, color='r', label='Ks', marker='x')
        plt.scatter(mags['x_pix_H'].values, mags['y_pix_H'].values, color='g', label='H', marker='x')
        plt.legend()
        plt.show()

    # (J-H) vs Ks, and (H-Ks) vs Ks
    plot_cmd(mags['mag_aper_J'].values, mags['mag_aper_H'].values, mags['mag_aper_Ks'].values, args.outfig)


def sort_by_coords_tree(table_j, table_ks, table_h, max_sep=0.0001):
    """
    For every star found by the photometry in the J band, look for the same star in the Ks, H bands
    """

    idxs = query_tree3(table_j, table_ks, table_h, max_sep)

    mags = {'mag_aper_J': table_j['mag_aper_J'].values[idxs[0]], 'magerr_aper_J': table_j['magerr_aper_J'].values[idxs[0]],
            'mag_aper_Ks': table_ks['mag_aper_Ks'].values[idxs[1]], 'magerr_aper_Ks': table_ks['magerr_aper_Ks'].values[idxs[1]],
            'mag_aper_H': table_h['mag_aper_H'].values[idxs[2]], 'magerr_aper_H': table_h['magerr_aper_H'].values[idxs[2]]}

    return pd.DataFrame(data=mags)

    # idxsj, idxsks, idxsh = query_tree3(table_j, table_ks, table_h, max_sep)
    #
    # spliced_table_j = table_j[idxsj]  # apply True/False array to splice out matched indexes
    # spliced_table_j.reset_index(inplace=True, drop=True)
    # spliced_table_ks = table_ks[idxsks]
    # spliced_table_ks.reset_index(inplace=True, drop=True)
    # spliced_table_h = table_h[idxsh]
    # spliced_table_h.reset_index(inplace=True, drop=True)
    #
    # sorted_table_j = spliced_table_j.sort_values('x_wcs_J')  # sort both tables by the x axis
    # sorted_table_ks = spliced_table_ks.sort_values('x_wcs_Ks')
    # sorted_table_h = spliced_table_h.sort_values('x_wcs_H')
    #
    # plt.plot(sorted_table_j['x_pix_J'].values, color='k', label='J')
    # plt.plot(sorted_table_ks['x_pix_Ks'].values, color='r', label='Ks')
    # plt.plot(sorted_table_h['x_pix_H'].values, color='g', label='H')
    # plt.legend()
    # plt.show()
    #
    # plt.plot(sorted_table_j['y_pix_J'].values, color='k', label='J')
    # plt.plot(sorted_table_ks['y_pix_Ks'].values, color='r', label='Ks')
    # plt.plot(sorted_table_h['y_pix_H'].values, color='g', label='H')
    # plt.legend()
    # plt.show()
    #
    # sorted_table_ks.drop(['kron_radius', 'fwhm_image', 'fwhm_world', 'flags'], axis=1, inplace=True)
    # sorted_table_h.drop(['kron_radius', 'fwhm_image', 'fwhm_world', 'flags'], axis=1, inplace=True)
    #
    # return pd.concat((sorted_table_j, sorted_table_ks, sorted_table_h), axis=1)


def query_tree3(table_j, table_ks, table_h, max_sep=0.0001):
    """
    For two pandas database files, with the ra/dec header names x_wcs/y_wcs as defined in read_sex
    """

    tree2 = KDTree(zip(table_ks['x_wcs_Ks'].values.ravel(), table_ks['y_wcs_Ks'].values.ravel()))
    tree3 = KDTree(zip(table_h['x_wcs_H'].values.ravel(), table_h['y_wcs_H'].values.ravel()))

    idxs1 = []
    idxs2 = []
    idxs3 = []
    dups2 = []
    dups3 = []

    for i in range(len(table_j.index)):
        d2, icoords2 = tree2.query((table_j['x_wcs_J'].values[i], table_j['y_wcs_J'].values[i]), distance_upper_bound=max_sep)
        d3, icoords3 = tree3.query((table_j['x_wcs_J'].values[i], table_j['y_wcs_J'].values[i]), distance_upper_bound=max_sep)
        if (d2 != np.inf) & (d3 != np.inf):
            if icoords2 in idxs2:
                dups2.append(icoords2)
            if icoords3 in idxs3:
                dups3.append(icoords3)

            idxs1.append(i)
            idxs2.append(icoords2)
            idxs3.append(icoords3)

    if (len(dups2) > 0) | (len(dups3) > 0):
        raise Exception, 'WARNING: duplicates found'

    return np.array((idxs1, idxs2, idxs3))

    # tree_ks = KDTree(zip(table_ks['x_wcs_Ks'].values.ravel(), table_ks['y_wcs_Ks'].values.ravel()))
    # tree_h = KDTree(zip(table_h['x_wcs_H'].values.ravel(), table_h['y_wcs_H'].values.ravel()))
    #
    # idxs_j = np.zeros(len(table_j), dtype=bool)
    # idxs_ks = np.zeros(len(table_ks), dtype=bool)
    # idxs_h = np.zeros(len(table_h), dtype=bool)
    #
    # for i in range(len(table_j)):
    #     d2, i2 = tree_ks.query((table_j['x_wcs_J'].values[i], table_j['y_wcs_J'].values[i]),
    #                            distance_upper_bound=max_sep)
    #     d3, i3 = tree_h.query((table_j['x_wcs_J'].values[i], table_j['y_wcs_J'].values[i]),
    #                           distance_upper_bound=max_sep)
    #     if (d2 != np.inf) & (d3 != np.inf):
    #         if (idxs_j[i] is True) | (idxs_ks[i2] is True) | (idxs_h[i3] is True):
    #             print '  WARNING: duplicates found for index {}'.format(i)
    #         idxs_j[i] = True
    #         idxs_ks[i2] = True
    #         idxs_h[i3] = True
    #     else:
    #         print 'no match index {}'.format(i)
    #
    # if (np.count_nonzero(idxs_j) != np.count_nonzero(idxs_ks)) | (np.count_nonzero(idxs_j) != np.count_nonzero(idxs_h)):
    #     print 'Number of indexes dont match: {}, {}, {}'.format(np.count_nonzero(idxs_j), np.count_nonzero(idxs_ks),
    #                                                             np.count_nonzero(idxs_h))
    #
    # return idxs_j, idxs_ks, idxs_h


def loop_remove_satur(image, sexfile, band=None):
    imgdata, imghead = phot_curves.read_fits(image)
    phottable = read_sex_band(sexfile, band=band)

    ix = phot_curves.detect_satur(imgdata, phottable['x_pix_{}'.format(band)].values,
                                  phottable['y_pix_{}'.format(band)].values, phottable['kron_radius'].values)
    unsat_phottable = phottable[ix]

    return unsat_phottable


def read_sex_band(sfile, skiplines=12, band=None):
    #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
    #   2 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
    #   3 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   4 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
    #   5 KRON_RADIUS            Kron apertures in units of A or B
    #   6 X_IMAGE                Object position along x                                    [pixel]
    #   7 Y_IMAGE                Object position along y                                    [pixel]
    #   8 X_WORLD                Barycenter position along world x axis                     [deg]
    #   9 Y_WORLD                Barycenter position along world y axis                     [deg]
    #  10 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
    #  11 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
    #  12 FLAGS                  Extraction flags

    if band is None:
        names = ['flux_aper', 'fluxerr_aper', 'mag_aper', 'magerr_aper', 'kron_radius', 'x_pix', 'y_pix', 'x_wcs',
                 'y_wcs', 'fwhm_image', 'fwhm_world', 'flags']
    else:
        names = ['flux_aper_{}'.format(band), 'fluxerr_aper_{}'.format(band), 'mag_aper_{}'.format(band),
                 'magerr_aper_{}'.format(band), 'kron_radius', 'x_pix_{}'.format(band), 'y_pix_{}'.format(band),
                 'x_wcs_{}'.format(band), 'y_wcs_{}'.format(band), 'fwhm_image', 'fwhm_world', 'flags']
    # read in file as pandas dataframe
    table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python', names=names)

    # drop all rows with saturated stars
    table = table[table[names[2]] != 99.]
    table.reset_index()

    return table


def plot_cmd(mag_j, mag_h, mag_ks, outfig):
    """
    Plot colour magnitude diagrams for (J-H) v. Ks, (H-Ks) v. Ks
    """

    fig, ax = plt.subplots(nrows=1, ncols=2)

    # ax[0] = plt.gca()
    ax[0].scatter(mag_j - mag_h, mag_ks, marker='.')
    # ax[0].errorbar(apertures, average_rel, yerr=yerr)

    ax[0].set_ylabel('Ks')
    ax[0].set_xlabel('J-H')

    # ax[1] = plt.gca()
    ax[1].scatter(mag_h - mag_ks, mag_ks, marker='.')
    # ax[1].errorbar(apertures, average_rel, yerr=yerr)

    ax[1].set_ylabel('Ks')
    ax[1].set_xlabel('H-Ks')

    ax[0].invert_yaxis()
    ax[1].invert_yaxis()

    if outfig is not None:
        plt.savefig(outfig)
    plt.show()


if __name__ == '__main__':
    main()
