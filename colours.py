__author__ = 'kawebb'

"""
Match stars from J,H,Ks bands, plot colour magnitude diagrams

python plot_cmd.py -is CB68/CB68_J_sub.fits CB68/CB68_Ks_sub.fits CB68/CB68_H_sub.fits -sfs phot_t3/CB68_J_sex_t3_ap30.txt phot_t3/CB68_Ks_sex_t3_ap24.txt phot_t3/CB68_H_sex_t3_ap30.txt -ofile phot_t3/CB68_allmags.txt -ofig CB68/CB68_cmd.png
python plot_cmd.py -is L429/L429_J_sub.fits L429/L429_Ks_sub.fits L429/L429_H_sub.fits -sfs phot_t3/L429_J_sex_t3_ap24.txt phot_t3/L429_Ks_sex_t3_ap24.txt phot_t3/L429_H_sex_t3_ap20.txt -ofile phot_t3/L429_allmags.txt -ofig L429/L429_cmd.png
python plot_cmd.py -is L1521E/L1521E_J_sub.fits L1521E/L1521E_Ks_sub.fits L1521E/L1521E_H_sub.fits -sfs phot_t3/L1521E_J_sex_t3_ap30.txt phot_t3/L1521E_Ks_sex_t3_ap28.txt phot_t3/L1521E_H_sex_t3_ap26.txt -ofile phot_t3/L1521E_allmags.txt -ofig L1521E/L1521E_cmd.png
python plot_cmd.py -is L1544/L1544_J_sub.fits L1544/L1544_Ks_sub.fits L1544/L1544_H_sub.fits -sfs phot_t3/L1544_J_sex_t3_ap28.txt phot_t3/L1544_Ks_sex_t3_ap28.txt phot_t3/L1544_H_sex_t3_ap28.txt -ofile phot_t3/L1544_allmags.txt -ofig L1544/L1544_cmd.png

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

    if args.outfig is None:
        print 'Warning: No output figure specified'

    if not os.path.exists(str(args.outfile)):
        if args.sexfiles is None:
            raise Exception, 'No source extractor files given'
        if args.images is None:
            raise Exception, 'No image files given'
        if args.outfile is None:
            print 'Warning: No output file specified'

        sexfile_j, sexfile_ks, sexfile_h = args.sexfiles
        fits_j, fits_ks, fits_h = args.images

        # remove stars saturated in the fits images, marked with 0 center
        phot_j = phot_curves.remove_saturated(fits_j, sexfile_j)
        phot_ks = phot_curves.remove_saturated(fits_ks, sexfile_ks)
        phot_h = phot_curves.remove_saturated(fits_h, sexfile_h)

        # match the stars in each source extractor file by nearest neighbour
        mags = sort_by_coords_tree(phot_j, phot_ks, phot_h, float(args.maxsep), args.toplot)
        mags.to_csv(args.outfile, index=False, delimater=' ')

    else:
        print 'Reading magnitude table from file {}'.format(args.outfile)
        mags = pd.read_csv(args.outfile)

    outfig_cc = args.outfig.split('.')[0] + '_cc.png'
    plot_cmd([mags['mag_aper_J'].values, mags['mag_aper_H'].values, mags['mag_aper_Ks'].values],
             [mags['magerr_aper_J'].values, mags['magerr_aper_H'].values, mags['magerr_aper_Ks'].values], args.outfig)
    plot_colourcolour([mags['mag_aper_J'].values, mags['mag_aper_H'].values, mags['mag_aper_Ks'].values], outfig_cc)


    # plot colour-colour plots for reddened and unreddened regions
    # centers = np.array([[2610.6899, 2868.1778], [2496.3232, 2158.1909], [2532.6025, 2753.7333], [2345.069, 2855.932],
    #                     [2710.337, 2593.2019]])  # in pixels
    # sizes = np.array([[1382.2207, 1227.3661], [1869.7259, 2345.7605], [1880.5105, 1777.3177], [1788.7782, 1570.9142],
    #                   [1639.7134, 177.3117]])  # in pixels
    # out_centers = np.array([[4242, 1782], [0, 0]])
    # out_sizes = np.array([[1176, 1024], [0, 0]])
    # compare_colourcolour(mags, centers[0], sizes[0], out_centers[0], out_sizes[0], 'CB68/CB68_cc_compare.png')


def sort_by_coords_tree(table_j, table_ks, table_h, max_sep=0.0001, toplot=False):
    """
    For every star found by the photometry in the J band, look for the same star in the Ks, H bands
    """

    idxs = query_tree3(table_j, table_ks, table_h, max_sep)

    if toplot:
        fig, ax = plt.subplots(1, 2)
        ax[0].plot(table_j['X_IMAGE'][idxs[0]], color='k', label='J')
        ax[0].plot(table_ks['X_IMAGE'][idxs[1]], color='r', label='Ks')
        ax[0].plot(table_h['X_IMAGE'][idxs[2]], color='g', label='H')

        ax[1].plot(table_j['Y_IMAGE'][idxs[0]], color='k', label='J')
        ax[1].plot(table_ks['Y_IMAGE'][idxs[1]], color='r', label='Ks')
        ax[1].plot(table_h['Y_IMAGE'][idxs[2]], color='g', label='H')
        plt.legend()
        plt.show()

    if toplot:
        # plot to make sure matching the same points
        plt.scatter(table_j['X_IMAGE'][idxs[0]], table_j['Y_IMAGE'][idxs[0]], color='k', label='J', marker='x')
        plt.scatter(table_ks['X_IMAGE'][idxs[1]], table_ks['Y_IMAGE'][idxs[1]], color='r', label='Ks', marker='x')
        plt.scatter(table_h['X_IMAGE'][idxs[2]], table_h['Y_IMAGE'][idxs[2]], color='g', label='H', marker='x')
        plt.legend()
        plt.show()

    mags = {'x_pix': table_j['X_IMAGE'][idxs[0]], 'y_pix': table_j['Y_IMAGE'][idxs[0]],
            'x_wcs': table_j['X_WORLD'][idxs[0]], 'y_wcs': table_j['Y_WORLD'][idxs[0]],
            'mag_aper_J': table_j['MAG_APER'][idxs[0]], 'magerr_aper_J': table_j['MAGERR_APER'][idxs[0]],
            'mag_aper_Ks': table_ks['MAG_APER'][idxs[1]], 'magerr_aper_Ks': table_ks['MAGERR_APER'][idxs[1]],
            'mag_aper_H': table_h['MAG_APER'][idxs[2]], 'magerr_aper_H': table_h['MAGERR_APER'][idxs[2]]}

    return pd.DataFrame(data=mags)


def query_tree3(table_j, table_ks, table_h, max_sep=0.0001):
    """
    For two pandas database files, with the ra/dec header names x_wcs/y_wcs as defined in read_sex
    """

    tree2 = KDTree(zip(table_ks['X_WORLD'].values.ravel(), table_ks['Y_WORLD'].values.ravel()))
    tree3 = KDTree(zip(table_h['X_WORLD'].values.ravel(), table_h['Y_WORLD'].values.ravel()))

    idxs1 = []
    idxs2 = []
    idxs3 = []
    dups2 = []
    dups3 = []

    for i in range(len(table_j.index)):
        d2, icoords2 = tree2.query((table_j['X_WORLD'][i], table_j['Y_WORLD'][i]), distance_upper_bound=max_sep)
        d3, icoords3 = tree3.query((table_j['X_WORLD'][i], table_j['Y_WORLD'][i]), distance_upper_bound=max_sep)
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


def plot_cmd(mags, magerrs, outfig=None):
    """
    Plot colour magnitude diagrams for (J-H) v. Ks, (H-Ks) v. Ks
    mags array of magnitudes J, H, Ks in that order
    same oder for magerrs
    """

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))

    # J - H vs Ks
    # xerr = np.sqrt((magerrs[0]/mags[0])**2 + (magerrs[1]/mags[1])**2)*(mags[0] - mags[1])
    ax[0].scatter(mags[0] - mags[1], mags[2], marker='.')
    # ax[0].errorbar(mags[0] - mags[1], mags[2], yerr=magerrs[2], xerr=xerr, ecolor='g')
    ax[0].set_ylabel('Ks')
    ax[0].set_xlabel('J-H')

    # H - Ks vs Ks
    # xerr = np.sqrt((magerrs[2]/mags[2])**2 + (magerrs[1]/mags[1])**2)*(mags[1] - mags[2])
    ax[1].scatter(mags[1] - mags[2], mags[2], marker='.')
    # ax[1].errorbar(mags[1] - mags[2], mags[2], yerr=magerrs[2], xerr=xerr, ecolor='g')
    ax[1].set_ylabel('Ks')
    ax[1].set_xlabel('H-Ks')

    ax[0].invert_yaxis()
    ax[1].invert_yaxis()

    if outfig is not None:
        plt.savefig(outfig)
    plt.show()


def plot_colourcolour(mags, outfig=None):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

    # J - H vs H - Ks
    ax.scatter(mags[1] - mags[2], mags[0] - mags[1], marker='.')
    # ax[0].errorbar(mags[0] - mags[1], mags[2], yerr=magerrs[2], xerr=xerr, ecolor='g')
    ax.set_xlabel('H - K')
    ax.set_ylabel('J - H')

    ax.invert_yaxis()
    ax.invert_xaxis()

    if outfig is not None:
        plt.savefig(outfig)
    plt.show()


def compare_colourcolour(mags, center, size, outcenter, outsize, outfig=None):
    in_x = mags.query('{} < x_pix < {}'.format(center[0] - size[0] / 2, center[0] + size[0] / 2))
    in_xy = in_x.query('{} < y_pix < {}'.format(center[1] - size[1] / 2, center[1] + size[1] / 2))

    inmags = [in_xy['mag_aper_J'].values, in_xy['mag_aper_H'].values, in_xy['mag_aper_Ks'].values]

    in_x = mags.query('{} < x_pix < {}'.format(outcenter[0] - outsize[0] / 2, outcenter[0] + outsize[0] / 2))
    in_xy = in_x.query('{} < y_pix < {}'.format(outcenter[1] - outsize[1] / 2, outcenter[1] + outsize[1] / 2))
    outmags = [in_xy['mag_aper_J'].values, in_xy['mag_aper_H'].values, in_xy['mag_aper_Ks'].values]

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5, 5), sharex=True, sharey=True)
    # J - H vs H - Ks
    incloud = ax[0].scatter(inmags[1] - inmags[2], inmags[0] - inmags[1], marker='.', color='r')
    outcloud = ax[1].scatter(outmags[1] - outmags[2], outmags[0] - outmags[1], marker='.', color='b')
    ax[0].set_xlabel('H - K')
    ax[1].set_xlabel('H - K')
    ax[0].set_ylabel('J - H')
    ax[0].set_xlim(-4, 4)
    ax[0].set_ylim(-3, 4)

    ax[0].invert_yaxis()
    ax[0].invert_xaxis()

    fig.subplots_adjust(hspace=0, wspace=0.1)
    lgd = fig.legend((incloud, outcloud), ('Reddened', 'Reference'), loc=(0.588, .878))
    plt.savefig(outfig, bbox_extra_artists=(lgd,), bbox_inches='tight')
    # plt.show()


if __name__ == '__main__':
    main()
