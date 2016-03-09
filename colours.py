__author__ = 'kawebb'

"""
Procedure:
Match stars from J,H,Ks bands, plot colour magnitude diagrams

Input consists of:
  a list of source extractor files corresponding to the J,H,K bands
  the fits images in the J,H,K bands
  the name of the photometry file to be output (this will be the input into NICEST)
  the name of the colour colour figure to be saved, by default not saved

  Note: if the output photometry file has already been produced, the script can be run without specifying
        the source extractor or images file, as long as the output photometry file is specified

Optional input:
  maximum separation for nearest neighbour search, default 0.0001 degrees

To produce the colour colour plots:
  1. remove saturated stars from each photometry catalogue using phot_curves function
  2. match the stars in each catalogue, only keep stars identified in all bands, uses nearest neighbour approach
  3. compare colours, and plot
  4. if a particular region is specifed (i.e. written directly into the script) the magnitude list can be cropped easily

Example commands to run:
python colours.py -is CB68/CB68_J_sub.fits CB68/CB68_Ks_sub.fits CB68/CB68_H_sub.fits -sfs phot_t3/CB68_J_sex_t3_ap30.txt phot_t3/CB68_Ks_sex_t3_ap24.txt phot_t3/CB68_H_sex_t3_ap30.txt -ofile phot_t3/CB68_allmags.txt -ofig CB68/CB68_cmd.png
python colours.py -is L429/L429_J_sub.fits L429/L429_Ks_sub.fits L429/L429_H_sub.fits -sfs phot_t3/L429_J_sex_t3_ap24.txt phot_t3/L429_Ks_sex_t3_ap24.txt phot_t3/L429_H_sex_t3_ap20.txt -ofile phot_t3/L429_allmags.txt -ofig L429/L429_cmd.png
python colours.py -is L1521E/L1521E_J_sub.fits L1521E/L1521E_Ks_sub.fits L1521E/L1521E_H_sub.fits -sfs phot_t3/L1521E_J_sex_t3_ap30.txt phot_t3/L1521E_Ks_sex_t3_ap28.txt phot_t3/L1521E_H_sex_t3_ap26.txt -ofile phot_t3/L1521E_allmags.txt -ofig L1521E/L1521E_cmd.png
python colours.py -is L1544/L1544_J_sub.fits L1544/L1544_Ks_sub.fits L1544/L1544_H_sub.fits -sfs phot_t3/L1544_J_sex_t3_ap28.txt phot_t3/L1544_Ks_sex_t3_ap28.txt phot_t3/L1544_H_sex_t3_ap28.txt -ofile phot_t3/L1544_allmags.txt -ofig L1544/L1544_cmd.png
python colours.py -is L1552/L1544_J_sub.fits L1552/L1552_Ks_sub.fits L1552/L1552_H_sub.fits -sfs phot_t3/L1552_J_sex_t3_ap28.txt phot_t3/L1552_Ks_sex_t3_ap28.txt phot_t3/L1552_H_sex_t3_ap28.txt -ofile phot_t3/L1552_allmags.txt -ofig L1552/L1552_cmd.png

"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import phot_curves
import os
from scipy.spatial import KDTree
from astropy.io import ascii


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

        i_j = phot_curves.half_light(phot_j['MAG_APER'], phot_j['FLUX_RADIUS'], toplot=False)
        phot_j = phot_j[i_j]
        i_k = phot_curves.half_light(phot_ks['MAG_APER'], phot_ks['FLUX_RADIUS'], toplot=False)
        phot_ks = phot_ks[i_k]
        i_h = phot_curves.half_light(phot_h['MAG_APER'], phot_h['FLUX_RADIUS'], toplot=False)
        phot_h = phot_h[i_h]

        # match the stars in each source extractor file by nearest neighbour
        mags = sort_by_coords_tree(phot_j, phot_ks, phot_h, float(args.maxsep), args.toplot)
        mags.to_csv(args.outfile, index=False, delimater=' ')
    else:
        print 'Reading magnitude table from file {}'.format(args.outfile)
        mags = ascii.read(args.outfile)

    outfig_cc = args.outfig.split('.')[0] + '_cc.png'
    mags_colours = [mags['mag_aper_J'].data, mags['mag_aper_H'].data, mags['mag_aper_Ks'].data]
    magerrs_colours = [mags['magerr_aper_J'].data, mags['magerr_aper_H'].data, mags['magerr_aper_Ks'].data]
    plot_cmd(mags_colours, magerrs_colours, args.outfig)
    plot_colourcolour(mags_colours, outfig_cc)


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

    tree2 = KDTree(zip(table_ks['X_WORLD'].ravel(), table_ks['Y_WORLD'].ravel()))
    tree3 = KDTree(zip(table_h['X_WORLD'].ravel(), table_h['Y_WORLD'].ravel()))

    idxs1 = []
    idxs2 = []
    idxs3 = []
    dups2 = []
    dups3 = []

    for i in range(len(table_j)):
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


def plot_colourcolour(mags, outfig=None, colour='k', xlim=(None, None), ylim=(None, None), alpha=1):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))

    # J - H vs H - Ks
    ax.scatter(mags[1] - mags[2], mags[0] - mags[1], marker='.', color=colour, alpha=alpha)
    # ax[0].errorbar(mags[0] - mags[1], mags[2], yerr=magerrs[2], xerr=xerr, ecolor='g')
    ax.set_xlabel('H - K')
    ax.set_ylabel('J - H')

    ax.invert_yaxis()
    ax.invert_xaxis()

    if xlim[0] is not None:
        ax.set_xlim(xlim[0], xlim[1])
    if ylim[0] is not None:
        ax.set_ylim(ylim[0], ylim[1])

    if outfig is not None:
        plt.savefig(outfig)
    plt.show()


def select_region(mags, center, size):
    in_x = mags[(center[0] - size[0] / 2. < mags['x_pix']) & (mags['x_pix'] < center[0] + size[0] / 2.)]
    in_xy = in_x[(center[1] - size[1] / 2. < in_x['y_pix']) & (in_x['y_pix'] < center[1] + size[1] / 2.)]
    return in_xy


def compare_colourcolour(mags, center, size, outcenter, outsize, outfig=None):
    inmags = select_region(mags, center, size)
    outmags = select_region(mags, outcenter, outsize)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5, 5), sharex=True, sharey=True)
    # J - H vs H - Ks
    incloud = ax[0].scatter(inmags['mag_aper_H'] - inmags['mag_aper_Ks'], inmags['mag_aper_J'] - inmags['mag_aper_H'],
                            marker='.', color='r')
    outcloud = ax[1].scatter(outmags['mag_aper_H'] - outmags['mag_aper_Ks'],
                             outmags['mag_aper_J'] - outmags['mag_aper_H'], marker='.', color='b')
    ax[0].set_xlabel('H - K')
    ax[1].set_xlabel('H - K')
    ax[0].set_ylabel('J - H')
    ax[0].set_xlim(-4, 4)
    ax[0].set_ylim(-3, 4)

    ax[0].invert_yaxis()
    ax[0].invert_xaxis()

    fig.subplots_adjust(hspace=0, wspace=0.1)
    lgd = fig.legend((incloud, outcloud), ('Reddened', 'Reference'), loc=(0.25, .818))
    if outfig is not None:
        plt.savefig(outfig, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
