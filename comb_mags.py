__author__ = 'kawebb'

import numpy as np
import argparse
import pandas as pd
import os

import matplotlib.pyplot as plt
from scipy.spatial import KDTree
import zeropoint


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

    if not os.path.exists(str(args.outfile)):


        # remove stars saturated in the fits images, marked with 0 center
        phot_j = zeropoint.remove_satur(fits_j, sexfile_j)
        phot_ks = zeropoint.remove_satur(fits_ks, sexfile_ks)
        phot_h = zeropoint.remove_satur(fits_h, sexfile_h)

        # match the stars in each source extractor file by nearest neighbour
        mags = sort_by_coords_tree(phot_j, phot_ks, phot_h, args.maxsep)

        # remove stars saturated in any of the source extractor files, marked with mag 99.0
        mags = mags[mags.mag_j != 99.]
        mags = mags[mags.mag_ks != 99.]
        mags = mags[mags.mag_h != 99.]

        mags.to_csv(args.outfile, index=False, sep=' ')
    else:
        mags = pd.read_csv(args.outfile, skiprows=1, sep=r"\s*", engine='python',
                           names=['mag_h', 'mag_j', 'mag_ks', 'magerr_h', 'magerr_j', 'magerr_ks'])

    # (J-H) vs Ks, and (H-Ks) vs Ks
    plot_cmd(mags.mag_j.values, mags.mag_h.values, mags.mag_ks.values, args.outfig)


def sort_by_coords_tree(table_j, table_ks, table_h, max_sep=0.0001):
    """
    For every star found by the photometry in the J band, look for the same star in the Ks, H bands
    """

    idxs, dups = query_separation_tree(table_j, table_ks, table_h, max_sep)

    mags = {'mag_j': table_j.mag_aper.values[idxs[0]], 'magerr_j': table_j.magerr_aper.values[idxs[0]],
            'mag_ks': table_ks.mag_aper.values[idxs[1]], 'magerr_ks': table_ks.magerr_aper.values[idxs[1]],
            'mag_h': table_h.mag_aper.values[idxs[2]], 'magerr_h': table_h.magerr_aper.values[idxs[2]]}

    return pd.DataFrame(data=mags)


def query_separation_tree(table1, table2, table3, max_sep=0.0001):
    """
    For two pandas database files, with the ra/dec header names x_wcs/y_wcs as defined in read_sex
    """

    tree2 = KDTree(zip(table2.x_wcs.values.ravel(), table2.y_wcs.values.ravel()))
    tree3 = KDTree(zip(table3.x_wcs.values.ravel(), table3.y_wcs.values.ravel()))

    idxs1 = []
    idxs2 = []
    idxs3 = []
    dups2 = []
    dups3 = []

    for i in range(len(table1.index)):
        d2, icoords2 = tree2.query((table1.x_wcs.values[i], table1.y_wcs.values[i]), distance_upper_bound=max_sep)
        d3, icoords3 = tree3.query((table1.x_wcs.values[i], table1.y_wcs.values[i]), distance_upper_bound=max_sep)
        if (d2 != np.inf) & (d3 != np.inf):
            if icoords2 in idxs2:
                dups2.append(icoords2)
            if icoords3 in idxs3:
                dups3.append(icoords3)

            idxs1.append(i)
            idxs2.append(icoords2)
            idxs3.append(icoords3)

    if (len(dups2) > 0) | (len(dups3) > 0):
        print '  WARNING: duplicates found'

    return np.array((idxs1, idxs2, idxs3)), np.array((dups2, dups3))


def plot_cmd(mag_j, mag_h, mag_ks, outfig):
    """
    Plot colour magnitude diagrams for (J-H) v. Ks, (H-Ks) v. Ks
    """

    fig, ax = plt.subplots(nrows=1, ncols=2)

    ax[0].scatter(mag_j - mag_h, mag_ks, marker='.')
    ax[0].set_ylabel('Ks')
    ax[0].set_xlabel('J-H')

    ax[1].scatter(mag_h - mag_ks, mag_ks, marker='.')
    ax[1].set_ylabel('Ks')
    ax[1].set_xlabel('H-Ks')

    if outfig is not None:
        plt.savefig(outfig)

    plt.show()


if __name__ == '__main__':
    main()
