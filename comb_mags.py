__author__ = 'kawebb'

import numpy as np
import argparse
import pandas as pd

import matplotlib.pyplot as plt
from scipy.spatial import KDTree


def main():
    parser = argparse.ArgumentParser(
        description='For photometry output from sextractor for 3 bands, compile magnitude arrays into singe file')
    parser.add_argument("--sexfile_j", '-sfj',
                        action='store',
                        default='CB68_J_sex_t3_r28.txt',
                        help='Object.')
    parser.add_argument("--sexfile_ks", '-sfks',
                        action='store',
                        default='CB68_Ks_sex_t3_r28.txt',
                        help='Object.')
    parser.add_argument("--sexfile_h", '-sfh',
                        action='store',
                        default='CB68_H_sex_t3_r28.txt',
                        help='Object.')
    parser.add_argument("--output", '-out',
                        action='store',
                        default='CB68_mags.txt',
                        help='Object.')
    parser.add_argument("--maxsep", '-ms',
                        action='store',
                        default=0.0001,
                        help='Object.')

    args = parser.parse_args()

    phot_j = read_sex(args.sexfile_j)
    phot_ks = read_sex(args.sexfile_ks)
    phot_h = read_sex(args.sexfile_h)

    sort_by_coords_tree(phot_j, phot_ks, phot_h)
    mags = sort_by_coords_tree(phot_j, phot_ks, phot_h, args.maxsep)
    np.savetxt(str(args.output), mags)


def sort_by_coords_tree(table_j, table_ks, table_h, max_sep=0.0001):
    """
    For every star found by the photometry in the J band, look for the same star in the Ks, H bands
    """

    idxs, dups = query_separation_tree(table_j, table_ks, table_h, 0.0001)

    mag_j = table_j.mag_aper.values[idxs[0]]
    mag_ks = table_ks.mag_aper.values[idxs[1]]
    mag_h = table_h.mag_aper.values[idxs[2]]

    return np.array((mag_j, mag_ks, mag_h))


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

    for i in table1.index:
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


def sort_by_coords(table_j, table_ks, table_h, sep=0.005):
    """
    For every star found by the photometry in the J band, look for the same star in the Ks, H bands
    """

    # make lists of indices for stars already matched, make sure not match a star twice
    imatch_ks = []
    imatch_h = []

    mag_table = np.zeros((len(table_j), 3))
    for i in table_j.index:
        ra = table_j['x_wcs'][i]
        dec = table_j['y_wcs'][i]

        i_ks, sep_ks = query_separation(table_ks, ra, dec, sep)
        i_h, sep_h = query_separation(table_h, ra, dec, sep)

        if i_ks in imatch_ks:
            print '  star in Ks band matched twice, {}'.format(i_ks)
        if i_h in imatch_h:
            print '  star in H band matched twice, {}'.format(i_h)

        mag_table[i, :] = table_j['mag_aper'][i], table_ks['mag_aper'][i_ks], table_h['mag_aper'][i_h]

        imatch_ks.append(i_ks)
        imatch_h.append(i_h)

    return mag_table


def query_separation(table, ra, dec, max_sep=0.01):
    sep = np.sqrt((ra - table.x_wcs.values) ** 2 + (dec - table.y_wcs.values) ** 2)
    if np.min(sep) > max_sep:
        print 'separation greater than expected for coordinates {},{}'.format(ra, dec)
    return np.argmin(sep), np.min(sep)


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


if __name__ == '__main__':
    main()
