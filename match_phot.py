__author__ = 'kawebb'

"""
Procedure:
Determine the best aperture for the photometry by calculating curve of growths for each
star, averaging curves, and fitting an optimal curve.

Example use:
python match_phot.py -img CB68/CB68_J_sub.fits -sf_f CB68/CB68_J_sex_2m_r{}.txt -aps 5 10 15
python match_phot.py -img CB68/CB68_Ks_sub.fits -sf_f CB68/curve_of_growths/CB68_Ks_sex_2m_r{}.txt -aps 5 10 15 18 20 22 24 28 30 32 34 36 38 40 --outfile CB68/CB68_Ks_mags.txt --outimage CB68/CB68_Ks_mags.png
python match_phot.py -img CB68/CB68_J_sub.fits -sf_f CB68/curve_of_growths/CB68_J_sex_2m_r{}.txt -aps 5 10 15 18 20 22 24 28 30 32 34 36 38 40 --outfile CB68/CB68_J_mags.txt --outimage CB68/CB68_J_mags.png
python match_phot.py -img CB68/CB68_H_sub.fits -sf_f CB68/curve_of_growths/CB68_H_sex_2m_r{}.txt -aps 5 10 15 18 20 22 24 28 30 32 34 36 38 40 --outfile CB68/CB68_H_mags.txt --outimage CB68/CB68_H_mags.png

python match_phot.py -img CB68/CB68_Ks_sub.fits -sf_f CB68/curve_of_growths/CB68_Ks_sex_2m_r{}_t7.txt -aps 5 10 15 18 20 22 24 28 30 32 34 36 38 40 42 44 46 48 50 --outfile CB68/CB68_Ks_mags_t7.txt --outimage CB68/CB68_Ks_mags_t7.png
python match_phot.py -img CB68/CB68_J_sub.fits -sf_f CB68/curve_of_growths/CB68_J_sex_2m_r{}_t7.txt -aps 5 10 15 18 20 22 24 28 30 32 34 36 38 40 42 44 46 48 50 --outfile CB68/CB68_J_mags_t7.txt --outimage CB68/CB68_J_mags_t7.png
python match_phot.py -img CB68/CB68_H_sub.fits -sf_f CB68/curve_of_growths/CB68_H_sex_2m_r{}_t7.txt -aps 5 10 15 18 20 22 24 28 30 32 34 36 38 40 42 44 46 48 50 --outfile CB68/CB68_H_mags_t7.txt --outimage CB68/CB68_H_mags_t7.png
"""

import os
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.spatial import KDTree
import zeropoint


def main():
    parser = argparse.ArgumentParser(
        description='For a given image, and list of star positions, return list of stars which are not saturated')
    parser.add_argument("--image", '-img',
                        action="store",
                        default=None,
                        help="Fits image of star field")
    parser.add_argument("--apertures", '-aps',
                        nargs='*',
                        action='store',
                        default=None,
                        help='Apertures of photometry corresponding to sexfils')
    parser.add_argument("--sexfiles", '-sfs',
                        nargs='*',
                        action='store',
                        default=None,
                        help='Source extractor files for J,H, or Ks bands.')
    parser.add_argument("--outfile", '-of',
                        action='store',
                        default=None,
                        help='Outfile of matched magnitudes.')
    parser.add_argument("--outimage",
                        action='store',
                        default=None,
                        help='Output plot of curve of growths.')
    parser.add_argument("--sexfile_format", '-sf_f',
                        action='store',
                        default=None,
                        help='Instead of list of sexfiles, a format.')
    parser.add_argument("--toplot",
                        action='store',
                        default=False,
                        help='Wheter or not to plot figures.')
    parser.add_argument("--max_sep", '-msep',
                        action='store',
                        default=0.0001,
                        help='Maximum separation for nearest neighbour search.')
    parser.add_argument("--r_fed",
                        action='store',
                        default=30,
                        help='Feducial radius.')

    args = parser.parse_args()

    if args.image is None:
        raise Exception, 'No image file given'
    if args.apertures is None:
        print 'Warning: No apertures specified'

    if args.sexfiles is None:
        print 'WARNING: No source extractor files given'

        if args.sexfile_format is None:
            raise Exception, 'No source extractor files or format given'
        else:
            print 'Using source extractor format {}'.format(args.sexfile_format)
            sexfiles = []
            for ap in args.apertures:
                sexfiles.append(args.sexfile_format.format(ap))
    else:
        sexfiles = args.sexfiles

    names = ['x_wcs', 'y_wcs']
    for ap in args.apertures:
        names.append('mag_aper_{}'.format(ap))
        names.append('magerr_aper_{}'.format(ap))

    # if not os.path.exists(args.outfile):
    phottable = match_phot(args.image, sexfiles, args.apertures, args.toplot, float(args.max_sep))

    if args.outfile is None:
        outfile = 'test.txt'
    else:
        outfile = args.outfile

    phottable[names].to_csv(outfile, index=False)
    # else:
    #     phottable = pd.read_csv(args.outfile, skiprows=1, names=names)

    curve_of_growths(phottable, np.array(args.apertures, dtype=np.int), int(args.r_fed), args.outimage)


def curve_of_growths(phottable, apertures, r_feducial, outimage=None):
    """
    Plot curve of growth, relative to a given feducial radius
    """
    magtable = np.ones((len(apertures), len(phottable)))
    magerrtable = np.ones((len(apertures), len(phottable)))

    for i, ap in enumerate(apertures):
        magtable[i, :] = phottable['mag_aper_{}'.format(ap)].values
        magerrtable[i, :] = phottable['magerr_aper_{}'.format(ap)].values

    iap = np.argmin(abs(r_feducial - np.array(apertures)))
    ap_fed = apertures[iap]
    print 'Given feducial radius {}, closest aperture in data set {}'.format(r_feducial, ap_fed)

    # magtable_rel = magtable*0.
    average = apertures * 0.

    for i in range(len(apertures)):
        # magtable_rel[i, :] = np.subtract(magtable[i, :], magtable[iap, :])
        average[i] = np.average(magtable[i, :], weights=(magerrtable[i, :]) ** (-2))
        # average[i] = np.average(magtable_rel[i, :], weights=(magerrtable[i, :])**(-2))

    average_rel = average * 0.
    for i in range(len(apertures)):
        average_rel[i] = average[i] - average[iap]

    ax = plt.gca()
    ax.plot(apertures, average_rel, '-b')
    ax.plot(apertures, apertures * 0., '-k')
    ax.invert_yaxis()
    ax.set_xlabel('apertures')
    ax.set_ylabel(r'$\Delta$ mag relative to aperture {}'.format(ap_fed))
    if outimage is not None:
        plt.savefig(outimage)
    plt.show()


def match_phot(image, sexfiles, apertures, toplot=False, max_sep=0.0001):
    """
    Read photometry files for several apertures, remove saturated stars, and match coordinate lists
    """

    imgdata, imghdr = read_fits(image)  # Access fits data to read locations of saturated stars

    sextable = read_sex(sexfiles[0], skiplines=13, aperture=apertures[0])  # Import source extractor file
    # detect saturated stars in fits image ('0' center), remove the from the source extractor phot. table
    # returns index of unsaturated stars as True, saturated as False
    unsat_idxs = zeropoint.detect_satur(imgdata, sextable['x'].values, sextable['y'].values, sextable['b'].values)
    sextable = sextable[unsat_idxs]  # splice phot. table according to index array
    sextable.reset_index(inplace=True, drop=True)  # drop index column (inherent to pandas dataframe) as out of order
    # drop unneeded data from photometry table
    sextable.drop(['flux_aper', 'fluxerr_aper', 'thresh', 'fwhm_image', 'fwhm_world', 'flag', 'x', 'y', 'b'], axis=1,
                  inplace=True)

    # iterate through all source extractor files input
    for i, sexfile in enumerate(sexfiles[1:]):
        i += 1  # because starting at index 1
        print '** {} {}'.format(sexfile, apertures[i])

        # do the same steps as above
        phottable = read_sex(sexfile, skiplines=13, aperture=apertures[i])
        unsat_idxs = zeropoint.detect_satur(imgdata, phottable['x'].values, phottable['y'].values,
                                            phottable['b'].values)
        unsat_stars = phottable[unsat_idxs]
        # match the two tables by their coordinates (nearest neighbour), return combined table
        sextable = match_tables(sextable, unsat_stars, max_sep, aper=apertures[i], toplot=toplot)


    # remove additional saturated stars marked with 99.
    for aper in apertures:
        sextable = sextable[sextable['mag_aper_{}'.format(aper)] != 99.]

    return sextable


def match_tables(reftable, newtable, max_sep, aper, toplot=False):
    """
    Match coordinates of table_new to table, return ordered list of indices
    """

    # match the stars in the reference list to the new table of stars
    refidxs, idxs = query_separation_tree(reftable.x_wcs.values, reftable.y_wcs.values, newtable.x_wcs.values,
                                          newtable.y_wcs.values, max_sep)

    spliced_reftable = reftable[refidxs]  # apply True/False array to splice out matched indexes
    spliced_reftable.reset_index(inplace=True, drop=True)
    spliced_newtable = newtable[idxs]
    spliced_newtable.reset_index(inplace=True, drop=True)

    assert len(spliced_reftable) == len(spliced_newtable), 'Matched tables not the same length'

    # sort both tables by the x axis
    sorted_reftable = spliced_reftable.sort_values('x_wcs')
    sorted_newtable = spliced_newtable.sort_values('x_wcs')

    if toplot:
        # plot to make sure matching the same points
        plt.plot(sorted_reftable.x_wcs.values, color='k')
        plt.plot(sorted_newtable.x_wcs.values, color='r')
        plt.show()
        plt.plot(sorted_reftable.y_wcs.values, color='k')
        plt.plot(sorted_newtable.y_wcs.values, color='r')
        plt.show()

    sorted_newtable.drop(
        ['flux_aper', 'fluxerr_aper', 'thresh', 'fwhm_image', 'fwhm_world', 'flag', 'x', 'y', 'b', 'x_wcs', 'y_wcs'],
        axis=1, inplace=True)

    return pd.concat(
        (sorted_reftable, sorted_newtable['mag_aper_{}'.format(aper)], sorted_newtable['magerr_aper_{}'.format(aper)]),
        axis=1)


def query_separation_tree(xref, yref, x, y, max_sep=0.0001):
    """
    For two pandas database files, with the ra/dec header names x_wcs/y_wcs as defined in read_sex
    """

    # establish 'tree' of coordinates
    tree = KDTree(zip(x.ravel(), y.ravel()))

    refidxs = np.zeros(len(xref), dtype=bool)  # splice pandas Dataframe with True/False array
    idxs = np.zeros(len(x), dtype=bool)  # By default False, when matched set to True

    for i in range(len(xref)):
        d, idx = tree.query((xref[i], yref[i]), distance_upper_bound=max_sep)
        if (d != np.inf):  # ignore results with no matches within max_sep
            if (refidxs[i] is True) | (idxs[idx] is True):
                print '  WARNING: duplicates found for refidxs[{}]'.format(i)

            refidxs[i] = True
            idxs[idx] = True

        else:
            print 'no match index {}'.format(i)

    if np.count_nonzero(refidxs) != np.count_nonzero(idxs):
        print 'Number of indexes dont match: ref {}, new {}'.format(np.count_nonzero(refidxs), np.count_nonzero(idxs))
        # print '  index where ref False: {}'.format([index for index,value in enumerate(refidxs) if value < 1])
        # print '  index where new False: {}'.format([index for index,value in enumerate(idxs) if value < 1])
        raise Exception

    return refidxs, idxs


def read_sex(sfile, skiplines=13, aperture=None):
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

    # c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13 = np.loadtxt(sfile, skiprows=skiplines, unpack=True)
    # dict = {'mag_aper_{}'.format(aperture): c3, 'magerr_aper_{}'.format(aperture): c4, 'x_wcs': c8, 'y_wcs': c9, 'b':c10}
    # return pd.DataFrame(data=dict)

    if aperture is not None:
        names = ['flux_aper', 'fluxerr_aper', 'mag_aper_{}'.format(aperture), 'magerr_aper_{}'.format(aperture),
                 'thresh', 'x', 'y', 'x_wcs', 'y_wcs', 'b', 'fwhm_image', 'fwhm_world', 'flag']
    else:
        names = ['flux_aper', 'fluxerr_aper', 'mag_aper', 'magerr_aper', 'thresh', 'x', 'y', 'x_wcs', 'y_wcs',
                 'b', 'fwhm_image', 'fwhm_world', 'flag']
    table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python', names=names)
    return table


def read_fits(ffile):
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
