__author__ = 'kawebb'

import argparse
from pyraf import iraf
# import pandas as pd
import numpy as np
from astropy.io import ascii

__author__ = 'kawebb'

"""
Create input for DAOPHOT ADDSTAR to generate artificial stars

python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -ap 30 -zp 30-0.192943311392
"""


def main():
    # clouds = ['CB68', 'L429', 'L1521E', 'L1544']  # , 'L1552']
    # sexfiles_j = ['phot_t3/CB68_J_sex_t3_ap30.txt', 'phot_t3/L429_J_sex_t3_ap24.txt',
    #               'phot_t3/L1521E_J_sex_t3_ap30.txt', 'phot_t3/L1544_J_sex_t3_ap28.txt']  # , 'phot_t3/']
    # sexfiles_ks = ['phot_t3/CB68_Ks_sex_t3_ap24.txt', 'phot_t3/L429_Ks_sex_t3_ap24.txt',
    #                'phot_t3/L1521E_Ks_sex_t3_ap28.txt', 'phot_t3/L1544_Ks_sex_t3_ap28.txt']  # , 'phot_t3/']
    # sexfiles_h = ['phot_t3/CB68_H_sex_t3_ap30.txt', 'phot_t3/L429_H_sex_t3_ap20.txt',
    #               'phot_t3/L1521E_H_sex_t3_ap26.txt', 'phot_t3/L1544_H_sex_t3_ap28.txt']  # , 'phot_t3/']
    #
    # sexfiles = (sexfiles_j, sexfiles_ks, sexfiles_h)
    #
    # aps = [[30, 24, 30, 28], [24, 24, 28, 28], [30, 20, 26, 28]]
    # zmags = [[30 - 0.192943311392, 30 - 0.249078133049, 30 - 0.0278771275042, 30 - 0.681144001009],
    #          [30 + 0.605746118244, 30 + 1.07900318805, 30 + 0.859804630004, 30 + 0.208893637398],
    #          [30 + 0.466909777698, 30 + 0.776624355589, 30 + 0.552873836659, 30 - 1.14661744067]]
    #
    # bands = ['J', 'Ks', 'H']
    # images = '{}/{}_{}_sub.fits'

    parser = argparse.ArgumentParser(
        description='For a given source extractor file, output a text file with coordinates in pixels')
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default=None,
                        help='Source extractor file.')
    parser.add_argument("--image", '-img',
                        action='store',
                        default=None,
                        help='Fits image file.')
    parser.add_argument("--aperture", '-ap',
                        action='store',
                        default=None,
                        help='Aperture for photometry.')
    parser.add_argument("--zeropoint", '-zp',
                        action='store',
                        default=None,
                        help='Zeropoint for photometry.')

    # optional
    parser.add_argument("--outcoords", '-coo',
                        action='store',
                        default=None,
                        help='Outfile of coordinates [pix].')
    parser.add_argument("--outphot", '-phot',
                        action='store',
                        default=None,
                        help='Outfile of daophot photometry.')
    parser.add_argument("--outpst", '-pst',
                        action='store',
                        default=None,
                        help='Outfile of daophot pst.')
    parser.add_argument("--outopst", '-opst',
                        action='store',
                        default=None,
                        help='Outfile of daophot pst.')
    parser.add_argument("--outgroup", '-group',
                        action='store',
                        default=None,
                        help='Outfile of daophot pst groupfile.')
    parser.add_argument("--outpsfimg", '-psfimg',
                        action='store',
                        default=None,
                        help='Outfile of daophot psf image.')
    args = parser.parse_args()

    run_daophot(args.sexfile, args.image, args.aperture, args.zeropoint, args.outcoords, args.outphot, args.outpst,
                args.outopst, args.group, args.outpsfimg)


def run_daophot(sexfile, image, aperture, zeropoint, coords=None, photfile=None, pstfile=None, opstfile=None,
                groupfile=None, psfimage=None,
                maxnpsf=100.):
    if sexfile is None:
        raise Exception, 'No input source extractor file given'
    if image is None:
        raise Exception, 'No input image file given'
    if aperture is None:
        raise Exception, 'No input aperture parameter given'
    if zeropoint is None:
        raise Exception, 'No input zeropoint value given'

    if coords is None:
        coords = sexfile.split('.')[0] + '.coo'
    if photfile is None:
        photfile = sexfile.split('.')[0] + '.phot'
    if pstfile is None:
        pstfile = sexfile.split('.')[0] + '.pst'
    if opstfile is None:
        opstfile = sexfile.split('.')[0] + '.opst'
    if groupfile is None:
        groupfile = sexfile.split('.')[0] + '.group'
    if psfimage is None:
        psfimage = sexfile.split('.')[0] + '.fits'

    data = ascii.read(sexfile, format='sextractor')
    np.savetxt(coords, np.transpose((np.array(data['X_IMAGE']), np.array(data['X_IMAGE']))), delimiter=' ')

    iraf.daophot()
    # phot -- do aperture photometry on a list of stars
    # phot image
    iraf.daophot.phot(image, coords=coords, output=photfile, wcsin="world", apertures=aperture,
                      annulus=aperture + 10, dannulus=10, zmag=zeropoint, verbose=False)

    # pstselect -- select candidate psf stars from a photometry file
    # pstselect image photfile pstfile maxnpsf
    iraf.daophot.pstselect(image, photfile, pstfile, maxnpsf, wcsin='world')

    # psf -- build the point spread function for an image
    # psf image photfile pstfile psfimage opstfile groupfile
    iraf.daophot.psf(image, photfile, pstfile, psfimage, opstfile, groupfile)

    # addstar -- add artificial stars to images
    # addstar image photfile psfimage addimage
    iraf.daophot.addstar(image, photfile, psfimage)  # 50% recovery


# def read_sex_band(sfile, skiplines=15, band=None):
#     #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
#     #   2 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
#     #   3 MAG_APER               Fixed aperture magnitude vector                            [mag]
#     #   4 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
#     #   5 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
#     #   6 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
#     #   7 KRON_RADIUS            Kron apertures in units of A or B
#     #   8 X_IMAGE                Object position along x                                    [pixel]
#     #   9 Y_IMAGE                Object position along y                                    [pixel]
#     #  10 X_WORLD                Barycenter position along world x axis                     [deg]
#     #  11 Y_WORLD                Barycenter position along world y axis                     [deg]
#     #  12 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
#     #  13 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
#     #  14 FLAGS                  Extraction flags
#
#     if band is None:
#         names = ['flux_aper', 'fluxerr_aper', 'mag_aper', 'magerr_aper', 'mag_auto', 'magerr_auto', 'kron_radius',
#                  'bkg', 'x_pix', 'y_pix', 'x_wcs', 'y_wcs', 'fwhm_image', 'fwhm_world', 'flags']
#     else:
#         names = ['flux_aper_{}'.format(band), 'fluxerr_aper_{}'.format(band), 'mag_aper_{}'.format(band),
#                  'magerr_aper_{}'.format(band), 'mag_auto_{}'.format(band), 'magerr_auto_{}'.format(band),
#                  'bkg', 'kron_radius', 'x_pix_{}'.format(band), 'y_pix_{}'.format(band), 'x_wcs_{}'.format(band),
#                  'y_wcs_{}'.format(band), 'fwhm_image', 'fwhm_world', 'flags']
#     # read in file as pandas dataframe
#     table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python', names=names)
#
#     # drop all rows with saturated stars
#     table = table[table[names[2]] != 99.]
#     table.reset_index()
#
#     return table


if __name__ == '__main__':
    main()
