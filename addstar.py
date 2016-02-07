__author__ = 'kawebb'

import numpy as np
# import plot_cmd
# from pyraf import iraf
import pandas as pd

__author__ = 'kawebb'

"""
Create input for DAOPHOT ADDSTAR to generate artificial stars

addstar -- add artificial stars to images
addstar image photfile psfimage addimage

psf -- build the point spread function for an image
psf image photfile pstfile psfimage opstfile groupfile

pstselect -- select candidate psf stars from a photometry file
pstselect image photfile pstfile maxnpsf
"""


def main():

    clouds = ['CB68', 'L429', 'L1521E', 'L1544']  # , 'L1552']
    sexfiles_j = ['phot_t3/CB68_J_sex_t3_ap30.txt', 'phot_t3/L429_J_sex_t3_ap24.txt',
                  'phot_t3/L1521E_J_sex_t3_ap30.txt', 'phot_t3/L1544_J_sex_t3_ap28.txt']  # , 'phot_t3/']
    sexfiles_ks = ['phot_t3/CB68_Ks_sex_t3_ap24.txt', 'phot_t3/L429_Ks_sex_t3_ap24.txt',
                   'phot_t3/L1521E_Ks_sex_t3_ap28.txt', 'phot_t3/L1544_Ks_sex_t3_ap28.txt']  # , 'phot_t3/']
    sexfiles_h = ['phot_t3/CB68_H_sex_t3_ap30.txt', 'phot_t3/L429_H_sex_t3_ap20.txt',
                  'phot_t3/L1521E_H_sex_t3_ap26.txt', 'phot_t3/L1544_H_sex_t3_ap28.txt']  # , 'phot_t3/']

    sexfiles = (sexfiles_j, sexfiles_ks, sexfiles_h)

    bands = ['J', 'Ks', 'H']
    images = '{}/{}_{}_sub.fits'

    output_photfile = 'addstars/{}_{}_phot.txt'

    output_pstfile = 'addstars/{}_{}_pst.txt'
    output_opstfile = 'addstars/{}_{}_opst.txt'
    output_groupfile = 'addstars/{}_{}_group.txt'

    output_psfimage = 'addstars/{}_{}_psf.fits'

    maxnpsf = 100

    for i, cloud in enumerate(clouds[:1]):
        for j, band in enumerate(bands[:1]):
            image = images.format(cloud, cloud, band)
            photfile = output_photfile.format(cloud, band)
            pstfile = output_pstfile.format(cloud, band)
            opstfile = output_opstfile.format(cloud, band)
            groupfile = output_groupfile.format(cloud, band)
            psfimage = output_psfimage.format(cloud, band)


            print sexfiles[j][i]

            # read in photometry from source extractor files, keep x,y,mag data
            # sex_phot = plot_cmd.read_sex_band(sexfiles[j][i])
            sex_phot = read_sex_band(sexfiles[j][i])

            phot = (range(len(sex_phot)), sex_phot.x_wcs.values, sex_phot.y_wcs.values, sex_phot.mag_aper.values)
            np.savetxt(photfile, phot, header='ID XCENTER YCENTER MAG MSKY', delimiter=' ')

            # NEEDS MSKY

            # PSTSELECT reads the input photometry file photfile , extracts the ID, XCENTER, YCENTER, MAG, and MSKY
            # fields for up to maxnpsf psf stars, and the results to pstfile . Pstfile automatically inherits the file
            # format of photfile .
            #
            # pstselect -- select candidate psf stars from a photometry file
            # pstselect image photfile pstfile maxnpsf
            iraf.daophot()
            iraf.daophot.pstselect(image, photfile, pstfile, maxnpsf, wcsin='world')

            # psf -- build the point spread function for an image
            # psf image photfile pstfile psfimage opstfile groupfile
            iraf.daophot.psf(image, photfile, pstfile, psfimage, opstfile, groupfile)

            # addstar -- add artificial stars to images
            # addstar image photfile psfimage addimage
            iraf.daophot.addstar(image, photfile, psfimage)


def read_sex_band(sfile, skiplines=15, band=None):
    #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
    #   2 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
    #   3 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   4 MAGERR_APER            RMS error vector for fixed aperture mag.                   [mag]
    #   5 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
    #   6 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
    #   7 KRON_RADIUS            Kron apertures in units of A or B
    #   8 X_IMAGE                Object position along x                                    [pixel]
    #   9 Y_IMAGE                Object position along y                                    [pixel]
    #  10 X_WORLD                Barycenter position along world x axis                     [deg]
    #  11 Y_WORLD                Barycenter position along world y axis                     [deg]
    #  12 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
    #  13 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
    #  14 FLAGS                  Extraction flags

    if band is None:
        names = ['flux_aper', 'fluxerr_aper', 'mag_aper', 'magerr_aper', 'mag_auto', 'magerr_auto', 'kron_radius',
                 'bkg', 'x_pix', 'y_pix', 'x_wcs', 'y_wcs', 'fwhm_image', 'fwhm_world', 'flags']
    else:
        names = ['flux_aper_{}'.format(band), 'fluxerr_aper_{}'.format(band), 'mag_aper_{}'.format(band),
                 'magerr_aper_{}'.format(band), 'mag_auto_{}'.format(band), 'magerr_auto_{}'.format(band),
                 'bkg', 'kron_radius', 'x_pix_{}'.format(band), 'y_pix_{}'.format(band), 'x_wcs_{}'.format(band),
                 'y_wcs_{}'.format(band), 'fwhm_image', 'fwhm_world', 'flags']
    # read in file as pandas dataframe
    table = pd.read_csv(sfile, skiprows=skiplines, sep=r"\s*", engine='python', names=names)

    # drop all rows with saturated stars
    table = table[table[names[2]] != 99.]
    table.reset_index()

    return table


if __name__ == '__main__':
    main()
