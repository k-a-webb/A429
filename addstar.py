__author__ = 'kawebb'

from pyraf import iraf
import pandas as pd
import numpy as np

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

    aps = [[30, 24, 30, 28], [24, 24, 28, 28], [30, 20, 26, 28]]

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
            # image = images.format(cloud, cloud, band)
            image = images.format('release', cloud, band)
            photfile = output_photfile.format(cloud, band)
            pstfile = output_pstfile.format(cloud, band)
            opstfile = output_opstfile.format(cloud, band)
            groupfile = output_groupfile.format(cloud, band)
            psfimage = output_psfimage.format(cloud, band)
            ap = aps[j][i]

            sex_phot = read_sex_band(sexfiles[j][i])
            np.savetext('addstars/coords.txt', np.transpose(sex_phot.x_wcs.values, sex_phot.y_wcs.values))

            # write_phot(sexfiles[j][i], photfile, image)

            iraf.daophot()
            # phot -- do aperture photometry on a list of stars
            # phot image
            print ap
            iraf.daophot.phot(image, coords="coords.txt", output=photfile, wcsin="world")

            # pstselect -- select candidate psf stars from a photometry file
            # pstselect image photfile pstfile maxnpsf

            iraf.daophot.pstselect(image, photfile, pstfile, maxnpsf, wcsin='world')

            # psf -- build the point spread function for an image
            # psf image photfile pstfile psfimage opstfile groupfile
            iraf.daophot.psf(image, photfile, pstfile, psfimage, opstfile, groupfile)

            # addstar -- add artificial stars to images
            # addstar image photfile psfimage addimage
            iraf.daophot.addstar(image, photfile, psfimage)  # 50% recovery


"""
def write_phot(sexfile, photfile, image):
    # read in photometry from source extractor files, keep x,y,mag data
    sex_phot = read_sex_band(sexfile)

    # IMAGE, XINIT, YINIT, ID, COORDS, LID, XCENTER, YCENTER, XSHIFT, YSHIFT, XERR, YERR, CIER CERROR, MSKY,
    # STDEV, SSKEW, NSKY, NSREJ, SIER SERROR, ITIME, XAIRMASS, IFILTER, OTIME, RAPERT, SUM, AREA, FLUX, MAG,
    # MERR, PIER PERROR

    with open(photfile, 'w') as outfile:
        outfile.write('#K IRAF       = NOAO/IRAFV2.16.1        version    %-23s \n')
        outfile.write('#K USER       = kawebb                  name       %-23s \n')
        outfile.write('#K HOST       = sangiovese.phys.uvic.ca computer   %-23s \n')
        outfile.write('#K DATE       = 2016-02-07              yyyy-mm-dd %-23s \n')
        outfile.write('# K TIME       = 18:11:19                hh:mm:ss   %-23s \n')
        outfile.write('# K PACKAGE    = apphot                  name       %-23s \n')
        outfile.write('# K TASK       = phot                    name       %-23s \n')
        outfile.write('# \n')
        outfile.write('# K SCALE      = 1.                      units      %-23.7g \n')
        outfile.write('# K FWHMPSF    = 2.5                     scaleunit  %-23.7g \n')
        outfile.write('# K EMISSION   = yes                     switch     %-23b   \n')
        outfile.write('# K DATAMIN    = INDEF                   counts     %-23.7g \n')
        outfile.write('# K DATAMAX    = INDEF                   counts     %-23.7g \n')
        outfile.write('# K EXPOSURE   = ""                      keyword    %-23s \n')
        outfile.write('# K AIRMASS    = ""                      keyword    %-23s \n')
        outfile.write('# K FILTER     = ""                      keyword    %-23s \n')
        outfile.write('# K OBSTIME    = ""                      keyword    %-23s \n')
        outfile.write('# \n')
        outfile.write('# K NOISE      = poisson                 model      %-23s \n')
        outfile.write('# K SIGMA      = 3.                      counts     %-23.7g \n')
        outfile.write('# K GAIN       = ""                      keyword    %-23s \n')
        outfile.write('# K EPADU      = 1.                      e-/adu     %-23.7g \n')
        outfile.write('# K CCDREAD    = ""                      keyword    %-23s \n')
        outfile.write('# K READNOISE  = 0.                      e-         %-23.7g \n')
        outfile.write('# \n')
        outfile.write('# K CALGORITHM = none                    algorithm  %-23s \n')
        outfile.write('# K CBOXWIDTH  = 5.                      scaleunit  %-23.7g \n')
        outfile.write('# K CTHRESHOLD = 0.                      sigma      %-23.7g \n')
        outfile.write('# K MINSNRATIO = 1.                      number     %-23.7g \n')
        outfile.write('# K CMAXITER   = 10                      number     %-23d \n')
        outfile.write('# K MAXSHIFT   = 1.                      scaleunit  %-23.7g \n')
        outfile.write('# K CLEAN      = no                      switch     %-23b \n')
        outfile.write('# K RCLEAN     = 1.                      scaleunit  %-23.7g \n')
        outfile.write('# K RCLIP      = 2.                      scaleunit  %-23.7g \n')
        outfile.write('# K KCLEAN     = 3.                      sigma      %-23.7g \n')
        outfile.write('# \n')
        outfile.write('# K SALGORITHM = mode                    algorithm  %-23s \n')
        outfile.write('# K ANNULUS    = 30.                     scaleunit  %-23.7g \n')
        outfile.write('# K DANNULUS   = 10.                     scaleunit  %-23.7g \n')
        outfile.write('# K SKYVALUE   = 0.                      counts     %-23.7g \n')
        outfile.write('# K KHIST      = 3.                      sigma      %-23.7g \n')
        outfile.write('# K BINSIZE    = 0.1                     sigma      %-23.7g \n')
        outfile.write('# K SMOOTH     = no                      switch     %-23b \n')
        outfile.write('# K SMAXITER   = 10                      number     %-23d \n')
        outfile.write('# K SLOCLIP    = 0.                      percent    %-23.7g \n')
        outfile.write('# K SHICLIP    = 0.                      percent    %-23.7g \n')
        outfile.write('# K SNREJECT   = 50                      number     %-23d   \n')
        outfile.write('# K SLOREJECT  = 3.                      sigma      %-23.7g \n')
        outfile.write('# K SHIREJECT  = 3.                      sigma      %-23.7g \n')
        outfile.write('# K RGROW      = 0.                      scaleunit  %-23.7g \n')
        outfile.write('# \n')
        outfile.write('# K WEIGHTING  = constant                model      %-23s \n')
        outfile.write('# K APERTURES  = 30                      scaleunit  %-23s \n')
        outfile.write('# K ZMAG       = 25.                     zeropoint  %-23.7g \n')
        outfile.write('# \n')
        outfile.write('# N IMAGE               XINIT     YINIT     ID    COORDS                 LID    \ \n')
        outfile.write('# U imagename           pixels    pixels    ##    filename               ##     \ \n')
        outfile.write('# F %-23s               %-10.3f   %-10.3f   %-6d  %-23s                  %-6d     \n')
        outfile.write('# \n')
        outfile.write('# N XCENTER    YCENTER    XSHIFT  YSHIFT  XERR    YERR            CIER CERROR   \ \n')
        outfile.write('# U pixels     pixels     pixels  pixels  pixels  pixels          ##   cerrors  \ \n')
        outfile.write('# F %-14.3f    %-11.3f    %-8.3f  %-8.3f  %-8.3f  %-15.3f         %-5d %-9s       \n')
        outfile.write('# \n')
        outfile.write('# N MSKY           STDEV          SSKEW          NSKY   NSREJ     SIER SERROR   \ \n')
        outfile.write('# U counts         counts         counts         npix   npix      ##   serrors  \ \n')
        outfile.write('# F %-18.7g        %-15.7g        %-15.7g        %-7d   %-9d      %-5d %-9s       \n')
        outfile.write('# \n')
        outfile.write('# N ITIME          XAIRMASS       IFILTER                OTIME                  \ \n')
        outfile.write('# U timeunit       number         name                   timeunit               \ \n')
        outfile.write('# F %-18.7g        %-15.7g        %-23s                  %-23s                    \n')
        outfile.write('# \n')
        outfile.write('# N RAPERT   SUM           AREA       FLUX          MAG    MERR   PIER PERROR   \ \n')
        outfile.write('# U scale    counts        pixels     counts        mag    mag    ##   perrors  \ \n')
        outfile.write('# F %-12.2f  %-14.7g       %-11.7g    %-14.7g       %-7.3f %-6.3f %-5d %-9s       \n')
        outfile.write('# \n')

        for i in range(len(sex_phot)):
            outfile.write('{} 0. 0. {} coords.txt {} \ \n'.format(image, i, i))
            outfile.write('{} {} 0. 0. INDEF INDEF 0 NoError \ \n'.format(sex_phot.x_wcs.values[i], sex_phot.y_wcs.values[i]))
            outfile.write('{} 0. 0. 0 0 0 NoError \ \n'.format(sex_phot.bkg.values[i]))
            outfile.write('1. INDEF INDEF INDEF \ \n')
            outfile.write('30. 0. 0. 0. {} {} 0 NoError \ \n'.format(sex_phot.mag_aper.values[i], sex_phot.magerr_aper.values[i]))
"""


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

"""
#K IRAF       = NOAO/IRAFV2.16.1        version    %-23s
#K USER       = kawebb                  name       %-23s
#K HOST       = sangiovese.phys.uvic.ca computer   %-23s
#K DATE       = 2016-02-07              yyyy-mm-dd %-23s
#K TIME       = 13:39:18                hh:mm:ss   %-23s
#K PACKAGE    = apphot                  name       %-23s
#K TASK       = phot                    name       %-23s
#
#K SCALE      = 1.                      units      %-23.7g
#K FWHMPSF    = 2.5                     scaleunit  %-23.7g
#K EMISSION   = yes                     switch     %-23b
#K DATAMIN    = INDEF                   counts     %-23.7g
#K DATAMAX    = INDEF                   counts     %-23.7g
#K EXPOSURE   = ""                      keyword    %-23s
#K AIRMASS    = ""                      keyword    %-23s
#K FILTER     = ""                      keyword    %-23s
#K OBSTIME    = ""                      keyword    %-23s
#
#K NOISE      = poisson                 model      %-23s
#K SIGMA      = 0.                      counts     %-23.7g
#K GAIN       = ""                      keyword    %-23s
#K EPADU      = 1.                      e-/adu     %-23.7g
#K CCDREAD    = ""                      keyword    %-23s
#K READNOISE  = 0.                      e-         %-23.7g
#
#K CALGORITHM = none                    algorithm  %-23s
#K CBOXWIDTH  = 5.                      scaleunit  %-23.7g
#K CTHRESHOLD = 0.                      sigma      %-23.7g
#K MINSNRATIO = 1.                      number     %-23.7g
#K CMAXITER   = 10                      number     %-23d
#K MAXSHIFT   = 1.                      scaleunit  %-23.7g
#K CLEAN      = no                      switch     %-23b
#K RCLEAN     = 1.                      scaleunit  %-23.7g
#K RCLIP      = 2.                      scaleunit  %-23.7g
#K KCLEAN     = 3.                      sigma      %-23.7g
#
#K SALGORITHM = mode                    algorithm  %-23s
#K ANNULUS    = 10.                     scaleunit  %-23.7g
#K DANNULUS   = 10.                     scaleunit  %-23.7g
#K SKYVALUE   = 0.                      counts     %-23.7g
#K KHIST      = 3.                      sigma      %-23.7g
#K BINSIZE    = 0.1                     sigma      %-23.7g
#K SMOOTH     = no                      switch     %-23b
#K SMAXITER   = 10                      number     %-23d
#K SLOCLIP    = 0.                      percent    %-23.7g
#K SHICLIP    = 0.                      percent    %-23.7g
#K SNREJECT   = 50                      number     %-23d
#K SLOREJECT  = 3.                      sigma      %-23.7g
#K SHIREJECT  = 3.                      sigma      %-23.7g
#K RGROW      = 0.                      scaleunit  %-23.7g
#
#K WEIGHTING  = constant                model      %-23s
#K APERTURES  = 3.                      scaleunit  %-23s
#K ZMAG       = 25.                     zeropoint  %-23.7g
#
#
#N IMAGE               XINIT     YINIT     ID    COORDS                 LID    \
#U imagename           pixels    pixels    ##    filename               ##     \
#F %-23s               %-10.3f   %-10.3f   %-6d  %-23s                  %-6d
#
#N XCENTER    YCENTER    XSHIFT  YSHIFT  XERR    YERR            CIER CERROR   \
#U pixels     pixels     pixels  pixels  pixels  pixels          ##   cerrors  \
#F %-14.3f    %-11.3f    %-8.3f  %-8.3f  %-8.3f  %-15.3f         %-5d %-9s
#
#N MSKY           STDEV          SSKEW          NSKY   NSREJ     SIER SERROR   \
#U counts         counts         counts         npix   npix      ##   serrors  \
#F %-18.7g        %-15.7g        %-15.7g        %-7d   %-9d      %-5d %-9s
#
#N ITIME          XAIRMASS       IFILTER                OTIME                  \
#U timeunit       number         name                   timeunit               \
#F %-18.7g        %-15.7g        %-23s                  %-23s
#
#N RAPERT   SUM           AREA       FLUX          MAG    MERR   PIER PERROR   \
#U scale    counts        pixels     counts        mag    mag    ##   perrors  \
#F %-12.2f  %-14.7g       %-11.7g    %-14.7g       %-7.3f %-6.3f %-5d %-9s
#
#
CB68_J_sub.fits        3427.207  148.638   1     coords.txt             1      \
   3427.207   148.638    0.000   0.000   INDEF   INDEF          0    NoError   \
   -9.509591      105.8633       -38.67295      881    62       0    NoError   \
   1.             INDEF          INDEF                  INDEF                  \
   3.00     15077.21      28.54634   15348.67      14.535 0.042 0    NoError
CB68_J_sub.fits        3172.455  148.505   2     coords.txt             2      \
   3172.455   148.505    0.000   0.000   INDEF   INDEF          0    NoError   \
   -1.919249      92.5009        -47.00073      871    78       0    NoError   \
   1.             INDEF          INDEF                  INDEF                  \
   3.00     5247.238      28.29497   5301.543      15.689 0.103 0    NoError
CB68_J_sub.fits        3221.972  150.446   3     coords.txt             3      \
   3221.972   150.446    0.000   0.000   INDEF   INDEF          0    NoError   \
   -14.67083      145.0668       68.19846       874    67       0    NoError   \
   1.             INDEF          INDEF                  INDEF                  \
   3.00     3313.528      28.77096   3735.622      16.069 0.231 0    NoError
CB68_J_sub.fits        3626.566  149.766   4     coords.txt             4      \
   3626.566   149.766    0.000   0.000   INDEF   INDEF          0    NoError   \
   -1.125241      58.79516       -32.04229      742    198      0    NoError   \
   1.             INDEF          INDEF                  INDEF                  \
   3.00     2237.937      28.50082   2270.007      16.610 0.155 0    NoError
"""
