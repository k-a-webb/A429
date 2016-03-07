import argparse
import os
import numpy as np
from astropy.io import ascii, fits
from scipy.spatial import KDTree
from scipy.stats import chisquare
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LogNorm

"""
Create artificial stars from a good candidate star in the image. Determine number of stars recovered.

Run addstars to generate model PSF and implant # of them into the image
Run source extractor with code output to measure photometry on image with artificial stars
Run addstars to counts number of artificial stars recovered
"""

_SEXFILE = 'phot_t3/CB68_J_sex_t3_ap30.txt'
# _IMAGE = 'CB68/CB68_J_sub.fits'
_IMAGE = 'release/CB68_J_sub.fits'
_ARTIMAGE = 'mag_lim/CB68_J_sub_{}_art.fits'

###### Default values
_MAGRANGE = [18.,24.]
_STAR_INDICES = [None, None, None, None, None, None, None, None, None, None, None, None]
_MAGS = [None, None, None, None, None, None, None, None, None, None, None, None]
_MAXMAGERR = 0.01

_ARTSEXFILE = 'mag_lim/CB68_J_sex_{}_art.txt'
_ARTPNG = 'mag_lim/CB68_J_sex_{}_art.png'

_ARTFILES = ['mag_lim/art1.coo', 'mag_lim/art2.coo', 'mag_lim/art3.coo', 'mag_lim/art4.coo', 'mag_lim/art5.coo',
            'mag_lim/art6.coo', 'mag_lim/art7.coo', 'mag_lim/art8.coo', 'mag_lim/art9.coo', 'mag_lim/art10.coo',
            'mag_lim/art11.coo', 'mag_lim/art12.coo']
#               1                2              3               4              5                6                 7
#               8                9              10              11             12
_XRANGE = [[700., 2000.], [2100., 3400.], [3600., 4800.], [700., 2000.], [2100., 3400.], [3600., 4800.], [700., 2000.],
           [2100., 3400.], [3600., 4800.], [700., 2000.], [2100., 3400.], [3600., 4800.]]
_YRANGE = [[3800., 4600.], [3800., 4600.], [3800., 4600.], [2800., 3600.], [2800., 3600.], [2800., 3600.],
           [1300., 2500.], [1300., 2500.], [1300., 2500.], [900., 1600.], [900., 1600.], [900., 1600.]]

# _ARTSEXFILE = _SEXFILE
# _ARTIMAGE = _IMAGE
# _MAGS = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]

# CB68 J
#
#_MAGRANGE = [21.0, 21.5]
#_STAR_INDICES = [2,3,4,7,5,7,7,1,2, 2,10,12]
#_MAGS = [21.1269,21.0681,21.1589,21.1441,21.0893,21.2558,21.0657,21.1308,21.1766,21.0787,21.2691,21.1766]
#_MAXMAGERR = 0.03

#_MAGRANGE = [21.5, 22.]
#_STAR_INDICES = [4, 2, 3, 4, 3, 7, 8, 3, 2, 5, 2, 1]
#_MAGS = [21.7131, 21.6848, 21.5841, 21.5804, 21.5654, 21.5035, 21.6421, 21.6441, 21.7821, 21.5449, 21.7918, 21.6072]
#_MAXMAGERR = 0.05

#_MAGRANGE = [22., 22.5]
#_STAR_INDICES = [0,0,0,2,1,0,1,1,0,1,1,0]
#_MAGS = [22.0183, 22.067, 22.0456, 22.0486, 22.0165, 22.0436, 22.0496, 22.0325, 22.1663, 22.0496, 22.2055, 22.0092]
#_MAXMAGERR = 0.06

#_MAGRANGE = [22.5, 23.]
#_STAR_INDICES = [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0]
#_MAGS = [22.6074, 22.5360, 22.6229, 22.5411, 22.5966, 22.6141, 22.56, 22.5042, 22.5526, 22.562, 22.5042, 22.969]
#_MAXMAGERR = 0.06

#_MAGRANGE = [23., 25.]
#_STAR_INDICES = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#_MAGS = [23.1053, 23.1185, 23.0592, 23.1472, 23.5124, 23.03, 23.0679, 23.016, 23.0441, 23.2377, 23.016, 23.3859]
#_MAXMAGERR = 0.06


# More default values
_DIMENSION = 25.
_NSTARS = 100.
_MAXCOUNTS = np.inf
_MINCOUNTS = -np.inf
_MIN_MAXCOUNTS = -np.inf
_TOPLOT = False
_CHI2MAX = 100.
_FUNCTION = 'gauss'


def main():
    # clouds = ['CB68', 'L429', 'L1521E', 'L1544']  # , 'L1552']
    # sexfiles_j = ['phot_t3/CB68_J_sex_t3_ap30.txt', 'phot_t3/L429_J_sex_t3_ap24.txt',
    #               'phot_t3/L1521E_J_sex_t3_ap30.txt', 'phot_t3/L1544_J_sex_t3_ap28.txt']  # , 'phot_t3/']
    # sexfiles_ks = ['phot_t3/CB68_Ks_sex_t3_ap24.txt', 'phot_t3/L429_Ks_sex_t3_ap24.txt',
    #                'phot_t3/L1521E_Ks_sex_t3_ap28.txt', 'phot_t3/L1544_Ks_sex_t3_ap28.txt']  # , 'phot_t3/']
    # sexfiles_h = ['phot_t3/CB68_H_sex_t3_ap30.txt', 'phot_t3/L429_H_sex_t3_ap20.txt',
    #               'phot_t3/L1521E_H_sex_t3_ap26.txt', 'phot_t3/L1544_H_sex_t3_ap28.txt']  # , 'phot_t3/']
    # aps = [[30, 24, 30, 28], [24, 24, 28, 28], [30, 20, 26, 28]]
    # zmags = [[30 - 0.192943311392, 30 - 0.249078133049, 30 - 0.0278771275042, 30 - 0.681144001009],
    #          [30 + 0.605746118244, 30 + 1.07900318805, 30 + 0.859804630004, 30 + 0.208893637398],
    #          [30 + 0.466909777698, 30 + 0.776624355589, 30 + 0.552873836659, 30 - 1.14661744067]]

    for i, artfile in enumerate(_ARTFILES):
        if not os.path.exists(artfile):
            generate_coords(artfile, int(_NSTARS), _XRANGE[i], _YRANGE[i])

    magstr = '{:<4.4f}'.format(np.mean(_MAGS)).replace('.', 'p')
    artsexfile = _ARTSEXFILE.format(magstr)

    if not os.path.exists(artsexfile):
        if _ARTIMAGE is None:
            raise Exception, 'No output image file format is given'
        if _SEXFILE is None:
            raise Exception, 'No input source extractor file given'
        if _IMAGE is None:
            raise Exception, 'No input image file given'
        if _MAGRANGE is None:
            raise Exception, 'No magnitude range given'
        assert len(_MAGRANGE) == 2, 'Magnitude input is two values, min and max'
        assert _FUNCTION == 'gauss' or _FUNCTION == 'lorentz', 'ERROR: function not "gauss" or "lorentz"'

        indices = []
        magnitudes = []
        x_art = np.array([])
        y_art = np.array([])

        print '*********************************************************************'
        print '**                           ADDSTAR                               **'
        print '*********************************************************************'

        with fits.open(_IMAGE) as hdu:
            imgdata = hdu[0].data

        for i in range(len(_XRANGE)):
            art_imgdata, x, y, idx, mag = addstar(_SEXFILE, _IMAGE, _ARTFILES[i], _MAGRANGE, _STAR_INDICES[i], _XRANGE[i],
                               _YRANGE[i], float(_MAXMAGERR), float(_CHI2MAX), float(_MINCOUNTS), float(_MAXCOUNTS),
                               float(_MIN_MAXCOUNTS), _DIMENSION, _FUNCTION, _TOPLOT)
            indices.append(idx)
            magnitudes.append(mag)
            imgdata += art_imgdata
            x_art = np.concatenate((x_art, x))
            y_art = np.concatenate((y_art, y))

        if _TOPLOT:
            plt.imshow(imgdata, origin='bottom', norm=LogNorm())
            plt.scatter(x_art, y_art, marker='x', color='k')
            plt.show()

        magstr = '{:<4.4f}'.format(np.mean(magnitudes)).replace('.', 'p')
        outimage = _ARTIMAGE.format(magstr).replace('temp', magstr)
        print '  Writing image with artificial stars: {}'.format(outimage)
        hdu = fits.PrimaryHDU(1)
        hdu.data = imgdata
        hdu.writeto(outimage, clobber=True)

        print indices
        print magnitudes

        outartsexfile = _ARTSEXFILE.format(magstr)
        print 'sex -c phot_t3.sex ../{} -CATALOG_NAME ../{} -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705 \n'.format(outimage, outartsexfile)


    if os.path.exists(artsexfile):
        fig_regions = _ARTPNG.format(magstr + '_each')
        fig_all = _ARTPNG.format(magstr + '_all')
        find_art_stars_allregions(artsexfile, _ARTFILES, _MAGS, _TOPLOT, fig_regions, fig_all)


def find_art_stars_allregions(artsexfile, artfiles, mags, toplot, fig_regions, fig_all):
    allartmags = np.array([])
    recovered = []
    artmagdiff = []

    fig = plt.figure(figsize=(10, 10))
    for i, artfile in enumerate(artfiles):
        ax = fig.add_subplot(4, 3, i + 1)

        artmags, recov, magdiff = find_art_stars(artsexfile, artfile, mags[i], toplot=toplot, writeto=toplot)

        if len(artmags) > 1:
            n, bins, patches = plt.hist(artmags, 50, color='b', alpha=0.5, label='{:<2.1f}%'.format(recov * 100.))
        ax.axvline(np.mean(artmags), color='k', label='mean {:<4.3}'.format(np.mean(artmags)))
        ax.axvline(mags[i], 0, 10, color='r', label='PSF mag {}'.format(mags[i]))
        ax.legend(loc=0, fontsize='small')
        ax.set_xlim(17, 24)

        allartmags = np.concatenate((allartmags, artmags))
        recovered.append(recov)
        artmagdiff.append(magdiff)

    plt.savefig(fig_regions)
    plt.show()

    plt.hist(allartmags, 50, color='b', alpha=0.5, label=str(np.sum(recovered * 100) / 12))
    for mag in mags:
        plt.axvline(mag, 0, 10, color='r')
    plt.xlim(17, 24)
    plt.savefig(fig_all)
    plt.legend()
    plt.show()

    print '  Recovered {}%'.format(np.sum(recovered * 100) / 12)
    print '  Mean magnitude offset {}'.format(np.nanmean(artmagdiff))


def find_art_stars(artsexfile, artfile, mag=None, max_sep=2., toplot=False, writeto=None):
    # fig = plt.figure(figsize=(16,8))
    n_art = 0.
    n_found = 0.

    # for j, artsexfile in enumerate(artsexfiles):
    sex = ascii.read(artsexfile, format='sextractor')
    x_pix, y_pix = np.loadtxt(artfile, unpack=True)
    n_art += len(x_pix)

    stars = sex[sex['MAG_APER'] < 99.]

    tree = KDTree(zip(stars['X_IMAGE'].ravel(), stars['Y_IMAGE'].ravel()))

    idx_found = []
    stars_idx = []
    for i in range(len(x_pix)):
        d, idx = tree.query((x_pix[i], y_pix[i]))
        if d < max_sep:
            idx_found.append(i)
            stars_idx.append(idx)
    if len(idx_found) == 0:
        print 'WARNING: No artificial stars found'
    # assert len(idx_found) > 0, 'No artificial stars found'
    n_found += len(idx_found)

    print '>> {}/{} recovered'.format(len(idx_found), len(x_pix))

    mags = stars['MAG_APER'][idx_found]

    mags_mean = np.mean(mags)
    print '  Mean magnitude measured, {:<2.3}'.format(mags_mean)
    print '  Magnitude shift, mag-mag_avg={:<2.3}'.format(float(mag) - mags_mean)

    if toplot or writeto:

        n, bins, patches = plt.hist(mags, 50, color='b', alpha=0.5, label='{}/{}'.format(len(idx_found), len(x_pix)*100))
        plt.axvline(mags_mean, color='k', label='mean {:<4.3}'.format(mags_mean))
        plt.axvline(mag, 0, 10, color='r', label='PSF mag {}'.format(mag))
        plt.legend(loc=0)
        plt.xlim(17, 24)
        # plt.tight_layout()

        if writeto is None:
            writeto = artsexfile.replace('.art1.txt', '.png').format('')
            print '  Writing to {}'.format(writeto)
            plt.savefig(writeto)
        if toplot:
            plt.show()

    return mags, n_found / n_art, float(mag) - mags_mean


def addstar(sexfile, image, artfile, (magmin, magmax), star_idx=None, (xmin, xmax)=(600., 4800.),
            (ymin, ymax)=(600., 4800.), mag_limit=0.01, chi2max=100., mincounts=-np.inf, maxcounts=np.inf,
            min_maxcounts=-np.inf, s=25., function='gauss', toplot=False):
    """
    Parse source extractor photometry file for stars within the ranges of magmin-magmax, xmin-xmax, ymin-ymax,
    mincounts-maxcounts. Plot 3D cutouts of each star which qualifies, prompt user for ID of selected star.
    Fit the selected star with a 2D-gaussian function, implant nstars # of artificial stars randomly into the image.
    Output coordinates of artificial stars into artfile, and image with artificial stars as artimage.
    """

    # open source extrator file as ascii table
    sextable = ascii.read(sexfile, format='sextractor')
    # retrieve fits data from image
    with fits.open(image) as hdu:
        imgdata = hdu[0].data

    # parse for stars in source extractor table with parameters in range of those specified
    stars_inx = sextable[(xmin < sextable['X_IMAGE']) & (sextable['X_IMAGE'] < xmax)]
    assert len(stars_inx) != 0, 'ERROR: No stars found within xrange'
    stars_inxy = stars_inx[(ymin < stars_inx['Y_IMAGE']) & (stars_inx['Y_IMAGE'] < ymax)]
    assert len(stars_inxy) != 0, 'ERROR: No stars found within, yrange'
    stars_mag = stars_inxy[(magmin < stars_inxy['MAG_APER']) & (stars_inxy['MAG_APER'] < magmax)]
    assert len(stars_mag) != 0, 'ERROR: No stars found within magrange'

    success = False
    while not success:
        stars = stars_mag[(stars_mag['MAGERR_APER'] / stars_mag['MAG_APER']) < mag_limit]
        if len(stars) > 0:
            success = True
        else:
            mag_limit += 0.01
    assert len(stars) != 0, 'ERROR: No stars found within magerr limit'

    if np.shape(stars)[0] == 1:
        star_idx = 0.
        print 'WARNING: Only one star found, turning plotting on'
        # toplot = True
    print '>>  {} Stars found'.format(np.shape(stars)[0])

    # If a particular star is not specified in the command line, plot potential stars until one chosen
    if star_idx is None:
        print stars
        print '  Finding star PSFs to model: type a number to select the corresponding indexed star, "enter" to continue, any string to kill'
        star_idx = False
        i = 0
        counter = 0  # not used for anything except having something in the exception line
        print '    idx      x [pix]    y [pix]    mag    chi2/dof'
        while not star_idx:
            x, y = stars['X_IMAGE'][i], stars['Y_IMAGE'][i]
            cutout = imgdata[int(y - s / 2.):int(y + s / 2.), int(x - s / 2.):int(x + s / 2.)]

            if (np.min(cutout) > mincounts) & (np.max(cutout) < maxcounts) & (np.max(cutout) > min_maxcounts):

                try:  # use a try/except statements to avoid errors if fit fails, just skip those unable to be fit
                    xx, yy = np.meshgrid(np.arange(s), np.arange(s))
                    if function == 'gauss':
                        popt, pcov = curve_fit(gaussian_2d, (xx, yy), cutout.ravel(),
                                               p0=[np.max(cutout), s / 2., s / 2., s / 2., s / 2., 0., 0.])
                        data_fitted = gaussian_2d((xx, yy), *popt)  # .reshape(s, s)
                    if function == 'lorentz':
                        popt, pcov = curve_fit(lorentz_2d, (xx, yy), cutout.ravel(),
                                               p0=[np.max(cutout), s / 2., s / 2., s / 50., s / 50., 0., 0.])
                        data_fitted = lorentz_2d((xx, yy), *popt)  # .reshape(s, s)

                    # determine "goodness of fit" by chi square statistic
                    dof = len(data_fitted) + len(popt)
                    chi2, p = chisquare(data_fitted, cutout.ravel())

                    # only consider those with chi squares less than user specifed value
                    if np.abs(chi2 / dof) < chi2max:
                        print '   {} {} {} {} {}'.format(i, stars['X_IMAGE'][i], stars['Y_IMAGE'][i],
                                                         stars['MAG_APER'][i], np.abs(chi2 / dof), np.max(cutout))

                        fig = plt.figure()
                        ax = fig.add_subplot(111, projection='3d')
                        xx, yy = np.meshgrid(np.arange(s), np.arange(s))
                        ax.plot_wireframe(xx, yy, cutout, rstride=1, cstride=1)
                        plt.show()

                        star_idx = raw_input('Star index: ')
                except:
                    counter += 1

            i += 1
            if i >= np.shape(stars)[0]:
                print "ERROR: exceeded number of stars, restarting loop"
                i = 0

    try:
        star = stars[int(star_idx)]
    except:
        raise Exception, 'Invalid star index specified'

    print '*******'
    print '  Star selected for PSF modeling:'
    print '  idx  x [pix]   y [pix]   mag'
    print '  {}   {} {} {}'.format(star_idx, star['X_IMAGE'], star['Y_IMAGE'], star['MAG_APER'])

    # cutout a region around the chosen star
    x, y = star['X_IMAGE'], star['Y_IMAGE']
    cutout = imgdata[int(y - s / 2.):int(y + s / 2.), int(x - s / 2.):int(x + s / 2.)]

    # fit a 2D gaussian to the data
    xx, yy = np.meshgrid(np.arange(s), np.arange(s))
    if function == 'gauss':
        popt, pcov = curve_fit(gaussian_2d, (xx, yy), cutout.ravel(),
                               p0=[np.max(cutout), s / 2., s / 2., s / 2., s / 2., 0., 0.])
        data_fitted = gaussian_2d((xx, yy), *popt).reshape(s, s)
        # create a 2D gaussian with (optimized) parameters returned from the fit, but centered
        star_centered = gaussian_2d((xx, yy), popt[0], s / 2., s / 2., popt[3], popt[4], popt[5], popt[6]).reshape(s, s)

    if function == 'lorentz':
        popt, pcov = curve_fit(lorentz_2d, (xx, yy), cutout.ravel(),
                               p0=[np.max(cutout), s / 2., s / 2., s / 50., s / 50., 0., 0.])
        data_fitted = lorentz_2d((xx, yy), *popt).reshape(s, s)
        star_centered = lorentz_2d((xx, yy), popt[0], s / 2., s / 2., popt[3], popt[4], popt[5], popt[6]).reshape(s, s)

    print '  Optimal parameters, A, x0, y0, sigma_x, sigma_y, theta, offset:  '
    print '    {} {} {} {} {} {} {}'.format(popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
    print '*******'

    if toplot:
        fig = plt.figure()
        ax = fig.add_subplot(131)
        ax.imshow(imgdata, cmap=plt.cm.jet, norm=LogNorm(), origin='bottom')
        ax.scatter(sextable['X_IMAGE'], sextable['Y_IMAGE'], marker='x', color='g')
        ax.scatter(x, y, marker='o', color='k')
        ax.set_xlim([int(x - s / 2.), int(x + s / 2.)])
        ax.set_ylim([int(y - s / 2.), int(y + s / 2.)])
        ax = fig.add_subplot(132)
        ax.hold(True)
        ax.imshow(cutout, cmap=plt.cm.jet, origin='bottom', extent=(xx.min(), xx.max(), yy.min(), yy.max()))
        ax.contour(xx, yy, data_fitted, 8, colors='w')
        ax = fig.add_subplot(133, projection='3d')
        ax.plot_wireframe(xx, yy, cutout, rstride=1, cstride=1, color='blue')
        ax.plot_wireframe(xx, yy, data_fitted, rstride=1, cstride=1, color='red')
        plt.show()

    # implant the gaussian "artificial star" at each location
    x_art, y_art = np.loadtxt(artfile, unpack=True, delimiter=' ')
    art_imgdata = imgdata*0.
    for n in range(len(x_art)):
        art_imgdata[y_art[n] - s / 2.:y_art[n] + s / 2., x_art[n] - s / 2.:x_art[n] + s / 2.] += star_centered

    return art_imgdata, x_art, y_art, int(star_idx), star['MAG_APER']


def generate_coords(artfile, nstars=100., (xmin, xmax)=(600., 4800.), (ymin, ymax)=(600., 4800.)):
    x = random_sample(nstars, min=xmin, max=xmax)
    y = random_sample(nstars, min=ymin, max=ymax)

    np.savetxt(artfile, np.transpose([x, y]), delimiter=' ')


def random_sample(n, min=None, max=None):
    if min is None:
        min = 0.
    if max is None:
        max = nstars
    idx = np.random.randint(min, max, n)
    return idx


def gaussian_2d((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
    g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                                       + c * ((y - yo) ** 2)))
    return g.ravel()


def lorentz_2d((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (sigma_x ** 2) + (np.sin(theta) ** 2) / (sigma_y ** 2)
    b = -(np.sin(2 * theta)) / (2 * sigma_x ** 2) + (np.sin(2 * theta)) / (2 * sigma_y ** 2)
    c = (np.sin(theta) ** 2) / (sigma_x ** 2) + (np.cos(theta) ** 2) / (sigma_y ** 2)
    g = offset + amplitude * (sigma_x ** 2 + sigma_y ** 2) / (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
                                                              + c * ((y - yo) ** 2) + (sigma_x ** 2 + sigma_y ** 2))
    return g.ravel()


if __name__ == '__main__':
    main()
