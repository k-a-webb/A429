import argparse
import os
import numpy as np
from astropy.io import ascii, fits
from scipy.spatial import KDTree
from scipy.stats import chisquare
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

"""
Create artificial stars from a good candidate star in the image. Determine number of stars recovered.

python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artfile mag_lim_psfex/CB68_J_sub_art{}.coo -artimg mag_lim_psfex/CB68_J_sub.art{}.fits -mags 18. 18.5 --chi2max 30.

sex -c phot_t3.sex ../mag_lim_psfex/CB68_J_sub.art18.254.fits -CATALOG_NAME ../mag_lim_psfex/CB68_J_sex_t3_ap30_add_18.254.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705

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
        description='Generate a PSF model from a ideal star, implant nstars number into the given image. Determine the '
                    'number of artificial stars recovered.')
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default=None,
                        help='Source extractor file.')
    parser.add_argument("--image", '-img',
                        action='store',
                        default=None,
                        help='Fits image file.')
    parser.add_argument("--artfile", '-artfile',
                        action='store',
                        default=None,
                        help='Output file of coordinates of artificial stars.')
    parser.add_argument("--artimage", '-artimg',
                        action='store',
                        default=None,
                        help='Output image with artificial stars added.')
    parser.add_argument("--magrange", '-mags',
                        nargs='*',
                        action='store',
                        type=float,
                        default=None,
                        help='Magnitude range of stars, must be a list of floats separated by a space, eg "18. 19.".')

    parser.add_argument("--artsexfile", '-artsf',
                        action='store',
                        default=None,
                        help='Source extractor file with artificial stars. Used in second function only.')

    # Optional, note default values
    parser.add_argument("--star_idx", '-si',
                        action='store',
                        type=float,
                        default=None,
                        help='Index of star with PSF chosen to build model from. Default: None')
    parser.add_argument("--xrange", '-xs',
                        nargs='*',
                        action='store',
                        type=float,
                        default=(600., 4800.),
                        help='x range [pixels], Default: 600. 4800.')
    parser.add_argument("--yrange", '-ys',
                        nargs='*',
                        action='store',
                        type=float,
                        default=(600., 4800.),
                        help='y range [pixels], Default: 600. 4800.')
    parser.add_argument("--dimension", '-d',
                        action='store',
                        default=25.,
                        type=float,
                        help='Dimension of PSF box, Default: 25.')
    parser.add_argument("--nstars", '-n',
                        action='store',
                        default=100,
                        type=int,
                        help='Number of artificial stars to add. Default: 100')
    parser.add_argument("--maxcounts",
                        action='store',
                        type=float,
                        default=np.inf,
                        help='Select only stars with counts less than this value [counts]. Default: inf')
    parser.add_argument("--mincounts",
                        action='store',
                        type=float,
                        default=-np.inf,
                        help='Select only stars with counts greater than this value [counts]. Default: -inf')
    parser.add_argument("--min_maxcounts",
                        action='store',
                        default=-np.inf,
                        help='Select only stars with maximum counts greater than this value [counts]. Default: -inf')
    parser.add_argument("--toplot",
                        action='store',
                        default=False,
                        help='True: plot all figures, False: plot only when finding PSF star. Defalut: False')
    parser.add_argument("--mag",
                        action='store',
                        default=None,
                        help='Magnitude of chosen PSF star. Used for plotting purposes only. Default: None')
    parser.add_argument("--chi2max",
                        action='store',
                        default=100.,
                        help='Only display PSFs which have a chi2 less than this value when fitting a 2d gaussian. Default: 100')
    args = parser.parse_args()

    if args.artfile is None:
        raise Exception, 'No coordinate file given for artificial stars'

    if not os.path.exists(args.artfile):
        if args.sexfile is None:
            raise Exception, 'No input source extractor file given'
        if args.image is None:
            raise Exception, 'No input image file given'
        if args.artimage is None:
            raise Exception, 'No output image file given'
        if args.magrange is None:
            raise Exception, 'No magnitude range given'
        assert len(args.magrange) == 2, 'Magnitude input is two values, min and max'

        addstar(args.sexfile, args.image, args.artfile, args.artimage, args.magrange, args.star_idx,
                args.dimension, int(args.nstars), args.xrange, args.yrange, float(args.mincounts),
                float(args.maxcounts), float(args.min_maxcounts), args.toplot, chi2max=float(args.chi2max))
    else:
        'File {} already exists'.format(args.artfile)

    if args.artsexfile is not None:
        find_art_stars(args.artsexfile, args.artfile, args.mag)


def find_art_stars(artsexfile, artfile, mag=None, max_sep=2.):
    sex = ascii.read(artsexfile, format='sextractor')
    x_pix, y_pix = np.loadtxt(artfile, unpack=True)

    tree = KDTree(zip(sex['X_IMAGE'].ravel(), sex['Y_IMAGE'].ravel()))

    idx_found = []
    sex_idx = []
    for i in range(len(x_pix)):
        d, idx = tree.query((x_pix[i], y_pix[i]))
        if d < max_sep:
            idx_found.append(i)
            sex_idx.append(idx)
    assert len(idx_found) > 0, 'No artificial stars found'
    print '{}/100 recovered'.format(len(idx_found))

    # plt.hist(np.array(mag), 100, color='g')
    plt.hist(sex['MAG_APER'][sex_idx], 100, color='b', alpha=0.5)
    if mag is not None:
        plt.axvline(mag, 0, 10, color='r')
    # plt.xlim(10, 30)
    plt.show()


def addstar(sexfile, image, artfile, artimage, (magmin, magmax), star_idx=None, s=25., nstars=100,
            (xmin, xmax)=(600., 4800.), (ymin, ymax)=(600., 4800.), mincounts=-np.inf,
            maxcounts=np.inf, min_maxcounts=-np.inf, toplot=False, mag_limit=0.01, chi2max=100.):
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
    stars_inxy = stars_inx[(ymin < stars_inx['Y_IMAGE']) & (stars_inx['Y_IMAGE'] < ymax)]
    stars_mag = stars_inxy[(magmin < stars_inxy['MAG_APER']) & (stars_inxy['MAG_APER'] < magmax)]
    stars = stars_mag[(stars_mag['MAGERR_APER'] / stars_mag['MAG_APER']) < mag_limit]

    # If a particular star is not specified in the command line, plot potential stars until one chosen
    print 'Finding star PSFs to model: type a number to select the corresponding indexed star, "enter" to continue, any string to kill'
    if star_idx is None:
        star_idx = False
        i = 0
        counter = 0  # not used for anything except having something in the exception line
        print '   idx     x [pix]   y [pix]   mag   chi2/dof'
        while not star_idx:
            x, y = stars['X_IMAGE'][i], stars['Y_IMAGE'][i]
            cutout = imgdata[int(y - s / 2.):int(y + s / 2.), int(x - s / 2.):int(x + s / 2.)]

            if (np.min(cutout) > mincounts) & (np.max(cutout) < maxcounts) & (np.max(cutout) > min_maxcounts):

                try:  # use a try/except statements to avoid errors if fit fails, just skip those unable to be fit
                    # fit the data in the cutout with a 2d gaussian
                    xx, yy = np.meshgrid(np.arange(s), np.arange(s))
                    popt, pcov = curve_fit(gaussian_2d, (xx, yy), cutout.ravel(),
                                           p0=[np.max(cutout), s / 2., s / 2., s / 2., s / 2., 0., 0.])
                    data_fitted = gaussian_2d((xx, yy), *popt)  # .reshape(s, s)

                    # determine "goodness of fit" by chi square statistic
                    dof = len(data_fitted) + len(popt)
                    chi2, p = chisquare(data_fitted, cutout.ravel(), ddof=dof)

                    # only consider those with chi squares less than user specifed value
                    if np.abs(chi2 / dof) < chi2max:

                        print '{} {} {} {} {}'.format(i, stars['X_IMAGE'][i], stars['Y_IMAGE'][i], stars['MAG_APER'][i],
                                                      np.abs(chi2 / dof))

                        if toplot:
                            fig = plt.figure()
                            ax = fig.add_subplot(121, projection='3d')
                            xx, yy = np.meshgrid(np.arange(s), np.arange(s))
                            ax.plot_wireframe(xx, yy, cutout, rstride=1, cstride=1)
                            ax = fig.add_subplot(122)
                            ax.imshow(np.log(imgdata), cmap=plt.cm.jet)
                            ax.scatter(sextable['X_IMAGE'], sextable['Y_IMAGE'], marker='x', color='g')
                            ax.scatter(x, y, marker='o', color='k')
                            ax.set_xlim([int(x - s / 2.), int(x + s / 2.)])
                            ax.set_ylim([int(y - s / 2.), int(y + s / 2.)])
                            plt.show()
                        else:
                            fig = plt.figure()
                            ax = fig.add_subplot(111, projection='3d')
                            xx, yy = np.meshgrid(np.arange(s), np.arange(s))
                            ax.plot_wireframe(xx, yy, cutout, rstride=1, cstride=1)
                            plt.show()

                        star_idx = raw_input('Star index: ')

                except:
                    counter += 1
            i += 1

    try:
        star = stars[int(star_idx)]
    except:
        raise Exception, 'Invalid star index specified'

    # cutout a region around the chosen star
    x, y = star['X_IMAGE'], star['Y_IMAGE']
    cutout = imgdata[int(y - s / 2.):int(y + s / 2.), int(x - s / 2.):int(x + s / 2.)]

    # fit a 2D gaussian to the data
    xx, yy = np.meshgrid(np.arange(s), np.arange(s))
    popt, pcov = curve_fit(gaussian_2d, (xx, yy), cutout.ravel(),
                           p0=[np.max(cutout), s / 2., s / 2., s / 2., s / 2., 0., 0.])

    if toplot:
        data_fitted = gaussian_2d((xx, yy), *popt).reshape(s, s)
        fig = plt.figure()
        ax = fig.add_subplot(131)
        ax.imshow(np.log(imgdata), cmap=plt.cm.jet,
                  origin='bottom')  # , extent=(xx.min(), xx.max(), yy.min(), yy.max()))
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

    # create a 2D gaussian with (optimized) parameters returned from the fit, but centered
    star_centered = gaussian_2d((xx, yy), popt[0], s / 2., s / 2., popt[3], popt[4], popt[5], popt[6]).reshape(s, s)

    if toplot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(xx, yy, star_centered, rstride=1, cstride=1, color='red')
        plt.show()

    # create random number of stars
    print xmin, xmax, ymin, ymax, len(sextable['X_IMAGE'])

    x = random_sample(nstars, min=xmin, max=xmax)
    y = random_sample(nstars, min=ymin, max=ymax)

    # implant the gaussian "artificial star" at each location
    art_imgdata = imgdata
    for n in range(nstars):
        art_imgdata[y[n] - s / 2.:y[n] + s / 2., x[n] - s / 2.:x[n] + s / 2.] += star_centered

    if toplot:
        plt.imshow(np.log(art_imgdata), origin='bottom', cmap='rainbow')
        plt.scatter(x, y, marker='x', color='k')
        plt.show()

    print '  Writing coordinates of artificial stars'
    np.savetxt(artfile.format(str(star['MAG_APER']).replace('.', '-')), np.transpose([x, y]), delimiter=' ')

    print '  Writing image with artificial stars'
    hdu = fits.PrimaryHDU(1)
    hdu.data = art_imgdata
    hdu.writeto(artimage.format(str(star['MAG_APER']).replace('.', '-')), clobber=True)


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


if __name__ == '__main__':
    main()
