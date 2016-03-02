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

python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_21p0074.art{}.fits -mags 21.006 21.1007
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_21p5337.art{}.fits -mags 21.5336 21.5338 --maxmagerr 0.1
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_22p0609.art{}.fits -mags 22.0608 22.061 --maxmagerr 0.1
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_22p0963.art{}.fits -mags 22.0962 22.0964 --maxmagerr 0.1
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_22p2896.art{}.fits -mags 22.2895 22.2897 --maxmagerr 0.1
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_22p6087.art{}.fits -mags 22.6086 22.6088 --maxmagerr 0.1
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_22p6291.art{}.fits -mags 22.6290 22.6292 --maxmagerr 0.1
python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -artimg mag_lim/CB68_J_sub_23p1608.art{}.fits -mags 23.1607 23.1609 --maxmagerr 0.2

sex -c phot_t3.sex ../mag_lim/CB68_J_sub_21p0074.art1.fits -CATALOG_NAME ../mag_lim/CB68_J_sex_21p0074.art1.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705
sex -c phot_t3.sex ../mag_lim/CB68_J_sub_21p0074.art2.fits -CATALOG_NAME ../mag_lim/CB68_J_sex_21p0074.art2.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705
sex -c phot_t3.sex ../mag_lim/CB68_J_sub_21p0074.art3.fits -CATALOG_NAME ../mag_lim/CB68_J_sex_21p0074.art3.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705
sex -c phot_t3.sex ../mag_lim/CB68_J_sub_21p0074.art4.fits -CATALOG_NAME ../mag_lim/CB68_J_sex_21p0074.art4.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705
sex -c phot_t3.sex ../mag_lim/CB68_J_sub_21p0074.art5.fits -CATALOG_NAME ../mag_lim/CB68_J_sex_21p0074.art5.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705

python addstar.py -artsf mag_lim/CB68_J_sex_21p0074.art{}.txt --mag 21.0074
python addstar.py -artsf mag_lim/CB68_J_sex_21p5337.art{}.txt --mag 21.5337
python addstar.py -artsf mag_lim/CB68_J_sex_22p0609.art{}.txt --mag 22.0609
python addstar.py -artsf mag_lim/CB68_J_sex_22p0963.art{}.txt --mag 22.0963
python addstar.py -artsf mag_lim/CB68_J_sex_22p2896.art{}.txt --mag 22.2896
python addstar.py -artsf mag_lim/CB68_J_sex_22p6087.art{}.txt --mag 22.6087
python addstar.py -artsf mag_lim/CB68_J_sex_22p6291.art{}.txt --mag 22.6291
python addstar.py -artsf mag_lim/CB68_J_sex_23p1608.art{}.txt --mag 23.1608

"""


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

    parser = argparse.ArgumentParser(
        description='Generate a PSF model from a ideal star, implant nstars number into the given image. Determine the number of artificial stars recovered.')
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default=None,
                        help='Source extractor file.')
    parser.add_argument("--image", '-img',
                        action='store',
                        default=None,
                        help="Fits image file.")
    parser.add_argument("--artimage", '-artimg',
                        action='store',
                        default=None,
                        help='Output images format with artificial stars added.')
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
    parser.add_argument("--artfiles", '-artfiles',
                        nargs='*',
                        action='store',
                        default=['mag_lim/art1.coo', 'mag_lim/art2.coo', 'mag_lim/art3.coo', 'mag_lim/art4.coo',
                                 'mag_lim/art5.coo'],
                        help='List of output file of coordinates of artificial stars. Must be same length as artimages')
    parser.add_argument("--star_idx", '-si',
                        action='store',
                        default=None,
                        help='Index of star with PSF chosen to build model from, WARNING: index will be from cropped table of stars in range of specified parameters. Default: None')
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
                        help='True: plot all figures, False: plot only when finding PSF star. Default: False')
    parser.add_argument("--mag",
                        action='store',
                        default=None,
                        help='Magnitude of chosen PSF star. Used for plotting purposes only. Default: None')
    parser.add_argument("--maxmagerr",
                        action='store',
                        default=0.05,
                        help='Maximum err level to continue, percent eg 0.05 is 5 percent. Default: 0.05')
    parser.add_argument("--function",
                        action='store',
                        default='gauss',
                        help='Function used to fit PSF. Options "gauss", "lorentz". Default: gauss')
    parser.add_argument("--chi2max",
                        action='store',
                        default=100.,
                        help='Only display PSFs which have a chi2 less than this value when fitting a 2d gaussian. Default: 100')
    args = parser.parse_args()

    if args.artsexfile is None:
        if args.artimage is None:
            raise Exception, 'No output image file format is given'
        artimages = []
        for i, artfile in enumerate(args.artfiles):
            if not os.path.exists(artfile):
                generate_coords(artfile, int(args.nstars), args.xrange, args.yrange)
            artimages.append(args.artimage.format(i + 1))
        if args.sexfile is None:
            raise Exception, 'No input source extractor file given'
        if args.image is None:
            raise Exception, 'No input image file given'
        if args.magrange is None:
            raise Exception, 'No magnitude range given'
        assert len(args.magrange) == 2, 'Magnitude input is two values, min and max'
        assert args.function == 'gauss' or args.function == 'lorentz', 'ERROR: function not "gauss" or "lorentz"'

        addstar(args.sexfile, args.image, args.artfiles, artimages, args.magrange, args.star_idx,
                args.dimension, args.xrange, args.yrange, float(args.mincounts),
                float(args.maxcounts), float(args.min_maxcounts), args.toplot, mag_limit=float(args.maxmagerr),
                chi2max=float(args.chi2max), function=args.function)

    if args.artsexfile is not None:
        artsexfiles = []
        for i, artfile in enumerate(args.artfiles):
            artsexfiles.append(args.artsexfile.format(i + 1))
        find_art_stars(artsexfiles, args.artfiles, args.mag, args.toplot)


def find_art_stars(artsexfiles, artfiles, mag=None, max_sep=2., toplot=False):
    fig = plt.figure(figsize=(16,8))
    n_art = 0.
    n_found = 0.
    N = 0
    for str in artsexfiles:
        N += 1

    for j, artsexfile in enumerate(artsexfiles):
        sex = ascii.read(artsexfile, format='sextractor')
        x_pix, y_pix = np.loadtxt(artfiles[j], unpack=True)
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
        assert len(idx_found) > 0, 'No artificial stars found'
        n_found += len(idx_found)

        print '{}/{} recovered'.format(len(idx_found),len(x_pix))

        ax = fig.add_subplot(1, N, j+1)
        ax.hist(stars['MAG_APER'][stars_idx], 100, color='b', alpha=0.5, label='{}/{} recovered'.format(len(idx_found),len(x_pix)))
        ax.axvline(mag, 0, 10, color='r', label='PSF mag {}'.format(mag))
        ax.legend()
        ax.set_xlim(17, 24)

    print '{}/{} recovered, {} percent'.format(n_found, n_art, n_found / n_art * 100.)

    plt.tight_layout()
    print '  Writing to {}'.format(artsexfile[0].replace('.coo', '.png').format(''))
    plt.savefig(artsexfile[0].replace('.coo', '.png').format(''))
    if toplot:
        plt.show()


def addstar(sexfile, image, artfiles, artimages, (magmin, magmax), star_idx=None, s=25.,
            (xmin, xmax)=(600., 4800.), (ymin, ymax)=(600., 4800.), mincounts=-np.inf,
            maxcounts=np.inf, min_maxcounts=-np.inf, toplot=False, mag_limit=0.01, chi2max=100., function='gauss'):
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
    stars = stars_mag[(stars_mag['MAGERR_APER'] / stars_mag['MAG_APER']) < mag_limit]
    assert len(stars) != 0, 'ERROR: No stars found within magerr limit'

    if np.shape(stars)[0] == 1:
        star_idx = 0.
        print 'WARNING: Only one star found, turning plotting on'
        # toplot = True

    print_again = False

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
                        print '{} {} {} {} {}'.format(i, stars['X_IMAGE'][i], stars['Y_IMAGE'][i], stars['MAG_APER'][i],
                                                      np.abs(chi2 / dof))
                        if toplot:
                            fig = plt.figure()
                            ax = fig.add_subplot(121, projection='3d')
                            xx, yy = np.meshgrid(np.arange(s), np.arange(s))
                            ax.plot_wireframe(xx, yy, cutout, rstride=1, cstride=1)
                            ax = fig.add_subplot(122)
                            ax.hold(True)
                            ax.imshow(cutout, cmap=plt.cm.jet, origin='bottom',
                                      extent=(xx.min(), xx.max(), yy.min(), yy.max()))
                            ax.contour(xx, yy, data_fitted, 8, colors='w')
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
            if i > np.shape(stars)[0]:
                print "ERROR: exceeded number of stars, restarting loop"
                i = 0.
    else:
        print_again = True

    try:
        star = stars[int(star_idx)]
    except:
        raise Exception, 'Invalid star index specified'

    if print_again:
        print 'Star selected for PSF modeling:'
        print '   idx     x [pix]   y [pix]   mag'
        print '{} {} {} {}'.format(star_idx, star['X_IMAGE'], star['Y_IMAGE'], star['MAG_APER'])

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
        star_centered = gaussian_2d((xx, yy), popt[0], s / 2., s / 2., popt[3], popt[4], popt[5], 0.).reshape(s, s)

    if function == 'lorentz':
        popt, pcov = curve_fit(lorentz_2d, (xx, yy), cutout.ravel(),
                               p0=[np.max(cutout), s / 2., s / 2., s / 50., s / 50., 0., 0.])
        data_fitted = lorentz_2d((xx, yy), *popt).reshape(s, s)
        star_centered = lorentz_2d((xx, yy), popt[0], s / 2., s / 2., popt[3], popt[4], popt[5], 0.).reshape(s, s)

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

    for i, artfile in enumerate(artfiles):
        # implant the gaussian "artificial star" at each location
        x, y = np.loadtxt(artfile, unpack=True, delimiter=' ')
        art_imgdata = imgdata
        for n in range(len(x)):
            art_imgdata[y[n] - s / 2.:y[n] + s / 2., x[n] - s / 2.:x[n] + s / 2.] += star_centered

        if toplot:
            plt.imshow(art_imgdata, origin='bottom', norm=LogNorm())
            plt.scatter(x, y, marker='x', color='k')
            plt.show()

        print '  Writing image with artificial stars'
        hdu = fits.PrimaryHDU(1)
        hdu.data = art_imgdata
        hdu.writeto(artimages[i].format(str(star['MAG_APER']).replace('.', 'p'), i), clobber=True)


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
