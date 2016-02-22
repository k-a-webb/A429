import argparse
from pyraf import iraf
import os
import numpy as np
from astropy.io import ascii, fits
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

"""
Create input for DAOPHOT ADDSTAR to generate artificial stars

python addstar.py -sf phot_t3/CB68_J_sex_t3_ap30.txt -img CB68/CB68_J_sub.fits -ap 30 -zp 29.80705
python addstar.py -img release/CB68_J_sub.fits -sf phot_t3/CB68_J_sex_t3_ap30.txt -path mag_lim/ -ap 30 -zp 29.80705
python ../addstar.py -img ../release/CB68_J_sub.fits -sf CB68_J_sex_t3_ap30.txt -path ../mag_lim/ -ap 30 -zp 29.80705

sex -c phot_t3.sex ../mag_lim/CB68_J_sub.add.fits -CATALOG_NAME ../mag_lim/CB68_J_sex_t3_ap30_add.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705

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
    parser.add_argument("--outpath", '-path',
                        action='store',
                        default=None,
                        help='File path for output files.')
    parser.add_argument("--artphotfile", '-artphot',
                        action='store',
                        default=None,
                        help='Input file of source extractor file with artificial stars added.')
    args = parser.parse_args()

    if args.sexfile is None:
        raise Exception, 'No input source extractor file given'
    if args.image is None:
        raise Exception, 'No input image file given'

    artphotfile = run_daophot(args.sexfile, args.image, args.outpath, args.artphotfile, args.aperture, args.zeropoint)

    print 'sex -c phot_t3.sex CB68_J_sub.add.fits -CATALOG_NAME CB68_J_sex_t3_ap30_add.txt -PHOT_APERTURES 30 -MAG_ZEROPOINT 29.80705'

    find_art_stars(args.addphot, artphotfile)


def find_art_stars(addphot, artphotfile, max_sep=0.1):
    sex = ascii.read(addphot, format='sextractor')
    x_wcs, y_wcs, mag = np.loadtxt(artphotfile, unpack=True)

    tree = KDTree(zip(sex['X_WORLD'].ravel(), sex['Y_WORLD'].ravel()))

    idx_found = []
    for i in range(len(x_wcs)):
        d, idx = tree.query((x_wcs[i], y_wcs[i]))
        if d < max_sep:
            idx_found.append(i)
    assert len(idx_found) > 0, 'No artificial stars found'

    plt.hist(np.array(mag), 100, color='g')
    plt.hist(np.array(mag)[idx_found], 100, color='b', alpha=0.5)
    plt.xlim(10, 30)
    plt.show()


def run_daophot(sexfile, image, outpath, artphotfile, aperture=None, zeropoint=None, maxnpsf=50., nstars=3000):
    try:
        prefix_img = (image.split('/')[-1]).split('.')[0]
    except IndexError:
        prefix_img = sexfile.split('.')[0]

    coofile = os.path.join(outpath, prefix_img + '.coo')
    coofile2 = os.path.join(outpath, prefix_img + '.coo.2')
    photfile = os.path.join(outpath, prefix_img + '.phot')
    sexphotfile = os.path.join(outpath, prefix_img + '.mag')
    pstfile = os.path.join(outpath, prefix_img + '.pst')
    opstfile = os.path.join(outpath, prefix_img + '.opst')
    groupfile = os.path.join(outpath, prefix_img + '.psg')
    psfimage = os.path.join(outpath, prefix_img + '.psf')
    addimage = os.path.join(outpath, prefix_img + '.add')

    iraf.daophot()

    # daofind - automatically detect objects in an image
    if not os.path.exists(coofile):
        iraf.daophot.daofind(image=image, output=coofile)
    print 'done with {}'.format(coofile)

    # remove saturated stars from coordinate file
    if not os.path.exists(coofile2):
        remove_saturated(image, sexfile, coofile, coofile2)
    print 'done with {}'.format(coofile2)

    # phot -- do aperture photometry on a list of stars
    if not os.path.exists(photfile):
        raise Exception, 'STOP HERE AND RUN DAOPHOT PHOT IN IRAF, bkgstd=100, fwhm=3'
        # if aperture is None:
        #     raise Exception, 'No input aperture parameter given'
        # if zeropoint is None:
        #     raise Exception, 'No input zeropoint value given'
        # iraf.daophot.phot(image, coords=coofile2, output=photfile, wcsin="world", apertures=float(aperture),
        #                   zmag=float(zeropoint), verbose=True)
    print 'done with {}'.format(photfile)

    # replace phot values with those from source extractor
    if not os.path.exists(sexphotfile):
        write_phot(photfile, sexphotfile, sexfile)
    print 'done with {}'.format(sexphotfile)

    # pstselect -- select candidate psf stars from a photometry file
    if not os.path.exists(pstfile):
        iraf.daophot.pstselect(image, sexphotfile, pstfile, maxnpsf)
    print 'done with {}'.format(pstfile)

    # psf -- build the point spread function for an image
    if not os.path.exists(psfimage):
        raise Exception, 'STOP HERE AND RUN DAOPHOT PSF IN IRAF, max=50000, rad=11'
        # iraf.daophot.psf(image, photfile=sexphotfile, pstfile=pstfile, psfimage=psfimage, opstfile=opstfile,
        #                  groupfile=groupfile, interactive=False)
    print 'done with {}'.format(psfimage)

    # define parameters of artificial stars
    if not os.path.exists(artphotfile):
        sex_data = ascii.read(sexfile, format='sextractor')
        x = random_sample(nstars, sex_data['X_WORLD'])
        y = random_sample(nstars, sex_data['Y_WORLD'])
        mags = random_sample(nstars, np.linspace(18., 28., nstars))
        np.savetxt(artphotfile, np.transpose([x, y, mags]), delimiter=' ')

    # addstar -- add artificial stars to images
    if not os.path.exists(addimage):
        iraf.daophot.addstar(image, artphotfile, psfimage, addimage)  # 50% recovery

    return artphotfile


def random_sample(n, arr, min=None, max=None):
    if min is None:
        min = 0.
    if max is None:
        max = len(arr)
    idx = np.random.randint(min, max, n)
    return arr[idx]


def write_phot(photfile, newphotfile, sexfile):
    phot = ascii.read(photfile, format='daophot')
    sex = ascii.read(sexfile, format='sextractor')

    tree = KDTree(zip(phot['XCENTER'].ravel(), phot['YCENTER'].ravel()))
    max_sep = 0.1

    sex_idxs = []
    phot_idxs = []
    for i in range(len(sex)):
        d, ii = tree.query((sex['X_IMAGE'][i], sex['Y_IMAGE'][i]), distance_upper_bound=max_sep)
        if (d != np.inf):
            sex_idxs.append(i)
            phot_idxs.append(ii)

    sex_matched = sex[sex_idxs]
    phot_matched = phot[phot_idxs]

    with open(newphotfile, 'w') as outfile:
        with open(photfile, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)

        for i in range(len(sex_idxs)):
            if (phot_matched['MAG'][i] > 0.) and (sex['MAG_APER'][i] < 99.):
                outfile.write("{:<23}{:<10.3f}{:<10.3f}{:<6d}{:<23}{:<6d} \\\n".format(
                    phot_matched['IMAGE'][i], phot_matched['XINIT'][i], phot_matched['YINIT'][i], phot_matched['ID'][i],
                    phot_matched['COORDS'][i], phot_matched['LID'][i]))
                outfile.write("   {:<11.3f}{:<11.3f}0.000   0.000   INDEF   INDEF          0    NoError   \\\n".format(
                    phot_matched['XCENTER'][i], phot_matched['YCENTER'][i]))
                outfile.write("   {:<15.7g}{:<15.7g}{:<15.7g}{:<7d}{:<9d}0    NoError   \\\n".format(
                    phot_matched['MSKY'][i], phot_matched['STDEV'][i], phot_matched['SSKEW'][i],
                    phot_matched['NSKY'][i], phot_matched['NSREJ'][i]))
                outfile.write("   1.             INDEF          INDEF                  INDEF                  \\\n")
                outfile.write("   8.00     {:<14.7g}{:<11.7g}{:<14.7g}{:<7.3f}{:<6.3f}0    NoError    \n".format(
                    phot_matched['SUM'][i], phot_matched['AREA'][i], sex_matched['FLUX_APER'][i],
                    sex_matched['MAG_APER'][i], sex_matched['MAGERR_APER'][i]))

    return


def remove_saturated(image, sexfile, coofile, newcoofile, max_zeros=1.):
    sex_data = ascii.read(sexfile, format='sextractor')
    coo_data = ascii.read(coofile, format='daophot')

    with fits.open(image) as hdu:
        img_data = hdu[0].data

    x = sex_data['X_IMAGE']
    y = sex_data['Y_IMAGE']
    r = sex_data['KRON_RADIUS']

    idx = []
    for i in range(len(x)):
        zeros = []
        for xx in range(int(x[i]) - int(r[i]), int(x[i]) + int(r[i])):
            yy_p = np.sqrt(r[i] ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            yy_m = -np.sqrt(r[i] ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            for yy in range(int(yy_m), int(yy_p)):
                if img_data[yy, xx] == 0.:
                    zeros.append(1.)
        if len(zeros) <= max_zeros:
            idx.append(i)

    with open(newcoofile, 'w') as outfile:
        with open(coofile, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
        for i in range(len(idx)):
            # N XCENTER   YCENTER   MAG      SHARPNESS   SROUND      GROUND      ID         \
            # U pixels    pixels    #        #           #           #           #          \
            # F %-13.3f   %-10.3f   %-9.3f   %-12.3f     %-12.3f     %-12.3f     %-6d       \
            outfile.write("{:<13.3f}{:<10.3f}{:<9.3f}{:<12.3f}{:<12.3f}{:<12.3f}{:<6d}\n".format(coo_data['XCENTER'][i],
                                                                                                 coo_data['YCENTER'][i],
                                                                                                 coo_data['MAG'][i],
                                                                                                 coo_data['SHARPNESS'][
                                                                                                     i],
                                                                                                 coo_data['SROUND'][i],
                                                                                                 coo_data['GROUND'][i],
                                                                                                 coo_data['ID'][i]))
    return


if __name__ == '__main__':
    main()
