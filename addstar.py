import argparse
from pyraf import iraf
import os
import numpy as np
from astropy.io import ascii
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
    parser.add_argument("--addphot", '-addphot',
                        action='store',
                        default=None,
                        help='Input file of source extractor file with artificial stars added.')

    # optional
    parser.add_argument("--outcoords", '-coo',
                        action='store',
                        default=None,
                        help='Outfile of coordinates [wcs].')
    parser.add_argument("--outphot", '-phot',
                        action='store',
                        default=None,
                        help='Outfile of daophot photometry.')
    parser.add_argument("--outsexphot", '-sexphot',
                        action='store',
                        default=None,
                        help='Outfile of daophot photometry with source extractor photometry.')
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
    parser.add_argument("--outartphot", '-artphot',
                        action='store',
                        default=None,
                        help='Outfile of artificial star positions and magnitudes.')
    parser.add_argument("--outaddimage", '-addimg',
                        action='store',
                        default=None,
                        help='Outfile image with artificial stars.')
    args = parser.parse_args()

    if args.sexfile is None:
        raise Exception, 'No input source extractor file given'
    if args.image is None:
        raise Exception, 'No input image file given'

    artphotfile = run_daophot(args.sexfile, args.image, args.outpath, args.aperture, args.zeropoint, args.outcoords,
                              args.outphot, args.outsexphot, args.outpst, args.outopst, args.outgroup, args.outpsfimg,
                              args.outartphot, args.outaddimage)

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
    plt.xlim(10,30)
    plt.show()


def run_daophot(sexfile, image, outpath, aperture=None, zeropoint=None, coords=None, photfile=None, sexphotfile=None,
                pstfile=None, opstfile=None, groupfile=None, psfimage=None, maxnpsf=25., artphotfile=None,
                addimage=None):
    try:
        prefix_sex = (sexfile.split('/')[-1]).split('.')[0]
    except IndexError:
        prefix_sex = sexfile.split('.')[0]
    try:
        prefix_img = (image.split('/')[-1]).split('.')[0]
    except IndexError:
        prefix_img = sexfile.split('.')[0]

    if coords is None:
        coords = os.path.join(outpath, prefix_sex + '.coo')
    if photfile is None:
        photfile = os.path.join(outpath, prefix_sex + '.phot')
    if sexphotfile is None:
        sexphotfile = os.path.join(outpath, prefix_sex + '.mag')
    if pstfile is None:
        pstfile = os.path.join(outpath, prefix_sex + '.pst')
    if opstfile is None:
        opstfile = os.path.join(outpath, prefix_sex + '.opst')
    if groupfile is None:
        groupfile = os.path.join(outpath, prefix_sex + '.psg')
    if psfimage is None:
        psfimage = os.path.join(outpath, prefix_img + '.psf.fits')
    if artphotfile is None:
        artphotfile = os.path.join(outpath, prefix_sex + '.art')
    if addimage is None:
        addimage = os.path.join(outpath, prefix_img + '.add.fits')

    sex_data = ascii.read(sexfile, format='sextractor')
    np.savetxt(coords, np.transpose((np.array(sex_data['X_WORLD']), np.array(sex_data['Y_WORLD']))), delimiter=' ')

    bkg = np.std(np.array(sex_data['BACKGROUND']))
    print bkg

    iraf.daophot()
    # phot -- do aperture photometry on a list of stars
    # phot image
    if not os.path.exists(photfile):
        if aperture is None:
            raise Exception, 'No input aperture parameter given'
        if zeropoint is None:
            raise Exception, 'No input zeropoint value given'
        iraf.daophot.phot(image, coords=coords, output=photfile, wcsin="world", apertures=float(aperture),
                          zmag=float(zeropoint), verbose=True)  # , annulus=aperture + 10, dannulus=10,
    else:
        print 'File {} already exists'.format(photfile)

    if not os.path.exists(sexphotfile):
        write_phot(photfile, sexphotfile, sexfile)

    # pstselect -- select candidate psf stars from a photometry file
    # pstselect image photfile pstfile maxnpsf
    if not os.path.exists(pstfile):
        iraf.daophot.pstselect(image, sexphotfile, pstfile, maxnpsf)
    else:
        print 'File {} already exists'.format(pstfile)

    # psf -- build the point spread function for an image
    # psf image photfile pstfile psfimage opstfile groupfile
    if not os.path.exists(psfimage):
        iraf.daophot.psf(image, photfile=sexphotfile, pstfile=pstfile, psfimage=psfimage, opstfile=opstfile,
                         groupfile=groupfile, interactive=False)
    else:
        print 'File {} already exists'.format(psfimage)

    # define parameters of artificial stars
    if not os.path.exists(artphotfile):
        nstars = 3000
        x = random_sample(nstars, sex_data['X_WORLD'])
        y = random_sample(nstars, sex_data['Y_WORLD'])
        mags = random_sample(nstars, np.linspace(18., 28., nstars))
        np.savetxt(artphotfile, np.transpose([x, y, mags]), delimiter=' ')
    else:
        print 'File {} already exists'.format(artphotfile)

    # addstar -- add artificial stars to images
    # addstar image photfile psfimage addimage
    if not os.path.exists(addimage):
        iraf.daophot.addstar(image, artphotfile, psfimage, addimage)  # 50% recovery
    else:
        print 'File {} already exists'.format(addimage)

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
    idxs = []
    for i in range(len(sex)):
        d, i = tree.query((sex['X_IMAGE'][i], sex['Y_IMAGE'][i]), distance_upper_bound=max_sep)
        if (d != np.inf):
            idxs.append(i)

    sex_matched = sex[idxs]

    with open(newphotfile, 'w') as outfile:
        with open(photfile, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)

        for i in range(len(phot)):
            if (phot['MAG'][i] > 0.) and (sex_matched['MAG_APER'][i] < 99.):
                outfile.write("{:<23}{:<10.3f}{:<10.3f}{:<6d}{:<23}{:<6d} \\\n".format(
                    phot['IMAGE'][i], phot['XINIT'][i], phot['YINIT'][i], phot['ID'][i], phot['COORDS'][i],
                    phot['LID'][i]))
                outfile.write("   {:<11.3f}{:<11.3f}0.000   0.000   INDEF   INDEF          0    NoError   \\\n".format(
                    phot['XCENTER'][i], phot['YCENTER'][i]))
                outfile.write("   {:<15.7g}{:<15.7g}{:<15.7g}{:<7d}{:<9d}0    NoError   \\\n".format(
                    phot['MSKY'][i], phot['STDEV'][i], phot['SSKEW'][i], phot['NSKY'][i], phot['NSREJ'][i]))
                outfile.write("   1.             INDEF          INDEF                  INDEF                  \\\n")
                outfile.write("   8.00     {:<14.7g}{:<11.7g}{:<14.7g}{:<7.3f}{:<6.3f}0    NoError    \n".format(
                    phot['SUM'][i], phot['AREA'][i], sex_matched['FLUX_APER'][i], sex_matched['MAG_APER'][i],
                    sex_matched['MAGERR_APER'][i]))

    return


if __name__ == '__main__':
    main()
