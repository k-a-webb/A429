__author__ = 'kawebb'

import os
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits


def main():
    parser = argparse.ArgumentParser(
        description='For a given image, and list of star positions, return list of stars which are not saturated')
    parser.add_argument("--image", '-img',
                        action="store",
                        default='/Users/kawebb/Dropbox/A429/L1544/L1544_J_sub.fits',
                        help="Fits image of star field")
    parser.add_argument("--sexfile", '-sf',
                        action='store',
                        default='/Users/kawebb/Dropbox/A429/L1544/L1544_sex.txt',
                        help='List of star positions as output by sextractor.')
    parser.add_argument("--outfile", '-of',
                        action='store',
                        default=None,
                        help='Out file, list of unsaturated stars')
    parser.add_argument("--tmass", '-tm',
                        action='store',
                        default='2mass_L1544.tbl',
                        help='File of 2Mass stars in region of interest')
    parser.add_argument("--toplot", '-plt',
                        action='store',
                        default=False,
                        help='Whether or not to plot saturated stars over top of the image')
    args = parser.parse_args()


    remove_satur(args.image, args.sexfile, args.outfile, args.toplot)
    find_2mass(args.image, args.tmass, args.outfile)


def find_2mass(ffile, tmfile, usfile):

    tmlist = read_2mass(tmfile)
    uslist = read_usfile(usfile)

    idx = np.zeros(len(uslist))
    for i in uslist.index:
        sig = 0.005
        stars = []
        while len(stars) == 0:
            ras = tmlist.query('({} < ra) & (ra < {})'.format(uslist.x_wcs[i]-sig, uslist.x_wcs[i]+sig))
            stars = ras.query('({} < dec) & (dec < {})'.format(uslist.y_wcs[i]-sig, uslist.y_wcs[i]+sig))
            sig += 0.001
        idx[i] = stars.index[0]

    tmstars = tmlist.loc[idx,:]

    tmstars.to_csv('outfile.txt', index=False, sep=' ')
    return tmstars

def remove_satur(ffile, sfile, ofile=None, toplot=False):

    if ofile is None:
        ofile = sfile.strip('.txt') + '_clean.txt'
    if os.path.exists(ofile):
        print '  List of unsaturated stars already exists: {}'.format(ofile)
        return
    else:
        print "  Creating list of unsaturated stars for image {}, writing to {}".format(ffile, ofile)

    stable = read_sex(sfile)
    fdata, fhdr = read_fits(ffile)

    ix = detect_satur(fdata, stable['x'].values, stable['y'].values, stable['b'].values)
    us_stable = stable[ix]

    us_stable.to_csv(ofile, index=False, sep=' ')

    if toplot:
        plt.imshow(fdata, origin='lower', norm=matplotlib.colors.LogNorm())
        s_stable = stable[np.invert(ix)]
        plt.plot(s_stable.x.values, s_stable.y.values, 'x')
        plt.title('Saturated stars')
        plt.show()

        plt.imshow(fdata, origin='lower', norm=matplotlib.colors.LogNorm())
        plt.plot(us_stable.x.values, us_stable.y.values, 'x')
        plt.title('Unsaturated stars')
        plt.show()

def detect_satur(fdata, x, y, r, nz=1.):

    idx = []
    n = np.arange(len(x))

    for i in n:
        zeros = []
        rr = int(r[i]) + 1
        # r**2 = (xx - x0)**2 + (yy - y0)**2
        for xx in range(int(x[i]) - rr, int(x[i]) + rr):
            yy_p = np.sqrt(rr ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            yy_m = -np.sqrt(rr ** 2 - (xx - int(x[i])) ** 2) + int(y[i])
            for yy in range(int(yy_m), int(yy_p)):
                if fdata[yy, xx] == 0.:
                    zeros.append(1.)
        if len(zeros) >= nz:
            idx.append(i)

    ix = np.invert(np.in1d(n.ravel(), idx).reshape(n.shape))
    return ix

def read_sex(sfile):
    #   1 NUMBER                 Running object number
    #   2 FLUXERR_ISO            RMS error for isophotal flux                               [count]
    #   3 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
    #   4 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
    #   5 X_IMAGE                Object position along x                                    [pixel]
    #   6 Y_IMAGE                Object position along y                                    [pixel]
    #   7 A_IMAGE                Profile RMS along major axis                               [pixel]
    #   8 B_IMAGE                Profile RMS along minor axis                               [pixel]
    #   9 FLAGS                  Extraction flags
    # table = pd.read_csv(sfile, skiprows=9, sep=r"\s*", engine='python',
    #                     names=['num', 'fluxerr_iso', 'flux_auto', 'fluxerr_auto', 'x', 'y', 'a', 'b', 'flag'])

    #   1 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
    #   2 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
    #   3 THRESHOLD              Detection threshold above background                       [count]
    #   4 X_IMAGE                Object position along x                                    [pixel]
    #   5 Y_IMAGE                Object position along y                                    [pixel]
    #   6 X_WORLD                Barycenter position along world x axis                     [deg]
    #   7 Y_WORLD                Barycenter position along world y axis                     [deg]
    #   8 B_IMAGE                Profile RMS along minor axis                               [pixel]
    #   9 FLAGS                  Extraction flags
    table = pd.read_csv(sfile, skiprows=9, sep=r"\s*", engine='python',
                        names=['flux_auto', 'fluxerr_auto', 'thresh', 'x', 'y', 'x_wcs', 'y_wcs', 'b', 'flag'])

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

def read_2mass(tmfile):
    # ra, dec, j_m, j_cmsig, j_msigcom, j_snr, h_m, h_cmsig, h_msigcom, h_snr, k_m, k_cmsig, k_msigcom, k_snr
    #   rd_flg, dist, angle, j_h, h_k, j_k
    table = pd.read_csv(tmfile, skiprows=55, sep=r"\s*", engine='python',
                        names=['ra', 'dec', 'j_m', 'j_cmsig', 'j_msigcom', 'j_snr', 'h_m', 'h_cmsig', 'h_msigcom',
                               'h_snr', 'k_m', 'k_cmsig', 'k_msigcom', 'k_snr', 'rd_flg', 'dist', 'angle', 'j_h',
                               'h_k', 'j_k'])
    return table

def read_usfile(usfile):
    # flux_auto fluxerr_auto thresh x y x_wcs y_wcs b flag
    table = pd.read_csv(usfile, skiprows=1, sep=r"\s*", engine='python',
                        names=['flux_auto', 'fluxerr_auto', 'thresh', 'x', 'y', 'x_wcs', 'y_wcs', 'b', 'flag'])
    return table

if __name__ == '__main__':
    main()
