__author__ = 'kawebb'

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='For photometry output from sextractor for 3 bands, compile magnitude arrays into singe file')
    parser.add_argument("--tmfile", '-tm',
                        action='store',
                        default='CB68/2mass_CB68.tbl',
                        help='Object.')
    parser.add_argument("--outfig", '-o',
                        action='store',
                        default=None,
                        help='Object.')

    args = parser.parse_args()
    tm_table = read_2mass(args.tmfile)

    plot_cmd(tm_table.j_m.values, tm_table.h_m.values, tm_table.k_m.values, args.outfig)


def plot_cmd(mag_j, mag_h, mag_ks, outfig=None):
    """
    Plot colour magnitude diagrams for (J-H) v. Ks, (H-Ks) v. Ks
    """

    fig, ax = plt.subplots(nrows=1, ncols=2)

    ax[0].scatter(mag_j - mag_h, mag_ks, marker='.')
    ax[0].set_ylabel('Ks')
    ax[0].set_xlabel('J-H')

    ax[1].scatter(mag_h - mag_ks, mag_ks, marker='.')
    ax[1].set_ylabel('Ks')
    ax[1].set_xlabel('H-Ks')

    if outfig is not None:
        plt.savefig(outfig)

    plt.show()

def read_2mass(tmfile, skiplines=55):
    # ra, dec, j_m, j_cmsig, j_msigcom, j_snr, h_m, h_cmsig, h_msigcom, h_snr, k_m, k_cmsig, k_msigcom, k_snr
    #   rd_flg, dist, angle, j_h, h_k, j_k
    table = pd.read_csv(tmfile, skiprows=skiplines, sep=r"\s*", engine='python',
                        names=['ra', 'dec', 'j_m', 'j_cmsig', 'j_msigcom', 'j_snr', 'h_m', 'h_cmsig', 'h_msigcom',
                               'h_snr', 'k_m', 'k_cmsig', 'k_msigcom', 'k_snr', 'rd_flg', 'dist', 'angle', 'j_h',
                               'h_k', 'j_k'])
    return table



if __name__ == '__main__':
    main()

