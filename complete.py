__author__ = 'kawebb'

import numpy as np
from scipy.optimize import curve_fit
from glob import glob
import plot_cmd
import matplotlib.pyplot as plt
from collections import Counter


def main():
    files_j = glob('phot_t3/*_J_sex*')
    files_ks = glob('phot_t3/*_Ks_sex*')
    files_h = glob('phot_t3/*_H_sex*')

    for file in files_j[:1]:
        print file
        phot = plot_cmd.read_sex_band(file)

        mags = phot['mag_aper'].values
        magerrs = phot['magerr_aper'].values

        i_1p = np.argmin(np.abs(magerrs - 0.01))
        mag_1p = mags[i_1p]

        plt.scatter(mags, magerrs, marker='.')
        plt.plot(mags, mags * 0. + 0.01, '-k')
        plt.scatter(mag_1p, magerrs[i_1p], marker='o', color='r', label='mag={}'.format(mag_1p))
        plt.xlabel('mag')
        plt.ylabel('mag err')
        plt.ylim(0., .3)
        plt.legend()
        plt.show()

        bins, n = np.unique(mags, return_counts=True)
        nbins = len(bins) / 100

        pow = lambda x, A, alpha: A * x ** alpha

        n, bins, patches = plt.hist(mags, bins=nbins)
        idx = np.argmax(n)

        popt, pcov = curve_fit(pow, bins[:idx], n[:idx], p0=[1., 1.])
        plt.plot(bins[:idx], pow(bins[:idx], popt[0], popt[1]), '-r', lw=3)
        plt.xlabel('mag')
        plt.ylabel('occurance')
        plt.show()


if __name__ == '__main__':
    main()
