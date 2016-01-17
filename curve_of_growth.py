__author__ = 'kawebb'

import matplotlib.pyplot as plt
import numpy as np
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='For a given object, searched for text file with curve of growth and plots')
    parser.add_argument("--object", '-o',
                        action="store",
                        default='L1544',
                        help="Object. Looks for file path 'object/curve_of_growths/curve_of_growth.txt'")
    args = parser.parse_args()

    INFILE = '{}/curve_of_growths/curve_of_growth.txt'.format(str(args.object))
    OUTFIG = '{}/curve_of_growths/curve_of_growth.png'.format(str(args.object))
    TITLE = (str(args.object))

    r_J = []
    mag_J = []
    r_Ks = []
    mag_Ks = []
    r_H = []
    mag_H = []

    with open(INFILE, 'r') as infile:
        copy = False
        for line in infile:
            if line.strip() == "Ks":
                copy = True
            elif line.strip() == "STOP_Ks":
                copy = False
            elif copy:
                values = line.split()
                r_Ks.append(float(values[0]))
                mag_Ks.append(float(values[4]))

    with open(INFILE, 'r') as infile:
        copy = False
        for line in infile:
            if line.strip() == "J":
                copy = True
            elif line.strip() == "STOP_J":
                copy = False
            elif copy:
                values = line.split()
                r_J.append(float(values[0]))
                mag_J.append(float(values[4]))

    with open(INFILE, 'r') as infile:
        copy = False
        for line in infile:
            if line.strip() == "H":
                copy = True
            elif line.strip() == "STOP_H":
                copy = False
            elif copy:
                values = line.split()
                r_H.append(float(values[0]))
                mag_H.append(float(values[4]))

    plt.plot(r_J, mag_J, label='J')
    plt.plot(r_Ks, mag_Ks, label='Ks')
    plt.plot(r_H, mag_H, label='H')
    plt.xlabel('aperture radius')
    plt.ylabel('average weighted magnitude')
    plt.title(TITLE)
    plt.legend()
    plt.savefig(OUTFIG)
    plt.show()

    plt.plot(r_J[1:], np.diff(mag_J), label='J')
    plt.plot(r_Ks[1:], np.diff(mag_Ks), label='Ks')
    plt.plot(r_H[1:], np.diff(mag_H), label='H')
    plt.xlabel('aperture radius')
    plt.ylabel('difference in average weighted magnitude')
    plt.title(TITLE + '_diff')
    plt.legend(loc=4)
    plt.ylim(-0.1,0.1)
    plt.savefig(OUTFIG.strip('.png') + '_diff.png')
    plt.show()


if __name__ == '__main__':
    main()
