#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def compare_extract(flat, trace, extract, title=''):
    """ compare flat and trace """
    flat = fits.open(flat)
    trace = fits.open(trace)
    extract = fits.open(extract)
    X = np.arange(2048)
    FIG = plt.figure(figsize=(10, 3.5))
    FIG.suptitle(title)

    for i in [1, 2, 3]:
        ax = FIG.add_subplot(2, 3, i)
        #ax.set_xticks([])
        #ax.set_ylabel('trace')
        #ax.set_yticks([-1, 0, 1, 2, 3,4 ,5, 6, 7, 8, 9])
        ax.set_title('CHIP%s' % i)

        edata = extract['CHIP%s' % i].data
        tdata = trace['CHIP%s' % i].data
        fdata = flat[i].data


        if edata is None or tdata is None or fdata is None:
            print('No data for CHIP%s, skipping.' % i)
            continue

        ax.imshow(fdata)

        for alla, upper, lower, order, tracenb, wave, curvA, curvB, curvC in tdata:
            pol = np.polyval(upper[::-1], X)
            ax.plot(X, pol, ':w')

            pol = np.polyval(lower[::-1], X)
            ax.plot(X, pol, ':w')

            pol = np.polyval(alla[::-1], X)
            ax.plot(X, pol, '--w')
            if np.isnan(pol[1024]):
                continue
            ax.text(1024, pol[1024], '%s-%s' % (order, tracenb), color='w', horizontalalignment='center',
                    verticalalignment='center', size=13)

        ax = FIG.add_subplot(2, 3, i+3)
        ax.set_xlabel('wavelength [nm]')
        ax.set_ylabel('intensity [a.u.]')

        data = np.empty((9, 2048))
        for j in range(2048):
            data[:len(edata[j][:9]), j] = edata[j][:9]

        for i, (alla, upper, lower, order, tracenb, wave) in enumerate(tdata[:9]):
            x = np.polyval(wave[::-1], X)
            y = data[i]
            ax.plot(x, y)
        #ax.imshow(data, aspect=50, interpolation='nearest')

    plt.show()


def compare_extract2(flat, trace, extract, title=''):
    """ compare flat and trace """
    flat = fits.open(flat)
    trace = fits.open(trace)
    extract = fits.open(extract)
    X = np.arange(2048)
    FIG = plt.figure(figsize=(10, 3.5))
    FIG.suptitle(title)

    for i in [1, 2, 3]:
        ax = FIG.add_subplot(2, 3, i)
        #ax.set_xticks([])
        #ax.set_ylabel('trace')
        #ax.set_yticks([-1, 0, 1, 2, 3,4 ,5, 6, 7, 8, 9])
        ax.set_title('CHIP%s' % i)

        edata = extract['CHIP%s' % i].data
        tdata = trace['CHIP%s' % i].data
        fdata = flat[i].data


        if edata is None or tdata is None or fdata is None:
            print('No data for CHIP%s, skipping.' % i)
            continue

        ax.imshow(fdata)

        for alla, upper, lower, order, tracenb, wave in tdata:
            pol = np.polyval(upper[::-1], X)
            ax.plot(X, pol, ':w')

            pol = np.polyval(lower[::-1], X)
            ax.plot(X, pol, ':w')

            pol = np.polyval(alla[::-1], X)
            ax.plot(X, pol, '--w')
            if np.isnan(pol[1024]):
                continue
            ax.text(1024, pol[1024], '%s-%s' % (order, tracenb), color='w', horizontalalignment='center',
                    verticalalignment='center', size=13)

        ax = FIG.add_subplot(2, 3, i+3)
        ax.set_ylabel('trace')
        ax.set_yticks([-1, 0, 1, 2, 3,4 ,5, 6, 7, 8, 9])

        data = edata.view('>f8', np.ndarray).reshape((2048, -1)).swapaxes(0, 1)

        ax.imshow(data, aspect=50, interpolation='nearest')

    plt.show()

def compare(flat, trace):
    """ compare flat and trace """
    flat = fits.open(flat)
    trace = fits.open(trace)
    X = np.arange(2048)
    FIG = plt.figure(figsize=(10, 3.5))

    for i in [1, 2, 3]:
        ax = FIG.add_subplot(1, 3, i)
        ax.set_xticks([])
        ax.set_yticks([])

        fdata = flat[i].data
        ax.imshow(fdata)
        axi = plt.axis()

        tdata = trace['CHIP%s' % i].data
        if tdata is None:
            print('No data for CHIP%s, skipping.' % i)
            continue
        for alla, upper, lower, order, tracenb, wave, curvA, curvB, curvC in tdata:
            pol = np.polyval(upper[::-1], X)
            ax.plot(X, pol, ':w')

            pol = np.polyval(lower[::-1], X)
            ax.plot(X, pol, ':w')

            pol = np.polyval(alla[::-1], X)
            ax.plot(X, pol, '--w')
            if np.isnan(pol[1024]):
                continue
            ax.text(1024, pol[1024], '%s-%s' % (order, tracenb), color='w', horizontalalignment='center',
                    verticalalignment='center', size=13)
            #print('order: %s, y: %s' % (tracenb, pol[1024]))

        plt.axis(axi)


    FIG.tight_layout(pad=0.02)
    plt.show()


if __name__ == '__main__':
    fname_flat = sys.argv[1]
    fname_trace = sys.argv[2]
    print('Comparing trace %s to image %s' % (fname_trace, fname_flat))
    compare(fname_flat, fname_trace)


