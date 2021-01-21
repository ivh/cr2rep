#!/usr/bin/env python3
import os
import sys
import numpy as np
np.seterr(all='raise')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit


EXTNAMES = ['CHIP%d.INT1'%i for i in [1,2,3]]
FIT_RAD=6 # x2 , to left and right of click
X = np.arange(2048)+1

def gauss(x, *p):
    A, mu, sigma, cont = p
    return cont + A* np.exp(-(x-mu)**2/(2.*sigma**2))

def onclick(event):
    if not event.button==2:
        return
    global fitter, ax, curspc, outf, alldiffs
    #print(event.x, event.y, event.xdata, event.ydata)
    peakx, peaky = event.xdata, event.ydata
    start,end=int(peakx)-FIT_RAD, int(peakx)+FIT_RAD
    centers = []
    for s in curspc:
        x,y = np.arange(end-start)+start, s[start:end]
        try:
            coeff, var_matrix = curve_fit(gauss, x, y, p0=(peaky,peakx,2,0))
        except:
            print('Bad fit! Try again...')
            return
        ax.plot(x,gauss(x,*coeff),'k-')
        centers.append(coeff[1])
    diffs = [np.abs(centers[i]-centers[i+1]) for i in np.arange(len(centers)-1)]
    print(diffs)
    np.savetxt(outf, [diffs])
    outf.flush()


def onkey(event):
    if event.key=='n':
        print('plotting next')
        plotnext()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def ev(p,x):
    return np.polyval(p[::-1],x)

def plotnext():
    global specs, ax, curspc
    ax.clear()
    try:
        curspc = [s - np.percentile(s,5) for s in specs.pop(0)]
    except IndexError:
        print("No more spectra! Exiting.")
        exit()
    for s in curspc:
        ax.plot(s,pickradius=5,picker=True)
    ax.figure.canvas.draw()


#globals
specs =[] # list of spec-lists
curspc = []
outf = open('diffs.dat','w')

def main(*fnames):
    global specs, ax
    FIG = plt.figure(figsize=(20, 7.5))
    ax = FIG.add_subplot(1, 1, 1)
    ax.set_yticks([])
    FIG.tight_layout(pad=0.02)

    files=[fits.open(fn) for fn in fnames]

    for i,ext in enumerate(EXTNAMES):
        for order in np.arange(9)+1:
            try:
                specs.append([fil[ext].data['%02d_01_SPEC'%order] for fil in files])
            except Exception as e:
                print('fail %s order %d, %s'%(ext,order,e))
                continue

    cid1 = FIG.canvas.mpl_connect('key_press_event', onkey)
    cid2 = FIG.canvas.mpl_connect('button_press_event', onclick)
    plotnext()
    plt.show()


if __name__ == '__main__':
    main(*sys.argv[1:])
