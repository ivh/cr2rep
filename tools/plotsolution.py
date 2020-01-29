#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
from glob import glob
import subprocess
from matplotlib import pyplot as plt

#C2RAWDIR = "/mnt.staging/DETDATA/PAE/201906_CR_cooldown13"
#C2RAWDIR = "/home/tom/pipes/tdata/CD13/"
#C2REDUCEDDIR = os.path.join(C2RAWDIR, "REDUCED")
C2RAWDIR = os.environ['C2RAWDIR']
C2REDUCEDDIR = os.environ['C2REDUCEDDIR']
os.chdir(C2REDUCEDDIR)

input_tracewave_name = 'cr2res_cal_wave_%s_tw.fits'

X = np.arange(2048)+1
fig1 = plt.figure(figsize=(14,4))
axs1=[fig1.add_subplot(131)]
axs1 += [fig1.add_subplot(132, sharey=axs1[0]),
         fig1.add_subplot(133, sharey=axs1[0])]

if len(sys.argv) > 1:
    dirlist = sys.argv[1:]
else:
    dirlist = glob('?_?_?')
dirlist.sort()

for dir in dirlist:
    os.chdir(dir)
    sett = dir[-5:]
    print('# ', sett)
    

    with fits.open(input_tracewave_name%sett) as twf:
        maindisp = twf[0].header['ESO INS GRAT1 DISP']
        for det in [1,2,3]:
            tw = twf['CHIP%d.INT1'%det].data
            h = twf['CHIP%d.INT1'%det].header
            #tw.sort(order='Order')
            for order in np.arange(9)+1:
                try:
                    guesstart = h.get('ESO INS WLEN BEGIN%d'%order)
                    guessend = h.get('ESO INS WLEN END%d'%order)
                    guessdisp = (guessend-guesstart)/2048.0
                    axs1[det-1].plot(1024,guessdisp,'ro')
                except: pass
                xc = h.get('ESO QC WAVE BESTXCORR-%02d-01'%order)
                print("o:%d, xc:%.1f"%(order,xc or 0.0))
                if xc ==None: continue
                if xc < 0.1: continue

                mask = tw['Order']==order
                try:
                    p = tw[mask]['Wavelength'][0][::-1]
                except: continue
                #print(det, order, p)
                poly = np.poly1d(p)
                deriv = poly.deriv()
                axs1[det-1].plot(deriv(X))
                axs1[det-1].text(X[100],deriv(X[100]), 'O:%d XC:%.2f'%(order,xc))
            axs1[det-1].plot(1024,maindisp,'kD')
            
    fig1.savefig('../solut_deriv_%s.png'%sett)
    [ax.clear() for ax in axs1]
    os.chdir(os.pardir)

