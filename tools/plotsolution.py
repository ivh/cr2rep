#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
from glob import glob
import subprocess
from matplotlib import pyplot as plt


X = np.arange(2048)+1
fig1 = plt.figure(figsize=(14,4))
axs1=[fig1.add_subplot(131)]
axs1 += [fig1.add_subplot(132, sharey=axs1[0]),
         fig1.add_subplot(133, sharey=axs1[0])]

for twn in sys.argv[1:]:
    #print(twn)
    with fits.open(twn) as twf:
        maindisp = twf[0].header['ESO INS GRAT1 DISP']
        for det in [1,2,3]:
            tw = twf['CHIP%d.INT1'%det].data
            if tw is None: continue
            h = twf['CHIP%d.INT1'%det].header
            ax=axs1[det-1]
            for order in np.arange(9)+1:
                try:
                    guesstart = h.get('ESO INS WLEN BEGIN%d'%order)
                    guessend = h.get('ESO INS WLEN END%d'%order)
                    guessdisp = (guessend-guesstart)/2048.0
                    ax.plot(1024,guessdisp,'ro')
                except:
                    continue
                
                xc = h.get('ESO QC WAVE BESTXCORR-%02d-01'%order)
                mask = tw['Order']==order
                try:
                    p = tw[mask]['Wavelength'][0][::-1]
                except:
                    continue
                poly = np.poly1d(p)
                deriv = poly.deriv()
                ax.plot(deriv(X))
                ax.text(X[100],deriv(X[100]), 'O:%d XC:%.2f'%(order,xc or -1))
            ax.plot(1024,maindisp,'kD')
            ax.set_xlim((1,2048))
            if det >1:
                for tk in ax.get_yticklabels():
                    tk.set_visible(False)
 

    fig1.subplots_adjust(wspace=0)
    outf = twn.replace('.fits','_disp.png')              
    fig1.savefig(outf)
    [ax.clear() for ax in axs1]

