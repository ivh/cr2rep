#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


EXTNAMES = ['CHIP%d.INT1'%i for i in [1,2,3]]
ZOOM = 1
YMAX = 5000
YMIN = -3000
CAT_FACTOR = 5
CAT_OFFSET = 200
SPEC_FACTOR=1

X = np.arange(2048)+1

def ev(p,x):
    return np.polyval(p[::-1],x)

def main(specname,catname=None,cat2name=None,tracename=None):

    FIG = plt.figure(figsize=(20, 7.5))
    ax = FIG.add_subplot(1, 1, 1)
    ax.set_yticks([])
    FIG.tight_layout(pad=0.02)

    spec_exts = fits.open(specname)
    sett = spec_exts[0].header['ESO INS WLEN ID']
    sett = sett.replace('/','_')
    if catname:
        cat_data = fits.open(catname)[1].data
        cat_wav, cat_ints = cat_data['Wavelength'], cat_data['Emission']
        ax.vlines(cat_wav, np.zeros_like(cat_ints)-CAT_OFFSET, -CAT_OFFSET-1*CAT_FACTOR*cat_ints, 'k',alpha=0.2)
    if cat2name:
        cat_data = fits.open(cat2name)[1].data
        cat_wav, cat_ints = cat_data['Wavelength'], cat_data['Emission']
        ax.vlines(cat_wav, np.zeros_like(cat_ints)-CAT_OFFSET, -CAT_OFFSET-1*CAT_FACTOR*cat_ints, 'g')

    for i,ext in enumerate(EXTNAMES):
        twd =None
        if tracename:
            try:
                tw=fits.open(tracename)
                twd = tw[ext].data

            except Exception as e:
                print('TRACEWAVE has no extension: %s'%ext)

        h=spec_exts[ext].header
        for order in np.arange(9)+1:
            try:
                wl = spec_exts[ext].data['%02d_01_WL'%order]
                spec = spec_exts[ext].data['%02d_01_SPEC'%order]
            except Exception as e:
                #print('fail %s order %d, %s'%(ext,order,e))
                continue
            if np.sum(np.isnan(spec)) == 2048:
                continue
            if twd is not None:
                p = twd[(twd['Order']==order) & (twd['TraceNb']==1)]\
                    ['Wavelength'][0]
                if not np.isnan(p).any():
                    wl = ev(p,X)
            #spec = np.ma.masked_where(np.isnan(spec),spec)
            spec *= SPEC_FACTOR
            spec = spec - np.median(spec) + 100
            ax.plot(wl,spec,label=str(order),color='tab:blue',alpha=0.8, linestyle='-')
            xcor = h.get('ESO QC WAVE BESTXCORR-%02d-01'%order)
            ax.text(wl.mean(),1000,'(Order %d, Detector %d, X-corr %.2f %%)'%(order,i+1,xcor or 0.0), fontsize=11,
                horizontalalignment='center')


            ax.axis((wl.min(),wl.max(),YMIN,YMAX))
            outname = specname.replace('.fits','_')
            plt.savefig(outname+'wavecal_o%d_d%d.png'%(order,i+1),dpi=120)



if __name__ == '__main__':
    main(*sys.argv[1:])
