#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

STEPS=900

EXTNAMES = ['CHIP%d.INT1'%i for i in [1,2,3]]
colors=['y','r','b']
lstyle=['-', '--','--','-']
lwidth=[1.0,3.0,1.0,3.0]
ZOOM = 1
YMAX = 5000
YMIN = -3000
CAT_FACTOR = 50
SPEC_FACTOR=1

X = np.arange(2048)+1

def onkey(event):
    if event.key=='w':
        print('%.2f '%event.xdata,end='',flush=True)
    if event.key=='e':
        print('%.2f'%event.xdata,flush=True)
    if event.key=='d':
        print('undo last',flush=True)


def ev(p,x):
    return np.polyval(p[::-1],x)

def main(specname,catname=None,cat2name=None,tracename=None):

    FIG = plt.figure(figsize=(20, 7.5))
    ax = FIG.add_subplot(1, 1, 1)
    ax.set_yticks([])
    FIG.tight_layout(pad=0.02)

    spec_exts = fits.open(specname)
    sett = spec_exts[0].header['ESO INS WLEN ID']
    shortsett = ''.join(sett.split('/')[:2])
    if catname:
        cat_data = fits.open(catname)[1].data
        cat_wav, cat_ints = cat_data['Wavelength'], cat_data['Emission']
        ax.vlines(cat_wav, np.zeros_like(cat_ints), -1*CAT_FACTOR*cat_ints, 'k',alpha=0.2)
    if cat2name:
        cat_data = fits.open(cat2name)[1].data
        cat_wav, cat_ints = cat_data['Wavelength'], cat_data['Emission']
        ax.vlines(cat_wav, np.zeros_like(cat_ints), -1*CAT_FACTOR*cat_ints, 'g')

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
                print('fail %s order %d, %s'%(ext,order,e))
                continue
            if np.sum(np.isnan(spec)) == 2048:
                continue
            #print('found %s order %d'%(ext,order))
            if twd is not None:
                p = twd[(twd['Order']==order) & (twd['TraceNb']==1)]\
                    ['Wavelength'][0]
                if not np.isnan(p).any():
                    wl = ev(p,X)
            xcor = h.get('ESO QC WAVE BESTXCORR-%02d-01'%order)
            #ax.hlines(np.percentile(spec,10),wl[0],wl[-1],'r')
            ax.plot(wl,spec*SPEC_FACTOR,label=str(order),color=colors[i],
                linestyle='-')
            ax.text(wl.mean(),3000,'(O:%d D:%d X:%.2f)'%(order,i+1,xcor or 0.0), fontsize=9,
                horizontalalignment='center')



    #x1,x2,y1,y2=ax.axis()
    if sett.startswith('Y'):
        x1,x2=950,1100
    elif sett.startswith('J'):
        x1,x2=1150,1350
    elif sett.startswith('H'):
        x1,x2=1430,1780
    elif sett.startswith('K'):
        x1,x2=1940,2400
    elif sett.startswith('L'):
        x1,x2=3500,4200
    elif sett.startswith('M'):
        x1,x2=4500,5000
    rang = x2-x1

    ax.axis((x1,x1+(rang/ZOOM),-YMAX,YMAX))

    cid = FIG.canvas.mpl_connect('key_press_event', onkey)

    plt.show()


if __name__ == '__main__':
    main(*sys.argv[1:])
