#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

STEPS=900

EXTNAMES = ['CHIP%d.INT1'%i for i in [1,2,3]]
ZOOM = 1
YMAX = 5000
YMIN = -3000
CAT_FACTOR = 5
CAT_OFFSET = 200
SPEC_FACTOR=1
FIT_NPIX=8 # x2 , to left and right of click

X = np.arange(2048)+1
catalog=None
PIX_is=[]
WL_should=[]
SPEC = []
LABEL = ''

def onkey(event):
    global WL_should, PIX_is, SPEC, LABEL
    if event.key=='w':
        print('%.2f '%event.xdata,end='',flush=True)
    if event.key=='e':
        print('%.2f'%event.xdata,flush=True)
    if event.key=='d':
        print('undo last',flush=True)
    if event.key=='x':
        cat = find_nearest(catalog,event.xdata)
        print('%.3f'%cat)
        WL_should.append(cat)
    if event.key=='F':
        try:
            p=np.polyfit(PIX_is,WL_should,deg=2)
            plt.plot(np.polyval(p,X),SPEC,'y')
            ext,order=LABEL.split()
            order=int(order)
            updatefits(p,ext,order)
        except Exception as E:
            print('fit failed: ', E)
        PIX_is=[]
        WL_should=[]

def onpick(event):
    line = event.artist
    xdata, ydata = line.get_data()
    peak = event.ind.mean()
    #print('on pick line:', np.array([xdata[ind], ydata[ind]]).T)
    #x=xdata[ind-FIT_NPIX-1:ind+FIT_NPIX]
    #y=ydata[ind-FIT_NPIX-1:ind+FIT_NPIX]
    #y-=y.min()
    #peak = np.sum(x*y)/np.sum(y)
    print(peak)
    global PIX_is, SPEC, LABEL
    LABEL = line.get_label()
    SPEC = ydata
    PIX_is.append(peak)

def updatefits(p,ext=None,order=None,trace=1):
    tw=fits.open(sys.argv[-1])
    twd=tw[ext].data
    for i,row in enumerate(twd):
        print(i,row['Order'])
        if row['Order'] == order:
            break
    print('Writing: ', ext,i, order, tw[ext].data[i]['Order'])
    tw[ext].data[i]['Wavelength']=p[::-1]
    tw.writeto(sys.argv[-1],overwrite=True)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

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
        global catalog
        catalog = cat_wav
        ax.vlines(cat_wav, np.zeros_like(cat_ints)-CAT_OFFSET, -CAT_OFFSET-CAT_FACTOR*cat_ints, 'tab:gray',alpha=0.5)
    if cat2name:
        cat_data = fits.open(cat2name)[1].data
        cat_wav, cat_ints = cat_data['Wavelength'], cat_data['Emission']
        ax.vlines(cat_wav, np.zeros_like(cat_ints)-CAT_OFFSET, -CAT_OFFSET-CAT_FACTOR*cat_ints, 'tab:green')

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
                tw.close()
            xcor = h.get('ESO QC WAVE BESTXCORR-%02d-01'%order)
            #ax.hlines(np.percentile(spec,10),wl[0],wl[-1],'r')
            spec *= SPEC_FACTOR
            spec = spec - np.median(spec) + 100
            ax.plot(wl,spec,label=' '.join((ext,str(order))),color='tab:blue',alpha=0.8,pickradius=5,picker=True)
            ax.text(wl.mean(),1000,'(O:%d D:%d X:%.2f)'%(order,i+1,xcor or 0.0), fontsize=11,
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

    cid1 = FIG.canvas.mpl_connect('key_press_event', onkey)
    cid2 = FIG.canvas.mpl_connect('pick_event', onpick)

    plt.show()


if __name__ == '__main__':
    main(*sys.argv[1:])
