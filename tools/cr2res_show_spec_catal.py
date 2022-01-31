#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

STEPS=900

EXTNAMES = ['CHIP%d.INT1'%i for i in [1,2,3]]
ZOOM = 1
YMAX = 15000
YMIN = -3000
CAT_FACTOR = 100
CAT_OFFSET = 500
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
    if event.key=='X':
        print('fit reset!')
        PIX_is=[]
        WL_should=[]
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
            try:
                coef,cfit = FitCon1(wl,spec,swin=700,deg=2,k2=.5)
                spec -= cfit
            except:
                pass
            ax.plot(wl,spec,label=' '.join((ext,str(order))),color='tab:blue',alpha=0.8,pickradius=5,picker=True)
            #ax.plot(wl,cfit,alpha=0.8,)
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
        x1,x2=1940,2500
    elif sett.startswith('L'):
        x1,x2=3500,4200
    elif sett.startswith('M'):
        x1,x2=4500,5000
    rang = x2-x1

    ax.axis((x1,x1+(rang/ZOOM),-YMAX,YMAX))

    cid1 = FIG.canvas.mpl_connect('key_press_event', onkey)
    cid2 = FIG.canvas.mpl_connect('pick_event', onpick)

    plt.show()

##############  END MAIN ########################

def FitCon1(
    wave,
    flux,
    deg=3,
    niter=10,
    snr=200,
    sig=None,
    swin=7,
    k1=1,
    k2=3,
    mask=None,
    plot=False,
):
    """
    Continuum normalize a spectrum
    based on IDL code by Oleg Kochukhov    Parameters
    ----------
    wave : array
        wavelength grid
    flux : arrau
        spectrum flux
    niter : int
        number of iterations
    deg : int
        polynomial degree of the continuum fit
    swin : int
        smoothing window
    sig : array, optional
        uncertainties on the spectral flux, by default None
    plot : bool, optional
        Whether to plot the results, by default False
    k1 : float, optional
        lower sigma clipping cutoff, by default 1
    k2 : float, optional
        upper sigma clipping cutoff, by default 3
    mask : array, optional
        mask for the points in the spectrum, mask == 1 is always used,
        mask == 2 is never used, mask == 0 is always used, by default None
    snr : float, optional
        signal to noise ratio, by default 200    Returns
    -------
    coeff: array
        polynomial coefficients of the continuum
    con : array
        Continuum points    Raises
    ------
    ValueError
        If poly_type is not a valid value
    """    # flag array
    if mask is None:
        mask = np.full(wave.shape, 0)    # choice of polynomial fitting routine
    if sig is not None:
        sig = 1 / sig
    else:
        sig = np.ones_like(wave)    # initial set of points to fit
    wave = wave - np.mean(wave)
    fmean = np.median(flux)
    flux = flux - fmean
    idx = np.where(mask != 2)
    #idx = np.where( ((flux > con - k1 * rms) & (flux < con + k2 * rms) & (mask != 2)) | (mask == 1))
    smooth = gaussian_filter1d(flux[idx], swin)
    coeff = np.polyfit(wave[idx], smooth, deg, w=sig[idx])
    con = np.polyval(coeff, wave)
    rms = np.sqrt(np.mean((con[idx] - flux[idx]) ** 2))    # iterate niter times
    for _ in range(niter):
        idx = np.where(
            ((flux > con - k1 * rms) & (flux < con + k2 * rms) & (mask != 2))
            | (mask == 1)
        )
        coeff = np.polyfit(wave[idx], flux[idx], deg, w=sig[idx])
        con = np.polyval(coeff, wave)
    rms = np.sqrt(np.mean((con[idx] - flux[idx]) ** 2))    # Re-add mean flux value
    con += fmean
    # optional plot
    if plot:
        pass

    return coeff, con


if __name__ == '__main__':
    main(*sys.argv[1:])

