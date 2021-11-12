#!/usr/bin/env python3
import os, sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy.constants import speed_of_light as sol

fig = plt.figure(figsize=(8,5))
ax = fig.subplots()
try:
    vmi,vma=[float(v) for v in sys.argv[1].split(',')]
    offs=1
except:
    vmi,vma=None,None
    offs=0

def dlam2vel(lam0):
    return lambda dlam: dlam/lam0*sol

def vel2dlam(lam0):
    return lambda vel: vel/sol*lam0

def dlam2pix(disp):
    return lambda dlam: dlam/disp

def pix2dlam(disp):
    return lambda pix: disp*pix

for filename in sys.argv[1+offs:]:
    f = fits.open(filename)
    sett = f[0].header['HIERARCH ESO INS WLEN ID']
    cwlen = f[0].header['HIERARCH ESO INS WLEN CWLEN']
    disp = f[0].header['HIERARCH ESO INS GRAT1 DISP']
    for d in [1,2,3]:
        dat = f['CHIP%d.INT1'%d].data
        xoff = (d-2)*0.2
        try:
            ax.plot(dat['Order']+xoff,dat['Delta_WL'],'.',label='D%d'%d)
        except:
            pass


    ax.set_title('Residual delta-lambda, %s'%sett)
    ax.legend(loc='upper right')
    ax.set_xlabel('Order #')
    ax.set_ylabel(r'$\Delta \lambda$ [nm]')
    ax.set_xticks(np.arange(9)+1)
    if vmi:
        ax.set(ylim=(vmi,vma))

    secax = ax.secondary_yaxis('right',functions=(dlam2vel(cwlen),vel2dlam(cwlen)))
    secax.set_ylabel(r'$\Delta v$ [m/s]')

    triax = ax.secondary_yaxis(-0.2,functions=(dlam2pix(disp),pix2dlam(disp)))
    triax.set_ylabel(r'$\Delta$ pixels ')

    fig.tight_layout()
    outf = filename.replace('.fits','_resid.png')
    fig.savefig(outf)
    ax.clear()
