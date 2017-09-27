#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

FIG = plt.figure(figsize=(20,8))
ax = FIG.add_subplot(111)

fname_spec = sys.argv[1]
fname_trace = sys.argv[2]
fname_cat = sys.argv[3]

spec= fits.open(fname_spec)
trace = fits.open(fname_trace)

detstyle=('-k','-b','-g')

for i in [1,2,3]:
    sdata = spec['CHIP%d'%i].data
    tdata = trace['CHIP%d'%i].data
    if tdata is None or sdata is None:
        print('No data for CHIP%d, skipping.'%i)
        continue

    for alla, upper, lower, order, tracenb, wave in tdata:
        if order==-1: continue
        try:
            sp = sdata['%02d_%02d_SPEC'%(order,tracenb)]
        except KeyError:
            print('%02d_%02d_SPEC not found'%(order,tracenb))
            continue


        nx = len(sp)
        X = np.arange(nx)
        wave = np.polyval(wave[::-1],X)
        ax.plot(wave, sp, detstyle[i-1])
        if np.isnan(wave[1024]): continue
        ax.text(wave[1024], -2500,'%s-%s'%(i,order),horizontalalignment='center',
            verticalalignment='center', size=9)



cat_w, cat_i = np.loadtxt(fname_cat,unpack=True)
ax.plot(cat_w,cat_i*10,'ro')


ax.set_yticks([])
FIG.tight_layout(pad=0.02)
ax.axis((1300,1360,-5E3,2E5))
plt.show()
