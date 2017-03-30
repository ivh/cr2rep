#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

FIG = plt.figure(figsize=(10,3.5))

fname_flat = sys.argv[1]
fname_trace = sys.argv[2]
print('Comparing trace %s to image %s'%(fname_trace, fname_flat))

flat = fits.open(fname_flat)
trace = fits.open(fname_trace)
X = np.arange(2048)

for i in [1,2,3]:
    ax = FIG.add_subplot(1,3,i)
    ax.set_xticks([])
    ax.set_yticks([])

    fdata = flat[i].data
    ax.imshow(fdata)
    axi = plt.axis()

    tdata = trace['CHIP%s'%i].data
    if tdata is None:
        print('No data for CHIP%s, skipping.'%i)
        continue
    for alla, upper, lower, order in tdata:
        pol = np.polyval(upper[::-1],X)
        ax.plot(X, pol, ':w')

        pol = np.polyval(lower[::-1],X)
        ax.plot(X, pol, ':w')

        pol = np.polyval(alla[::-1],X)
        ax.plot(X, pol, '--w')
        if np.isnan(pol[1024]): continue
        ax.text(1024,pol[1024] ,'%s'%order,color='w',horizontalalignment='center',
            verticalalignment='center', size=18)
        print('order: %s, y: %s'%(order,pol[1024]))

    plt.axis(axi)


FIG.tight_layout(pad=0.02)
plt.show()
