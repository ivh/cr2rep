#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

HALFDET = 1024


FIG = plt.figure(figsize=(10,3.5))

fname_img = sys.argv[1]
fname_tilt = sys.argv[2]
print('Showing tilt %s on image %s'%(fname_tilt, fname_img))

img = fits.open(fname_img)
tilt = fits.open(fname_tilt)
X = np.arange(2*HALFDET)

for i in [1,2,3]:
    ax = FIG.add_subplot(1,3,i)
    ax.set_xticks([])
    ax.set_yticks([])

    idata = img[i].data
    ax.imshow(idata)
    axi = plt.axis()

    tdata = tilt['CHIP%s'%i].data
    if tdata is None:
        print('No data for CHIP%s, skipping.'%i)
        continue
    for alla, upper, lower, order, tiltnb, wave in tdata:
        pol = np.polyval(upper[::-1],X)
        ax.plot(X, pol, ':w')

        pol = np.polyval(lower[::-1],X)
        ax.plot(X, pol, ':w')

        pol = np.polyval(alla[::-1],X)
        ax.plot(X, pol, '--w')
        if np.isnan(pol[HALFDET]): continue
        ax.text(HALFDET,pol[HALFDET] ,'%s-%s'%(order,tiltnb),color='w',horizontalalignment='center',
            verticalalignment='center', size=19)
        print('order: %s, y: %s'%(tiltnb,pol[HALFDET]))

    plt.axis(axi)


FIG.tight_layout(pad=0.02)
plt.show()
