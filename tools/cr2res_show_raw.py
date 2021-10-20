#!/usr/bin/env python3
import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

FIG = plt.figure(figsize=(10, 3.5))

vmi,vma=None,None
try:
    vmi,vma=[float(v) for v in sys.argv[1].split(',')]
    offs=1
except:
    offs=0

for fname_img in sys.argv[1+offs:]:
    if not fname_img.endswith('.fits'): continue
    print(fname_img)
    FIG.clf()
    axs = FIG.subplots(1,3)
    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])

    img = fits.open(fname_img)

    head = img[0].header
    dit = head['HIERARCH ESO DET SEQ1 DIT']
    ndit = head['HIERARCH ESO DET NDIT']

    for i in [1, 2, 3]:
        imgdata = img['CHIP%d.INT1'%i].data
        imgdata = np.nan_to_num(imgdata)
        #imgdata /= float(dit)
        vmin = vmi or np.percentile(imgdata,5)
        vmax = vma or np.percentile(imgdata,98)
        axs[i-1].imshow(imgdata, origin='lower', cmap='plasma',
            vmin=vmin, vmax=vmax)

        #axs[i-1].axis((1,2048,1,2048))


    FIG.suptitle(f'{fname_img} - cuts [{vmin:.1f} : {vmax:.1f}], DIT:%.01f NDIT:%d'%(dit,ndit), fontsize=8)
    FIG.tight_layout(pad=0.02)
    figname = fname_img.replace('.fits','.png')
    plt.savefig(figname,dpi=240)


