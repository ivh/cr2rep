#!/usr/bin/env python3
import os, sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

fig = plt.figure()
ax = fig.subplots()

for filename in sys.argv[1:]:
    f = fits.open(filename)
    sett = f[0].header['HIERARCH ESO INS WLEN ID']
    for d in [1,2,3]:
        dat = f['CHIP%d.INT1'%d].data
        m = dat['FPET_Order']
        lamb = dat['Catalog_WL']
        try:
            ax.plot(lamb,m*lamb,'.',label='D%d'%d)
        except:
            pass


    ax.set_title('Etalon diagnostics, %s'%sett)
    ax.legend()
    ax.set_xlabel('nm')
    ax.set_ylabel('λ * m')
    outf = filename.replace('.fits','_fpet.png')
    fig.savefig(outf)
    ax.clear()
