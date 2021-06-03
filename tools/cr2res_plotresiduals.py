#!/usr/bin/env python3
import os, sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from scipy.constants import speed_of_light as c

fig = plt.figure()
ax = fig.subplots()

for filename in sys.argv[1:]:
    f = fits.open(filename)
    sett = f[0].header['HIERARCH ESO INS WLEN ID']
    for d in [1,2,3]:
        dat = f['CHIP%d.INT1'%d].data
        xoff = (d-2)*0.2
        ax.plot(dat['Order']+xoff,dat['Delta_WL'],'.',label='D%d'%d)



    ax.set_title('Residual delta-lambda, %s'%sett)
    ax.legend(loc='upper right')
    ax.set_xlabel('Order #')
    ax.set_ylabel('nm')
    ax.set_xticks(np.arange(9)+1)
    outf = filename.replace('.fits','_resid.png')
    fig.savefig(outf)
    ax.clear()
