#!/usr/bin/env python3
import sys
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("svg")
import matplotlib.pyplot as plt

fig = plt.figure()
axs=[fig.add_subplot(131)]
axs += [fig.add_subplot(132, sharey=axs[0]),
         fig.add_subplot(133, sharey=axs[0])]


for filename in sys.argv[1:]:
    hdu = fits.open(filename)

    for i in range(1, len(hdu)):
        data = hdu[i].data
        slit_image = np.zeros((2048, 2048))

        for j in range(10):
            try:
                x = data[f"{j:02}_01_POSITIONX"]
                y = data[f"{j:02}_01_POSITIONY"]
                slit_frac = data[f"{j:02}_01_SLIT_FRACTION"]
            except KeyError:
                continue
            mask = (np.isfinite(x)) & (np.isfinite(y)) & (np.isfinite(slit_frac))
            x = x[mask].astype(int) - 1
            y = y[mask].astype(int) - 1
            slit_frac = slit_frac[mask]
            slit_image[y, x] = slit_frac

        ax = axs[i-1]
        ax.imshow(slit_image, origin="lower")
        ax.set_xlabel("x [px]")
        ax.set_ylabel("y [px]")

    outf = filename.replace('.fits','_slitfrac.png')
    fig.savefig(outf)
