#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from turbo_colormap import turbo_colormap_data
Turbo = ListedColormap(turbo_colormap_data)

def compare(fname_trace, fname_img=None):
    """ compare img and trace """
    trace = fits.open(fname_trace)
    if fname_img:
        img = fits.open(fname_img)
        linecol = "w"
    else:
        linecol = "k"
    
    X = np.arange(2048)
    FIG = plt.figure(figsize=(15, 5))

    for det in [1, 2, 3]:
        ax = FIG.add_subplot(1, 3, det)
        ax.set_xticks([])
        ax.set_yticks([])

        try:
            tdata = trace['CHIP%d.INT1'%det].data
        except:
            print("extension %s is missing, skipping." % i)
            continue
        if tdata is None:
            print("Data for CHIP%s is empty, skipping." % i)
            continue

        if fname_img:
            imgdata = img['CHIP%d.INT1'%det].data
            imgdata = np.ma.masked_where(np.isnan(imgdata), imgdata)
            vmin, vmax = np.percentile(imgdata.compressed(), (5, 95))
            vmax += (vmax-vmin)*0.4
            ax.imshow(imgdata, origin="lower", vmin = vmin, vmax=vmax,
                cmap='viridis')

        for t in tdata:
            upper = t["Upper"]
            lower = t["Lower"]
            alla = t["All"]
            slitfrac = t["SlitFraction"]
            order = t["Order"]
            wave = t["Wavelength"]
            ca = t["SlitPolyA"]
            cb = t["SlitPolyB"]
            cc = t["SlitPolyC"]

            upper = np.polyval(upper[::-1], X)
            ax.plot(X, upper, ":" + linecol)

            lower = np.polyval(lower[::-1], X)
            ax.plot(X, lower, ":" + linecol)

            middle = np.polyval(alla[::-1], X)
            ax.plot(X, middle, "--" + linecol)

            i1 = tdata[tdata["order"] == order]["Slitfraction"][:, 1]
            i2 = tdata[tdata["order"] == order]["All"]
            coeff = [np.interp(0.5, i1, i2[:, k]) for k in range(i2.shape[1])]


            for i in range(30, 2048, 200):
                ew = [int(middle[i] - lower[i]), int(upper[i] - middle[i])]
                x = np.zeros(ew[0] + ew[1] + 1)
                y = np.arange(-ew[0], ew[1] + 1).astype(float)

                # Evaluate the curvature polynomials to get coefficients
                a = np.polyval(ca[::-1], i)
                b = np.polyval(cb[::-1], i)
                c = np.polyval(cc[::-1], i)
                yc = np.polyval(coeff[::-1], i)

                # Shift polynomials to the local frame of reference
                a = a - i + yc * b + yc * yc * c
                b += 2 * yc * c

                for j, yt in enumerate(y):
                    x[j] = i + yt * b + yt ** 2 * c
                y += middle[i]
                plt.plot(x, y, "-"+linecol)

            if np.isnan(middle[1024]):
                continue
            ax.text(
                500,
                middle[1024],
                "order: %s\ntrace: %s" % (t["order"], t["TraceNb"]),
                color=linecol,
                horizontalalignment="center",
                verticalalignment="center",
                size=9,
            )
            ax.axis((1,2048,1,2048))

    FIG.tight_layout(pad=0.02)
    #plt.show()
    figname = fname_trace.replace(".fits", ".png")
    plt.savefig(figname, dpi=180)


if __name__ == "__main__":
    compare(*sys.argv[1:])
