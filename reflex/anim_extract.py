#!/usr/bin/env python3
import sys
import time
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def getflatimg(img, axis=0):
    idx = np.indices(img.shape)[axis]
    return idx.flatten(), img.flat


def getspecvar(img):
    ny, nx = img.shape
    nimg = img / np.nansum(img, axis=1)[:, None]
    x = np.indices(img.shape)[1]
    return x.flatten(), nimg.flat


def getslitvar(img, xoff, osample=1):
    x = np.indices(img.shape)[0]
    x = x - xoff[None, :]
    return x.flatten() * osample, img.flat


def plot_slitfunction(img, spec, slitf, model, ycen, onum, left, right, osample):

    ny, nx = img.shape
    ny_orig = ny
    size = img.size
    ny_os = len(slitf)
    os = (ny_os - 1) / (ny + 1)

    mask = np.ma.getmaskarray(img)
    idx = np.indices(img.shape)

    di = img.data / np.sum(img)
    ds = spec.data / np.sum(spec)

    df = np.copy(slitf)
    dm = model / np.sum(model)

    max_bad_pixels = 100

    xbad = np.full(max_bad_pixels, np.nan)
    ybad = np.full(max_bad_pixels, np.nan)
    ibad = np.zeros(max_bad_pixels)
    jbad = np.zeros(max_bad_pixels)
    n = np.ma.count_masked(img)
    xbad[:n] = (di / np.nansum(di, axis=1)[:, None])[mask][:max_bad_pixels]
    ybad[:n] = (di[mask] * nx)[:max_bad_pixels]
    ibad[:n] = idx[0][mask][:max_bad_pixels]
    jbad[:n] = idx[1][mask][:max_bad_pixels]

    # on first execution of plot create a new figure
    if not hasattr(plot_slitfunction, "fig"):
        # If no figure exists, create a new one
        plt.ion()
        FIG = plt.figure(figsize=(12, 4))
        FIG.tight_layout(pad=0.05)

        AX1 = FIG.add_subplot(231)
        AX1.set_title("Swath")
        AX2 = FIG.add_subplot(132)
        AX2.set_title("Spectrum")
        AX3 = FIG.add_subplot(133)
        AX3.set_title("Slit")
        AX4 = FIG.add_subplot(234)
        AX4.set_title("Model")

        im1 = AX1.imshow(di, aspect="auto", origin="lower")
        line1, = AX1.plot(ny_orig / 2 + ycen, "-r")
        im4 = AX4.imshow(dm, aspect="auto", origin="lower")

        specvar, = AX2.plot(*getspecvar(di), ".r", ms=2, alpha=0.6)
        slitvar, = AX3.plot(*getslitvar(di * nx, ycen), ".r", ms=2, alpha=0.6)
        slitfu, = AX3.plot(
            np.linspace(0, di.shape[0] + 1, len(df)), df, "-k", lw=3
        ) #TODO which limits for the x axis?

        masked, = AX3.plot(ibad, ybad, "+g")
        masked2, = AX2.plot(jbad, xbad, "+g")

        line2, = AX2.plot(ds, "-k")

        # Save plots to this function
        setattr(plot_slitfunction, "fig", FIG)
        setattr(plot_slitfunction, "ax1", AX1)
        setattr(plot_slitfunction, "ax2", AX2)
        setattr(plot_slitfunction, "ax3", AX3)
        setattr(plot_slitfunction, "ax4", AX4)
        setattr(plot_slitfunction, "ny", ny)
        setattr(plot_slitfunction, "nx", nx)

        setattr(plot_slitfunction, "im1", im1)
        setattr(plot_slitfunction, "line1", line1)
        setattr(plot_slitfunction, "line2", line2)
        setattr(plot_slitfunction, "masked", masked)
        setattr(plot_slitfunction, "masked2", masked2)
        setattr(plot_slitfunction, "specvar", specvar)
        setattr(plot_slitfunction, "slitvar", slitvar)
        setattr(plot_slitfunction, "slitfu", slitfu)
        setattr(plot_slitfunction, "im4", im4)
    else:
        FIG = plot_slitfunction.fig
        AX1 = plot_slitfunction.ax1
        AX2 = plot_slitfunction.ax2
        AX3 = plot_slitfunction.ax3
        AX4 = plot_slitfunction.ax4
        im1 = plot_slitfunction.im1
        line2 = plot_slitfunction.line2
        im4 = plot_slitfunction.im4
        line1 = plot_slitfunction.line1
        masked = plot_slitfunction.masked
        masked2 = plot_slitfunction.masked2
        specvar = plot_slitfunction.specvar
        slitvar = plot_slitfunction.slitvar
        slitfu = plot_slitfunction.slitfu

        ny = plot_slitfunction.ny
        nx = plot_slitfunction.nx

    # Fix size of array
    if di.shape[0] > ny:
        di = di[:ny, :]
        df = df[: ny + 2]
        dm = dm[:ny, :]
    elif di.shape[0] < ny:
        ypad = 0, ny - di.shape[0]
        di = np.pad(di, (ypad, (0, 0)), "constant", constant_values=np.nan)
        df = np.pad(df, (0, ny + 2 - df.shape[0]), "constant", constant_values=np.nan)
        dm = np.pad(dm, (ypad, (0, 0)), "constant", constant_values=np.nan)

    if di.shape[1] > nx:
        di = di[:, :nx]
        ds = ds[:nx]
        dm = dm[:, :nx]
        ycen = ycen[:nx]
    elif di.shape[1] < nx:
        xpad = 0, nx - di.shape[1]
        di = np.pad(di, ((0, 0), xpad), "constant", constant_values=np.nan)
        ds = np.pad(ds, xpad, "constant", constant_values=np.nan)
        dm = np.pad(dm, ((0, 0), xpad), "constant", constant_values=np.nan)
        ycen = np.pad(ycen, xpad, "constant", constant_values=np.nan)

    # Update data
    FIG.suptitle("Order %i, Columns %i - %i" % (onum, left, right))

    im1.set_norm(mcolors.Normalize(vmin=np.nanmin(di), vmax=np.nanmax(di)))
    im1.set_data(di)
    im4.set_norm(mcolors.Normalize(vmin=np.nanmin(dm), vmax=np.nanmax(dm)))
    im4.set_data(dm)

    line1.set_ydata(ny_orig / 2 + ycen)
    line2.set_ydata(ds)

    slitvar.set_data(*getslitvar(di * nx, ycen))
    slitfu.set_ydata(df)
    specvar.set_data(*getspecvar(di))

    masked.set_xdata(ibad)
    masked.set_ydata(ybad)

    masked2.set_xdata(jbad)
    masked2.set_ydata(xbad)

    # Set new limits
    AX1.set_xlim((0, img.shape[1]))
    AX1.set_ylim((0, img.shape[0]))
    AX4.set_xlim((0, img.shape[1]))
    AX4.set_ylim((0, img.shape[0]))

    AX2.set_xlim((0, len(spec)))
    AX3.set_xlim((0, img.shape[0]))
    AX2.set_ylim((0, np.nanmax(di / np.nansum(di, axis=1)[:, None]) * 1.1))
    AX3.set_ylim((0, np.nanmax(di) * nx * 1.1))

    FIG.canvas.draw()
    FIG.canvas.flush_events()
    # plt.show()



class MyHandler(PatternMatchingEventHandler):
    patterns=["*.fits"]

    def process(self, event):
        """
        event.event_type
            'modified' | 'created' | 'moved' | 'deleted'
        event.is_directory
            True | False
        event.src_path
            path/to/observed/file
        """
        if event.src_path == FILE1:

            img = fits.open(FILE1)[0].data
            spec = fits.open(FILE2)[0].data
            slitf = fits.open(FILE3)[0].data
            model = fits.open(FILE4)[0].data
            ycen = fits.open(FILE5)[0].data

            plot_slitfunction(img, spec, slitf, model, ycen, onum=0, left=0, right=2048, osample=2)
        else:
            return

    def on_modified(self, event):
        self.process(event)

    def on_created(self, event):
        self.process(event)


if __name__ == '__main__':

    FILE1 = './debug_img_sw.fits'
    FILE2 = './debug_spc.fits'
    FILE3 = './debug_slitfu.fits'
    FILE4 = './debug_model_sw.fits'
    FILE5 = './debug_ycen.fits'

    img = fits.open(FILE1)[0].data
    spec = fits.open(FILE2)[0].data
    slitf = fits.open(FILE3)[0].data
    model = fits.open(FILE4)[0].data
    ycen = fits.open(FILE5)[0].data

    plot_slitfunction(img, spec, slitf, model, ycen, onum=0, left=0, right=2048, osample=2)



    args = sys.argv[1:]
    observer = Observer()
    observer.schedule(MyHandler(), path=args[0] if args else '.')
    observer.start()

    try:
        while True:
            plt.pause(0.01)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()
