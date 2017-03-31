#!/usr/bin/env python3
import sys
import time
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
FIG = plt.figure(figsize=(12,4))

AX1 = FIG.add_subplot(231)
AX2 = FIG.add_subplot(132)
AX3 = FIG.add_subplot(133)
AX4 = FIG.add_subplot(234)

FILE1 = './debug_img_sw.fits'
FILE2 = './debug_spc.fits'
FILE3 = './debug_slitfu.fits'
FILE4 = './debug_model_sw.fits'
FILE5 = './debug_ycen.fits'

def getflatimg(img,axis=0):
    idx = np.indices(img.shape)[axis]
    return idx.flatten(), img.flat

def getspecvar(img):
    ny,nx=img.shape
    nimg = np.transpose(np.transpose(img) / img.sum(axis=1))
    x = np.indices(img.shape)[1]
    return x.flatten(), nimg.flat

def getslitvar(img):
    x = np.indices(img.shape)[0]
    x = x.astype('f')
    with fits.open(FILE5) as f:
        xoff = f[0].data
        for i in range(x.shape[0]):
            x[i,:] -= xoff-1
    return x.flatten(), img.flat

with fits.open(FILE1) as f:
    di = f[0].data
    di /= di.sum()
    ny,nx = di.shape
    im1 = AX1.imshow(di)
    specvar, = AX2.plot(*getspecvar(di),'.r',ms=2,alpha=0.6)
with fits.open(FILE2) as f:
    d = f[0].data
    d /= d.sum()
    spec, = AX2.plot(d,'-k')
with fits.open(FILE3) as f:
    d = f[0].data
    ny_os = len(d)
    os = (ny_os-1) / (ny+1)
    slitvar, = AX3.plot(*getslitvar(di*nx), '.r',ms=2,alpha=0.6)
    X = np.arange(len(d)) / len(d) * (ny)
    slitfu, = AX3.plot(X, d,'-k',lw=3)
with fits.open(FILE4) as f:
    dm = f[0].data
    dm /= dm.sum()
    im2 = AX4.imshow(dm)

FIG.tight_layout(pad=0.05)


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
            with fits.open(FILE1) as f:
                d = f[0].data
                d /= d.sum()
                im1.set_data(d)
                slitvar.set_data(*getslitvar(d*nx))
                specvar.set_data(*getspecvar(d))
            with fits.open(FILE2) as f:
                d2=f[0].data
                d2 /= d2.sum()
                spec.set_ydata(d2)
            with fits.open(FILE3) as f:
                slitfu.set_ydata(f[0].data)
            with fits.open(FILE4) as f:
                d = f[0].data
                d /= d.sum()
                im2.set_data(d)

        else:
            return

    def on_modified(self, event):
        self.process(event)

    def on_created(self, event):
        self.process(event)


if __name__ == '__main__':

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
