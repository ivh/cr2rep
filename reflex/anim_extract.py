#!/usr/bin/env python3
import sys
import time
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def getosidx(i,os):
    return os*(i+1)+1
def getflatimg(img,axis=0,os=1):
    idx = np.indices(img.shape)[axis]
    return getosidx(idx.flatten(),os), img.flat

plt.ion()
FIG = plt.figure(figsize=(10,6))
AX1 = FIG.add_subplot(121)
AX2 = FIG.add_subplot(222)
AX3 = FIG.add_subplot(224)

FILE1 = './img_sw.fits'
FILE2 = './spc.fits'
FILE3 = './slitfu.fits'

with fits.open(FILE1) as f:
    di = f[0].data
    ny,nx = di.shape
    im1 = AX1.imshow(di)
    specvar, = AX2.plot(*getflatimg(di,1),'.r',ms=2,alpha=0.6)
with fits.open(FILE2) as f:
    spec, = AX2.plot(f[0].data,'-k')
with fits.open(FILE3) as f:
    d = f[0].data
    ny_os = len(d)
    os = (ny_os-1) / (ny+1)
    print('os: %s %s %s'%(os,ny_os,ny))
    slitvar, = AX3.plot(*getflatimg(di,0,os),'.r',ms=2,alpha=0.6)
    slitfu, = AX3.plot(d,'-k',lw=3)



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
                im1.set_data(d)
                d1 = getflatimg(d,0)
                slitvar.set_data(*d1)
                d1 = getflatimg(d,1)
                specvar.set_data(*d1)
            with fits.open(FILE2) as f:
                spec.set_ydata(f[0].data)
            with fits.open(FILE3) as f:
                slitfu.set_ydata(f[0].data)
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
