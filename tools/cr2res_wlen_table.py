#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np


X = np.arange(2048)+1

for twn in sys.argv[1:]:
    print('#File: %s'%twn)
    print('#Order Detector XC StartWL EndWL')
    with fits.open(twn) as twf:
        for order in [9,8,7,6,5,4,3,2,1]:      
            for det in [1,2,3]:
                tw = twf['CHIP%d.INT1'%det].data
                if tw is None: continue
                h = twf['CHIP%d.INT1'%det].header
                xc = h.get('ESO QC WAVE BESTXCORR-%02d-01'%order, 0.0)
                mask = tw['Order']==order
                try:
                    p = tw[mask]['Wavelength'][0][::-1]
                except:
                    continue
                poly = np.poly1d(p)
                #deriv = poly.deriv()
                print(order, det, '%.1f%%'%(xc*100), '%.3f'%poly(1), '%.3f'%poly(2048))
                