#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np
from glob import glob
import subprocess
from matplotlib import pyplot as plt

C2RAWDIR = os.environ['C2RAWDIR']
C2REDUCEDDIR = os.environ['C2REDUCEDDIR']
os.chdir(C2REDUCEDDIR)

input_cat_name = sys.argv[1]
linefile = fits.open(input_cat_name)
lines = np.sort(linefile[1].data, order='Wavelength')[:-1]

input_selection_name = sys.argv[2]
GOOD = open(input_selection_name).readlines()

mask = lines['Wavelength'] < 0.0

for goo in GOOD:
    if goo.startswith('#'): continue
    start,end=[ float(x) for x in goo.split(',') ]
    m = (lines['Wavelength'] < end) & (lines['Wavelength'] > start)
    mask = mask | m

selbasename = os.path.split(input_selection_name)[-1]
outname = '.'.join(selbasename.split('.')[:-1]) + '.fits'
linefile[1].data=np.sort(lines[mask], order='Wavelength')
linefile.writeto(outname,overwrite=True)