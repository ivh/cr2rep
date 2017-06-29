#!/usr/bin/env python3
import os, sys
from astropy.io import fits as F
from pylab import *

filename = sys.argv[1]
Tab = F.open(filename)
all_diffs = [array(ext.data).astype('f') for ext in Tab[1:]]

stds = [diff.std() for diff in all_diffs]
meds = [median(diff) for diff in all_diffs]

for i,diff in enumerate(all_diffs):
    med = median(diff)
    std = diff.std()
    bins=linspace(med/3,med*2.2,100)
    bins.sort()
    clf()
    hist(diff,bins=bins)
    axvline(med,ls='--',color='k')
    axvline(med/2,ls='--',color='k')
    axvline(med*2,ls='--',color='k')
    axvline(diff.mean(),ls='-',color='k')
    axvline(med-std,ls='-',color='k')
    axvline(med+std,ls='-',color='k')
    savefig('hist_%s.png'%i,dpi=100)

