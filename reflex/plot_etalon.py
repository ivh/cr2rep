#!/usr/bin/env python3
import os, sys
import matplotlib.pyplot as plt
plt.style.use('seaborn-talk')
import numpy as np
from astropy.io import fits as F

filename = sys.argv[1]
Tab = F.open(filename)
all_diffs = [np.array(ext.data).astype('f') for ext in Tab[1:]]

stds = [diff.std() for diff in all_diffs]
meds = [np.median(diff) for diff in all_diffs]

for i,diff in enumerate(all_diffs):
    med = np.median(diff)
    std = diff.std()
    bins=np.linspace(med/2,med*1.5,100)
    bins.sort()
    plt.clf()
    plt.hist(diff,bins=bins)
    plt.axvline(med,ls='--',color='k')
    plt.axvline(med/2,ls='--',color='k')
    plt.axvline(med*2,ls='--',color='k')
    plt.axvline(diff.mean(),ls='-',color='k')
    plt.axvline(med-std,ls='-',color='k')
    plt.axvline(med+std,ls='-',color='k')
    plt.xlabel('Δλ')
    plt.ylabel('N')
    plt.title('Detector %d, order #%s'%(i//9+1,i%9))
    plt.savefig('hist_%02d.png'%i,dpi=60)

