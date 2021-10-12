import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

# Usage:
# python3 cr2rem_drs2molecfit.py [your-cr2res-extracted-spectrum.fits]
# where the argument is a cr2res DRS extracted spectrum

# Output:
# (1) SCIENCE.fits - SCIENCE file for Molecfit
# (2) WAVE_INCLUDE.fits - WAVE_INCLUDE file for Molecfit
# (3) ATM_PARAMETERS_EXT.fits - ATM_PARAMETERS_EXT file for Molecfit


fi = sys.argv[1] # input file you want to convert
output = "SCIENCE.fits" # name for output file

hd = fits.open(fi) # open data fits file
primary_hdu = hd[0]  # keep primary HDU with header

hdul_output = fits.HDUList([primary_hdu]) # initialize output HDU list with primary HDU
hdul_winc = fits.HDUList([primary_hdu]) # initialise WAVE_INCLUDE HDU list
hdul_atm =  fits.HDUList([primary_hdu]) # initialise ATM_PARAMETERS_EXT HDU list

wmin,wmax,map2chip, atm_parms_ext = [], [], [], [] # initialize empty arrays
jjj = 1 # initialize counter for order/detector

for idet in range(3): # loop on detectors
    chip = 'CHIP'+str(idet+1)+'.INT1' # set up chip name
    data = hd[chip].data # isolate data for given detector
    cp = np.sort([i[0:6] for i in data.dtype.names if "WL" in i]) # find all orders
    for iord in range(len(cp)): # loop on orders
        cpw,cps,cpe = cp[iord]+'WL', cp[iord]+'SPEC', cp[iord]+'ERR'  # set up column name for WAVE SPEC and ERR
        w,s,e = hd[chip].data[cpw], hd[chip].data[cps],hd[1].data[cpe] # isolate WAVE SPEC and ERR for given order/detector
        col1 = fits.Column(name='WAVE', format='D', array=w) # create fits column for WAVE
        col2 = fits.Column(name='SPEC', format='D', array=s) # create fits column for SPEC
        col3 = fits.Column(name='ERR', format='D', array=e) # create fits column for ERR
        table_hdu = fits.BinTableHDU.from_columns([col1, col2, col3]) # create fits table HDU with WAVE SPEC ERR
        hdul_output.append(table_hdu) # append table HDU to output HDU list
        wmin.append(np.min(w)*0.001) # append wmin for given order to wmin array and convert to micron
        wmax.append(np.max(w)*0.001) # append wmax for giver order to wmax array and convert to micron
        map2chip.append(jjj) # append order counter to map2chip
        if jjj == 1:
            atm_parms_ext.append(0)
        else:
            atm_parms_ext.append(1)
        jjj+=1

# setup the columns for the WAVE_INCLUDE output file
wmin_winc = fits.Column(name='LOWER_LIMIT', format='D', array=wmin)
wmax_winc = fits.Column(name='UPPER_LIMIT', format='D', array=wmax)
map2chip_winc = fits.Column(name='MAPPED_TO_CHIP', format='I', array=map2chip)
# creates WAVE_INCLUDE table HDU and appends to the HDU list
table_hdu_winc = fits.BinTableHDU.from_columns([wmin_winc, wmax_winc,map2chip_winc])
hdul_winc.append(table_hdu_winc)

# Write to the output file
hdul_output.writeto(output, overwrite=True)

# Write to WAVE_INCLUDE file
hdul_winc.writeto('WAVE_INCLUDE.fits', overwrite=True)

col1_atm = fits.Column(name='ATM_PARAMETERS_EXT', format='I', array=atm_parms_ext)
table_hdu_atm = fits.BinTableHDU.from_columns([col1_atm])
hdul_atm.append(table_hdu_atm)
hdul_atm.writeto('MAPPING_ATMOSPHERIC.fits', overwrite=True)
hd.close()
