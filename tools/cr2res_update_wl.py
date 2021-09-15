#!/usr/bin/env python3
from astropy.io import fits
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser("CRIRES+ Update Wavelength Tool", description="Use the wavelength polynomials in the trace wave to replace the wavelength columns in the extracted spectrum")
    parser.add_argument("extr1d", help="FITS file with the extracted spectrum")
    parser.add_argument("trace_wave", help="FITS file with the trace wave table")
    parser.add_argument("--out", help="Filename of the new FITS file, if not given will overwrite the existing extr1d file")
    args = parser.parse_args()

    extr1d_fname = args.extr1d
    tw_fname = args.trace_wave
    out_fname = args.out

    # Load the files
    extr1d_hdu = fits.open(extr1d_fname)
    tw_hdu = fits.open(tw_fname)
    if out_fname is None:
        out_fname = extr1d_fname

    for chip in ["CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1"]:
        tw_chip = tw_hdu[chip]
        extr1d_chip = extr1d_hdu[chip]

        tw_order = tw_chip.data["Order"]
        tw_traceNb = tw_chip.data["TraceNb"]
        tw_wave = tw_chip.data["Wavelength"]

        # We only care about the WL columns
        extr1d_colnames = extr1d_chip.data.columns.names
        extr1d_wave_order_trace = [(c, *c.split("_")[:2]) for c in extr1d_colnames if c.endswith("WL")]

        # The number of pixels is equal to the number of columns
        x = np.arange(len(extr1d_chip.data))

        for colname, order, trace in extr1d_wave_order_trace:
            order = int(order)
            trace = int(trace)

            # Find the correct polynomials
            tw_idx = (tw_order == order) & (tw_traceNb == trace)
            try:
                tw_coef = tw_wave[tw_idx][0]
            except:
                print("Could not find data for Order {order} TraceNb {trace} in the trace wave")
                continue
            # and calculate the new wavelength points
            new_wave = np.polyval(tw_coef[::-1], x)
            # So we can replace them in the extr1d file
            extr1d_chip.data[colname] = new_wave

    extr1d_hdu.writeto(out_fname, overwrite=True, output_verify="fix+warn")
