#!/usr/bin/env python3
import os, sys
from astropy.io import fits
import numpy as np

from numbers import Number, Integral


def get_angles(fname_trace):
    """ compare img and trace """
    trace = fits.open(fname_trace)
    header = trace[0].header
    setting = header["ESO INS WLEN ID"]
    angles = {}
    for det in [1, 2, 3]:
        try:
            tdata = trace["CHIP%d.INT1" % det].data
            thead = trace["CHIP%d.INT1" % det].header
        except:
            print("extension %s is missing, skipping." % i)
            continue
        if tdata is None:
            print("Data for CHIP%s is empty, skipping." % i)
            continue

        for t in tdata:
            order = t["Order"]
            alla = t["All"]
            poly = np.polynomial.Polynomial(alla)
            deriv = poly.deriv(1)

            if order not in angles:
                angles[order] = {"setting": setting, "order": order}

            for i in [100, 1024, 1948]:
                # Evaluate the curvature polynomials to get coefficients
                b = deriv(i)
                # angle:
                #     /|
                #    / | b
                # --/--|--
                # A/  1
                angle = np.arctan(b) * 180 / np.pi
                angles[order][f"D{det},x={i}"] = angle

    angles = list(angles.values())
    return angles


if __name__ == "__main__":
    files = sys.argv[1:]

    angles = []
    for file in files:
        angles += get_angles(file)

    keys = ["setting", "order"]
    for d in [1, 2, 3]:
        for x in [100, 1024, 1948]:
            keys += [f"D{d},x={x}"]

    header = "\t".join(keys)
    text = [header]
    for angle in angles:
        values = []
        for key in keys:
            value = angle.get(key, "-")
            if isinstance(value, Integral):
                value = str(value)
            elif isinstance(value, Number):
                value = f"{value:.2f}"                
            values += [value]
        line = "\t".join(values)
        text += [line]

    text = "\n".join(text)
    print(text)
