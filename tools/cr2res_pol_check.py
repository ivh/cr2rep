#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["numpy", "astropy", "matplotlib"]
# ///
"""Quick-look assessment of cr2res_obs_pol products.

Usage:
    cr2res_pol_check.py [-v] [--plot] [--trend [PNG]] cr2res_obs_pol_pol_spec[AB].fits ...

For each file it prints the recipe's QC parameters, the noise level of the
Stokes and null spectra relative to their formal errors, continuum offsets,
and how much polarisation signal there is above the null.

It also measures how much of each spectrum follows dI/dlambda, expressed as a
velocity shift. In the null that reveals a wavelength shift between the beams
or exposures. In Stokes V its sign tracks that of the longitudinal field, but
the convention is not calibrated and neither is the size: the weak-field relation does not hold for Zeeman-resolved lines and
telluric lines dilute it.

The QC parameters SNR, POL NULL SCATTER, POL NULL MED and POL BEAMRATIO are read
from the headers. Products of cr2res versions that do not write them yet get
the same numbers computed here, in the same way as the recipe. The beam ratio
needs the extracted beams (cr2res_obs_pol_extracted[AB]_[1-4][ud].fits, written
with --save_group=1) next to the pol_spec file.

With --trend the per-file medians of the QC parameters are plotted against
the time of observation, to follow them over a night or over months. The
colour tells the four polarimeters apart: circular (V) or linear (Q, U), in
YJ or HK.

When both nodding positions of a set are given (specA and specB from the same
directory), the two are correlated against each other: a real signal repeats
in A and B, noise and most systematics do not.
"""

import argparse
import os
import re
import sys

import numpy as np
from astropy.io import fits

EDGE = 20  # pixels dropped at both ends of each order
HIPASS = 101  # running-median width that removes continuum offsets [pixels]
C_KMS = 299792.458
# Stokes rms over null rms above which an order counts as showing signal
SIGNAL_RATIO = 1.3

TREND_DEFAULT = "cr2res_pol_check_trend.png"

# Limits behind the flags in the verdict line
MAX_NULL_SCATTER = 2.0
MAX_NULL_OFFSET = 5e-3
BEAMRATIO_RANGE = (0.8, 1.25)
MIN_SNR = 100.0


def robust_std(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def running_median(x, width):
    half = width // 2
    padded = np.pad(x, half, mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(padded, width)
    return np.median(windows, axis=1)


def highpass(x):
    return x - running_median(np.where(np.isfinite(x), x, 0.0), HIPASS)


class Order:
    """One spectral order of one detector, with its statistics."""

    def __init__(self, chip, idx, table, header, recompute=False):
        self.chip, self.idx = chip, idx
        self.qc, self.computed = {}, set()
        self._recipe_qc(table, header, recompute)
        sl = slice(EDGE, -EDGE)
        col = lambda name: np.asarray(table[f"{idx}_{name}"], dtype=float)[sl]
        self.wl = col("WL")
        self.stokes, self.stokes_err = col("STOKES"), col("STOKES_ERR")
        self.null, self.null_err = col("NULL"), col("NULL_ERR")
        self.intens = col("INTENS")
        self.ok = np.isfinite(self.stokes) & np.isfinite(self.null) & np.isfinite(self.intens) & (self.intens > 0)
        self.usable = self.ok.sum() > 200
        if not self.usable:
            return
        self.stokes_hp = highpass(self.stokes)
        self.null_hp = highpass(self.null)
        ok = self.ok
        self.wl_mid = np.median(self.wl[ok])
        self.stokes_med = np.median(self.stokes[ok])
        self.null_med = np.median(self.null[ok])
        # cr2res before the fix of cr2res_pol_resample writes NaN errors
        self.stokes_scatter = self._scatter(self.stokes_hp, self.stokes_err)
        self.null_scatter = self._scatter(self.null_hp, self.null_err)
        self.stokes_rms = np.std(self.stokes_hp[ok])
        self.null_rms = np.std(self.null_hp[ok])
        self.excess = np.sqrt(max(self.stokes_rms**2 - self.null_rms**2, 0.0))
        self.peak = np.percentile(np.abs(self.stokes_hp[ok]), 99.5)

    def _scatter(self, spec, err):
        ok = self.ok & np.isfinite(err) & (err > 0)
        return robust_std(spec[ok] / err[ok]) if ok.sum() > 200 else np.nan

    def _recipe_qc(self, table, header, recompute):
        """The recipe's per-order QC, from the header or computed as the recipe does."""
        n = int(self.idx)
        full = lambda name: np.asarray(table[f"{self.idx}_{name}"], dtype=float)

        def resolve(key, name, compute):
            value = None if recompute else header.get(f"ESO QC {name}{n}")
            if value is None:
                value = compute()
                self.computed.add(key)
            # the recipe writes -1 where it could not compute a value
            self.qc[key] = None if key != "null_med" and value is not None and value < 0 else value

        def snr():
            # cr2res_qc_compute_snr: median of spec/err, unusable errors set to 1
            spec, err = full("INTENS"), full("INTENS_ERR").copy()
            err[~np.isfinite(err) | (np.abs(err) < 1e-3)] = 1.0
            value = np.nanmedian(spec / err) if np.isfinite(spec).any() else np.nan
            return float(value) if np.isfinite(value) else -1.0

        def null_scatter():
            # cr2res_qc_pol_null_scatter: 1.4826 * MAD of null/err
            null, err = full("NULL"), full("NULL_ERR")
            ok = np.isfinite(null) & np.isfinite(err) & (err > 0)
            if ok.sum() < 10:
                return None if not np.isfinite(err).any() else -1.0
            ratio = null[ok] / err[ok]
            return float(1.4826 * np.median(np.abs(ratio - np.median(ratio))))

        def null_med():
            # cr2res_qc_pol_null_offset: median of the finite null pixels
            null = full("NULL")
            return float(np.median(null[np.isfinite(null)])) if np.isfinite(null).any() else None

        resolve("snr", "SNR", snr)
        resolve("null_scatter", "POL NULL SCATTER", null_scatter)
        resolve("null_med", "POL NULL MED", null_med)
        self.qc["realorder"] = header.get(f"ESO QC REALORDER{n}")
        value = None if recompute else header.get(f"ESO QC POL BEAMRATIO{n}")
        self.qc["beamratio"] = value if value is None or value > 0 else None
        if value is None:
            self.computed.add("beamratio")  # filled in by PolSpec if the beams are there

    def derivative_shift(self, which):
        """Slope of spectrum against dlnI/dlambda, as a velocity [m/s].

        Only pixels with a significant intensity gradient are used: elsewhere
        the gradient is noise, which would pull the slope towards zero.
        """
        spec = self.stokes_hp if which == "stokes" else self.null_hp
        with np.errstate(all="ignore"):
            x = np.gradient(self.intens, self.wl) / self.intens
        ok = self.ok & np.isfinite(x)
        ok &= np.abs(x) > 4 * robust_std(x[ok])
        if ok.sum() < 20:
            return None
        slope_nm = np.sum(x[ok] * spec[ok]) / np.sum(x[ok] ** 2)
        return 1e3 * C_KMS * slope_nm / self.wl_mid


def fit_shift(orders, which):
    """Median over orders of the dI/dlambda slope [m/s], its error, the significance."""
    shifts = np.array([v for o in orders if (v := o.derivative_shift(which)) is not None])
    if shifts.size < 3:
        return np.nan, np.nan, np.nan
    shift = np.median(shifts)
    err = 1.25 * robust_std(shifts) / np.sqrt(shifts.size)
    return shift, err, shift / err if err > 0 else np.nan


class PolSpec:
    def __init__(self, path, recompute=False):
        self.path = path
        with fits.open(path) as hdul:
            self.header = hdul[0].header.copy()
            self.orders = []
            for hdu in hdul[1:]:
                if hdu.data is None:
                    continue
                indices = sorted(
                    m.group(1) for name in hdu.data.names if (m := re.match(r"(\d+)_STOKES$", name))
                )
                for idx in indices:
                    self.orders.append(Order(hdu.name, idx, hdu.data, hdu.header, recompute))
        self.good = [o for o in self.orders if o.usable]
        h = self.header
        self.stokes = str(h.get("ESO INS POL TYPE", "?")).strip()
        catg = str(h.get("ESO PRO CATG", ""))
        self.nod = catg[-1] if catg.endswith(("SPECA", "SPECB")) else "?"
        band = str(h.get("ESO INS WLEN ID", "?")).strip()[:1]
        arm = "YJ" if band in "YJ" else "HK" if band in "HK" else "?"
        kind = "circ" if self.stokes == "V" else "lin" if self.stokes in ("Q", "U") else "?"
        self.polarimeter = f"{kind} {arm}"
        self.key = (os.path.dirname(os.path.abspath(path)), h.get("DATE-OBS"), self.stokes)
        self.beams_found = None
        if any("beamratio" in o.computed for o in self.orders):
            self.beams_found = self._beam_ratios()

    def _beam_ratios(self):
        """cr2res_qc_pol_beam_ratio from the extracted beams, if they are next to the file.

        Geometric mean over the four exposures of median(up) / median(down).
        """
        folder, name = os.path.split(os.path.abspath(self.path))
        match = re.match(r"(.*)pol_spec([AB])(\.fits.*)$", name)
        if not match:
            return False
        prefix, nod, suffix = match.groups()
        paths = [os.path.join(folder, f"{prefix}extracted{nod}_{k}{b}{suffix}") for k in "1234" for b in "ud"]
        if not all(os.path.exists(p) for p in paths):
            return False
        beams = [fits.open(p) for p in paths]
        try:
            for o in self.orders:
                if "beamratio" not in o.computed:
                    continue
                logsum = 0.0
                for k in range(4):
                    medians = []
                    for hdul in beams[2 * k: 2 * k + 2]:
                        try:
                            spec = np.asarray(hdul[o.chip].data[f"{o.idx}_01_SPEC"], dtype=float)
                        except KeyError:
                            spec = np.array([np.nan])
                        finite = spec[np.isfinite(spec)]
                        medians.append(np.median(finite) if finite.size else np.nan)
                    if not (medians[0] > 0 and medians[1] > 0):
                        logsum = np.nan
                        break
                    logsum += np.log(medians[0] / medians[1])
                o.qc["beamratio"] = float(np.exp(logsum / 4)) if np.isfinite(logsum) else None
        finally:
            for hdul in beams:
                hdul.close()
        return True

    def title(self):
        h = self.header
        target = h.get("ESO OBS TARG NAME") or h.get("OBJECT")
        return (
            f"{target}  {h.get('ESO INS WLEN ID')}  Stokes {self.stokes}  nod {self.nod}  "
            f"DIT {h.get('ESO DET SEQ1 DIT')} s x NDIT {h.get('ESO DET NDIT')}  {h.get('DATE-OBS')}  "
            f"[{h.get('ESO PRO REC1 PIPE ID')}]"
        )


def med(values):
    values = [v for v in values if v is not None and np.isfinite(v)]
    return np.median(values) if values else np.nan


def span(values, fmt):
    values = [v for v in values if v is not None and np.isfinite(v)]
    if not values:
        return "n/a"
    return f"{format(np.median(values), fmt)} ({format(min(values), fmt)} .. {format(max(values), fmt)})"


def report(spec, verbose):
    print("=" * 100)
    print(spec.path)
    print(spec.title())
    good = spec.good
    if not good:
        print("  no usable orders")
        return
    print(f"  {len(good)} of {len(spec.orders)} orders usable")

    if verbose:
        print("\n  chip        ord  real   wl[nm]    SNR  beamratio  null_scat  S_scat  "
              "null_med[%]  S_med[%]  excess[%]  peak[%]")
        for o in good:
            q = o.qc
            num = lambda v, f: format(v, f) if v is not None else format("n/a", ">" + f.split(".")[0])
            print(f"  {o.chip:10s}  {o.idx}   {num(q['realorder'], '3d')}  {o.wl_mid:8.2f}  "
                  f"{num(q['snr'], '5.0f')}  {num(q['beamratio'], '9.3f')}  {o.null_scatter:9.2f}  "
                  f"{o.stokes_scatter:6.2f}  {o.null_med*100:11.3f}  {o.stokes_med*100:8.3f}  "
                  f"{o.excess*100:9.3f}  {o.peak*100:7.3f}")

    def source(key):
        n = sum(key in o.computed for o in good)
        return "" if n == 0 else "   [computed here]" if n == len(good) else f"   [computed here for {n} orders]"

    print("\n  Recipe QC, median over orders (min .. max)")
    print(f"    SNR per pixel                 {span([o.qc['snr'] for o in good], '.0f')}{source('snr')}")
    if spec.beams_found is False:
        print("    beam ratio up/down            n/a   [not in the header; needs the extracted beams next to "
              "this file, from --save_group=1]")
    else:
        print(f"    beam ratio up/down            {span([o.qc['beamratio'] for o in good], '.3f')}{source('beamratio')}")
    print(f"    null scatter / formal error   {span([o.qc['null_scatter'] for o in good], '.2f')}{source('null_scatter')}")
    print(f"    null median [%]               {span([None if o.qc['null_med'] is None else 100*o.qc['null_med'] for o in good], '.3f')}{source('null_med')}")

    print("\n  Measured here, continuum offsets removed with a running median")
    print(f"    null scatter / formal error   {span([o.null_scatter for o in good], '.2f')}    (1 = errors correct)")
    print(f"    Stokes scatter / formal error {span([o.stokes_scatter for o in good], '.2f')}")
    print(f"    continuum offset, null [%]    {span([100*o.null_med for o in good], '.3f')}")
    print(f"    continuum offset, Stokes [%]  {span([100*o.stokes_med for o in good], '.3f')}")

    print(f"\n  Polarisation signal, Stokes {spec.stokes}")
    print(f"    rms of Stokes / null [%]      {100*med([o.stokes_rms for o in good]):.3f} / "
          f"{100*med([o.null_rms for o in good]):.3f}")
    print(f"    rms above the null [%]        {span([100*o.excess for o in good], '.3f')}")
    print(f"    99.5th percentile |Stokes| [%] {span([100*o.peak for o in good], '.3f')}")
    ratio = med([o.stokes_rms / o.null_rms for o in good])
    n_sig = sum(o.stokes_rms > SIGNAL_RATIO * o.null_rms for o in good)
    print(f"    rms ratio Stokes / null       {ratio:.2f}, above {SIGNAL_RATIO} in {n_sig} of {len(good)} orders")

    ss, ss_err, ss_sig = fit_shift(good, "stokes")
    ns, ns_err, ns_sig = fit_shift(good, "null")
    print("\n  Part that follows dI/dlambda, as a velocity shift (median over orders)")
    print(f"    Stokes                        {ss:+.1f} +- {ss_err:.1f} m/s  ({abs(ss_sig):.1f} sigma)"
          + ("   sign tracks the sign of <Bz>; neither sign convention nor size is calibrated" if spec.stokes == "V" else ""))
    print(f"    null                          {ns:+.1f} +- {ns_err:.1f} m/s  ({abs(ns_sig):.1f} sigma, expect 0)")

    flags = []
    snr = med([o.qc["snr"] for o in good])
    if np.isfinite(snr) and snr < MIN_SNR:
        flags.append(f"low SNR ({snr:.0f})")
    null_scatter = med([o.null_scatter for o in good])
    if not np.isfinite(null_scatter):
        flags.append("no formal errors for Stokes and null (NaN, fixed in later cr2res), so no scatter check")
    elif null_scatter > MAX_NULL_SCATTER:
        flags.append(f"null scatter {null_scatter:.1f} x formal error")
    worst_null = max(abs(o.null_med) for o in good)
    if worst_null > MAX_NULL_OFFSET:
        flags.append(f"null continuum offset up to {100*worst_null:.2f} %")
    ratios = [o.qc["beamratio"] for o in good if o.qc["beamratio"] is not None]
    if ratios and (min(ratios) < BEAMRATIO_RANGE[0] or max(ratios) > BEAMRATIO_RANGE[1]):
        flags.append(f"beam ratio {min(ratios):.2f} .. {max(ratios):.2f}")
    if np.isfinite(ns_sig) and abs(ns_sig) > 5:
        flags.append(f"null follows dI/dlambda at {abs(ns_sig):.0f} sigma, {ns:+.0f} m/s (beams shifted in wavelength?)")
    if len(good) < len(spec.orders):
        flags.append(f"{len(spec.orders) - len(good)} orders unusable")

    if n_sig >= max(3, len(good) // 3):
        signal = f"Stokes {spec.stokes} signal in {n_sig} of {len(good)} orders"
    elif n_sig:
        signal = f"weak Stokes {spec.stokes} signal, in {n_sig} of {len(good)} orders only"
    else:
        signal = f"Stokes {spec.stokes} not above the null"
    print(f"\n  VERDICT: {signal}; " + ("data look fine" if not flags else "CHECK: " + "; ".join(flags)))


def correlate(a, b):
    """Correlation of the nod A and nod B spectra, order by order."""
    by_key = {(o.chip, o.idx): o for o in b.good}
    result = {"stokes": [], "null": []}
    for oa in a.good:
        ob = by_key.get((oa.chip, oa.idx))
        if ob is None:
            continue
        for which in result:
            xa = getattr(oa, which + "_hp")
            xb = np.interp(oa.wl, ob.wl, np.where(ob.ok, getattr(ob, which + "_hp"), 0.0))
            ok = oa.ok
            if ok.sum() > 200:
                result[which].append(np.corrcoef(xa[ok], xb[ok])[0, 1])
    return result


def report_pair(a, b):
    result = correlate(a, b)
    print("=" * 100)
    print(f"Nod A against nod B: {os.path.dirname(a.path) or '.'}  Stokes {a.stokes}")
    lines = {}
    for which, r in result.items():
        r = np.array(r)
        if r.size < 2:
            print(f"  {which}: too few common orders")
            return
        mean, err = r.mean(), r.std(ddof=1) / np.sqrt(r.size)
        lines[which] = (mean, err)
        print(f"  {which:6s} correlation  {mean:+.3f} +- {err:.3f}  ({mean/err:+.1f} sigma, {r.size} orders, "
              f"best order {r.max():+.2f})")
    s, n = lines["stokes"], lines["null"]
    if s[0] / s[1] > 5:
        text = "the Stokes signal repeats in A and B, so it is real"
    elif s[0] / s[1] > 3:
        text = "marginal common Stokes signal"
    else:
        text = "no common Stokes signal in A and B"
    if abs(n[0] / n[1]) > 3:
        text += "; the NULL also repeats, which means a systematic effect or signal leaking into the null"
    print(f"  VERDICT: {text}")


def plot(spec):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)
    for o in spec.good:
        ok = o.ok
        axes[0].plot(o.wl[ok], o.intens[ok] / np.percentile(o.intens[ok], 98), lw=0.5)
        axes[1].plot(o.wl[ok], 100 * o.stokes_hp[ok], lw=0.5)
        axes[2].plot(o.wl[ok], 100 * o.null_hp[ok], lw=0.5)
    lim = 100 * max(max(o.peak for o in spec.good) * 1.5, 1e-3)
    axes[0].set_ylabel("I / I(98%)")
    axes[1].set_ylabel(f"{spec.stokes}/I [%]")
    axes[2].set_ylabel("null/I [%]")
    for ax in axes[1:]:
        ax.set_ylim(-lim, lim)
        ax.axhline(0, color="k", lw=0.3)
    axes[2].set_xlabel("wavelength [nm]")
    axes[0].set_title(spec.title(), fontsize=9)
    out = os.path.splitext(spec.path)[0] + "_check.png"
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"  plot: {out}")


def plot_trend(specs, out):
    """Time series of the per-file medians over orders of the QC parameters."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from astropy.time import Time

    panels = (
        ("SNR per pixel", lambda s: med([o.qc["snr"] for o in s.good]), None),
        ("beam ratio up/down", lambda s: med([o.qc["beamratio"] for o in s.good]), 1.0),
        ("null scatter / error", lambda s: med([o.qc["null_scatter"] for o in s.good]), 1.0),
        ("null median [%]", lambda s: 100 * med([o.qc["null_med"] for o in s.good]), 0.0),
        ("rms Stokes / null", lambda s: med([o.stokes_rms / o.null_rms for o in s.good]), 1.0),
    )
    colours = {"circ YJ": "tab:blue", "circ HK": "tab:red", "lin YJ": "tab:green", "lin HK": "tab:orange"}
    markers = {"A": "o", "B": "s"}

    points = []
    for spec in specs:
        date = spec.header.get("DATE-OBS")
        if not spec.good or not date:
            continue
        try:
            points.append((Time(date, format="isot", scale="utc").to_datetime(), spec))
        except ValueError:
            print(f"{spec.path}: cannot read DATE-OBS '{date}', left out of the trend plot", file=sys.stderr)
    if not points:
        print("trend plot: no usable files", file=sys.stderr)
        return
    points.sort(key=lambda p: p[0])

    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 2.3 * len(panels)), sharex=True)
    for ax, (label, value, reference) in zip(axes, panels):
        if reference is not None:
            ax.axhline(reference, color="0.6", lw=0.8)
        for when, spec in points:
            ax.plot(when, value(spec), ls="none", marker=markers.get(spec.nod, "x"), ms=6,
                    color=colours.get(spec.polarimeter, "k"), mfc="none" if spec.nod == "B" else None)
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
    handles = [plt.Line2D([], [], ls="none", marker="o", color=c, label=k)
               for k, c in colours.items() if any(s.polarimeter == k for _, s in points)]
    handles += [plt.Line2D([], [], ls="none", marker=m, color="k", mfc="none" if k == "B" else None, label=f"nod {k}")
                for k, m in markers.items() if any(s.nod == k for _, s in points)]
    axes[0].legend(handles=handles, ncol=len(handles), fontsize=8, loc="best")
    axes[0].set_title(f"cr2res_obs_pol QC, median over orders, {len(points)} files", fontsize=10)
    locator = mdates.AutoDateLocator()
    axes[-1].xaxis.set_major_locator(locator)
    axes[-1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    axes[-1].set_xlabel("DATE-OBS (UTC)")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print(f"trend plot: {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="*", help="cr2res_obs_pol_pol_spec[AB].fits")
    parser.add_argument("-v", "--verbose", action="store_true", help="also print a line per order")
    parser.add_argument("--plot", action="store_true", help="save <file>_check.png next to each input")
    parser.add_argument("--trend", nargs="?", const=TREND_DEFAULT, metavar="PNG",
                        help="plot the per-file QC medians against time [cr2res_pol_check_trend.png]")
    parser.add_argument("--recompute", action="store_true",
                        help="ignore the QC parameters in the headers and compute them here")
    args = parser.parse_args()
    # "--trend a.fits b.fits" hands the first input to --trend as its file name
    if args.trend and re.search(r"\.fits?(\.\w+)?$", args.trend, re.I):
        args.files.insert(0, args.trend)
        args.trend = TREND_DEFAULT

    if not args.files:
        parser.error("no input files")

    specs = []
    for path in args.files:
        try:
            spec = PolSpec(path, args.recompute)
        except (OSError, KeyError, ValueError) as exc:
            print(f"{path}: cannot read as a pol_spec file: {exc}", file=sys.stderr)
            continue
        specs.append(spec)
        report(spec, args.verbose)
        if args.plot and spec.good:
            plot(spec)

    if args.trend:
        plot_trend(specs, args.trend)

    pairs = {}
    for spec in specs:
        pairs.setdefault(spec.key, {})[spec.nod] = spec
    for pair in pairs.values():
        if "A" in pair and "B" in pair:
            report_pair(pair["A"], pair["B"])
    return 0 if specs else 1


if __name__ == "__main__":
    sys.exit(main())
