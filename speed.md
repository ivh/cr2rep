# Speed branch: faster slit-decomposition extraction

Status notes for the port of the fast extraction algorithm from CharSlit
into cr2rep. Last updated 2026-06-11.

## Status

- Branch `speed`, commit `632dd74` (based on `26f0671`, v1.6.12 paranal release).
- The algorithm from `~/CharSlit.git` (branch `speed`, commit `4681cbe`) is
  plugged into `cr2res/cr2res_extract.c`, replacing the old
  `cr2res_extract_slit_func_curved` / `cr2res_extract_xi_zeta_tensors`.
  Net -276 lines.
- Builds clean with `-Wall -Wextra`; all 9 unit test suites pass
  (incl. `cr2res_extract-test`, no memory leaks per CPL diagnostics).
- Benchmarked on real data (2026-06-11): `cr2res_obs_nodding` on a
  CD-33_7795 M4368 combined 200-frame SOF (extract_swath_width=2048,
  extract_height=45, extract_oversample=10, no flat),
  old (26f0671) vs new. Wall: 106 s -> 51 s; user CPU: 78 s -> 25 s
  (frame loading/combination overhead is shared, so the extraction-only
  speedup is larger than the 2x wall ratio).
- Spectra agree within numerical errors (mostly < 0.01 sigma, many orders
  bit-identical) after restoring the sum-of-|sL| normalization (see below).
  Exceptions: order 05 (4324-4410 nm, inside the opaque CO2 band, pure
  noise / negative flux -> ill-conditioned decomposition differs) and a
  few isolated pixels where kappa-clipping decisions flip.

## What was ported (where the 3-5x comes from)

- `cr2res_extract_zeta_tensors`: builds only the zeta tensor. The xi tensor
  (subpixel -> pixel mapping) and its 4-corner (LL/LR/UL/UR) bookkeeping are
  gone; the zeta insertions were identical for all corner cases, factored
  into `cr2res_extract_zeta_add`.
- New memory layout: zeta entries of one detector pixel are contiguous
  (`zeta_index`/`mzeta_index` macros changed accordingly).
- Pixel-centric SLE fills: both band matrices are sums over detector pixels
  of all pairs of zeta entries of that pixel. Entries sharing the same
  subpixel (sL system) or column (sP system) are merged first into small
  scratch buffers `zw`/`zk` (size 3*(osample+1)); masked pixels are skipped
  entirely. Matrices are symmetric: only upper bands are accumulated, then
  mirrored.
- Band matrices are now row-major (band entries of one row contiguous),
  solved by a new static `cr2res_extract_bandsol_rowmajor`.
- Convergence: reduced chi-square based, stops when cost improvement
  < ftol=1e-7; maxiter is an unconditional bound (also for non-finite cost).
  Replaces the old sP_change/sP_stop criterion (the 5e-5 literal at the
  call site is gone, `sP_old` no longer needed).
- Diagonal regularization (max_diag * 1e-10 floor) on both matrices
  prevents singular systems from fully masked rows/columns.

## Interface adaptations (cr2res vs CharSlit)

- Public API unchanged: `cr2res_extract_slitdec_curved` signature identical,
  so recipes and tests needed no changes.
- Slit curvature: caller flattens the trace-table A/B/C polynomials into
  `double slitcurve_sw[swath*6]` (CharSlit layout, degrees up to 5; cr2res
  fills only c1, c2). Local-frame shift done analytically:
  `c0 = 0, c1 = B(x) + 2*yc*C(x), c2 = C(x)`. No more `cpl_polynomial`
  objects / eval calls in the hot path.
- Kept caller-side: workspace pre-allocation reused across swaths
  (l_Aij/p_Aij/l_bj/p_bj/zeta/m_zeta, plus new zw/zk), the conservative
  whole-order `delta_x` computation, `int *mask` (CharSlit uses uchar).
- CharSlit's `slitdeltas` input dropped (cr2res has no such data; if ever
  added, it must also enter the delta_x bound - see CharSlit CLAUDE.md).
- Kept cr2res-specific blocks the CharSlit code lacks: `pclip` pre-clipping,
  extraction with fixed input slit function (`slit_func_in` skips the sL
  solve), sign-flip on negative convergence, both `error_factor` uncertainty
  modes (Horne 1986 for -1, flux-based otherwise). CharSlit's own unc
  estimate and `info[5]` output were not taken.
- CharSlit's edge zeroing of sP/unc within delta_x of the swath borders was
  NOT ported; the caller's overlap weights already discard those columns
  (matches old cr2res behavior).

## Gotchas / things learned

- `cr2res_extract_slitdec_bandsol` (public, column-major) MUST stay: used by
  `cr2res_utils.c` (polynomial fitting etc.). The new row-major bandsol is a
  separate static function.
- `kappa` outlier rejection is now gated by `kappa > 0` (CharSlit behavior);
  old code always applied it. Identical for the usual kappa>0.
- Slit function normalization MUST stay sum of |sL| (old cr2res), not the
  plain sum (CharSlit). The model sP*sL is invariant either way, but the
  flux scale of sP changes by sum/sum|sL| when sL has negative parts —
  which it does in nodding A-B images (background residuals). Observed up
  to 23% per-order flux offsets with the plain sum before reverting
  (commit after 632dd74). Where sL >= 0 both conventions are identical.
- Removed along with the old code: `debug_output()` dump-on-failure helper,
  `img_mad` debug image, `debug_img_mad_*.fits` output.

## Next steps

1. Compare QC parameters old vs new on a reference dataset (obs_nodding
   benchmark done, see Status; cal_flat not yet).
2. Consider syncing future CharSlit improvements; the algorithm core was
   kept structurally close to CharSlit's `slitdec.c` to ease diffing.
   NOTE: do not sync CharSlit's plain-sum sL normalization (see above).
