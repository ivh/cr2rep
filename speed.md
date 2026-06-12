# Speed branch: faster slit-decomposition extraction

Status notes for the port of the fast extraction algorithm from CharSlit
into cr2rep. Last updated 2026-06-12.

## Round 2 (2026-06-12): dense-window fills, in-order model merge, OpenMP

A second optimization pass on top of the CharSlit port. All changes keep
the products BIT-IDENTICAL to the round-1 code (verified on a synthetic
full-frame benchmark: spectra, errors, slit functions, models and the
extract_traces tables compared with cmp over several configurations,
incl. curvature, bad pixels, cosmics, pclip, both error_factor modes and
the input-slit-function path; valgrind clean; all 9 unit suites pass).

What changed, in decreasing order of measured impact:

1. Model computation loop in cr2res_extract_slit_func_curved walks x
   outermost so the zeta tensor (the largest array, ~16 B/entry) is read
   sequentially; it was striding ~23 kB per step and was ~30% of the
   runtime. Per-pixel sums accumulate in the same order as before.
2. The SLE fill loops no longer search a list of unique keys to merge a
   pixel's zeta entries: the zeta build now records per pixel the key
   ranges (zeta_rng: min/max of iy and of x), and the merge scatters into
   a dense window zw[key - min]. The iy range of one pixel is at most
   2*osample (the band width the matrix assumes anyway; checked, with the
   old search path kept as fallback), the x range at most 2*delta_x (by
   construction). The band accumulation walks the window row-wise so the
   inner loop is contiguous (vectorizable) in both operands. Bit-exact
   because per-slot merge order is unchanged and each band element gets
   at most one contribution per pixel, so pair enumeration order cannot
   matter; window gaps contribute exact zeros.
3. The per-trace model merge in cr2res_extract_traces used
   hdrl_image_get/set_pixel per pixel (4.2M hdrl calls per trace ~ as
   expensive as the decomposition itself); replaced by direct buffer
   access with identical semantics (write data+error and accept the
   pixel where the trace model is non-zero and not rejected).
4. OpenMP: the trace loop in cr2res_extract_traces runs in parallel
   (schedule(dynamic,1)); each trace is computed entirely within one
   thread, results are stored per trace and merged sequentially in trace
   order, so products do not depend on thread count or scheduling.
   configure.ac gained the same ESO_ENABLE_OPENMP([yes]) macro hdrl uses
   (--disable-openmp to turn off). CPL's error state is omp-threadprivate
   and the system/ESO CPL builds link libgomp; like hdrl's own parallel
   regions this assumes the default CPL memory mode (no thread-unsafe
   xmemory tracking, i.e. don't combine with CPL_MEMORY_MODE=1/2).
   The two fill inner loops carry "#pragma omp simd" (inert without
   -fopenmp); element-wise independent, so no FP reordering.
5. cr2res_extract_zeta_add is static inline and maintains the zeta_rng
   ranges; the zeta build remains once per swath.

Measured on the synthetic benchmark (2048x2048, height 45, 4-px shear,
noise+cosmics+bad pixels, niter 30, kappa 10), best-of-N wall time of the
extraction call against round 1:

- slitdec_curved, swath 800 / osample 7 / pclip 0.1: 0.321 -> 0.212 s
- slitdec_curved, swath 2048 / osample 10 / Horne unc: 0.327 -> 0.214 s
- extract_traces, 4 orders, defaults: 1.127 -> 0.635 s (1 thread),
  0.264 s (4 threads). Real recipes (cal_flat, obs_nodding) extract
  6-9 orders per detector, so the parallel win applies there.
  (These wall times still include the per-call rectify/median-collapse
  setup, so the core-decomposition speedup is larger than the ratios.)

Profile after round 2 is dominated by the three unavoidable sequential
sweeps over zeta per iteration (two SLE fills + model) and the zeta build
itself - i.e. memory bandwidth on the tensor. Further gains would need
slimming the zeta entries (e.g. packed 12-B entries, ~25% traffic) or
swath-level parallelism; both judged not worth the complexity now.

## Status

- Branch `speed`, commit `632dd74` (based on `26f0671`, v1.6.12 paranal release).
- 2026-06-12: restored master's spectrum-change stop criterion (see
  Convergence below). With it, old vs new products match at machine
  precision across the board: SPEC relmax ~5-7e-14, ERR ~3-5e-11,
  SLIT_FUNC ~1e-11, WL identical — including the CO2 order 05 and the
  flipped clipping pixels listed below, which no longer differ at all.
  Speed unchanged: 55 s wall / 26 s user on the benchmark.
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
- (2026-06-11, before the stop-criterion restore:) spectra agreed within
  numerical errors (mostly < 0.01 sigma, many orders bit-identical) after
  restoring the sum-of-|sL| normalization (see below). Exceptions then:
  order 05 (4324-4410 nm, inside the opaque CO2 band, pure noise /
  negative flux -> ill-conditioned decomposition differed) and a few
  isolated pixels where kappa-clipping decisions flipped.

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
- Convergence: CharSlit's reduced-chi-square criterion was initially taken,
  but reverted (2026-06-12) to master's historic criterion: stop when the
  largest per-pixel spectrum change between iterations drops below
  sP_stop * |median(sP)| (sP_stop=5e-5 literal at the call site,
  caller-allocated `sP_old` buffer, as in master). The cost criterion left
  iteration-path differences (1e-4 ripple, divergent solutions in
  noise-only orders); the restored one reproduces master exactly. Do not
  sync CharSlit's cost-based criterion. One fix kept vs master: the
  no-convergence warning checks `iter > maxiter` (master's `== maxiter`
  could never fire).
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
   NOTE: do not sync CharSlit's plain-sum sL normalization or its
   cost-based convergence criterion (see above).
