# EsoRex Installation and Configuration

## Installation Summary

EsoRex (ESO Recipe Execution Tool) version 3.13.10 has been successfully installed from the cr2re-kit-1.6.10.

### Installed Components

- **EsoRex binary**: `/opt/cr2res/bin/esorex`
- **Symlink**: `/usr/local/bin/esorex` (for easy access)
- **Configuration**: `/opt/cr2res/etc/esorex.rc`
- **CR2RES recipes**: `/usr/local/lib/esopipes-plugins/cr2re-1.6.11/`

### CR2RES Recipes Available (21 total)

Calibration recipes:
- `cr2res_cal_dark` - Dark recipe
- `cr2res_cal_flat` - Flat recipe
- `cr2res_cal_detlin` - Detector Linearity recipe
- `cr2res_cal_wave` - Wavelength Calibration

Observation recipes:
- `cr2res_obs_nodding` - Nodding Observation recipe
- `cr2res_obs_2d` - 2D Observation recipe
- `cr2res_obs_staring` - Staring Observation recipe
- `cr2res_obs_pol` - Polarimetry Observation recipe

Utility recipes:
- `cr2res_util_bpm_merge` - BPM merging utility
- `cr2res_util_bpm_split` - BPM splitting utility
- `cr2res_util_calib` - Calibration utility
- `cr2res_util_extract` - Optimal Extraction utility
- `cr2res_util_genlines` - Generate spectrum calibration FITS tables
- `cr2res_util_genstd` - Generate standard star FITS tables
- `cr2res_util_normflat` - Flat Normalization utility
- `cr2res_util_plot` - Plotting utility
- `cr2res_util_slit_curv` - Slit Curvature utility
- `cr2res_util_splice` - Splicing utility
- `cr2res_util_trace` - Trace utility
- `cr2res_util_trace_map` - TRACE_WAVE maps creation
- `cr2res_util_wave` - Wavelength Calibration

## Usage

### List all available recipes:
```bash
esorex --recipes
```

### Get help for a specific recipe:
```bash
esorex --help cr2res_cal_flat
```

### Run a recipe:
```bash
esorex cr2res_cal_flat sof_file.sof
```

### Get detailed manual page:
```bash
esorex --man-page cr2res_cal_flat
```

## Libraries Used

- CPL = 7.3.2
- CFITSIO = 4.3.1
- WCSLIB
- FFTW (normal precision) = 3.3.10-sse2
- FFTW (single precision) = 3.3.10-sse2
- OPENMP = 201511

## Configuration

The esorex configuration file is located at `/opt/cr2res/etc/esorex.rc` and has been configured to automatically find CR2RES recipes at:
```
esorex.caller.recipe-dir=/usr/local/lib/esopipes-plugins/cr2re-1.6.11
```

This means you don't need to specify `--recipe-dir` when running esorex.
