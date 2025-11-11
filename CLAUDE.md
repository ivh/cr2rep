# CR2RES Pipeline - Quick Setup Guide for Claude

This guide provides the fastest path to get the CR2RES pipeline compiled and running with esorex.

## Prerequisites

This repository requires dependencies that are bundled in ESO's official distribution kit. The git repository alone is incomplete.

## Quick Setup (20 minutes)

### 1. Download ESO CR2RES Kit

```bash
cd /tmp
wget https://ftp.eso.org/pub/dfs/pipelines/instruments/cr2res/cr2re-kit-1.6.10.tar.gz
tar -xzf cr2re-kit-1.6.10.tar.gz
cd cr2re-kit-1.6.10
```

### 2. Install System Dependencies

```bash
apt-get update
apt-get install -y build-essential autoconf automake libtool \
    libcpl-dev libgsl-dev liberfa-dev libcurl4-openssl-dev
```

### 3. Extract Missing Dependencies from Kit

The git repository is missing `hdrl/`, `irplib/`, `regtests/`, and m4 macros. Extract them:

```bash
# Extract CR2RES source from kit
tar -xzf cr2re-1.6.10.tar.gz

# Copy missing dependencies to git repository
cd /path/to/cr2rep
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/hdrl ./
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/irplib ./
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/regtests ./
cp /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/m4macros/*.m4 ./m4macros/
```

### 4. Build CR2RES Pipeline

```bash
./autogen.sh
./configure
make -j$(nproc)
sudo make install
```

**Installation location**: `/usr/local/lib/esopipes-plugins/cr2re-1.6.11/`

### 5. Install and Configure EsoRex

```bash
cd /tmp/cr2re-kit-1.6.10

# Extract and build esorex
tar -xzf esorex-3.13.10.tar.gz
cd esorex-3.13.10
./configure --prefix=/opt/cr2res
make
sudo make install

# Create symlink for easy access
sudo ln -sf /opt/cr2res/bin/esorex /usr/local/bin/esorex

# Configure recipe directory
sudo sed -i 's|^esorex.caller.recipe-dir=.*|esorex.caller.recipe-dir=/usr/local/lib/esopipes-plugins/cr2re-1.6.11|' \
    /opt/cr2res/etc/esorex.rc
```

### 6. Verify Installation

```bash
# Check esorex version
esorex --version

# List available recipes (should show 21 cr2res_* recipes)
esorex --recipes

# Get help for a specific recipe
esorex --help cr2res_cal_flat
```

## What Gets Built

**Libraries:**
- `libcr2res.so` (3.0MB) - Main CR2RES library
- `libhdrl.la` - High-level Data Reduction Library
- `libirplib.la` - Instrument Reduction Pipeline Library

**Recipes (21 total):**
- 4 Calibration recipes (dark, flat, detlin, wave)
- 4 Observation recipes (nodding, 2d, staring, pol)
- 13 Utility recipes (extraction, tracing, calibration, etc.)

## Key Files Modified

When setting up from scratch, these files need the dependencies integrated:
- `hdrl/` - Complete HDRL library source
- `irplib/` - Complete irplib library source
- `regtests/` - Regression test framework
- `m4macros/hdrl.m4`, `erfa_pkg.m4`, `gsl.m4`, `libcurl_pkg.m4` - Build macros

## Configuration Files

- **EsoRex config**: `/opt/cr2res/etc/esorex.rc`
- **Recipe directory**: `/usr/local/lib/esopipes-plugins/cr2re-1.6.11/`
- **EsoRex binary**: `/usr/local/bin/esorex` (symlink to `/opt/cr2res/bin/esorex`)

## Common Issues

**Issue**: `configure: error: No ERFA/GSL/LIBCURL available`
- **Fix**: Install missing system library: `apt-get install liberfa-dev libgsl-dev libcurl4-openssl-dev`

**Issue**: `make: *** No rule to make target '../irplib/libirplib.la'`
- **Fix**: Missing dependencies. Extract `hdrl/`, `irplib/`, and `regtests/` from kit.

**Issue**: `./configure: line XXXX: syntax error near unexpected token 'hdrl'`
- **Fix**: Missing m4 macros. Copy all `.m4` files from kit's m4macros directory.

**Issue**: `esorex --recipes` shows no cr2res recipes
- **Fix**: Update `/opt/cr2res/etc/esorex.rc` with correct recipe-dir path, or use `--recipe-dir` flag.

## Environment Variables (Optional)

If you don't want to modify esorex.rc:

```bash
export ESOREX_RECIPES=/usr/local/lib/esopipes-plugins/cr2re-1.6.11
```

## Running Recipes

Basic usage:

```bash
# Create a SOF (Set of Frames) file
cat > example.sof << EOF
/path/to/raw_flat_1.fits FLAT
/path/to/raw_flat_2.fits FLAT
/path/to/dark.fits DARK
EOF

# Run the recipe
esorex cr2res_cal_flat example.sof
```

## Development Notes

- The git repository is a **development checkout** and requires ESO distribution kit dependencies
- For production use, ESO recommends using their complete kit with `./install_pipeline` script
- This setup allows development while having all necessary dependencies

## Quick Rebuild After Code Changes

```bash
make -j$(nproc)
sudo make install
# No need to reinstall esorex or reconfigure
```

## Links

- **ESO CR2RES Page**: https://www.eso.org/sci/software/pipelines/cr2res/
- **Kit Download**: https://ftp.eso.org/pub/dfs/pipelines/instruments/cr2res/
- **Pipeline Manual**: Included in kit as `cr2re-pipeline-manual-1.6.10.pdf`
