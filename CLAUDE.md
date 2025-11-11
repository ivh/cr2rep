# CR2RES Pipeline - Quick Setup Guide for Claude

This guide provides the fastest path to get the CR2RES pipeline compiled and running with esorex.

## Prerequisites

This repository requires dependencies that are bundled in ESO's official distribution kit. The git repository alone is incomplete.

## Quick Start (Copy-Paste for Claude)

**For automated installation, run these commands in sequence:**

```bash
# Navigate to repository
cd /home/user/cr2rep

# Step 1: Download kit (90 seconds, 83MB)
cd /tmp && wget -q https://ftp.eso.org/pub/dfs/pipelines/instruments/cr2res/cr2re-kit-1.6.10.tar.gz && tar -xzf cr2re-kit-1.6.10.tar.gz

# Step 2: Install system dependencies (30 seconds)
apt-get update && apt-get install -y build-essential autoconf automake libtool libcpl-dev libgsl-dev liberfa-dev libcurl4-openssl-dev

# Step 3: Extract and copy dependencies (5 seconds)
cd /tmp/cr2re-kit-1.6.10 && tar -xzf cr2re-1.6.10.tar.gz && \
cd /home/user/cr2rep && \
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/{hdrl,irplib,regtests} ./ && \
cp /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/m4macros/*.m4 ./m4macros/

# Step 4: Build CR2RES pipeline (3-5 minutes)
cd /home/user/cr2rep && ./autogen.sh && ./configure && make -j$(nproc) && make install

# Step 5: Build and install esorex (2 minutes)
cd /tmp/cr2re-kit-1.6.10 && tar -xzf esorex-3.13.10.tar.gz && \
cd esorex-3.13.10 && ./configure --prefix=/opt/cr2res && make && make install

# Step 6: Configure esorex
ln -sf /opt/cr2res/bin/esorex /usr/local/bin/esorex && \
sed -i 's|^esorex.caller.recipe-dir=.*|esorex.caller.recipe-dir=/usr/local/lib/esopipes-plugins/cr2re-1.6.11|' /opt/cr2res/etc/esorex.rc

# Step 7: Verify (should show version 3.13.10 and 21 recipes)
esorex --version && echo "---" && esorex --recipes | grep -c cr2res
```

**Total time: ~10 minutes**

## Detailed Setup (Step-by-Step)

### 1. Download ESO CR2RES Kit (~90 seconds)

```bash
cd /tmp
wget https://ftp.eso.org/pub/dfs/pipelines/instruments/cr2res/cr2re-kit-1.6.10.tar.gz
tar -xzf cr2re-kit-1.6.10.tar.gz
cd cr2re-kit-1.6.10
```

**Success indicator**: Directory `/tmp/cr2re-kit-1.6.10` contains 16 files including `cr2re-1.6.10.tar.gz` and `esorex-3.13.10.tar.gz`

### 2. Install System Dependencies (~30 seconds)

```bash
apt-get update
apt-get install -y build-essential autoconf automake libtool \
    libcpl-dev libgsl-dev liberfa-dev libcurl4-openssl-dev
```

**Success indicator**: All packages show "Setting up..." messages without errors

### 3. Extract Missing Dependencies from Kit (~5 seconds)

The git repository is missing `hdrl/`, `irplib/`, `regtests/`, and m4 macros. Extract them:

```bash
# Extract CR2RES source from kit (if still in /tmp/cr2re-kit-1.6.10)
tar -xzf cr2re-1.6.10.tar.gz

# Copy missing dependencies to git repository
cd /home/user/cr2rep
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/hdrl ./
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/irplib ./
cp -r /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/regtests ./
cp /tmp/cr2re-kit-1.6.10/cr2re-1.6.10/m4macros/*.m4 ./m4macros/
```

**Success indicator**:
- `ls hdrl irplib regtests` shows all three directories
- `ls m4macros/*.m4` shows `erfa_pkg.m4`, `gsl.m4`, `hdrl.m4`, `libcurl_pkg.m4`

### 4. Build CR2RES Pipeline (~3-5 minutes)

```bash
cd /home/user/cr2rep
./autogen.sh
./configure
make -j$(nproc)
make install
```

**Success indicators**:
- `autogen.sh` outputs: "Don't forget to run ./configure"
- `configure` ends with: "config.status: executing libtool commands"
- `make` completes without errors (warnings are OK)
- `make install` installs to `/usr/local/lib/esopipes-plugins/cr2re-1.6.11/`

**Note**: `autogen.sh` may take 30-60 seconds with no output - this is normal.

### 5. Install and Configure EsoRex (~2 minutes)

```bash
cd /tmp/cr2re-kit-1.6.10

# Extract and build esorex
tar -xzf esorex-3.13.10.tar.gz
cd esorex-3.13.10
./configure --prefix=/opt/cr2res
make
make install

# Create symlink for easy access
ln -sf /opt/cr2res/bin/esorex /usr/local/bin/esorex

# Configure recipe directory
sed -i 's|^esorex.caller.recipe-dir=.*|esorex.caller.recipe-dir=/usr/local/lib/esopipes-plugins/cr2re-1.6.11|' \
    /opt/cr2res/etc/esorex.rc
```

**Success indicators**:
- `configure` finds CPL libraries and creates config files
- `make` and `make install` complete without errors
- `/usr/local/bin/esorex` symlink exists

**Note**: `sudo` removed - not needed if running as root

### 6. Verify Installation

```bash
# Check esorex version (should show 3.13.10)
esorex --version

# List available recipes (should show 21 cr2res_* recipes)
esorex --recipes

# Count recipes
esorex --recipes | grep -c cr2res

# Get help for a specific recipe
esorex --help cr2res_cal_flat
```

**Expected output**:
- Version: `ESO Recipe Execution Tool, version 3.13.10`
- Libraries: `CPL = 7.3.2, CFITSIO = 4.3.1`
- Recipe count: `21`
- All recipes start with `cr2res_`

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

## Troubleshooting

### Build Issues

**Issue**: `configure: error: No ERFA/GSL/LIBCURL available`
- **Cause**: Missing system libraries
- **Fix**: `apt-get install liberfa-dev libgsl-dev libcurl4-openssl-dev`

**Issue**: `make: *** No rule to make target '../irplib/libirplib.la'`
- **Cause**: Missing dependencies from kit
- **Fix**: Extract and copy `hdrl/`, `irplib/`, and `regtests/` from kit (see Step 3)
- **Verify**: `ls hdrl irplib regtests` should show all three directories

**Issue**: `./configure: line XXXX: syntax error near unexpected token 'hdrl'`
- **Cause**: Missing m4 macros
- **Fix**: Copy all `.m4` files from kit's m4macros directory
- **Verify**: `ls m4macros/*.m4` should show 4+ files

**Issue**: `autogen.sh` appears to hang
- **Cause**: Normal behavior - autotools can take 30-60 seconds
- **Fix**: Wait patiently. If >5 minutes, Ctrl+C and restart

**Issue**: `./bootstrap` fails with timeout
- **Cause**: Equivalent to `autogen.sh`, can timeout in some environments
- **Fix**: Use `./autogen.sh` instead (it's what bootstrap calls anyway)

### Runtime Issues

**Issue**: `esorex --recipes` shows no cr2res recipes
- **Cause**: Recipe directory not configured
- **Fix**: Update `/opt/cr2res/etc/esorex.rc`:
  ```bash
  sed -i 's|^esorex.caller.recipe-dir=.*|esorex.caller.recipe-dir=/usr/local/lib/esopipes-plugins/cr2re-1.6.11|' \
      /opt/cr2res/etc/esorex.rc
  ```
- **Alternative**: Set environment variable: `export ESOREX_RECIPES=/usr/local/lib/esopipes-plugins/cr2re-1.6.11`

**Issue**: `esorex: command not found`
- **Cause**: Symlink not created or not in PATH
- **Fix**: `ln -sf /opt/cr2res/bin/esorex /usr/local/bin/esorex`
- **Verify**: `which esorex` should show `/usr/local/bin/esorex`

**Issue**: Recipe count is less than 21
- **Cause**: Incomplete installation
- **Fix**: Check `/usr/local/lib/esopipes-plugins/cr2re-1.6.11/` for `.so` files
- **Verify**: `ls /usr/local/lib/esopipes-plugins/cr2re-1.6.11/*.so | wc -l` should show 21

## Clean Restart

If installation fails and you need to start over:

```bash
# Remove built artifacts from repository
cd /home/user/cr2rep
make distclean 2>/dev/null || true
rm -rf hdrl irplib regtests
git checkout m4macros/  # Restore to original state

# Remove installed files
rm -rf /usr/local/lib/esopipes-plugins/cr2re-1.6.11
rm -rf /usr/local/lib/cr2re-1.6.11
rm -f /usr/local/bin/esorex
rm -rf /opt/cr2res

# Clear temporary files
rm -rf /tmp/cr2re-kit-1.6.10

# Now start from Step 1
```

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

## Installation Verification Checklist

After installation, verify everything is working:

```bash
# 1. Check esorex is accessible
which esorex
# Expected: /usr/local/bin/esorex

# 2. Verify version
esorex --version | head -1
# Expected: ESO Recipe Execution Tool, version 3.13.10

# 3. Count recipes
esorex --recipes | grep cr2res | wc -l
# Expected: 21

# 4. List all CR2RES recipes
esorex --recipes | grep cr2res
# Expected: List of 21 recipes starting with cr2res_

# 5. Check library files exist
ls -lh /usr/local/lib/esopipes-plugins/cr2re-1.6.11/*.so | wc -l
# Expected: 21

# 6. Check main library
ls -lh /usr/local/lib/cr2re-1.6.11/libcr2res.so
# Expected: ~3MB shared library

# 7. Verify config file
grep recipe-dir /opt/cr2res/etc/esorex.rc
# Expected: esorex.caller.recipe-dir=/usr/local/lib/esopipes-plugins/cr2re-1.6.11

# 8. Test recipe help
esorex --help cr2res_cal_flat | head -5
# Expected: Recipe description and usage
```

**All checks passed?** Installation complete! ✓

## Tips for Claude

1. **Use the Quick Start**: The command sequence at the top can be copied and executed step-by-step
2. **Watch for timing**: If a step takes >2x the estimated time, something may be wrong
3. **Check success indicators**: After each step, verify the expected output
4. **Don't rush autogen.sh**: It's normal for it to appear "stuck" for 30-60 seconds
5. **Use the clean restart**: If anything fails midway, use the clean restart procedure
6. **Parallel builds**: Use `make -j$(nproc)` for faster compilation
7. **No sudo needed**: If running as root (typical in containers), omit sudo from all commands

## Links

- **ESO CR2RES Page**: https://www.eso.org/sci/software/pipelines/cr2res/
- **Kit Download**: https://ftp.eso.org/pub/dfs/pipelines/instruments/cr2res/
- **Pipeline Manual**: Included in kit as `cr2re-pipeline-manual-1.6.10.pdf`
- **CPL Documentation**: https://www.eso.org/sci/software/cpl/

---

**Last Updated**: Based on successful installation on 2025-11-11
**Tested On**: Ubuntu 24.04 (Noble) with CPL 7.3.2
**Total Install Time**: ~10 minutes on modern hardware
