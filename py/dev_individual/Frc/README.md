# ERA5 Forcing Code Improvements (for future enhancement)

Author: shjo  
Last Updated: 2025-05-26  

---

## ✅ TODO List (unordered priority)

### 0. Time index based coordinates slice

### 1. Time Resampling Feature
- Purpose: Convert 3-hourly input to daily or other intervals
- Method: Use `xarray.resample()` or manual averaging
- Variable handling:
  - `tair`, `qair`, `Pair` → mean
  - `rain`, `srf`, `lrf` → sum

### 2. Flux Calculation Options
- Purpose: Support multiple physical formulas or methods
- Suggested configuration:
  - `dqdsst_method`: `default`, `bulk`, `none`
  - `qair_method`: `saturation_formula`, `empirical`
- Implement by branching with a `method` argument in each module

### 3. Dry-run Mode
- Purpose: Test structure and I/O paths without writing files
- Implementation:
  - Add `dry_run: true` to `config.yml`
  - Wrap file-writing functions in conditionals
- Example output: `[DRY-RUN] Would create ERA5_3hourly_qair.nc (shape: [2920, 150, 180])`

### 4. NetCDF Format Selection
- Purpose: Use NetCDF-4 during development, convert to NetCDF-3 for delivery
- Implementation:
  - `config.yml`: `format: NETCDF4` or `NETCDF3_CLASSIC`
  - Pass format to `Dataset(..., format=config['format'])`
  - Optionally post-process with `nccopy` to ensure compatibility

### 5. Input and Variable Validation
- Detect missing keys, wrong shapes, or unexpected types
- Use `warnings.warn()` or `assert` for soft or hard failure modes

### 6. Execution Logging
- Replace silent logic with informative `print()` messages
- Optional: replace with logging module later
- Example: `print("[INFO] Qair calculated. shape=(..., ...)")`

---

## Optional Enhancements

- Save intermediate outputs using `np.save()` or `pickle`
- Implement multiprocessing over the time dimension
- Monitor file sizes and automatically verify variable integrity

