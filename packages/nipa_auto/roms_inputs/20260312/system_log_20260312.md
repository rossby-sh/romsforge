[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] --- START: Download HYCOM frcst ---
[INFO] HYCOM FMRC url: https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_all/FMRC_ESPC-D-V02_all_best.ncd
[INFO] window: 2026-03-11 21:00:00 to 2026-03-13 03:00:00
[INFO] save step (hours): 3
[INFO] region lat: (5, 60) lon: (100, 176)
[INFO] output: /home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[INFO] opening dataset (decode_times=False)...
[INFO] lat/lon coords: {'lat': 'lat', 'lon': 'lon'}
[INFO] decoding all time-like coords (time/time1/...) ...
[INFO] resampling each variable to target_times and unifying time dim ...
  - surf_el: time1 -> time, nt=11
  - water_temp: time -> time, nt=11
  - salinity: time -> time, nt=11
  - water_u: time -> time, nt=11
  - water_v: time -> time, nt=11
[INFO] converting unified time coord to days since 2000-01-01 ...
[INFO] writing netcdf ...
Traceback (most recent call last):
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc.tmp',), 'a', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '75eafb36-491c-47a9-ad7a-34427233151f']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/get_hycom_frcst.py", line 319, in <module>
    raise SystemExit(main())
                     ^^^^^^
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/get_hycom_frcst.py", line 311, in main
    out.to_netcdf(tmp, encoding=enc, unlimited_dims=["time"])
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/core/dataset.py", line 2102, in to_netcdf
    return to_netcdf(  # type: ignore[return-value]  # mypy cannot resolve the overloads:(
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/api.py", line 2088, in to_netcdf
    store = get_writable_netcdf_store(
            ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/api.py", line 1898, in get_writable_netcdf_store
    return store_open(target, mode=mode, format=format, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/netCDF4_.py", line 468, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/netCDF4_.py", line 398, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/netCDF4_.py", line 477, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/netCDF4_.py", line 471, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/usr/local/miniconda3/envs/fetch/lib/python3.11/site-packages/xarray/backends/file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src/netCDF4/_netCDF4.pyx", line 2521, in netCDF4._netCDF4.Dataset.__init__
  File "src/netCDF4/_netCDF4.pyx", line 2158, in netCDF4._netCDF4._ensure_nc_success
PermissionError: [Errno 13] Permission denied: '/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc.tmp'
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/get_hycom_frcst.py` failed. (See above for error)
[ERR] --- FAILED: Download HYCOM frcst ---
[INFO] --- START: CMEMS download ---
[OK] created /home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK] raw_dir  = /home/shjo/github/romsforge/packages/nipa_auto/data/cmems_raw
[OK] cmems_dir= /home/shjo/github/romsforge/packages/nipa_auto/data/cmems
[OK] --- DONE: CMEMS download ---
[INFO] --- START: INI build ---
== Initial Condition Build (ini) ===============================================
[FAIL] [01] Load configuration and input metadata | dur=0.083s | FileNotFoundError: [Errno 2] No such file or directory: '/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc'
       at _netCDF4.pyx:2158 in netCDF4._netCDF4._ensure_nc_success
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_ini_single_new.py", line 31, in <module>
    ogcm = tl.load_ogcm_metadata(cfg.ogcm_name, cfg.ogcm_var_name)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/shjo/github/romsforge/libs/utils.py", line 227, in load_ogcm_metadata
    with Dataset(ogcm_name) as nc:
         ^^^^^^^^^^^^^^^^^^
  File "src/netCDF4/_netCDF4.pyx", line 2521, in netCDF4._netCDF4.Dataset.__init__
  File "src/netCDF4/_netCDF4.pyx", line 2158, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: '/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc'
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_ini_single_new.py` failed. (See above for error)
[ERR] --- FAILED: INI build ---
[INFO] --- START: INI bio add ---
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_ini3.py", line 84
    bio_yaml = getattr(cfg, "bio_yaml", bio_vars_file")
                                                     ^
SyntaxError: unterminated string literal (detected at line 84)
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_ini3.py` failed. (See above for error)
[ERR] --- FAILED: INI bio add ---
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
[FAIL] [01] Load configuration and input metadata | dur=0.036s | FileNotFoundError: [Errno 2] No such file or directory: '/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc'
       at _netCDF4.pyx:2158 in netCDF4._netCDF4._ensure_nc_success
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_bry_single_new2.py", line 34, in <module>
    ogcm = tl.load_ogcm_metadata(cfg.ogcm_name, cfg.ogcm_var_name)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/shjo/github/romsforge/libs/utils.py", line 227, in load_ogcm_metadata
    with Dataset(ogcm_name) as nc:
         ^^^^^^^^^^^^^^^^^^
  File "src/netCDF4/_netCDF4.pyx", line 2521, in netCDF4._netCDF4.Dataset.__init__
  File "src/netCDF4/_netCDF4.pyx", line 2158, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: '/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc'
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_bry_single_new2.py` failed. (See above for error)
[ERR] --- FAILED: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.035s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.020s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
[OK]   [03] Create initial NetCDF | dur=0.001s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.036s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.114s
[OK]   [08] Remap (weights) | dur=0.397s
[OK]   [09] Flood: horizontal | dur=0.771s
[OK]   [10] Flood: vertical | dur=0.615s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.486s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.061s
[OK]   [08] Remap (weights) | dur=0.388s
[OK]   [09] Flood: horizontal | dur=0.769s
[OK]   [10] Flood: vertical | dur=0.082s
[OK]   [11] Mask | dur=0.013s
[OK]   [12] z→σ & save | dur=0.074s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.036s
[OK]   [08] Remap (weights) | dur=0.396s
[OK]   [09] Flood: horizontal | dur=0.775s
[OK]   [10] Flood: vertical | dur=0.082s
[OK]   [11] Mask | dur=0.013s
[OK]   [12] z→σ & save | dur=0.072s
[WARN] skip write (not in bry): NO3_north
[WARN] skip write (not in bry): NO3_south
[WARN] skip write (not in bry): NO3_east
[WARN] skip write (not in bry): NO3_west
[WARN] skip write (not in bry): phyt_north
[WARN] skip write (not in bry): phyt_south
[WARN] skip write (not in bry): phyt_east
[WARN] skip write (not in bry): phyt_west
[WARN] skip write (not in bry): zoop_north
[WARN] skip write (not in bry): zoop_south
[WARN] skip write (not in bry): zoop_east
[WARN] skip write (not in bry): zoop_west
[WARN] skip write (not in bry): detritus_north
[WARN] skip write (not in bry): detritus_south
[WARN] skip write (not in bry): detritus_east
[WARN] skip write (not in bry): detritus_west
[OK]   [13] Write variables | dur=0.001s
== Summary =====================================================================
--- Time elapsed: 5.267s ---
[OK] --- DONE: BRY bio add ---
[INFO] PBS submitted
[INFO] JOBID=1170.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] pipeline_status exists → resume
[INFO] --- START: Download HYCOM frcst ---
[INFO] HYCOM FMRC url: https://tds.hycom.org/thredds/dodsC/FMRC_ESPC-D-V02_all/FMRC_ESPC-D-V02_all_best.ncd
[INFO] window: 2026-03-11 21:00:00 to 2026-03-13 03:00:00
[INFO] save step (hours): 3
[INFO] region lat: (5, 60) lon: (100, 176)
[INFO] output: /home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[INFO] opening dataset (decode_times=False)...
[INFO] lat/lon coords: {'lat': 'lat', 'lon': 'lon'}
[INFO] decoding all time-like coords (time/time1/...) ...
[INFO] resampling each variable to target_times and unifying time dim ...
  - surf_el: time1 -> time, nt=11
  - water_temp: time -> time, nt=11
  - salinity: time -> time, nt=11
  - water_u: time -> time, nt=11
  - water_v: time -> time, nt=11
[INFO] converting unified time coord to days since 2000-01-01 ...
[INFO] writing netcdf ...
[OK] saved: /home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[OK] --- DONE: Download HYCOM frcst ---
[SKIP] CMEMS download already done
[INFO] --- START: INI build ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/roms_grd_15km_nearest.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_15km_nrst.nc
· ogcm=/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.062s
[('/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc', 1, cftime.DatetimeGregorian(2026, 3, 12, 0, 0, 0, 0, has_year_zero=False), np.float64(9567.0))]
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Deactivate initiating biological variables ---
--- [+] Initial file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.315s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.002s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[OK]   [05] Load OGCM fields | dur=4.061s
[OK]   [06] Remap (weights) | dur=0.327s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=13.839s
[OK]   [08] Flood: vertical | dur=0.425s
[OK]   [09] Mask & clean | dur=0.009s
[OK]   [10] Rotate & stagger (u,v) | dur=0.113s
[OK]   [11] z→σ interpolation | dur=1.000s
[OK]   [12] Conserve volume & fix barotropic | dur=0.173s
[OK]   [13] Write variables | dur=0.036s
== Summary =====================================================================
Total elapsed: 20.368s
[OK] --- DONE: INI build ---
[INFO] --- START: INI bio add ---
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_ini3.py", line 84
    bio_yaml = getattr(cfg, "bio_yaml", bio_vars_file")
                                                     ^
SyntaxError: unterminated string literal (detected at line 84)
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_ini3.py` failed. (See above for error)
[ERR] --- FAILED: INI bio add ---
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
· src=OGCM grid=roms_grd_15km_nearest.nc
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/hycom_wght_15km_nrst.nc
[OK]   [01] Load configuration and input metadata | dur=0.035s
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Initiating biological variables: None type ---
--- [+] Boundary file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.063s
· Use existing wght file /.../fixed/hycom_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.011s
· Biological variables: npzd
· Flood method (boundary): edt
[OK]   [05] List & group OGCM files | dur=0.000s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.912s
[OK]   [08] Remap (weights) | dur=0.334s
[OK]   [09] Flood H/V | dur=1.301s
[OK]   [10] Mask & rotate | dur=0.227s
[OK]   [11] z→σ & save bry | dur=0.652s
[DONE] 2026-03-11 21 |  8.431s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.683s
[OK]   [08] Remap (weights) | dur=0.310s
[OK]   [09] Flood H/V | dur=0.748s
[OK]   [10] Mask & rotate | dur=0.198s
[OK]   [11] z→σ & save bry | dur=0.207s
[DONE] 2026-03-12 00 | 15.578s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.685s
[OK]   [08] Remap (weights) | dur=0.311s
[OK]   [09] Flood H/V | dur=0.704s
[OK]   [10] Mask & rotate | dur=0.157s
[OK]   [11] z→σ & save bry | dur=0.164s
[DONE] 2026-03-12 03 | 22.599s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.621s
[OK]   [08] Remap (weights) | dur=0.307s
[OK]   [09] Flood H/V | dur=0.701s
[OK]   [10] Mask & rotate | dur=0.137s
[OK]   [11] z→σ & save bry | dur=0.163s
[DONE] 2026-03-12 06 | 29.529s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.371s
[OK]   [08] Remap (weights) | dur=0.309s
[OK]   [09] Flood H/V | dur=0.698s
[OK]   [10] Mask & rotate | dur=0.129s
[OK]   [11] z→σ & save bry | dur=0.165s
[DONE] 2026-03-12 09 | 36.201s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.334s
[OK]   [08] Remap (weights) | dur=0.307s
[OK]   [09] Flood H/V | dur=0.705s
[OK]   [10] Mask & rotate | dur=0.121s
[OK]   [11] z→σ & save bry | dur=0.166s
[DONE] 2026-03-12 12 | 42.834s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.398s
[OK]   [08] Remap (weights) | dur=0.303s
[OK]   [09] Flood H/V | dur=0.702s
[OK]   [10] Mask & rotate | dur=0.127s
[OK]   [11] z→σ & save bry | dur=0.169s
[DONE] 2026-03-12 15 | 48.533s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.673s
[OK]   [08] Remap (weights) | dur=0.308s
[OK]   [09] Flood H/V | dur=0.703s
[OK]   [10] Mask & rotate | dur=0.121s
[OK]   [11] z→σ & save bry | dur=0.173s
[DONE] 2026-03-12 18 | 54.511s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.397s
[OK]   [08] Remap (weights) | dur=0.305s
[OK]   [09] Flood H/V | dur=0.695s
[OK]   [10] Mask & rotate | dur=0.126s
[OK]   [11] z→σ & save bry | dur=0.167s
[DONE] 2026-03-12 21 | 60.202s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.601s
[OK]   [08] Remap (weights) | dur=0.311s
[OK]   [09] Flood H/V | dur=0.703s
[OK]   [10] Mask & rotate | dur=0.122s
[OK]   [11] z→σ & save bry | dur=0.163s
[DONE] 2026-03-13 00 | 66.102s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.308s
[OK]   [08] Remap (weights) | dur=0.306s
[OK]   [09] Flood H/V | dur=0.703s
[OK]   [10] Mask & rotate | dur=0.127s
[OK]   [11] z→σ & save bry | dur=0.169s
[DONE] 2026-03-13 03 | 71.715s
[OK]   [06] Open source file | dur=71.715s
[OK]   [12] Write variables | dur=0.009s
== Summary =====================================================================
Total elapsed: 71.840s
[OK] --- DONE: BRY build ---
[SKIP] BRY bio add already done
[INFO] PBS submitted
[INFO] JOBID=1171.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[INFO] --- START: INI bio add ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/roms_grd_15km_nearest.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
· ogcm(daily)=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
· bio_model=npzd
[OK]   [01] Load configuration and input metadata | dur=0.136s
[OK]   [02] Time index matching & relative time | dur=0.006s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· ini_vars=['NO3', 'phytoplankton', 'zooplankton', 'detritus']
[OK]   [03b] Load bio vars (YAML) | dur=0.021s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.001s
!!! SHJO: FIX ME !!!
[OK]   [05] Load OGCM fields | dur=0.112s
[OK]   [06] Remap (weights) | dur=0.385s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=13.417s
[OK]   [08] Flood: vertical | dur=0.460s
[OK]   [10] Mask & rotate | dur=0.011s
[OK]   [11] z→σ interpolation | dur=1.119s
== Summary =====================================================================
Total elapsed: 15.673s
[WARN] skip write (not in ini): NO3
[WARN] skip write (not in ini): phytoplankton
[WARN] skip write (not in ini): zooplankton
[WARN] skip write (not in ini): detritus
[OK]   [12] Write variables | dur=0.002s
== Summary =====================================================================
[OK] --- DONE: INI bio add ---
[SKIP] BRY build already done
[SKIP] BRY bio add already done
[INFO] PBS submitted
[INFO] JOBID=1172.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
· src=OGCM grid=roms_grd_15km_nearest.nc
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/hycom_wght_15km_nrst.nc
[OK]   [01] Load configuration and input metadata | dur=0.036s
[OK]   [02] Time index matching & relative time | dur=0.006s
--- [NOTE] Initiating biological variables: None type ---
--- [+] Boundary file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.066s
· Use existing wght file /.../fixed/hycom_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.012s
· Biological variables: npzd
· Flood method (boundary): edt
[OK]   [05] List & group OGCM files | dur=0.000s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.110s
[OK]   [08] Remap (weights) | dur=0.307s
[OK]   [09] Flood H/V | dur=1.201s
[OK]   [10] Mask & rotate | dur=0.212s
[OK]   [11] z→σ & save bry | dur=0.650s
[DONE] 2026-03-11 21 |  6.484s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.828s
[OK]   [08] Remap (weights) | dur=0.304s
[OK]   [09] Flood H/V | dur=0.714s
[OK]   [10] Mask & rotate | dur=0.121s
[OK]   [11] z→σ & save bry | dur=0.166s
[DONE] 2026-03-12 00 | 12.619s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.992s
[OK]   [08] Remap (weights) | dur=0.311s
[OK]   [09] Flood H/V | dur=0.721s
[OK]   [10] Mask & rotate | dur=0.127s
[OK]   [11] z→σ & save bry | dur=0.167s
[DONE] 2026-03-12 03 | 18.937s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.474s
[OK]   [08] Remap (weights) | dur=0.307s
[OK]   [09] Flood H/V | dur=0.698s
[OK]   [10] Mask & rotate | dur=0.122s
[OK]   [11] z→σ & save bry | dur=0.173s
[DONE] 2026-03-12 06 | 24.711s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.357s
[OK]   [08] Remap (weights) | dur=0.308s
[OK]   [09] Flood H/V | dur=0.703s
[OK]   [10] Mask & rotate | dur=0.120s
[OK]   [11] z→σ & save bry | dur=0.164s
[DONE] 2026-03-12 09 | 30.364s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.473s
[OK]   [08] Remap (weights) | dur=0.307s
[OK]   [09] Flood H/V | dur=0.707s
[OK]   [10] Mask & rotate | dur=0.121s
[OK]   [11] z→σ & save bry | dur=0.177s
[DONE] 2026-03-12 12 | 36.149s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.399s
[OK]   [08] Remap (weights) | dur=0.309s
[OK]   [09] Flood H/V | dur=0.698s
[OK]   [10] Mask & rotate | dur=0.116s
[OK]   [11] z→σ & save bry | dur=0.168s
[DONE] 2026-03-12 15 | 41.839s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.349s
[OK]   [08] Remap (weights) | dur=0.304s
[OK]   [09] Flood H/V | dur=0.701s
[OK]   [10] Mask & rotate | dur=0.125s
[OK]   [11] z→σ & save bry | dur=0.169s
[DONE] 2026-03-12 18 | 47.488s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.485s
[OK]   [08] Remap (weights) | dur=0.303s
[OK]   [09] Flood H/V | dur=0.700s
[OK]   [10] Mask & rotate | dur=0.124s
[OK]   [11] z→σ & save bry | dur=0.165s
[DONE] 2026-03-12 21 | 53.266s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.398s
[OK]   [08] Remap (weights) | dur=0.307s
[OK]   [09] Flood H/V | dur=0.701s
[OK]   [10] Mask & rotate | dur=0.120s
[OK]   [11] z→σ & save bry | dur=0.159s
[DONE] 2026-03-13 00 | 58.951s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.287s
[OK]   [08] Remap (weights) | dur=0.303s
[OK]   [09] Flood H/V | dur=0.701s
[OK]   [10] Mask & rotate | dur=0.132s
[OK]   [11] z→σ & save bry | dur=0.172s
[DONE] 2026-03-13 03 | 64.546s
[OK]   [06] Open source file | dur=64.547s
[OK]   [12] Write variables | dur=0.008s
== Summary =====================================================================
Total elapsed: 64.675s
[OK] --- DONE: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.036s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
[OK]   [03] Create initial NetCDF | dur=0.014s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.014s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.037s
[OK]   [08] Remap (weights) | dur=0.355s
[OK]   [09] Flood: horizontal | dur=0.620s
[OK]   [10] Flood: vertical | dur=0.511s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.468s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.032s
[OK]   [08] Remap (weights) | dur=0.356s
[OK]   [09] Flood: horizontal | dur=0.620s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.065s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.028s
[OK]   [08] Remap (weights) | dur=0.355s
[OK]   [09] Flood: horizontal | dur=0.620s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.066s
[FAIL] [13] Write variables | dur=0.002s | RuntimeError: [ERR] variable missing in bry file: NO3_north
available variables: ['spherical', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'Tcline', 'hc', 'sc_r', 'Cs_r', 'sc_w', 'Cs_w', 'bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time', 'zeta_north', 'zeta_south', 'zeta_east', 'zeta_west', 'ubar_north', 'ubar_south', 'ubar_east', 'ubar_west', 'vbar_north', 'vbar_south', 'vbar_east', 'vbar_west', 'u_north', 'u_south', 'u_east', 'u_west', 'v_north', 'v_south', 'v_east', 'v_west', 'temp_north', 'temp_south', 'temp_east', 'temp_west', 'salt_north', 'salt_south', 'salt_east', 'salt_west', 'bio_time']
       at add_bio_bry.py:237 in <module>
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py", line 237, in <module>
    raise RuntimeError(
RuntimeError: [ERR] variable missing in bry file: NO3_north
available variables: ['spherical', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'Tcline', 'hc', 'sc_r', 'Cs_r', 'sc_w', 'Cs_w', 'bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time', 'zeta_north', 'zeta_south', 'zeta_east', 'zeta_west', 'ubar_north', 'ubar_south', 'ubar_east', 'ubar_west', 'vbar_north', 'vbar_south', 'vbar_east', 'vbar_west', 'u_north', 'u_south', 'u_east', 'u_west', 'v_north', 'v_south', 'v_east', 'v_west', 'temp_north', 'temp_south', 'temp_east', 'temp_west', 'salt_north', 'salt_south', 'salt_east', 'salt_west', 'bio_time']
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py` failed. (See above for error)
[ERR] --- FAILED: BRY bio add ---
[INFO] PBS submitted
[INFO] JOBID=1173.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.034s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[OK]   [03] Create initial NetCDF | dur=0.003s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.014s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.037s
[OK]   [08] Remap (weights) | dur=0.354s
[OK]   [09] Flood: horizontal | dur=0.624s
[OK]   [10] Flood: vertical | dur=0.510s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.468s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.031s
[OK]   [08] Remap (weights) | dur=0.355s
[OK]   [09] Flood: horizontal | dur=0.615s
[OK]   [10] Flood: vertical | dur=0.041s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.065s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.028s
[OK]   [08] Remap (weights) | dur=0.355s
[OK]   [09] Flood: horizontal | dur=0.622s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.065s
[FAIL] [13] Write variables | dur=0.002s | RuntimeError: [ERR] variable missing in bry file: NO3_north
available variables: ['spherical', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'Tcline', 'hc', 'sc_r', 'Cs_r', 'sc_w', 'Cs_w', 'bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time', 'zeta_north', 'zeta_south', 'zeta_east', 'zeta_west', 'ubar_north', 'ubar_south', 'ubar_east', 'ubar_west', 'vbar_north', 'vbar_south', 'vbar_east', 'vbar_west', 'u_north', 'u_south', 'u_east', 'u_west', 'v_north', 'v_south', 'v_east', 'v_west', 'temp_north', 'temp_south', 'temp_east', 'temp_west', 'salt_north', 'salt_south', 'salt_east', 'salt_west', 'bio_time']
       at add_bio_bry.py:237 in <module>
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py", line 237, in <module>
    raise RuntimeError(
RuntimeError: [ERR] variable missing in bry file: NO3_north
available variables: ['spherical', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'Tcline', 'hc', 'sc_r', 'Cs_r', 'sc_w', 'Cs_w', 'bry_time', 'zeta_time', 'temp_time', 'salt_time', 'v2d_time', 'v3d_time', 'zeta_north', 'zeta_south', 'zeta_east', 'zeta_west', 'ubar_north', 'ubar_south', 'ubar_east', 'ubar_west', 'vbar_north', 'vbar_south', 'vbar_east', 'vbar_west', 'u_north', 'u_south', 'u_east', 'u_west', 'v_north', 'v_south', 'v_east', 'v_west', 'temp_north', 'temp_south', 'temp_east', 'temp_west', 'salt_north', 'salt_south', 'salt_east', 'salt_west', 'bio_time']
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py` failed. (See above for error)
[ERR] --- FAILED: BRY bio add ---
[INFO] PBS submitted
[INFO] JOBID=1174.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.034s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[OK]   [03] Create initial NetCDF | dur=0.003s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.013s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.036s
[OK]   [08] Remap (weights) | dur=0.356s
[OK]   [09] Flood: horizontal | dur=0.634s
[OK]   [10] Flood: vertical | dur=0.512s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.470s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.031s
[OK]   [08] Remap (weights) | dur=0.353s
[OK]   [09] Flood: horizontal | dur=0.626s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.064s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.027s
[OK]   [08] Remap (weights) | dur=0.349s
[OK]   [09] Flood: horizontal | dur=0.625s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.065s
[OK]   [13] Write variables | dur=0.005s
== Summary =====================================================================
--- Time elapsed: 4.340s ---
[OK] --- DONE: BRY bio add ---
[INFO] PBS submitted
[INFO] JOBID=1175.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
· src=OGCM grid=roms_grd_15km_nearest.nc
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/hycom_wght_15km_nrst.nc
[OK]   [01] Load configuration and input metadata | dur=0.034s
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Initiating biological variables: None type ---
--- [+] Boundary file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.062s
· Use existing wght file /.../fixed/hycom_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.011s
· Biological variables: npzd
· Flood method (boundary): edt
[OK]   [05] List & group OGCM files | dur=0.000s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.047s
[OK]   [08] Remap (weights) | dur=0.297s
[OK]   [09] Flood H/V | dur=1.198s
[OK]   [10] Mask & rotate | dur=0.123s
[OK]   [11] z→σ & save bry | dur=0.567s
[DONE] 2026-03-11 21 |  6.237s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.018s
[OK]   [08] Remap (weights) | dur=0.299s
[OK]   [09] Flood H/V | dur=0.714s
[OK]   [10] Mask & rotate | dur=0.121s
[OK]   [11] z→σ & save bry | dur=0.163s
[DONE] 2026-03-12 00 | 11.552s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.031s
[OK]   [08] Remap (weights) | dur=0.297s
[OK]   [09] Flood H/V | dur=0.691s
[OK]   [10] Mask & rotate | dur=0.118s
[OK]   [11] z→σ & save bry | dur=0.160s
[DONE] 2026-03-12 03 | 16.849s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.021s
[OK]   [08] Remap (weights) | dur=0.300s
[OK]   [09] Flood H/V | dur=0.693s
[OK]   [10] Mask & rotate | dur=0.115s
[OK]   [11] z→σ & save bry | dur=0.157s
[DONE] 2026-03-12 06 | 22.135s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.036s
[OK]   [08] Remap (weights) | dur=0.294s
[OK]   [09] Flood H/V | dur=0.703s
[OK]   [10] Mask & rotate | dur=0.120s
[OK]   [11] z→σ & save bry | dur=0.164s
[DONE] 2026-03-12 09 | 27.453s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.004s
[OK]   [08] Remap (weights) | dur=0.299s
[OK]   [09] Flood H/V | dur=0.699s
[OK]   [10] Mask & rotate | dur=0.115s
[OK]   [11] z→σ & save bry | dur=0.160s
[DONE] 2026-03-12 12 | 32.730s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.019s
[OK]   [08] Remap (weights) | dur=0.299s
[OK]   [09] Flood H/V | dur=0.699s
[OK]   [10] Mask & rotate | dur=0.122s
[OK]   [11] z→σ & save bry | dur=0.161s
[DONE] 2026-03-12 15 | 38.029s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.031s
[OK]   [08] Remap (weights) | dur=0.300s
[OK]   [09] Flood H/V | dur=0.697s
[OK]   [10] Mask & rotate | dur=0.119s
[OK]   [11] z→σ & save bry | dur=0.156s
[DONE] 2026-03-12 18 | 43.333s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.001s
[OK]   [08] Remap (weights) | dur=0.298s
[OK]   [09] Flood H/V | dur=0.704s
[OK]   [10] Mask & rotate | dur=0.122s
[OK]   [11] z→σ & save bry | dur=0.162s
[DONE] 2026-03-12 21 | 48.622s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.041s
[OK]   [08] Remap (weights) | dur=0.299s
[OK]   [09] Flood H/V | dur=0.696s
[OK]   [10] Mask & rotate | dur=0.118s
[OK]   [11] z→σ & save bry | dur=0.161s
[DONE] 2026-03-13 00 | 53.937s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.027s
[OK]   [08] Remap (weights) | dur=0.300s
[OK]   [09] Flood H/V | dur=0.709s
[OK]   [10] Mask & rotate | dur=0.114s
[OK]   [11] z→σ & save bry | dur=0.159s
[DONE] 2026-03-13 03 | 59.245s
[OK]   [06] Open source file | dur=59.246s
[OK]   [12] Write variables | dur=0.007s
== Summary =====================================================================
Total elapsed: 59.365s
[OK] --- DONE: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.036s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[OK]   [03] Create initial NetCDF | dur=0.107s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.012s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.036s
[OK]   [08] Remap (weights) | dur=0.357s
[OK]   [09] Flood: horizontal | dur=0.625s
[OK]   [10] Flood: vertical | dur=0.512s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.468s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.032s
[OK]   [08] Remap (weights) | dur=0.354s
[OK]   [09] Flood: horizontal | dur=0.629s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.065s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.028s
[OK]   [08] Remap (weights) | dur=0.350s
[OK]   [09] Flood: horizontal | dur=0.627s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.066s
[OK]   [13] Write variables | dur=0.004s
== Summary =====================================================================
--- Time elapsed: 4.444s ---
[OK] --- DONE: BRY bio add ---
[INFO] PBS submitted
[INFO] JOBID=1176.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
9567
./pipe01.sh: line 84: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 94: DT_NL: unbound variable
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 136: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 136: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 136: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 136: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 136: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 170: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe01.sh: line 170: dsa: command not found
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[INFO] --- START: INI build ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/roms_grd_15km_nearest.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_15km_nrst.nc
· ogcm=/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.129s
[('/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc', 1, cftime.DatetimeGregorian(2026, 3, 12, 0, 0, 0, 0, has_year_zero=False), np.float64(9567.0))]
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Deactivate initiating biological variables ---
--- [+] Initial file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.426s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.004s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[OK]   [05] Load OGCM fields | dur=4.269s
[OK]   [06] Remap (weights) | dur=0.327s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=13.998s
[OK]   [08] Flood: vertical | dur=0.449s
[OK]   [09] Mask & clean | dur=0.009s
[OK]   [10] Rotate & stagger (u,v) | dur=0.111s
[OK]   [11] z→σ interpolation | dur=0.982s
[OK]   [12] Conserve volume & fix barotropic | dur=0.201s
[OK]   [13] Write variables | dur=0.039s
== Summary =====================================================================
Total elapsed: 20.951s
[OK] --- DONE: INI build ---
[INFO] --- START: INI bio add ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/roms_grd_15km_nearest.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
· ogcm(daily)=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
· bio_model=npzd
[OK]   [01] Load configuration and input metadata | dur=0.141s
[OK]   [02] Time index matching & relative time | dur=0.006s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· ini_vars=['NO3', 'phytoplankton', 'zooplankton', 'detritus']
[OK]   [03b] Load bio vars (YAML) | dur=0.020s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.001s
!!! SHJO: FIX ME !!!
[OK]   [05] Load OGCM fields | dur=0.102s
[OK]   [06] Remap (weights) | dur=0.373s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=12.818s
[OK]   [08] Flood: vertical | dur=0.459s
[OK]   [10] Mask & rotate | dur=0.011s
[OK]   [11] z→σ interpolation | dur=1.113s
== Summary =====================================================================
Total elapsed: 15.048s
[WARN] skip write (not in ini): NO3
[WARN] skip write (not in ini): phytoplankton
[WARN] skip write (not in ini): zooplankton
[WARN] skip write (not in ini): detritus
[OK]   [12] Write variables | dur=0.002s
== Summary =====================================================================
[OK] --- DONE: INI bio add ---
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
[FAIL] [01] Load configuration and input metadata | dur=0.030s | AssertionError: No OGCM files found in ogcm_path!
       at mk_bry_single_new2.py:33 in <module>
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_bry_single_new2.py", line 33, in <module>
    assert len(filelist) > 0, "No OGCM files found in ogcm_path!"
           ^^^^^^^^^^^^^^^^^
AssertionError: No OGCM files found in ogcm_path!
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_bry_single_new2.py` failed. (See above for error)
[ERR] --- FAILED: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.035s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[FAIL] [03] Create initial NetCDF | dur=0.023s | RuntimeError: NetCDF: String match to name in use: (variable 'NO3_north', group '/')
       at _netCDF4.pyx:2164 in netCDF4._netCDF4._ensure_nc_success
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py", line 77, in <module>
    status = cn.create_bry_bio(
             ^^^^^^^^^^^^^^^^^^
  File "/home/shjo/github/romsforge/libs/create.py", line 363, in create_bry_bio
    var = ncfile.createVariable(varname, "f4", dims_d)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src/netCDF4/_netCDF4.pyx", line 3007, in netCDF4._netCDF4.Dataset.createVariable
  File "src/netCDF4/_netCDF4.pyx", line 4231, in netCDF4._netCDF4.Variable.__init__
  File "src/netCDF4/_netCDF4.pyx", line 2164, in netCDF4._netCDF4._ensure_nc_success
RuntimeError: NetCDF: String match to name in use: (variable 'NO3_north', group '/')
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py` failed. (See above for error)
[ERR] --- FAILED: BRY bio add ---
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
[FAIL] [01] Load configuration and input metadata | dur=0.030s | AssertionError: No OGCM files found in ogcm_path!
       at mk_bry_single_new2.py:33 in <module>
Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_bry_single_new2.py", line 33, in <module>
    assert len(filelist) > 0, "No OGCM files found in ogcm_path!"
           ^^^^^^^^^^^^^^^^^
AssertionError: No OGCM files found in ogcm_path!
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/mk_bry_single_new2.py` failed. (See above for error)
[ERR] --- FAILED: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.036s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[FAIL] [03] Create initial NetCDF | dur=0.009s | ValueError: cannot find dimension s_rho in this group or parent groups
       at utils.py:45 in _find_dim
Traceback (most recent call last):
  File "/usr/local/miniconda3/envs/romsforge/lib/python3.11/site-packages/netCDF4/utils.py", line 39, in _find_dim
    dim = group.dimensions[dimname]
          ^^^^^^^^^^^^^^^^
AttributeError: 'NoneType' object has no attribute 'dimensions'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/usr/local/miniconda3/envs/romsforge/lib/python3.11/site-packages/netCDF4/utils.py", line 43, in _find_dim
    group = group.parent
            ^^^^^^^^^^^^
AttributeError: 'NoneType' object has no attribute 'parent'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py", line 77, in <module>
    status = cn.create_bry_bio(
             ^^^^^^^^^^^^^^^^^^
  File "/home/shjo/github/romsforge/libs/create.py", line 363, in create_bry_bio
    var = ncfile.createVariable(varname, "f4", dims_d)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src/netCDF4/_netCDF4.pyx", line 3005, in genexpr
  File "src/netCDF4/_netCDF4.pyx", line 3005, in genexpr
  File "/usr/local/miniconda3/envs/romsforge/lib/python3.11/site-packages/netCDF4/utils.py", line 45, in _find_dim
    raise ValueError("cannot find dimension %s in this group or parent groups" % dimname)
ValueError: cannot find dimension s_rho in this group or parent groups
ERROR conda.cli.main_run:execute(127): `conda run python /home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/add_bio_bry.py` failed. (See above for error)
[ERR] --- FAILED: BRY bio add ---
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
· src=OGCM grid=roms_grd_15km_nearest.nc
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/hycom_wght_15km_nrst.nc
[OK]   [01] Load configuration and input metadata | dur=0.034s
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Initiating biological variables: None type ---
--- [+] Boundary file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.063s
· Use existing wght file /.../fixed/hycom_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.011s
· Biological variables: npzd
· Flood method (boundary): edt
[OK]   [05] List & group OGCM files | dur=0.000s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=5.004s
[OK]   [08] Remap (weights) | dur=0.295s
[OK]   [09] Flood H/V | dur=1.215s
[OK]   [10] Mask & rotate | dur=0.206s
[OK]   [11] z→σ & save bry | dur=0.606s
[DONE] 2026-03-11 21 |  7.330s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.325s
[OK]   [08] Remap (weights) | dur=0.297s
[OK]   [09] Flood H/V | dur=0.715s
[OK]   [10] Mask & rotate | dur=0.122s
[OK]   [11] z→σ & save bry | dur=0.159s
[DONE] 2026-03-12 00 | 12.948s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.284s
[OK]   [08] Remap (weights) | dur=0.294s
[OK]   [09] Flood H/V | dur=0.699s
[OK]   [10] Mask & rotate | dur=0.120s
[OK]   [11] z→σ & save bry | dur=0.158s
[DONE] 2026-03-12 03 | 18.505s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.165s
[OK]   [08] Remap (weights) | dur=0.298s
[OK]   [09] Flood H/V | dur=0.693s
[OK]   [10] Mask & rotate | dur=0.116s
[OK]   [11] z→σ & save bry | dur=0.153s
[DONE] 2026-03-12 06 | 23.928s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.095s
[OK]   [08] Remap (weights) | dur=0.292s
[OK]   [09] Flood H/V | dur=0.702s
[OK]   [10] Mask & rotate | dur=0.125s
[OK]   [11] z→σ & save bry | dur=0.162s
[DONE] 2026-03-12 09 | 29.304s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.126s
[OK]   [08] Remap (weights) | dur=0.293s
[OK]   [09] Flood H/V | dur=0.695s
[OK]   [10] Mask & rotate | dur=0.110s
[OK]   [11] z→σ & save bry | dur=0.159s
[DONE] 2026-03-12 12 | 34.687s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.167s
[OK]   [08] Remap (weights) | dur=0.298s
[OK]   [09] Flood H/V | dur=0.694s
[OK]   [10] Mask & rotate | dur=0.114s
[OK]   [11] z→σ & save bry | dur=0.156s
[DONE] 2026-03-12 15 | 40.117s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.102s
[OK]   [08] Remap (weights) | dur=0.293s
[OK]   [09] Flood H/V | dur=0.695s
[OK]   [10] Mask & rotate | dur=0.123s
[OK]   [11] z→σ & save bry | dur=0.161s
[DONE] 2026-03-12 18 | 45.490s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.192s
[OK]   [08] Remap (weights) | dur=0.292s
[OK]   [09] Flood H/V | dur=0.697s
[OK]   [10] Mask & rotate | dur=0.122s
[OK]   [11] z→σ & save bry | dur=0.157s
[DONE] 2026-03-12 21 | 50.951s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.116s
[OK]   [08] Remap (weights) | dur=0.297s
[OK]   [09] Flood H/V | dur=0.701s
[OK]   [10] Mask & rotate | dur=0.112s
[OK]   [11] z→σ & save bry | dur=0.155s
[DONE] 2026-03-13 00 | 56.332s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.110s
[OK]   [08] Remap (weights) | dur=0.294s
[OK]   [09] Flood H/V | dur=0.694s
[OK]   [10] Mask & rotate | dur=0.116s
[OK]   [11] z→σ & save bry | dur=0.154s
[DONE] 2026-03-13 03 | 61.701s
[OK]   [06] Open source file | dur=61.701s
[OK]   [12] Write variables | dur=0.029s
== Summary =====================================================================
Total elapsed: 61.845s
[OK] --- DONE: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_15km_nrst.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.045s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.006s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[OK]   [03] Create initial NetCDF | dur=0.107s
· Use existing wght file /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_15km_nrst.nc
[OK]   [04] Prepare weights | dur=0.013s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.131s
[OK]   [08] Remap (weights) | dur=0.355s
[OK]   [09] Flood: horizontal | dur=0.622s
[OK]   [10] Flood: vertical | dur=0.531s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.470s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.033s
[OK]   [08] Remap (weights) | dur=0.357s
[OK]   [09] Flood: horizontal | dur=0.622s
[OK]   [10] Flood: vertical | dur=0.041s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.065s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.056s
[OK]   [08] Remap (weights) | dur=0.353s
[OK]   [09] Flood: horizontal | dur=0.625s
[OK]   [10] Flood: vertical | dur=0.040s
[OK]   [11] Mask | dur=0.011s
[OK]   [12] z→σ & save | dur=0.066s
[OK]   [13] Write variables | dur=0.004s
== Summary =====================================================================
--- Time elapsed: 4.596s ---
[OK] --- DONE: BRY bio add ---
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
./pipe04.sh: line 48: JOBID: unbound variable
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[INFO] --- START: INI build ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/NWP4_grd_314_10m.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc
· ogcm=/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.099s
[('/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc', 1, cftime.DatetimeGregorian(2026, 3, 12, 0, 0, 0, 0, has_year_zero=False), np.float64(9567.0))]
[OK]   [02] Time index matching & relative time | dur=0.007s
--- [NOTE] Deactivate initiating biological variables ---
--- [+] Initial file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.136s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc ---
[OK]   [04] Prepare weights | dur=5.146s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[OK]   [05] Load OGCM fields | dur=4.809s
[OK]   [06] Remap (weights) | dur=0.103s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=3.432s
[OK]   [08] Flood: vertical | dur=0.136s
[OK]   [09] Mask & clean | dur=0.002s
[OK]   [10] Rotate & stagger (u,v) | dur=0.028s
[OK]   [11] z→σ interpolation | dur=0.438s
[OK]   [12] Conserve volume & fix barotropic | dur=0.049s
[OK]   [13] Write variables | dur=0.015s
== Summary =====================================================================
Total elapsed: 14.401s
[OK] --- DONE: INI build ---
[INFO] --- START: INI bio add ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/NWP4_grd_314_10m.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
· ogcm(daily)=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
· bio_model=npzd
[OK]   [01] Load configuration and input metadata | dur=0.115s
[OK]   [02] Time index matching & relative time | dur=0.006s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· ini_vars=['NO3', 'phytoplankton', 'zooplankton', 'detritus']
[OK]   [03b] Load bio vars (YAML) | dur=0.008s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc ---
[+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
[OK]   [04] Prepare weights | dur=1.074s
!!! SHJO: FIX ME !!!
[OK]   [05] Load OGCM fields | dur=0.079s
[OK]   [06] Remap (weights) | dur=0.118s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=3.140s
[OK]   [08] Flood: vertical | dur=0.152s
[OK]   [10] Mask & rotate | dur=0.003s
[OK]   [11] z→σ interpolation | dur=0.437s
== Summary =====================================================================
Total elapsed: 5.135s
[WARN] skip write (not in ini): NO3
[WARN] skip write (not in ini): phytoplankton
[WARN] skip write (not in ini): zooplankton
[WARN] skip write (not in ini): detritus
[OK]   [12] Write variables | dur=0.002s
== Summary =====================================================================
[OK] --- DONE: INI bio add ---
[INFO] --- START: BRY build ---
== Boundary Build ==============================================================
· src=OGCM grid=NWP4_grd_314_10m.nc
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/hycom_wght_24km_314.nc
[OK]   [01] Load configuration and input metadata | dur=0.015s
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Initiating biological variables: None type ---
--- [+] Boundary file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.047s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc ---
[OK]   [04] Prepare weights | dur=4.154s
· Biological variables: npzd
· Flood method (boundary): edt
[OK]   [05] List & group OGCM files | dur=0.000s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-11 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.845s
[OK]   [08] Remap (weights) | dur=0.102s
[OK]   [09] Flood H/V | dur=0.718s
[OK]   [10] Mask & rotate | dur=0.033s
[OK]   [11] z→σ & save bry | dur=0.484s
[DONE] 2026-03-11 21 |  5.185s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.030s
[OK]   [08] Remap (weights) | dur=0.103s
[OK]   [09] Flood H/V | dur=0.233s
[OK]   [10] Mask & rotate | dur=0.037s
[OK]   [11] z→σ & save bry | dur=0.052s
[DONE] 2026-03-12 00 |  9.640s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=4.020s
[OK]   [08] Remap (weights) | dur=0.103s
[OK]   [09] Flood H/V | dur=0.232s
[OK]   [10] Mask & rotate | dur=0.032s
[OK]   [11] z→σ & save bry | dur=0.051s
[DONE] 2026-03-12 03 | 14.079s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 06 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.845s
[OK]   [08] Remap (weights) | dur=0.101s
[OK]   [09] Flood H/V | dur=0.235s
[OK]   [10] Mask & rotate | dur=0.033s
[OK]   [11] z→σ & save bry | dur=0.053s
[DONE] 2026-03-12 06 | 18.347s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 09 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.784s
[OK]   [08] Remap (weights) | dur=0.101s
[OK]   [09] Flood H/V | dur=0.236s
[OK]   [10] Mask & rotate | dur=0.033s
[OK]   [11] z→σ & save bry | dur=0.052s
[DONE] 2026-03-12 09 | 22.553s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 12 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.810s
[OK]   [08] Remap (weights) | dur=0.103s
[OK]   [09] Flood H/V | dur=0.233s
[OK]   [10] Mask & rotate | dur=0.033s
[OK]   [11] z→σ & save bry | dur=0.054s
[DONE] 2026-03-12 12 | 26.785s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 15 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.818s
[OK]   [08] Remap (weights) | dur=0.101s
[OK]   [09] Flood H/V | dur=0.234s
[OK]   [10] Mask & rotate | dur=0.034s
[OK]   [11] z→σ & save bry | dur=0.055s
[DONE] 2026-03-12 15 | 31.027s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 18 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.810s
[OK]   [08] Remap (weights) | dur=0.101s
[OK]   [09] Flood H/V | dur=0.235s
[OK]   [10] Mask & rotate | dur=0.031s
[OK]   [11] z→σ & save bry | dur=0.054s
[DONE] 2026-03-12 18 | 35.259s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 21 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.811s
[OK]   [08] Remap (weights) | dur=0.101s
[OK]   [09] Flood H/V | dur=0.235s
[OK]   [10] Mask & rotate | dur=0.036s
[OK]   [11] z→σ & save bry | dur=0.053s
[DONE] 2026-03-12 21 | 39.496s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 00 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.786s
[OK]   [08] Remap (weights) | dur=0.100s
[OK]   [09] Flood H/V | dur=0.233s
[OK]   [10] Mask & rotate | dur=0.033s
[OK]   [11] z→σ & save bry | dur=0.052s
[DONE] 2026-03-13 00 | 43.700s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-13 03 file=hycom_korea_20260312.nc
[OK]   [07] Load OGCM fields | dur=3.812s
[OK]   [08] Remap (weights) | dur=0.102s
[OK]   [09] Flood H/V | dur=0.234s
[OK]   [10] Mask & rotate | dur=0.035s
[OK]   [11] z→σ & save bry | dur=0.052s
[DONE] 2026-03-13 03 | 47.936s
[OK]   [06] Open source file | dur=47.936s
[OK]   [12] Write variables | dur=0.007s
== Summary =====================================================================
Total elapsed: 52.165s
[OK] --- DONE: BRY build ---
[INFO] --- START: BRY bio add ---
== Boundary Build (bio) ========================================================
· out=/.../20260312/roms_nipa_bry_20260312.nc wght=/.../fixed/cmems_wght_24km_314.nc
· daily_files=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.015s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· bio_model=npzd
· bry_vars=['NO3', 'phyt', 'zoop', 'detritus']
[OK]   [02] Load bio vars (YAML) | dur=0.008s
· bry steps=3
[OK]   [03] Time index matching & relative time | dur=0.005s
OPEN BRY FILE: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
ABS PATH: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_bry_20260312.nc
[OK]   [03] Create initial NetCDF | dur=0.070s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc ---
[+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
[OK]   [04] Prepare weights | dur=1.101s
[OK]   [05] Allocate bry buffers | dur=0.000s
[OK]   [06] Group time entries by file | dur=0.000s
[NOTE] Flood method: edt
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.026s
[OK]   [08] Remap (weights) | dur=0.116s
[OK]   [09] Flood: horizontal | dur=0.210s
[OK]   [10] Flood: vertical | dur=0.522s
[OK]   [11] Mask | dur=0.003s
[OK]   [12] z→σ & save | dur=0.448s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.077s
[OK]   [08] Remap (weights) | dur=0.115s
[OK]   [09] Flood: horizontal | dur=0.206s
[OK]   [10] Flood: vertical | dur=0.012s
[OK]   [11] Mask | dur=0.003s
[OK]   [12] z→σ & save | dur=0.023s
!!! SHJO: FIXME phyto=chl/(0.02*6.625*12) !!!
[OK]   [07] Load OGCM fields | dur=0.023s
[OK]   [08] Remap (weights) | dur=0.114s
[OK]   [09] Flood: horizontal | dur=0.201s
[OK]   [10] Flood: vertical | dur=0.010s
[OK]   [11] Mask | dur=0.003s
[OK]   [12] z→σ & save | dur=0.023s
[OK]   [13] Write variables | dur=0.004s
== Summary =====================================================================
--- Time elapsed: 3.345s ---
[OK] --- DONE: BRY bio add ---
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1178.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[INFO] --- START: INI build ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/NWP4_grd_314_10m.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc
· ogcm=/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.014s
[('/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc', 1, cftime.DatetimeGregorian(2026, 3, 12, 0, 0, 0, 0, has_year_zero=False), np.float64(9567.0))]
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Deactivate initiating biological variables ---
--- [+] Initial file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.096s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc ---
[OK]   [04] Prepare weights | dur=4.127s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[OK]   [05] Load OGCM fields | dur=3.726s
[OK]   [06] Remap (weights) | dur=0.101s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=3.466s
[OK]   [08] Flood: vertical | dur=0.134s
[OK]   [09] Mask & clean | dur=0.002s
[OK]   [10] Rotate & stagger (u,v) | dur=0.029s
[OK]   [11] z→σ interpolation | dur=0.376s
[OK]   [12] Conserve volume & fix barotropic | dur=0.050s
[OK]   [13] Write variables | dur=0.014s
== Summary =====================================================================
Total elapsed: 12.141s
[OK] --- DONE: INI build ---
[INFO] --- START: INI bio add ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/NWP4_grd_314_10m.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
· ogcm(daily)=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
· bio_model=npzd
[OK]   [01] Load configuration and input metadata | dur=0.015s
[OK]   [02] Time index matching & relative time | dur=0.006s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· ini_vars=['NO3', 'phytoplankton', 'zooplankton', 'detritus']
[OK]   [03b] Load bio vars (YAML) | dur=0.008s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc ---
[+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
[OK]   [04] Prepare weights | dur=0.988s
!!! SHJO: FIX ME !!!
[OK]   [05] Load OGCM fields | dur=0.025s
[OK]   [06] Remap (weights) | dur=0.120s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=3.161s
[OK]   [08] Flood: vertical | dur=0.175s
[OK]   [10] Mask & rotate | dur=0.003s
[OK]   [11] z→σ interpolation | dur=0.438s
== Summary =====================================================================
Total elapsed: 4.943s
[WARN] skip write (not in ini): NO3
[WARN] skip write (not in ini): phytoplankton
[WARN] skip write (not in ini): zooplankton
[WARN] skip write (not in ini): detritus
[OK]   [12] Write variables | dur=0.002s
== Summary =====================================================================
[OK] --- DONE: INI bio add ---
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1179.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[INFO] --- START: INI build ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/NWP4_grd_314_10m.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc
· ogcm=/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc
[OK]   [01] Load configuration and input metadata | dur=0.014s
[('/home/shjo/github/romsforge/packages/nipa_auto/data/hycom/hycom_korea_20260312.nc', 1, cftime.DatetimeGregorian(2026, 3, 12, 0, 0, 0, 0, has_year_zero=False), np.float64(9567.0))]
[OK]   [02] Time index matching & relative time | dur=0.005s
--- [NOTE] Deactivate initiating biological variables ---
--- [+] Initial file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc ---
[OK]   [03] Create initial NetCDF | dur=0.100s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/hycom_wght_24km_314.nc ---
[OK]   [04] Prepare weights | dur=4.146s
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[WARN] utils.py:291 — invalid value encountered in divide | ts=2026-03-12 00:00:00
[OK]   [05] Load OGCM fields | dur=3.736s
[OK]   [06] Remap (weights) | dur=0.105s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=3.434s
[OK]   [08] Flood: vertical | dur=0.136s
[OK]   [09] Mask & clean | dur=0.002s
[OK]   [10] Rotate & stagger (u,v) | dur=0.029s
[OK]   [11] z→σ interpolation | dur=0.372s
[OK]   [12] Conserve volume & fix barotropic | dur=0.050s
[OK]   [13] Write variables | dur=0.014s
== Summary =====================================================================
Total elapsed: 12.144s
[OK] --- DONE: INI build ---
[INFO] --- START: INI bio add ---
== Initial Condition Build (ini) ===============================================
· grid=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/NWP4_grd_314_10m.nc
· ini_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312/roms_nipa_ini_20260312.nc
· wght_out=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
· ogcm(daily)=/home/shjo/github/romsforge/packages/nipa_auto/data/cmems/cmems_bio_20260312.nc
· bio_model=npzd
[OK]   [01] Load configuration and input metadata | dur=0.015s
[OK]   [02] Time index matching & relative time | dur=0.005s
· bio_yaml=/home/shjo/github/romsforge/packages/nipa_auto/src/preprocess/bio_vars.yml
· ini_vars=['NO3', 'phytoplankton', 'zooplankton', 'detritus']
[OK]   [03b] Load bio vars (YAML) | dur=0.009s
[OK]   [03] Create initial NetCDF | dur=0.098s
--- [+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc ---
[+] Weight file created: /home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/fixed/cmems_wght_24km_314.nc
[OK]   [04] Prepare weights | dur=0.979s
!!! SHJO: FIX ME !!!
[OK]   [05] Load OGCM fields | dur=0.025s
[OK]   [06] Remap (weights) | dur=0.119s
[NOTE] Flood method: griddata
[OK]   [07] Flood: horizontal | dur=3.151s
[OK]   [08] Flood: vertical | dur=0.152s
[OK]   [10] Mask & rotate | dur=0.003s
[OK]   [11] z→σ interpolation | dur=0.443s
== Summary =====================================================================
Total elapsed: 5.002s
[OK]   [12] Write variables | dur=0.013s
== Summary =====================================================================
[OK] --- DONE: INI bio add ---
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1180.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1181.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
./pipe04.sh: line 48: JOBID: unbound variable
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1183.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1185.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1187.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1188.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1190.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
./pipe04.sh: line 11: PBS_NODEFILE: unbound variable
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1191.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1192.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1193.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1194.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1195.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1197.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1198.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1202.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04.sh
[INFO] PBS submitted
[INFO] JOBID=1208.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04_mvp2.sh
./pipeline.sh: line 12: pipe04.sh: No such file or directory
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
./pipeline.sh: line 9: ./pipe02.sh: Permission denied
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
./pipe02.sh: line 15: /pipeline_status.md: Permission denied
grep: /pipeline_status.md: No such file or directory
[INFO] --- START: Download HYCOM frcst ---
./pipe02.sh: line 51: /bin/conda: No such file or directory
[ERR] --- FAILED: Download HYCOM frcst ---
sed: can't read /pipeline_status.md: No such file or directory
grep: /pipeline_status.md: No such file or directory
[INFO] --- START: CMEMS download ---
./pipe02.sh: line 51: /bin/conda: No such file or directory
[ERR] --- FAILED: CMEMS download ---
sed: can't read /pipeline_status.md: No such file or directory
grep: /pipeline_status.md: No such file or directory
[INFO] --- START: INI build ---
./pipe02.sh: line 51: /bin/conda: No such file or directory
[ERR] --- FAILED: INI build ---
sed: can't read /pipeline_status.md: No such file or directory
grep: /pipeline_status.md: No such file or directory
[INFO] --- START: INI bio add ---
./pipe02.sh: line 51: /bin/conda: No such file or directory
[ERR] --- FAILED: INI bio add ---
sed: can't read /pipeline_status.md: No such file or directory
grep: /pipeline_status.md: No such file or directory
[INFO] --- START: BRY build ---
./pipe02.sh: line 51: /bin/conda: No such file or directory
[ERR] --- FAILED: BRY build ---
sed: can't read /pipeline_status.md: No such file or directory
grep: /pipeline_status.md: No such file or directory
[INFO] --- START: BRY bio add ---
./pipe02.sh: line 51: /bin/conda: No such file or directory
[ERR] --- FAILED: BRY bio add ---
sed: can't read /pipeline_status.md: No such file or directory
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04_mvp2.sh
#!/bin/bash
set -eu

#PBS -q workq
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -N nipa_20260312
#PBS -j oe
#PBS -o /home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312/pbs_20260312.log

module load intel21/compiler-21
module load intel21/openmpi-5.0.0
module load intel21/netcdf-4.6.1

export HDF5_USE_FILE_LOCKING=FALSE
export HDF5_COLL_METADATA_WRITE=0

NP=$(wc -l < "$PBS_NODEFILE")
cd /home/shjo/github/romsforge/packages/nipa_auto/roms_application/roms_820_nemuscs_latest

mpirun -np 16   -x OMP_NUM_THREADS   -x HDF5_USE_FILE_LOCKING   -x HDF5_COLL_METADATA_WRITE   --map-by slot   --bind-to core   --report-bindings   --mca coll ^hcoll,han   /home/shjo/github/romsforge/packages/nipa_auto/roms_application/roms_820_nemuscs_latest/execM_nl /home/shjo/github/romsforge/packages/nipa_auto/roms_application/roms_820_nemuscs_latest/nl_nipa.in > /home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312/mdl_log_20260312.md 2>&1

[mpiexec@nft00] set_default_values (ui/mpich/utils.c:1488): no executable specified
[mpiexec@nft00] HYD_uii_mpx_get_parameters (ui/mpich/utils.c:1739): setting default values failed
[mpiexec@nft00] main (ui/mpich/mpiexec.c:149): error parsing parameters
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04_mvp2.sh
[INFO] PBS submitted
[INFO] JOBID=1209.nft00
[INFO] RUN_DATE=20260312
[INFO] INPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_inputs/20260312
[INFO] OUTPUTS_DIR=/home/shjo/github/romsforge/packages/nipa_auto/roms_outputs/20260312
[INFO] INI_DATE=2026-03-12 00:00:00
./pipe02.sh
[INFO] pipeline_status exists → resume
[SKIP] Download HYCOM frcst already done
[SKIP] CMEMS download already done
[SKIP] INI build already done
[SKIP] INI bio add already done
[SKIP] BRY build already done
[SKIP] BRY bio add already done
./pipe04_mvp2.sh
[INFO] PBS submitted
[INFO] JOBID=1210.nft00
