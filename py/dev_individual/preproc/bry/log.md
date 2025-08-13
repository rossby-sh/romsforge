[90m== Boundary Build [0m[90m==============================================================[0m
[90m· src=OGCM grid=roms_grd_fennel_15km_smooth_v2.nc[0m
[90m· out=/.../bry/test.nc wght=/.../wght/NWP_weight.nc[0m
[OK]   [01] Load configuration and input metadata | dur=0.016s
[OK]   [02] Time index matching & relative time | dur=4.758s
--- [NOTE] Initiating biological variables: npzd type ---
--- [+] Boundary file created: /data/share/DATA/ROMS_INPUTS/bry/test.nc ---
[OK]   [03] Create initial NetCDF | dur=0.031s
[90m· Use existing wght file /.../wght/NWP_weight.nc[0m
[OK]   [04] Prepare weights | dur=0.012s
[90m· Biological variables: npzd[0m
[90m· Flood method (boundary): edt[0m
[OK]   [05] List & group OGCM files | dur=0.000s
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-04-30 00 file=HYCOM_20250430_00UTC.nc
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-04-30 00 file=HYCOM_20250430_00UTC.nc
[OK]   [07] Load OGCM fields | dur=1.964s
[OK]   [08] Remap (weights) | dur=0.333s
[OK]   [09] Flood H/V | dur=1.392s
[OK]   [10] Mask & rotate | dur=0.151s
[OK]   [11] z→σ & save bry | dur=0.640s
[94m[DONE][0m 2025-04-30 00 |  4.486s
[OK]   [06] Open source file | dur=4.487s
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-05-01 00 file=HYCOM_20250501_00UTC.nc
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-05-01 00 file=HYCOM_20250501_00UTC.nc
[OK]   [07] Load OGCM fields | dur=2.428s
[OK]   [08] Remap (weights) | dur=0.356s
[OK]   [09] Flood H/V | dur=0.840s
[OK]   [10] Mask & rotate | dur=0.152s
[OK]   [11] z→σ & save bry | dur=0.209s
[94m[DONE][0m 2025-05-01 00 |  3.991s
[OK]   [06] Open source file | dur=3.991s
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-05-02 00 file=HYCOM_20250502_00UTC.nc
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-05-02 00 file=HYCOM_20250502_00UTC.nc
[OK]   [07] Load OGCM fields | dur=4.753s
[OK]   [08] Remap (weights) | dur=0.357s
[OK]   [09] Flood H/V | dur=0.817s
[OK]   [10] Mask & rotate | dur=0.146s
[OK]   [11] z→σ & save bry | dur=0.204s
[94m[DONE][0m 2025-05-02 00 |  6.282s
[OK]   [06] Open source file | dur=6.283s
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-05-03 00 file=HYCOM_20250503_00UTC.nc
[93m[WARN][0m utils.py:277 — invalid value encountered in divide | ts=2025-05-03 00 file=HYCOM_20250503_00UTC.nc
[OK]   [07] Load OGCM fields | dur=4.717s
[OK]   [08] Remap (weights) | dur=0.347s
[OK]   [09] Flood H/V | dur=0.820s
[OK]   [10] Mask & rotate | dur=0.142s
[OK]   [11] z→σ & save bry | dur=0.207s
[94m[DONE][0m 2025-05-03 00 |  6.238s
[OK]   [06] Open source file | dur=6.239s
[OK]   [12] Write variables | dur=0.015s
[90m== Summary [0m[90m=====================================================================[0m
Total elapsed: 25.833s
