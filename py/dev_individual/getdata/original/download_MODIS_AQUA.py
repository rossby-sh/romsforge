from datetime import datetime, timedelta
import xarray as xr

start_date = datetime(2025, 6, 2)
end_date = datetime(2025, 7, 1)

delta = timedelta(days=1)

while start_date <= end_date:
    yyyy = start_date.strftime('%Y')
    mmdd = start_date.strftime('%m%d')

    yyyymmdd = start_date.strftime('%Y%m%d')
    
    url = f'https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/{yyyy}/{mmdd}/AQUA_MODIS.{yyyymmdd}.L3m.DAY.CHL.chlor_a.9km.NRT.nc'
           #'https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2025/0101/AQUA_MODIS.20250101.L3m.DAY.CHL.chlor_a.9km.nc.dmr.html'
    try:
        ds = xr.open_dataset(url)
        chl_sub = ds['chlor_a'].sel(lon=slice(100, 170), lat=slice(60, 5))
        out_file = f'/data/share/DATA/RAW/MODIS_AQUA_CHL_NRT/AQUA_MODIS.{yyyymmdd}.NW_Pacific.chlor_a_NRT.nc'
        chl_sub.to_netcdf(out_file)
        print(f'Saved: {out_file}')
    except Exception as e:
        print(f'Failed for {yyyymmdd}: {e}')

    start_date += delta
