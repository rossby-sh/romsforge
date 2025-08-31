from netCDF4 import Dataset
import numpy as np

def add_global_attrs_to_obs(ncfiles):
    if isinstance(ncfiles, str):
        ncfiles = [ncfiles]

    for ncfile in ncfiles:
        print(f"🔧 Updating: {ncfile}")
        with Dataset(ncfile, 'a') as ds:
            # ===== 글로벌 속성 =====
            ds.title = 'North Western Pacific 25 km'
            ds.grd_file = 'roms_grd_fennel_15km_smooth_v2'
            ds.grid_Lm_Mm_N = np.array([354, 276, 36], dtype='i4')
            ds.variance_units = 'squared state variable units'

            ds.state_variables = "\n".join([
                "1: free-surface (m)",
                "2: vertically integrated u-momentum component (m/s)",
                "3: vertically integrated v-momentum component (m/s)",
                "4: u-momentum component (m/s)",
                "5: v-momentum component (m/s)",
                "6: temperature (Celsius)",
                "7: salinity",
                "8: NO3",
                "9: phytoplankton",
                "10: zooplankton",
                "11: detritus"
            ])

            ds.obs_provenance = "\n".join([
                "1: AVISO SLA",
                "2: blended SST",
                "3: XBT temperature",
                "4: CTD temperature",
                "5: CTD salinity",
                "6: ARGO temperature",
                "7: ARGO salinity",
                "8: CalCOFI temperature",
                "9: CalCOFI salinity",
                "10: GLOBEC temperature",
                "11: GLOBEC salinity",
                "12: buoy temp",
                "13: CTD Chlorophyll",
                "14: Satellite Chlorophyll"
            ])

            ds.obs_sources = "\n".join([
                "http://opendap.aviso.oceanobs.com/",
                "http://hadobs.metoffice.com/",
                "http://oceandata.sci.gsfc.nasa.gov/"
            ])

            # ===== 변수 속성 (flag_values, flag_meanings) 추가 =====
            if 'obs_type' in ds.variables:
                var = ds.variables['obs_type']
                var.setncattr('flag_values', np.arange(1, 12, dtype='i4'))
                var.setncattr('flag_meanings',
                              "zeta ubar vbar u v temperature salinity NO3 phytoplankton zooplankton detritus")

            if 'obs_provenance' in ds.variables:
                var = ds.variables['obs_provenance']
                var.setncattr('flag_values', np.arange(1, 15, dtype='i4'))
                var.setncattr('flag_meanings',
                              "AVISO_SLA blended_SST XBT_temperature CTD_temperature CTD_salinity "
                              "ARGO_temperature ARGO_salinity CalCOFI_temperature CalCOFI_salinity "
                              "GLOBEC_temperature GLOBEC_salinity buoy_temp CTD_Chlorophyll Satellite_Chlorophyll")

        print(f"✅ Done: {ncfile}\n")
        

pth='D:/shjo/ROMS_inputs/obs/pro/'
        
NC=Dataset('D:/shjo/ROMS_inputs/roms_grd_fennel_15km_smooth_v2.nc')



add_global_attrs_to_obs([pth+'obs_SST_OSTIA_30km.nc',pth+'obs_phyt_27km.nc',pth+'obs_KODC.nc'])



add_global_attrs_to_obs([pth+'obs_SST_OSTIA_APR_30km_re.nc',pth+'obs_phyt_27km.nc',pth+'obs_KODC.nc'])














