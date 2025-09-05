

# Import modules
import copernicusmarine

copernicusmarine.subset(
  dataset_id="METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2",
  variables=["analysed_sst", "analysis_error", "mask"],
  minimum_longitude=110,
  maximum_longitude=140,
  minimum_latitude=30,
  maximum_latitude=55,
  start_datetime="2025-06-01T00:00:00",
  end_datetime="2025-06-30T00:00:00",
)

# Print loaded dataset information
#print(sst_l3s)
