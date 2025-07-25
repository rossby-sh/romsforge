from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import copernicusmarine

# Product and dataset IDs
DATASET_ID = "cmems_mod_glo_bgc_my_0.25deg_P1M-m"
VARIABLE = ["chl","fe","no3","nppv","o2","ph",
            "phyc","po4","si","spco2"]

# Coordinates
longitude = (100, 170)
latitude = (5, 60)

# Depth range (if applicable)
depth_min = 0.5057600140571594 # Update with your depth range
depth_max = 5902.05810546875   # Update with your depth range

# Output directory
output_directory = "data"

# Boundary dates
start_date = datetime(1993, 1, 1, 00)
end_date = datetime(2022, 12, 1, 00)

# Delta time
# delta_t = timedelta(days=30)  # Interval to download by request
delta_t = relativedelta(months=1)  #

# Loop over the time range with the specified delta time
current_date = start_date
while current_date <= end_date:
    start_datetime = current_date
    # end_datetime = min(current_date + delta_t, end_date)
    end_datetime = (current_date + relativedelta(months=1)) - timedelta(days=1)

    # Output filename
    # output_filename = "CMEMS_data_{}_{}_{}_{}_{}_{}.nc".format(
    # latitude[0], latitude[1], longitude[0], longitude[1],
    # start_datetime.strftime("%Y-%m-%d"),
    # end_datetime.strftime("%Y-%m-%d"))
    output_filename = "CMEMS_data_{}.nc".format(    
    start_datetime.strftime("%Y-%m-%d"))    

    # Call the subset function
    
    copernicusmarine.subset(
        dataset_id=DATASET_ID,
        dataset_version="202406",
        variables=VARIABLE,
        minimum_longitude=longitude[0],
        maximum_longitude=longitude[1],
        minimum_latitude=latitude[0],
        maximum_latitude=latitude[1],
        start_datetime=start_datetime,
        end_datetime=end_datetime,
        minimum_depth=depth_min,
        maximum_depth=depth_max,
        output_filename=output_filename,
        output_directory=output_directory,
        force_download = True   
    )
    print(f"Requesting data from {start_datetime} to {end_datetime}")

    # Move to the next time interval
    current_date += delta_t
