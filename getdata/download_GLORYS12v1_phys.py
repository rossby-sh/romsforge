from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import copernicusmarine

# Product and dataset IDs
# DATASET_ID = "cmems_mod_glo_phy_my_0.083deg_P1M-m"   # Monthly Mean
DATASET_ID = "cmems_mod_glo_phy_my_0.083deg_P1D-m"   # Daily Mean
VARIABLE = ["zos","thetao","so","uo","vo"]

# Coordinates
longitude = (100, 160)
latitude = (10, 55)

# Depth range (if applicable)
depth_min = 0.5057600140571594 # Update with your depth range
depth_max = 5902.05810546875   # Update with your depth range

# Output directory
output_directory = "Phy/Daily_mean/2021/"

# Boundary dates
start_date = datetime(2021, 1, 1, 00)
end_date = datetime(2021, 12, 31, 00)

# Delta time
delta_t = timedelta(days=1)  # Interval to download by request
#delta_t = relativedelta(months=1)  #

# Loop over the time range with the specified delta time
current_date = start_date
while current_date <= end_date:
    start_datetime = current_date
    # end_datetime = min(current_date + delta_t, end_date)

    # end_datetime = (current_date + relativedelta(months=1)) - timedelta(days=1)  # Monthly Mean
    # end_datetime = (current_date + timedelta(days=1))                              # Daily Mean
    end_datetime = (current_date)                              # Daily Mean

    if current_date >= datetime(2021, 7, 1, 00):
        # DATASET_ID = "cmems_mod_glo_phy_myint_0.083deg_P1M-m"  # Monthly Mean
        DATASET_ID = "cmems_mod_glo_phy_myint_0.083deg_P1M-m"  # Daily Mean

    # Output filename
    # output_filename = "CMEMS_data_{}_{}_{}_{}_{}_{}.nc".format(
    # latitude[0], latitude[1], longitude[0], longitude[1],
    # start_datetime.strftime("%Y-%m-%d"),
    # end_datetime.strftime("%Y-%m-%d"))
    
    # output_filename = "CMEMS_Phy_data_{}.nc".format(     # Monthly Mean
    # start_datetime.strftime("%Y-%m"))

    output_filename = "CMEMS_Phy_data_{}.nc".format(       # Daily Mean
    start_datetime.strftime("%Y-%m-%d"))    

    # Call the subset function
    
    copernicusmarine.subset(
        dataset_id=DATASET_ID,        
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