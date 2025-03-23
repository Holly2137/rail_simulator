import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from shapely.geometry import LineString
from datetime import datetime
import contextily as ctx

# Define base directory dynamically from the script location
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Define relative paths
RAW_DATA_DIR = os.path.join(BASE_DIR, "rawdata")
SHAPEFILE_DIR = os.path.join(BASE_DIR, "shapefiles")

# Paths to data files
SEGMENTS_FILE = os.path.join(RAW_DATA_DIR, "Segments.xlsx")
TIMETABLE_FILE = os.path.join(RAW_DATA_DIR, "Timetable.xlsx")
LINE_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Line.shp")
STOPS_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Stops2.shp")

# Load timetable and segments
segments_df = pd.read_excel(SEGMENTS_FILE)
timetable_df = pd.read_excel(TIMETABLE_FILE)

# Load shapefiles
railline = gpd.read_file(LINE_SHP)
stations = gpd.read_file(STOPS_SHP)

# Reproject to EPSG:3857 if necessary
if railline.crs.to_string() != "EPSG:3857":
    railline = railline.to_crs(epsg=3857)
if stations.crs.to_string() != "EPSG:3857":
    stations = stations.to_crs(epsg=3857)

print("Data successfully loaded.")

 # Clean station names
stations["StopName"] = stations["StopName"].str.replace(" Train Station", "", regex=False)
stations["StopName"] = stations["StopName"].str.strip()

# Ensure the rail line is a single LineString
if railline.geometry.iloc[0].geom_type == 'MultiLineString':
    track = railline.geometry.iloc[0].unary_union  # Merge into one LineString
else:
    track = railline.geometry.iloc[0]  # Already a LineString
    
# Use milepost values instead of shapefile-projected distances
station_mileposts = timetable_df.drop_duplicates(subset=["ID", "Station"])[["Station", "Milepost(m)"]].copy()
station_mileposts = station_mileposts.dropna().drop_duplicates(subset="Station").set_index("Station")["Milepost(m)"].to_dict()

stations["milepost"] = stations["StopName"].map(station_mileposts)

# Interpolation mapper: map milepost to geometry along the track
min_milepost = stations["milepost"].min()
max_milepost = stations["milepost"].max()

# Normalize mileposts to range 0.0 - 1.0 and interpolate geometry
stations["milepost_position"] = stations["milepost"].apply(lambda m: (m - min_milepost) / (max_milepost - min_milepost) if pd.notnull(m) else np.nan)
stations["geometry_interpolated"] = stations["milepost_position"].apply(lambda ratio: track.interpolate(ratio * track.length) if pd.notnull(ratio) else None)
    
# Merge interpolated milepost positions into timetable for later use
timetable_df["From_milepost"] = timetable_df["From"].map(station_mileposts)
timetable_df["To_milepost"] = timetable_df["To"].map(station_mileposts)



 # Create a shared plot
fig, ax = plt.subplots(figsize=(10, 10))

 # Plot the rail lines and stations
railline.plot(ax=ax, color='red', label='Rail Line')
stations.plot(ax=ax, color='green', marker='o',zorder=3, markersize=10, label='Stations')

# # Add station labels for verification
# for idx, row in stations.iterrows():
#     label = f"{row['StopName']}\n{row['milepost']}m"
#     ax.text(row.geometry.x, row.geometry.y, label, fontsize=7, ha='right', color='black')
  
# Add a basemap
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

# # Merge interpolated milepost positions into timetable for later use
# segments_df["From_milepost"] = segments_df["From"].map(station_mileposts)
# segments_df["To_milepost"] = segments_df["To"].map(station_mileposts)

    
# Visual test: plot one train path (e.g., ID 101)
sample_train = timetable_df[timetable_df["ID"] == 101].copy()

for _, row in sample_train.iterrows():
    from_mile = row["From_milepost"]
    to_mile = row["To_milepost"]

    if pd.notnull(from_mile) and pd.notnull(to_mile):
        start_ratio = (from_mile - min_milepost) / (max_milepost - min_milepost)
        end_ratio = (to_mile - min_milepost) / (max_milepost - min_milepost)

        start_geom = track.interpolate(start_ratio * track.length)
        end_geom = track.interpolate(end_ratio * track.length)

        ax.plot([start_geom.x, end_geom.x], [start_geom.y, end_geom.y], color='blue', linewidth=2, alpha=0.5)

plt.show()



        
