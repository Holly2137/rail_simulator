import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from shapely.geometry import LineString
from shapely.ops import substring
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
    track = railline.geometry.iloc[0].unary_union
else:
    track = railline.geometry.iloc[0]

# Use existing milepost columns from segments
segments_df["From_milepost"] = segments_df["From Milepost"]
segments_df["To_milepost"] = segments_df["To Milepost"]

# Compute station mileposts from segments
station_mileposts = pd.concat([
    segments_df[['From', 'From_milepost']].rename(columns={'From': 'Station', 'From_milepost': 'Milepost'}),
    segments_df[['To', 'To_milepost']].rename(columns={'To': 'Station', 'To_milepost': 'Milepost'})
]).dropna().drop_duplicates(subset='Station').set_index('Station')['Milepost'].to_dict()

# Map mileposts to station geometries
stations["milepost"] = stations["StopName"].map(station_mileposts)

# Interpolation mapper: map milepost to geometry along the track
min_milepost = stations["milepost"].min()
max_milepost = stations["milepost"].max()

# Normalize mileposts to range 0.0 - 1.0 and interpolate geometry
stations["milepost_position"] = stations["milepost"].apply(lambda m: (m - min_milepost) / (max_milepost - min_milepost) if pd.notnull(m) else np.nan)
stations["geometry"] = stations["milepost_position"].apply(lambda ratio: track.interpolate(ratio * track.length) if pd.notnull(ratio) else None)

# Convert departure time to seconds since midnight
def time_to_seconds(tstr):
    try:
        t = datetime.strptime(tstr, "%H:%M")
        return t.hour * 3600 + t.minute * 60
    except:
        return np.nan

timetable_df["Departs_s"] = timetable_df["Depart time"].apply(time_to_seconds)

# Clean and validate timetable
timetable_df["Station"] = timetable_df["Station"].str.strip()
timetable_df = timetable_df.dropna(subset=["Station", "Departs_s"])
timetable_df = timetable_df[timetable_df["Station"].isin(station_mileposts.keys())]

# Create a shared plot
fig, ax = plt.subplots(figsize=(10, 10))

# Plot the rail lines and stations
railline.plot(ax=ax, color='red', label='Rail Line')
stations.plot(ax=ax, color='green', marker='o', zorder=3, markersize=10, label='Stations')

# Add station labels for verification
for idx, row in stations.iterrows():
    label = f"{row['StopName']}\n{row['milepost']}m"
    ax.text(row.geometry.x, row.geometry.y, label, fontsize=7, ha='right', color='black')

# Add a basemap
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

# Visual test: plot segments as true rail geometry using shapely substring
for _, row in segments_df.iterrows():
    from_mile = row["From_milepost"]
    to_mile = row["To_milepost"]

    if pd.notnull(from_mile) and pd.notnull(to_mile):
        start_ratio = (from_mile - min_milepost) / (max_milepost - min_milepost)
        end_ratio = (to_mile - min_milepost) / (max_milepost - min_milepost)

        start_distance = start_ratio * track.length
        end_distance = end_ratio * track.length

        segment_geom = substring(track, start_distance, end_distance)
        x, y = segment_geom.xy
        ax.plot(x, y, color='blue', linewidth=2, alpha=0.7)

# 🚆 Animate trains
train_ids = timetable_df['ID'].unique()
train_dots = {}

for train_id in train_ids:
    train_data = timetable_df[timetable_df['ID'] == train_id].sort_values("Departs_s")
    mileposts = train_data["Station"].map(station_mileposts)
    ratios = (mileposts - min_milepost) / (max_milepost - min_milepost)
    positions = ratios * track.length
    times = train_data["Departs_s"].values

    if len(times) < 2:
        continue

    dot, = ax.plot([], [], 'o', label=f"Train {train_id}", markersize=6)
    train_dots[train_id] = {
        'dot': dot,
        'positions': positions.values,
        'times': times
    }

max_time = timetable_df["Departs_s"].max()

clock_text = ax.text(0.5, 0.95, "Time: 00:00", transform=ax.transAxes, ha='center', fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

FRAMES = 1000
INTERVAL = 20

def update(frame):
    t = frame * (max_time / FRAMES)
    artists = []

    for train_id, data in train_dots.items():
        times = data['times']
        pos = data['positions']

        idx = np.searchsorted(times, t, side='right') - 1
        if idx < 0 or idx >= len(pos) - 1:
            continue

        t0, t1 = times[idx], times[idx+1]
        p0, p1 = pos[idx], pos[idx+1]
        alpha = (t - t0) / (t1 - t0) if t1 > t0 else 0
        dist = (1 - alpha) * p0 + alpha * p1
        pt = track.interpolate(dist)

        data['dot'].set_data(pt.x, pt.y)
        artists.append(data['dot'])

    clock_text.set_text(f"Time: {datetime.utcfromtimestamp(t).strftime('%H:%M')}")
    artists.append(clock_text)
    return artists

ani = FuncAnimation(fig, update, frames=FRAMES, interval=INTERVAL, blit=True)
plt.legend()
plt.show()


















