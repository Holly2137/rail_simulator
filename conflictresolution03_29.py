import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
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

sligo_geom = stations.loc[stations["StopName"] == "Sligo", "geometry"].values[0]
dublin_geom = stations.loc[stations["StopName"] == "Connolly", "geometry"].values[0]

sligo_pos = track.project(sligo_geom)
dublin_pos = track.project(dublin_geom)

if sligo_pos < dublin_pos:
    track = LineString(list(track.coords)[::-1])

sligo_pos = track.project(sligo_geom)
dublin_pos = track.project(dublin_geom)

segments_df["From_milepost"] = segments_df["From Milepost"]
segments_df["To_milepost"] = segments_df["To Milepost"]

station_mileposts = pd.concat([
    segments_df[['From', 'From_milepost']].rename(columns={'From': 'StopName', 'From_milepost': 'Milepost'}),
    segments_df[['To', 'To_milepost']].rename(columns={'To': 'StopName', 'To_milepost': 'Milepost'})
]).dropna().drop_duplicates(subset='StopName').set_index('StopName')['Milepost'].to_dict()

stations["milepost"] = stations["StopName"].map(station_mileposts)

min_milepost = stations["milepost"].min()
max_milepost = stations["milepost"].max()
stations["milepost_ratio"] = stations["milepost"].apply(lambda m: (m - min_milepost) / (max_milepost - min_milepost) if pd.notnull(m) else np.nan)
stations["geometry"] = stations["milepost_ratio"].apply(lambda r: track.interpolate(r * track.length) if pd.notnull(r) else None)

segments = []

for _, row in segments_df.iterrows():
    from_stop = row["From"].strip()
    to_stop = row["To"].strip()

    segment = {
        "from": from_stop,
        "to": to_stop,
        "track_type": row["Track"].strip(),
        "tight_run": row["Tight Run"],
        "free_pass": row["Track"].strip().lower() == "double"
    }
    segments.append(segment)

timetable_df["Passing"] = timetable_df["Passing"].astype(str).str.strip().str.upper()
passing_points = set(timetable_df.loc[timetable_df["Passing"] == "Y", "Station"].unique())

timetable_df["Stop"] = timetable_df["Stop"].astype(str).str.strip().str.upper()
timetable_df["IsStop"] = timetable_df["Stop"].isin(["START", "Y", "END"])

timetable_df["Minimum Dwell"] = pd.to_numeric(timetable_df["Minimum Dwell"], errors="coerce")
timetable_df.loc[timetable_df["IsStop"], "Minimum Dwell"] = timetable_df.loc[timetable_df["IsStop"], "Minimum Dwell"].fillna(60).astype(int)
timetable_df["Minimum Dwell"] = timetable_df["Minimum Dwell"].fillna(0).astype(int)

def time_to_seconds(tstr):
    for fmt in ("%H:%M", "%H:%M:%S"):
        try:
            t = datetime.strptime(str(tstr).strip(), fmt)
            return t.hour * 3600 + t.minute * 60 + t.second
        except:
            continue
    return np.nan

timetable_df["Departs_s"] = timetable_df["Depart time"].apply(time_to_seconds)
timetable_df["Station"] = timetable_df["Station"].str.strip()
timetable_df = timetable_df.dropna(subset=["Station", "Departs_s"])
timetable_df = timetable_df[timetable_df["Station"].isin(station_mileposts.keys())]

fig, ax = plt.subplots(figsize=(10, 10))
plt.subplots_adjust(bottom=0.15)

railline.plot(ax=ax, color='red', label='Rail Line')
stations.plot(ax=ax, color='green', marker='o', zorder=3, markersize=10, label='Stations')

passing_points = ["Killucan", "Carrick loop"]
for point in passing_points:
    if point in station_mileposts:
        milepost = station_mileposts[point]
        ratio = (milepost - min_milepost) / (max_milepost - min_milepost)
        geom = track.interpolate(ratio * track.length)
        ax.plot(geom.x, geom.y, 'o', color='yellow', markersize=3, zorder=4, label="Passing Point" if point == passing_points[0] else "")
        ax.text(geom.x, geom.y, point, fontsize=7, ha='right', color='darkgoldenrod')

for idx, row in stations.iterrows():
    label = f"{row['StopName']}\n{row['milepost']}m"
    ax.text(row.geometry.x, row.geometry.y, label, fontsize=7, ha='right', color='darkgoldenrod')

ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

ax_slider = plt.axes([0.4, 0.05, 0.2, 0.07])
speed_slider = Slider(ax_slider, 'Speed Factor', valmin=1, valmax=10, valinit=3, valstep=1)

clock_text = ax.text(0.5, 0.95, "Time: 00:00", transform=ax.transAxes, ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.7))

train_ids = timetable_df['ID'].unique()
train_dots = {}

for train_id in train_ids:
    train_data = timetable_df[timetable_df['ID'] == train_id].sort_values("Departs_s")
    mileposts = train_data["Station"].map(station_mileposts)
    times = train_data["Departs_s"].values
    ratios = (mileposts - min_milepost) / (max_milepost - min_milepost)
    positions = ratios * track.length
    min_dwell = train_data["Minimum Dwell"].values
    actual_dwell = min_dwell.copy()  # ✅ Currently same as min dwell, can be adjusted later
    arrivals = times - actual_dwell 

    dot, = ax.plot([], [], 'o', markersize=6, label=f"Train {train_id}")

    train_dots[train_id] = {
        'dot': dot,
        'positions': positions.values,
        'times': times,
        'arrivals': arrivals,
        'departs': times,
        'min_dwell': min_dwell.tolist(),
        'actual_dwell': actual_dwell.tolist()
    }

print("First 10 entries for Train", train_id)
for i in range(min(10, len(times))):
    print(f"Station {i+1}: Depart={times[i]}, Arrival={arrivals[i]}, Dwell={actual_dwell[i]}")
start_time = timetable_df["Departs_s"].min() - 2000
end_time = timetable_df["Departs_s"].max() + 3000

FRAME_STEP = 5
FRAMES = 12000
INTERVAL = 15

sim_time = start_time









# === Conflict Detection Preparation (before animation setup) ===

from collections import defaultdict, deque

def find_last_scheduled_departure(train_id, before_station):
    """Returns the last station in that train's list before `before_station` that has a fixed departure time."""
    train_rows = timetable_df[timetable_df['ID'] == train_id].sort_values("Departs_s")
    prev_depart = None
    for i, row in enumerate(train_rows.itertuples()):
        if row.Station == before_station:
            break
        if row.Stop in ('START', 'Y'):  # Fixed scheduled stop
            prev_depart = row.Station
    return prev_depart

 #(before animation setup) ===

# This dictionary will hold per-train segment movement info
from collections import defaultdict, deque

detected_segments = {}
segment_graph = defaultdict(list)

for _, row in segments_df.iterrows():
    a = row['From']
    b = row['To']
    segment_graph[a].append(b)
    segment_graph[b].append(a)

# Function to find physical segments between any two stations (BFS path search)
def find_path_segments(from_station, to_station):
    visited = set()
    queue = deque([(from_station, [])])

    while queue:
        current, path = queue.popleft()
        if current == to_station:
            return path
        if current in visited:
            continue
        visited.add(current)
        for neighbor in segment_graph[current]:
            segment_row = segments_df[((segments_df['From'] == current) & (segments_df['To'] == neighbor)) |
                                      ((segments_df['From'] == neighbor) & (segments_df['To'] == current))]
            if not segment_row.empty:
                queue.append((neighbor, path + [segment_row.iloc[0]]))
    return []

# Build per-train segment structure
for train_id in timetable_df['ID'].unique():
    train_data = timetable_df[timetable_df['ID'] == train_id].sort_values("Departs_s")
    segment_info = []

    for i in range(len(train_data) - 1):
        row_from = train_data.iloc[i]
        row_to = train_data.iloc[i + 1]

        from_station = row_from['Station']
        to_station = row_to['Station']
        departs = row_from['Departs_s']
        arrival = row_to['Departs_s'] - row_to['Minimum Dwell']

        if pd.isna(departs) or pd.isna(arrival):
            departs = departs if not pd.isna(departs) else None
            arrival = arrival if not pd.isna(arrival) else None

        segment_path = find_path_segments(from_station, to_station)

        for seg in segment_path:
            intermediate_from = seg['From']
            intermediate_to = seg['To']
            direction = 'forward' if intermediate_from == from_station else 'reverse'
            segment_info.append({
                'from': intermediate_from,
                'to': intermediate_to,
                'enter': departs,
                'exit': arrival,
                'direction': direction,
                'track': seg['Track'],
                'tight_runtime': seg['Tight Run'],
                'train_id': train_id
            })

    detected_segments[train_id] = segment_info

# === Conflict Detection & Resolution (Single Track Only, pre-Connolly) ===
conflicts = []
train_ids = list(detected_segments.keys())

for i in range(len(train_ids)):
    train_a = train_ids[i]
    for seg_a in detected_segments[train_a]:
        if seg_a['track'].lower() == 'double':
            continue
        if seg_a['to'] == 'Connolly':
            continue

        for j in range(i + 1, len(train_ids)):
            train_b = train_ids[j]
            for seg_b in detected_segments[train_b]:
                if seg_b['track'].lower() == 'double':
                    continue
                if seg_b['to'] == 'Connolly':
                    continue

                same_segment = (
                    (seg_a['from'] == seg_b['from'] and seg_a['to'] == seg_b['to']) or
                    (seg_a['from'] == seg_b['to'] and seg_a['to'] == seg_b['from'])
                )

                if not same_segment:
                    continue

                # Check for overlapping time windows
                if seg_a['enter'] is None or seg_a['exit'] is None:
                    continue
                if seg_b['enter'] is None or seg_b['exit'] is None:
                    continue

                latest_start = max(seg_a['enter'], seg_b['enter'])
                earliest_end = min(seg_a['exit'], seg_b['exit'])

                if latest_start < earliest_end:
                    # Determine which train arrived first
                    if seg_a['enter'] < seg_b['enter']:
                        later_train_id = train_b
                        later_seg = seg_b
                        earlier_seg = seg_a
                    else:
                        later_train_id = train_a
                        later_seg = seg_a
                        earlier_seg = seg_b

                    # Adjust speed of earlier train to exit 60s before later train's entry
                    for seg in detected_segments[earlier_seg['train_id']]:
                        if seg['from'] == earlier_seg['from'] and seg['to'] == earlier_seg['to']:
                            seg['exit'] = min(seg['exit'], later_seg['enter'] - 60) if seg['exit'] else later_seg['enter'] - 60

                    conflicts.append({
                        'segment': f"{seg_a['from']} → {seg_a['to']}",
                        'train_a': train_a,
                        'train_b': train_b,
                        'overlap_start': latest_start,
                        'overlap_end': earliest_end,
                        'conflict_time': latest_start,
                        'milepost_from': station_mileposts.get(seg_a['from'], 'N/A'),
                        'milepost_to': station_mileposts.get(seg_a['to'], 'N/A')
                    })

# Print detected conflicts
print("Detected Conflicts (with timing and milepost position):")
for conflict in conflicts:
    print(f"Conflict on {conflict['segment']} between Train {conflict['train_a']} and {conflict['train_b']} "
      f"from {conflict['overlap_start']} to {conflict['overlap_end']} (at {conflict['conflict_time']}s), "
      f"mileposts: {conflict['milepost_from']} → {conflict['milepost_to']}")



