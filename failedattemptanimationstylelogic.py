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
from collections import defaultdict

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



# === Build Per-Train Segment Timing ===

detected_segments = {}

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

        track_row = segments_df[((segments_df['From'] == from_station) & (segments_df['To'] == to_station)) |
                                ((segments_df['From'] == to_station) & (segments_df['To'] == from_station))]
        track_type = track_row.iloc[0]['Track'] if not track_row.empty else 'Single'

        segment_info.append({
            'from': from_station,
            'to': to_station,
            'enter': departs,
            'exit': arrival,
            'track': track_type,
            'train_id': train_id
        })

    detected_segments[train_id] = segment_info

# === Conflict Detection (Second-by-Second Simulation Style) ===

seen_conflicts = set()

# Occupation map for each time step: { time: { segment: train_id } }
segment_occupancy = defaultdict(dict)
conflicts = []

# Determine simulation bounds
sim_start = timetable_df["Departs_s"].min() - 300  # start 5 min early
sim_end = timetable_df["Departs_s"].max() + 600    # buffer at end

for t in range(int(sim_start), int(sim_end) + 1):
    for train_id, segments in detected_segments.items():
        for seg in segments:
            if seg['track'].lower() == 'double':
                continue
            if seg['to'] == 'Connolly':
                continue

            if seg['enter'] is None or seg['exit'] is None:
                continue
            if not (seg['enter'] <= t <= seg['exit']):
                continue

            segment_key = (min(seg['from'], seg['to']), max(seg['from'], seg['to']))
            occupying_train = segment_occupancy[t].get(segment_key)

            if occupying_train and occupying_train != train_id:
                # Conflict: another train is already here at this time
                conflict_key = (segment_key, occupying_train, train_id)
                if conflict_key not in seen_conflicts:
                    seen_conflicts.add(conflict_key)
                    conflicts.append({
                        'segment': f"{seg['from']} → {seg['to']}",
                        'train_a': occupying_train,
                        'train_b': train_id,
                        'conflict_time': t,
                        'milepost_from': station_mileposts.get(seg['from'], 'N/A'),
                        'milepost_to': station_mileposts.get(seg['to'], 'N/A')
                    })
            else:
                segment_occupancy[t][segment_key] = train_id

# === Conflict Resolution (Speed Adjustment) ===

for conflict in conflicts:
    segment_key = tuple(sorted(conflict['segment'].split(' → ')))
    conflict_time = conflict['conflict_time']

    # Determine train that enters segment first
    seg_a = next(seg for seg in detected_segments[conflict['train_a']] if tuple(sorted((seg['from'], seg['to']))) == segment_key)
    seg_b = next(seg for seg in detected_segments[conflict['train_b']] if tuple(sorted((seg['from'], seg['to']))) == segment_key)

    if seg_a['enter'] < seg_b['enter']:
        earlier_train = conflict['train_a']
        later_train_enter = seg_b['enter']
    else:
        earlier_train = conflict['train_b']
        later_train_enter = seg_a['enter']

    # Find how far back this train goes before the conflict segment
    segments = detected_segments[earlier_train]
    idx = next((i for i, seg in enumerate(segments)
                if tuple(sorted((seg['from'], seg['to']))) == segment_key), None)

    if idx is None:
        continue

    # Step back to find last fixed departure
    adjust_start_idx = idx
    for j in range(idx - 1, -1, -1):
        if segments[j]['enter'] is not None:
            adjust_start_idx = j
            break

    # Total original time span
    original_enter = segments[adjust_start_idx]['enter']
    target_exit = later_train_enter - 60  # must clear 60s before next train
    total_time = target_exit - original_enter

    # Total tight runtime
    tight_sum = 0
    for seg in segments[adjust_start_idx:idx + 1]:
        runtime = seg.get('tight_runtime')
        if pd.notna(runtime):
            tight_sum += runtime

    if tight_sum == 0:
        continue  # avoid div by 0

    # Adjust segment times proportionally
    elapsed = 0
    for seg in segments[adjust_start_idx:idx + 1]:
        runtime = seg.get('tight_runtime')
        if pd.notna(runtime):
            seg_duration = total_time * (runtime / tight_sum)
            seg['enter'] = original_enter + elapsed
            seg['exit'] = seg['enter'] + seg_duration
            elapsed += seg_duration

# Print detected conflicts (BEFORE adjustment)
print("Detected Conflicts (second-by-second style):")
for conflict in conflicts:
    print(f"Conflict on {conflict['segment']} between Train {conflict['train_a']} and {conflict['train_b']} "
          f"at {conflict['conflict_time']}s, mileposts: {conflict['milepost_from']} → {conflict['milepost_to']}")

# Print adjustments made

# === Recheck Conflicts After Adjustment ===
print("=== Rechecking for Remaining Conflicts ===")
segment_occupancy = defaultdict(dict)
remaining_conflicts = []
seen_conflicts = set()

for t in range(int(sim_start), int(sim_end) + 1):
    for train_id, segments in detected_segments.items():
        for seg in segments:
            if seg['track'].lower() == 'double':
                continue
            if seg['to'] == 'Connolly':
                continue
            if seg['enter'] is None or seg['exit'] is None:
                continue
            if not (seg['enter'] <= t <= seg['exit']):
                continue

            segment_key = (min(seg['from'], seg['to']), max(seg['from'], seg['to']))
            occupying_train = segment_occupancy[t].get(segment_key)

            if occupying_train and occupying_train != train_id:
                conflict_key = (segment_key, occupying_train, train_id)
                if conflict_key not in seen_conflicts:
                    seen_conflicts.add(conflict_key)
                    remaining_conflicts.append({
                        'segment': f"{seg['from']} → {seg['to']}",
                        'train_a': occupying_train,
                        'train_b': train_id,
                        'conflict_time': t
                    })
            else:
                segment_occupancy[t][segment_key] = train_id

if remaining_conflicts:
    for conflict in remaining_conflicts:
        print(f"Remaining conflict on {conflict['segment']} between Train {conflict['train_a']} and {conflict['train_b']} "
              f"at {conflict['conflict_time']}s")
else:
    print("No remaining conflicts found. ✅")
print("Adjusted Segment Timings:")
for train_id, segments in detected_segments.items():
    for seg in segments:
        print(f"Train {train_id}: {seg['from']} → {seg['to']} | enter: {seg['enter']}, exit: {seg['exit']}")








# def animate_trains(frame):
#     global sim_time
#     SPEED_FACTOR = speed_slider.val
#     sim_time += FRAME_STEP * SPEED_FACTOR
#     t = (sim_time - start_time) % (end_time - start_time) + start_time

#     updated_artists = []

#     for train_id, data in train_dots.items():
#         arrivals = data['arrivals']
#         departs = data['departs']
#         positions = data['positions']
#         dwell = data['actual_dwell']
#         dot = data['dot']

#         first_depart = departs[0]
#         if first_depart - 300 <= t < first_depart:
#             pos_pt = track.interpolate(positions[0])
#             dot.set_data([pos_pt.x], [pos_pt.y])
#             updated_artists.append(dot)
#             continue

#         # Final hold at last station
#         final_depart = departs[-1]
#         final_pos = positions[-1]
#         if final_depart <= t <= final_depart + 300:
#             pos_pt = track.interpolate(final_pos)
#             dot.set_data([pos_pt.x], [pos_pt.y])
#             updated_artists.append(dot)
#             continue

#         if t > final_depart + 300:
#             dot.set_data([], [])
#             continue

#         idx = np.searchsorted(departs, t, side='right') - 1
#         if idx < 0 or idx >= len(positions) - 1:
#             dot.set_data([], [])
#             continue

#         t_depart = departs[idx]
#         t_arrive_next = arrivals[idx + 1]
#         t_depart_next = departs[idx + 1]
#         p0 = positions[idx]
#         p1 = positions[idx + 1]

#         # Before departure: hold at platform
#         if t < t_depart:
#             dist = p0

#         # Between departure and next arrival: move based on computed speed
#         elif t_depart <= t < t_arrive_next:
#             segment_duration = t_arrive_next - t_depart
#             if segment_duration <= 0:
#                 dist = p1
#             else:
#                 progress = (t - t_depart) / segment_duration
#                 dist = (1 - progress) * p0 + progress * p1

#         # At arrival, hold until next departure
#         elif t_arrive_next <= t < t_depart_next:
#             dist = p1

#         else:
#             dot.set_data([], [])
#             continue

#         pos_pt = track.interpolate(dist)
#         dot.set_data([pos_pt.x], [pos_pt.y])
#         updated_artists.append(dot)

#     clock_text.set_text(f"Time: {datetime.utcfromtimestamp(t).strftime('%H:%M')}")
#     updated_artists.append(clock_text)

#     if frame == 0:
#         ax.legend(loc='upper right')

#     return updated_artists
# ani = FuncAnimation(fig, animate_trains, frames=FRAMES, interval=INTERVAL, blit=True)


plt.legend()
plt.show()

