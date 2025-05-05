import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import streamlit as st
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from matplotlib.widgets import Button
from matplotlib.animation import FFMpegWriter
from shapely.geometry import LineString
from shapely.ops import substring
from datetime import datetime
import contextily as ctx
import matplotlib as mpl


# Define base directory 
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Define relative paths
RAW_DATA_DIR = os.path.join(BASE_DIR, "rawdata")
SHAPEFILE_DIR = os.path.join(BASE_DIR, "shapefiles")
ANIMATIONS_DIR = os.path.join(BASE_DIR, "animations")

#For MP3 Videos
mpl.rcParams['animation.ffmpeg_path'] = os.path.join(BASE_DIR, "ffmpeg.exe")

# Paths to data files
SEGMENTS_FILE = os.path.join(RAW_DATA_DIR, "Segments.xlsx")
TIMETABLE_FILE = os.path.join(RAW_DATA_DIR, "Timetable2manual.xlsx")
LINE_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Line.shp")
STOPS_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Stops.shp")

# Load timetable and segments
segments_df = pd.read_excel(SEGMENTS_FILE)
timetable_df = pd.read_excel(TIMETABLE_FILE)

# Load shapefiles
railline = gpd.read_file(LINE_SHP)
stations = gpd.read_file(STOPS_SHP)

# Reproject to EPSG:3857 
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

fig, ax = plt.subplots(figsize=(14, 14))
plt.subplots_adjust(bottom=0.15)

railline.plot(ax=ax, color='orange', linewidth=2.5, label='Rail Line')
stations.plot(ax=ax, color='green', marker='o', zorder=4, markersize=15, label='Stations')

passing_points = ["Killucan", "Carrick loop"]
for point in passing_points:
    if point in station_mileposts:
        milepost = station_mileposts[point]
        ratio = (milepost - min_milepost) / (max_milepost - min_milepost)
        geom = track.interpolate(ratio * track.length)
        ax.plot(geom.x, geom.y, 'o', color='green', markersize=4, zorder=4, label="Pass Point" if point == passing_points[0] else "")
        ax.text(geom.x+10000, geom.y+4000, point, fontsize=7, ha='right', color='darkgoldenrod')

for idx, row in stations.iterrows():
        label = f"{row['StopName']}"
        ax.text(row.geometry.x -4000, row.geometry.y-4000, label, fontsize=7, ha='right', color='darkgoldenrod')

ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)
               #ctx.providers.CartoDB.Voyager)
                #ctx.providers.CartoDB.DarkMatter)
                #ctx.providers.OpenStreetMap.Mapnik)
                

#speed slider
ax_slider = plt.axes([0.5, 0.05, 0.2, 0.07])
speed_slider = Slider(ax_slider, 'Speed Factor', valmin=1, valmax=10, valinit=3, valstep=1)

clock_text = ax.text(0.5, 0.95, "Time: 00:00", transform=ax.transAxes, ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.7))

# Global pause state
is_paused = False

# Add pause button
ax_pause = plt.axes([0.3, 0.06, 0.1, 0.03])  # [left, bottom, width, height]
pause_button = Button(ax_pause, 'Pause', color='aqua', hovercolor='red')

# Pause button callback
def toggle_pause(event):
    global is_paused
    is_paused = not is_paused
    pause_button.label.set_text('Resume' if is_paused else 'Pause')
    fig.canvas.draw_idle()

pause_button.on_clicked(toggle_pause)

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

    # Color and label
    is_test_train = str(train_id).startswith("T")
    is_evening_test = str(train_id).startswith("TE")
    is_dublin_bound = str(train_id).startswith("1")
    is_sligo_bound = str(train_id).startswith("2")
    if is_test_train:
        color = 'red'
    elif is_dublin_bound:
        color = 'black'
    elif is_sligo_bound:
        color = 'black'
    else:
        color = 'black'
    # if is_evening_test:
    #     color = 'yellow'
    #     marker_shape = '^'
    #     marker_size = 9
            
    marker_size = 8 if is_test_train else 5
    marker_shape = '^' if is_test_train else 'o'
    label_color = 'darkred' if is_test_train else 'dimgray'

    dot, = ax.plot([], [], marker_shape, markersize=marker_size, label=f"Train {train_id}", color=color)
    label = ax.text(0, 0, str(train_id), fontsize=7, color=label_color, ha='left', va='center', visible=not is_test_train)
    
        
    train_dots[train_id] = {
        'dot': dot,
        'label': label,
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
start_time = timetable_df["Departs_s"].min() - 1000
end_time = timetable_df["Departs_s"].max() + 1000

FRAME_STEP = 5
FRAMES = 4500
INTERVAL = 15

sim_time = start_time

# Animation function
def animate_trains(frame):
   
    global sim_time, last_artists
    if is_paused:
        return last_artists  # return previous frame's state without change


    global sim_time
    SPEED_FACTOR = speed_slider.val
    sim_time += FRAME_STEP * SPEED_FACTOR
  #  time_slider.set_val(sim_time)
    t = (sim_time - start_time) % (end_time - start_time) + start_time

    updated_artists = []

    for train_id, data in train_dots.items():
        arrivals = data['arrivals']
        departs = data['departs']
        positions = data['positions']
        dwell = data['actual_dwell']
        dot = data['dot']
        label = data['label']

        first_depart = departs[0]
        if first_depart - 300 <= t < first_depart:
            pos_pt = track.interpolate(positions[0])
            dot.set_data([pos_pt.x], [pos_pt.y])
            label.set_position((pos_pt.x + 2000, pos_pt.y))
            label.set_visible(True)
            updated_artists.extend([dot, label])
            continue

        # Final hold at last station
        final_depart = departs[-1]
        final_pos = positions[-1]
        if final_depart <= t <= final_depart + 300:
            pos_pt = track.interpolate(final_pos)
            dot.set_data([pos_pt.x], [pos_pt.y])
            label.set_position((pos_pt.x + 2000, pos_pt.y))
            label.set_visible(True)
            updated_artists.extend([dot, label])
            continue

        if t > final_depart + 300:
            dot.set_data([], [])
            label.set_visible(False)
            continue

        idx = np.searchsorted(departs, t, side='right') - 1
        if idx < 0 or idx >= len(positions) - 1:
            dot.set_data([], [])
            label.set_visible(False)
            continue

        t_depart = departs[idx]
        t_arrive_next = arrivals[idx + 1]
        t_depart_next = departs[idx + 1]
        p0 = positions[idx]
        p1 = positions[idx + 1]

        # Before departure: hold at platform
        if t < t_depart:
            dist = p0

        # Between departure and next arrival: move based on computed speed
        elif t_depart <= t < t_arrive_next:
            segment_duration = t_arrive_next - t_depart
            if segment_duration <= 0:
                dist = p1
            else:
                progress = (t - t_depart) / segment_duration
                dist = (1 - progress) * p0 + progress * p1

        # At arrival, hold until next departure
        elif t_arrive_next <= t < t_depart_next:
            dist = p1

        else:
            dot.set_data([], [])
            label.set_visible(False)
            continue

        pos_pt = track.interpolate(dist)
        dot.set_data([pos_pt.x], [pos_pt.y])
        label.set_position((pos_pt.x + 2000, pos_pt.y))
        label.set_visible(True)
        updated_artists.extend([dot, label])

    clock_text.set_text(f"Time: {datetime.utcfromtimestamp(t).strftime('%H:%M')}")
    updated_artists.append(clock_text)

    if frame == 0:
        
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    last_artists = updated_artists  # store for paused state
    return updated_artists


ani = FuncAnimation(fig, animate_trains, frames=FRAMES, interval=INTERVAL, blit=True, repeat=False)
plt.legend()
plt.show()

# writer = FFMpegWriter(fps=30)
# ani.save(os.path.join(ANIMATIONS_DIR, "Ani_2.mp4"), writer=writer)
# print("🎬 Saving animation to Ani_2.mp4...")


