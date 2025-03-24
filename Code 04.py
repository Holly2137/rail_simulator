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
    
# ✅ Auto-fix track direction: Connolly should be milepost 0, Sligo high
sligo_geom = stations.loc[stations["StopName"] == "Sligo", "geometry"].values[0]
dublin_geom = stations.loc[stations["StopName"] == "Connolly", "geometry"].values[0]

sligo_pos = track.project(sligo_geom)
dublin_pos = track.project(dublin_geom)

print(f"Projected: Connolly={dublin_pos:.2f}, Sligo={sligo_pos:.2f}")

if sligo_pos < dublin_pos:
    track = LineString(list(track.coords)[::-1])
    print("🔁 Track reversed to match milepost direction")

# Confirm direction
sligo_pos = track.project(sligo_geom)
dublin_pos = track.project(dublin_geom)
print(f"After reverse: Connolly={dublin_pos:.2f}, Sligo={sligo_pos:.2f}")


# 🚩 Interpolate accurate geometry using milepost data from Segments
segments_df["From_milepost"] = segments_df["From Milepost"]
segments_df["To_milepost"] = segments_df["To Milepost"]

# Build milepost lookup
station_mileposts = pd.concat([
    segments_df[['From', 'From_milepost']].rename(columns={'From': 'StopName', 'From_milepost': 'Milepost'}),
    segments_df[['To', 'To_milepost']].rename(columns={'To': 'StopName', 'To_milepost': 'Milepost'})
]).dropna().drop_duplicates(subset='StopName').set_index('StopName')['Milepost'].to_dict()

# Assign mileposts to station geometries
stations["milepost"] = stations["StopName"].map(station_mileposts)

# Normalize mileposts and interpolate onto the track
min_milepost = stations["milepost"].min()
max_milepost = stations["milepost"].max()
stations["milepost_ratio"] = stations["milepost"].apply(lambda m: (m - min_milepost) / (max_milepost - min_milepost) if pd.notnull(m) else np.nan)
stations["geometry"] = stations["milepost_ratio"].apply(lambda r: track.interpolate(r * track.length) if pd.notnull(r) else None)




# Convert departure time to seconds since midnight
def time_to_seconds(tstr):
    for fmt in ("%H:%M", "%H:%M:%S"):
        try:
            t = datetime.strptime(str(tstr).strip(), fmt)
            return t.hour * 3600 + t.minute * 60 + t.second
        except:
            continue
    return np.nan

timetable_df["Departs_s"] = timetable_df["Depart time"].apply(time_to_seconds)

# Clean and validate timetable
timetable_df["Station"] = timetable_df["Station"].str.strip()
timetable_df = timetable_df.dropna(subset=["Station", "Departs_s"])
timetable_df = timetable_df[timetable_df["Station"].isin(station_mileposts.keys())]

# Filter to known station list
timetable_df = timetable_df[timetable_df["Station"].isin(station_mileposts.keys())]

# Create a shared plot
fig, ax = plt.subplots(figsize=(10, 10))
plt.subplots_adjust(bottom=0.15)  # Add space for slider

# Plot the rail lines and stations
railline.plot(ax=ax, color='red', label='Rail Line')
stations.plot(ax=ax, color='green', marker='o', zorder=3, markersize=10, label='Stations')

passing_points = ["Killucan", "Carrick loop"]

for point in passing_points:
    if point in station_mileposts:
        milepost = station_mileposts[point]
        ratio = (milepost - min_milepost) / (max_milepost - min_milepost)
        geom = track.interpolate(ratio * track.length)

        # Plot yellow dot
        ax.plot(geom.x, geom.y, 'o', color='yellow', markersize=3, zorder=4, label="Passing Point" if point == passing_points[0] else "")

        # Optional: label it
        ax.text(geom.x, geom.y, point, fontsize=7, ha='right', color='darkgoldenrod')


# Add station labels for verification
for idx, row in stations.iterrows():
    label = f"{row['StopName']}\n{row['milepost']}m"
    ax.text(row.geometry.x, row.geometry.y, label, fontsize=7, ha='right', color='darkgoldenrod')

# Add a basemap
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)


# Add slider for SPEED_FACTOR
ax_slider = plt.axes([0.4, 0.02, 0.2, 0.05])  # [left, bottom, width, height]
speed_slider = Slider(
    ax_slider,
    'Speed Factor',
    valmin=1,
    valmax=10,
    valinit=5,
    valstep=1,
)



# Create train markers and store their data
train_ids = timetable_df['ID'].unique()
train_dots = {}  # Dictionary to store train markers and their data

for train_id in train_ids:
    # Filter timetable data for this train
    train_data = timetable_df[timetable_df['ID'] == train_id].sort_values("Departs_s")
    
    # Get mileposts and times for this train
    mileposts = train_data["Station"].map(station_mileposts)
    times = train_data["Departs_s"].values
    
    # Calculate positions along the track
    ratios = (mileposts - min_milepost) / (max_milepost - min_milepost)
    positions = ratios * track.length
    
    # Create a dot for this train
    dot, = ax.plot([], [], 'o', markersize=6, label=f"Train {train_id}")
    
    # Store the train's data
    train_dots[train_id] = {
        'dot': dot,  # The train marker
        'positions': positions.values,  # Positions along the track
        'times': times  # Departure times
    }

# Add a clock to display the current simulation time
clock_text = ax.text(
    0.5, 0.95, "Time: 00:00",  # Position and initial text
    transform=ax.transAxes,  # Use axes coordinates
    ha='center', fontsize=12,  # Text alignment and size
    bbox=dict(facecolor='white', alpha=0.7)  # Background box
)

# Define animation parameters
start_time = timetable_df["Departs_s"].min() - 2000  #  buffer
end_time = timetable_df["Departs_s"].max() + 3000  #  buffer
#end_time = 83500    # 11:12 PM in seconds#18000  # 05:00 AM in seconds
total_simulation_seconds = end_time - start_time
# FRAME_STEP = 10  # Simulate 30 seconds per frame
# FRAMES = total_simulation_seconds // FRAME_STEP  # Total number of frames
# INTERVAL = 20  # Delay between frames in milliseconds

FRAME_STEP = 5  # constant base time step
FRAMES = 12000  # Or however long you want the loop to run


#SPEED_FACTOR = speed_slider.val
# FRAME_STEP = 10 * SPEED_FACTOR
INTERVAL = 15  # Fast and smooth frame rendering

#INTERVAL = max(1, int(1000 / SPEED_FACTOR))
# FRAMES = 1000  # or: int(total_simulation_seconds // FRAME_STEP)


# # Debugging: Print simulation parameters
# print(f"Start Time: {start_time} (05:00 AM)")
# print(f"End Time: {end_time} (11:12 PM)")
# print(f"Total Simulation Seconds: {total_simulation_seconds}")


sim_time = start_time

# # Animation function
# def animate_trains(frame):
#     # Calculate the current simulation time
#     global sim_time
#     SPEED_FACTOR = speed_slider.val
#     sim_range = end_time - start_time
#     t = start_time + (frame * FRAME_STEP * SPEED_FACTOR) % sim_range
    

def animate_trains(frame):
    global sim_time
    SPEED_FACTOR = speed_slider.val
    sim_time += FRAME_STEP * SPEED_FACTOR  # ✅ Accumulate time forward
    t = (sim_time - start_time) % (end_time - start_time) + start_time



    # t = start_time + frame * FRAME_STEP  # Simulated time in seconds
    # print(f"Frame: {frame}, Time: {t}")  # Debug: Print the current frame and time

    # Clear the previous train positions
    for train_id, data in train_dots.items():
        data['dot'].set_data([], [])  # Clear the train marker

    # Update train positions
    for train_id, data in train_dots.items():
        times = data['times']
        pos = data['positions']

        # Find the segment of the schedule the train is currently on
        idx = np.searchsorted(times, t, side='right') - 1
        if idx < 0 or idx >= len(pos) - 1:
            continue  # Skip if the train is not active at this time

        # Interpolate the train's position
        t0, t1 = times[idx], times[idx + 1]
        p0, p1 = pos[idx], pos[idx + 1]
        alpha = (t - t0) / (t1 - t0) if t1 > t0 else 0
        dist = (1 - alpha) * p0 + alpha * p1
        pt = track.interpolate(dist)

        # Update the train marker position
        data['dot'].set_data([pt.x], [pt.y])  # Pass as lists
        
        


    # Update the clock
    clock_text.set_text(f"Time: {datetime.utcfromtimestamp(t).strftime('%H:%M')}")
   
    if frame == 0:
        ax.legend(loc='upper right')  # or your preferred position

    # Return the artists to update
    return [data['dot'] for data in train_dots.values()] + [clock_text]

# Create the animation
ani = FuncAnimation(
    fig,  # The figure to animate
    animate_trains,  # The update function
    frames=FRAMES,  # Total number of frames
    interval=INTERVAL,  # Delay between frames in milliseconds
    blit=True  # Optimize animation by only redrawing changed parts
)

# Show the plot
plt.legend()
plt.show()


# # Inspect full stop list and times for Train 223
# train_223 = timetable_df[timetable_df['ID'] == 223].sort_values("Departs_s")
# print("\n🚂 Train 223 timetable:")
# print(train_223[["Station", "Departs_s"]])


# print("\n🛠 Stations with mileposts:")
# print(sorted(station_mileposts.keys()))

# raw_223 = pd.read_excel(TIMETABLE_FILE)
# print("\n📋 Raw entries for Train 223:")
# print(raw_223[raw_223["ID"] == 223][["Station", "Depart time"]])





























