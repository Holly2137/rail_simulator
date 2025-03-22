import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from shapely.geometry import LineString
from datetime import datetime
import contextily as ctx



# Load shapefiles
railline = gpd.read_file(r'shapefiles\Sligo_Dublin_Line.shp')
stations = gpd.read_file(r'shapefiles\Sligo_Dublin_Stops2.shp')

# Reproject to EPSG:3857 if necessary
if railline.crs.to_string() != "EPSG:3857":
    railline = railline.to_crs(epsg=3857)
if stations.crs.to_string() != "EPSG:3857":
    stations = stations.to_crs(epsg=3857)

# # Ensure geometry types are correct    Checked Fine
# print(railline.geom_type.unique())  # Should be ['LineString']
# print(stations.geom_type.unique())  # Should be ['Point']

# Clean station names
stations["StopName"] = stations["StopName"].str.replace(" Train Station", "", regex=False)

# Create a shared plot
fig, ax = plt.subplots(figsize=(10, 10))

# Plot the rail lines and stations
railline.plot(ax=ax, color='red', label='Rail Line')
stations.plot(ax=ax, color='green', marker='o',zorder=3, markersize=10, label='Stations')
# Add a basemap
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)
# Add station names as labels
#for idx, row in stations.iterrows():
  #  station_name = row["StopName"]
  #  ax.text(row.geometry.x, row.geometry.y, station_name, fontsize=8, ha='right', color='blue')

# Ensure the rail line is a single LineString
if railline.geometry.iloc[0].geom_type == 'MultiLineString':
    track = railline.geometry.iloc[0].unary_union  # Merge into one LineString
else:
    track = railline.geometry.iloc[0]  # Already a LineString

# Sort station points along the track
stations["distance_along_track"] = stations.geometry.apply(lambda point: track.project(point))
stations = stations.sort_values(by="distance_along_track")

# Load timetable data
timetable_csv = r'rawdata\SligoDublinTimetableMerged02.csv'
timetable_df = pd.read_csv(timetable_csv)

# Convert time to seconds
def time_to_seconds(time_str):
    try:
        time_obj = datetime.strptime(time_str, "%H:%M")
        return time_obj.hour * 3600 + time_obj.minute * 60
    except ValueError:
        return None

# Apply time conversion
timetable_df["Departs_s"] = timetable_df["Departs"].apply(time_to_seconds)
timetable_df["Arrives_s"] = timetable_df["Arrives"].apply(time_to_seconds)

# Standardize station names
timetable_df["From"] = timetable_df["From"].str.strip()
timetable_df["To"] = timetable_df["To"].str.strip()
stations["StopName"] = stations["StopName"].str.strip()

# Fix station naming mismatches
station_name_corrections = {
    "Carrick-on-Shannon": "Carrick on Shannon",
    "Dublin Connolly": "Connolly"
}
timetable_df["From"] = timetable_df["From"].replace(station_name_corrections)
timetable_df["To"] = timetable_df["To"].replace(station_name_corrections)
stations["StopName"] = stations["StopName"].replace(station_name_corrections)

# Merge track positions
timetable_df = timetable_df.merge(
    stations[["StopName", "distance_along_track"]],
    left_on="From",
    right_on="StopName",
    how="left"
).rename(columns={"distance_along_track": "From_Track_Position"})

timetable_df = timetable_df.merge(
    stations[["StopName", "distance_along_track"]],
    left_on="To",
    right_on="StopName",
    how="left"
).rename(columns={"distance_along_track": "To_Track_Position"})

# Sort timetable by departure time
timetable_df = timetable_df.sort_values(by="Departs_s")

pd.set_option('display.max_columns', None)
#print(timetable_df.head())


# # Select a train journey to animate (Example: Train 101 - Sligo to Dublin)
# train_id = 101
# train_journey = timetable_df[timetable_df["ID"] == train_id].sort_values(by="Departs_s")


# # Assign track positions & departure times
# track_positions = train_journey["From_Track_Position"].tolist() + [train_journey["To_Track_Position"].iloc[-1]]
# time_stamps = train_journey["Departs_s"].tolist() + [train_journey["Arrives_s"].iloc[-1]]

# print(track_positions)
# #print(time_stamps)
# print(train_journey[["From", "To", "Departs_s", "Arrives_s", "From_Track_Position", "To_Track_Position"]])

# # Normalize time values (start from zero)
# time_stamps = [t - time_stamps[0] for t in time_stamps]

# # Initialize train dot
# start_geom = track.interpolate(track_positions[0])  # First station's position
# train_dot, = ax.plot([start_geom.x], [start_geom.y], 'go', markersize=8, zorder=10, label="Train")

# # Convert track positions to kilometers
# track_positions_km = [pos / 1000 for pos in track_positions]
# #print("Track Positions (km):", track_positions_km)
# # Calculate total track length in kilometers
# total_length_km = (track_positions[-1] - track_positions[0]) / 1000
# print(f"Total Track Length: {total_length_km} km")

# Speed multiplier (Adjustable)
SPEED_FACTOR = 50  # Increase for faster animation

def update(frame):
    global timetable_df, track, ax
    
    # Scale frame time to move faster
    t = frame * (max_time / FRAMES)
    
    animated_objects = []  # Store train dots for animation
    
    for train_id, train_data in train_dict.items():
        time_stamps = train_data['time_stamps']
        track_positions = train_data['track_positions']
        dwell_times = train_data['dwell_times']
        train_dot = train_data['dot']

        # Find the closest station index
        index = np.searchsorted(time_stamps, t, side="right") - 1
        index = max(0, min(index, len(track_positions) - 2))  # Clamp index

        # Get station times and positions
        t_start, t_end = time_stamps[index], time_stamps[index + 1]
        pos_start, pos_end = track_positions[index], track_positions[index + 1]
        dwell_time = dwell_times[index]  # Get dwell time at the station
        
        # Calculate when train starts moving after dwell time
        t_movement_start = max(t_start, min(t_start + dwell_time, t_end - 1))
        
        if t < t_movement_start:
            # Train is dwelling at station, do not move
            interp_position = pos_start
        else:
            # Train is moving, interpolate its position along the track
            if t_end - t_movement_start > 0:
                alpha = (t - t_movement_start) / (t_end - t_movement_start)
            else:
                alpha = 1  # If arrival time == departure time, stay at position
            
            interp_position = (1 - alpha) * pos_start + alpha * pos_end
        
        pos = track.interpolate(max(0, min(interp_position, track.length - 1)))  # Ensure within track range

        if pos and not np.isnan(interp_position):  # Ensure valid geometry
            train_dot.set_data([pos.x], [pos.y])
        else:
            print(f"🚨 Warning: Invalid position at frame {frame} for train {train_id}")
        
        animated_objects.append(train_dot)  # Collect train dot for animation
    
    # 🕰️ Update clock display
    journey_start_time = datetime.strptime("05:00", "%H:%M")
    real_time_seconds = t  # Start simulation at 5 AM (18000 seconds after midnight)  # Current animation time in seconds
    current_time = (journey_start_time + pd.to_timedelta(real_time_seconds, unit='s')).strftime("%H:%M")
    clock_text.set_text(f"Time: {current_time}")
    animated_objects.append(clock_text)
    
    return animated_objects  # Return all train dots + clock for animation

# Initialize trains from the CSV
train_dict = {}
max_time = timetable_df["Arrives_s"].max() - timetable_df["Departs_s"].min() + 4000  # Extra buffer to prevent early reset

for train_id in timetable_df['ID'].unique():
    train_journey = timetable_df[timetable_df['ID'] == train_id].sort_values(by="Departs_s")
    
    track_positions = train_journey["From_Track_Position"].tolist() + [train_journey["To_Track_Position"].iloc[-1]]
    time_stamps = [t - 18000 for t in train_journey["Departs_s"].tolist()] + [train_journey["Arrives_s"].iloc[-1] - 18000]
    dwell_times = train_journey["Dwell"].tolist() + [0]  # Assume no dwell at the last stop
    
    start_geom = track.interpolate(track_positions[0])
    train_dot, = ax.plot([start_geom.x], [start_geom.y], 'o', markersize=8, label=f"Train {train_id}")
    
    train_dict[train_id] = {
        'time_stamps': time_stamps,
        'track_positions': track_positions,
        'dwell_times': dwell_times,
        'dot': train_dot
    }

# Add clock display to plot
clock_text = ax.text(0.5, 0.95, "Time: 00:00", transform=ax.transAxes, ha='center', fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

# 🎬 Create the animation
FRAMES = 2000  # More frames = smoother animation
INTERVAL = 20   # Lower interval = faster movement
ani = FuncAnimation(fig, update, frames=FRAMES, interval=INTERVAL, blit=True, repeat=True)

# 🖥 Show the plot
plt.legend()
plt.show()

