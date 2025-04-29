import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

# Define base directory dynamically from the script location
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Define relative paths
RAW_DATA_DIR = os.path.join(BASE_DIR, "rawdata")
TIMETABLE_FILE = os.path.join(RAW_DATA_DIR, "Timetable2manual.xlsx")

# Load timetable
timetable_df = pd.read_excel(TIMETABLE_FILE)

# Clean station names and times
timetable_df["Station"] = timetable_df["Station"].astype(str).str.strip()

# Time conversion function
def time_to_seconds(tstr):
    for fmt in ("%H:%M", "%H:%M:%S"):
        try:
            t = datetime.strptime(str(tstr).strip(), fmt)
            return t.hour * 3600 + t.minute * 60 + t.second
        except:
            continue
    return np.nan

# Apply time conversion
timetable_df["Departs_s"] = timetable_df["Depart time"].apply(time_to_seconds)

# Filter necessary columns and drop invalid entries
timetable_df = timetable_df.dropna(subset=["Station", "Departs_s", "Milepost(m)"])
timetable_df = timetable_df[["ID", "Station", "Departs_s", "Milepost(m)", "Minimum Dwell"]]

# Calculate arrival times based on dwell
# Ensure dwell is numeric and fill NaNs
timetable_df["Minimum Dwell"] = pd.to_numeric(timetable_df["Minimum Dwell"], errors="coerce").fillna(0)
timetable_df["Arrives_s"] = timetable_df["Departs_s"] - timetable_df["Minimum Dwell"]

# Get station order and mapping to km for y-axis labels
station_order = timetable_df.drop_duplicates("Station").sort_values("Milepost(m)", ascending=True)
station_to_km = dict(zip(station_order["Station"], station_order["Milepost(m)"] / 1000))

station_ticks = list(station_to_km.values())
station_labels = list(station_to_km.keys())

# Manual label offsets (labels only, not data)
manual_offsets = {
    "Carrick loop": 3.0,
    "Connolly": -1.0,
    "Broombridge": 1.0
}



# --- PLOTTING TRAIN GRAPH ---
fig, ax = plt.subplots(figsize=(14, 8))

# Group by train ID
grouped = timetable_df.groupby("ID")

for train_id, group in grouped:
    group = group.sort_values(by="Departs_s").reset_index(drop=True)
    # colour = "red" if str(train_id).startswith("T") else "black"

    train_id_str = str(train_id)
    if train_id_str.startswith("TA"):
        colour = "red"
    elif train_id_str.startswith("T"):
        colour = "darkgoldenrod"
    else:
        colour = "black"






    for i in range(1, len(group)):
        prev_row = group.iloc[i - 1]
        curr_row = group.iloc[i]

        # First segment: from previous departure to current arrival
        ax.plot([
            prev_row["Departs_s"] / 3600,
            curr_row["Arrives_s"] / 3600
        ], [
            prev_row["Milepost(m)"] / 1000,
            curr_row["Milepost(m)"] / 1000
        ], color=colour, linewidth=1)

        # Second segment: dwell (horizontal)
        if curr_row["Minimum Dwell"] > 0:
            ax.plot([
                curr_row["Arrives_s"] / 3600,
                curr_row["Departs_s"] / 3600
            ], [
                curr_row["Milepost(m)"] / 1000,
                curr_row["Milepost(m)"] / 1000
            ], color=colour, linewidth=1)

    # Label at start and end
    if not group.empty:
        start_label = group.iloc[0]
        end_label = group.iloc[-1]
        ax.text(start_label["Departs_s"] / 3600, start_label["Milepost(m)"] / 1000, str(train_id),
                fontsize=7, verticalalignment='bottom', horizontalalignment='right', color=colour)
        ax.text(end_label["Departs_s"] / 3600, end_label["Milepost(m)"] / 1000, str(train_id),
                fontsize=7, verticalalignment='top', horizontalalignment='left', color=colour)

# X-axis time formatting
ax.set_xlim(5, 23)  # 5 AM to 11 PM
ax.set_xticks(range(5, 24))
ax.set_xticklabels([f"{h}:00" for h in range(5, 24)])

# Y-axis ticks and adjusted labels
ax.set_yticks(station_ticks)
ax.set_yticklabels(["" for _ in station_labels])
for label, tick in zip(station_labels, station_ticks):
    offset = manual_offsets.get(label, 0)
    ax.text(4.95, tick + offset, label, ha='right', va='center', fontsize=8)

# Secondary y-axis for km labels with same offsets
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(station_ticks)
ax2.set_yticklabels(["" for _ in station_labels])
for label, tick in zip(station_labels, station_ticks):
    offset = manual_offsets.get(label, 0)
    ax2.text(ax.get_xlim()[1] + 0.1, tick + offset, f"{station_to_km[label]:.1f} km", ha='left', va='center', fontsize=8)

ax2.set_ylabel("Distance (km from Dublin)")
ax.set_xlabel("Time of Day (hrs)")
ax.set_ylabel("Station")
ax.set_title("Route 7 Operational and Test Services (Stringline Diagram)")
ax.yaxis.set_label_coords(-0.09, 0.5)     # Pushes "Station" to the left
ax2.yaxis.set_label_coords(1.07, 0.5)    # Pushes "Distance (km from Dublin)" further right
ax.grid(True)
plt.tight_layout()
plt.show()
