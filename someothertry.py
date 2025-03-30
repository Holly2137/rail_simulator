import os
import pandas as pd
import numpy as np
from datetime import timedelta
from collections import deque

# --- File Paths ---
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(BASE_DIR, "rawdata")

SEGMENTS_FILE = os.path.join(DATA_DIR, "Segments3.xlsx")
TIMETABLE_FILE = os.path.join(DATA_DIR, "Timetable4.xlsx")

# --- Load Data ---
segments_df = pd.read_excel(SEGMENTS_FILE)
timetable_df = pd.read_excel(TIMETABLE_FILE)

# Convert Depart time to seconds since midnight
from datetime import datetime as dt

def time_to_seconds(t):
    try:
        return int(dt.strptime(str(t).strip(), "%H:%M").hour * 3600 + dt.strptime(str(t).strip(), "%H:%M").minute * 60)
    except:
        return np.nan

timetable_df['secafter'] = timetable_df['Depart time'].apply(time_to_seconds)

print("Timetable converted to seconds.")

# Clean station names (strip whitespace)
timetable_df["Station"] = timetable_df["Station"].astype(str).str.strip()
segments_df["From"] = segments_df["From"].astype(str).str.strip()
segments_df["To"] = segments_df["To"].astype(str).str.strip()

# --- Build milepost dictionary from segments ---
segments_df["From_milepost"] = segments_df["From Milepost"]
segments_df["To_milepost"] = segments_df["To Milepost"]

station_mileposts = pd.concat([
    segments_df[["From", "From_milepost"]].rename(columns={"From": "Station", "From_milepost": "Milepost"}),
    segments_df[["To", "To_milepost"]].rename(columns={"To": "Station", "To_milepost": "Milepost"})
]).dropna().drop_duplicates(subset='Station').set_index("Station")["Milepost"].to_dict()

print("Station mileposts loaded.")

# --- Build Passing Blocks ---
def get_passing_blocks(timetable_df):
    passing_stations = timetable_df[timetable_df["Passing"].str.upper() == "Y"]["Station"].drop_duplicates().tolist()
    ordered_stations = timetable_df.drop_duplicates(subset="Station")["Station"].tolist()
    passing_blocks = []
    last_pass = None
    for station in ordered_stations:
        if station in passing_stations:
            if last_pass is not None:
                passing_blocks.append((last_pass, station))
            last_pass = station
    return passing_blocks

# --- Build segment usage based on passing blocks ---
def build_segment_usage_blocks(timetable_df, train_ids, passing_blocks):
    usage = []
    for train_id in train_ids:
        train_df = timetable_df[timetable_df["ID"] == train_id].sort_values("secafter").reset_index(drop=True)
        for from_station, to_station in passing_blocks:
            try:
                t_dep = train_df.loc[train_df["Station"] == from_station, "secafter"].values[0]
                dwell = train_df.loc[train_df["Station"] == to_station, "Minimum Dwell"].values[0]
                t_arr = train_df.loc[train_df["Station"] == to_station, "secafter"].values[0] - dwell
                usage.append({
                    "TrainID": train_id,
                    "BlockFrom": from_station,
                    "BlockTo": to_station,
                    "EnterTime": int(t_dep),
                    "ExitTime": int(t_arr)
                })
            except:
                continue
    return pd.DataFrame(usage)

# Build and show usage for trains 101 and 102
passing_blocks = get_passing_blocks(timetable_df)
segment_usage_df = build_segment_usage_blocks(timetable_df, train_ids=[101, 102], passing_blocks=passing_blocks)

print("\n--- Segment Usage by Passing Blocks ---")
print(segment_usage_df.head(20))


# --- Simulate test train ---
def simulate_test_train(start_time, t_template, passing_blocks, segment_usage_df, timetable_df):
    current_time = start_time
    delay_log = []
    print(f"\nSimulating test train starting at {start_time} seconds")

    for from_station, to_station in passing_blocks:
        block_dwell = t_template.loc[t_template["Station"] == to_station, "Minimum Dwell"].values[0]
        time_offset = t_template.loc[t_template["Station"] == to_station, "secafter"].values[0] - \
                      t_template.loc[t_template["Station"] == from_station, "secafter"].values[0]

        planned_arrival = current_time + time_offset
        planned_departure = planned_arrival + block_dwell

        while True:
            conflict = segment_usage_df[
                ((segment_usage_df["BlockFrom"] == from_station) & (segment_usage_df["BlockTo"] == to_station)) &
                (segment_usage_df["EnterTime"] < planned_departure) &
                (segment_usage_df["ExitTime"] > current_time)
            ]
            if conflict.empty:
                print(f"  ✓ {from_station} → {to_station}: clear [{current_time} → {planned_departure}]")
                break
            print(f"  ⏸ {from_station} → {to_station} blocked at {current_time}, waiting...")
            delay_log.append((from_station, to_station, current_time, "WAIT"))
            current_time += 60  # Wait 1 minute

        current_time = planned_departure

    print(f"Test train completed in {current_time - start_time} seconds")
    return current_time - start_time, delay_log






























































