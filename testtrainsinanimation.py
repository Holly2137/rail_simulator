import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import csv
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from matplotlib.widgets import Button
from shapely.geometry import LineString
from shapely.ops import substring
from datetime import datetime
from datetime import datetime, timedelta
from datetime import timezone
import contextily as ctx

# Define base directory dynamically from the script location
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Define relative paths
RAW_DATA_DIR = os.path.join(BASE_DIR, "rawdata")
SHAPEFILE_DIR = os.path.join(BASE_DIR, "shapefiles")

# Paths to data files
SEGMENTS_FILE = os.path.join(RAW_DATA_DIR, "Segments.xlsx")
TIMETABLE_FILE = os.path.join(RAW_DATA_DIR, "Timetable2.xlsx")
TEST_TRAIN_FILE = os.path.join(RAW_DATA_DIR, "testtrains2.csv")
LINE_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Line.shp")
STOPS_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Stops2.shp")

# Load timetable and segments
segments_df = pd.read_excel(SEGMENTS_FILE)
timetable_df = pd.read_excel(TIMETABLE_FILE)
test_trains_df = pd.read_csv(TEST_TRAIN_FILE)


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





def build_segment_occupancy(timetable_df: pd.DataFrame) -> dict:
    """
    Build occupancy mapping for each train.
    Mapping format:
      { (start_station, end_station): [ { "train": train_id, "start_time": <seconds>, "end_time": <seconds> }, ... ] }
    Occupancy is only added between stations with Passing flag "Y" or "END" and valid times.
    """
    segment_occupancy = {}

    # Ensure Passing column is standardized
    timetable_df["Passing"] = timetable_df["Passing"].astype(str).str.strip().str.upper()

    # Compute arrival time as departure minus dwell
    timetable_df["Minimum Dwell"] = pd.to_numeric(timetable_df["Minimum Dwell"], errors="coerce").fillna(60).astype(int)
    timetable_df["arrivals"] = timetable_df["Departs_s"] - timetable_df["Minimum Dwell"]

    # Group by train ID
    for train_id, group in timetable_df.groupby("ID"):
        group = group.sort_values("Departs_s")
        last_passing = None

        for _, row in group.iterrows():
            if row["Passing"] in ["Y", "END"]:
                if last_passing is not None:
                    start_time = last_passing.get("Departs_s")
                    end_time = row.get("arrivals")
                    if pd.notnull(start_time) and pd.notnull(end_time):
                        seg_key = (last_passing["Station"].strip(), row["Station"].strip())
                        segment_occupancy.setdefault(seg_key, []).append({
                            "train": train_id,
                            "start_time": start_time,
                            "end_time": end_time
                        })
                last_passing = row
    return segment_occupancy

# ✅ Now build the occupancy mapping
segment_occupancy = build_segment_occupancy(timetable_df)

# ✅ Optional: Preview a few entries
print("Sample segment occupancy entries:")
for i, (seg, occs) in enumerate(segment_occupancy.items()):
    print(f"{seg} → {len(occs)} trains")
    for occ in occs:
        start_str = datetime.utcfromtimestamp(occ['start_time']).strftime('%H:%M')
        end_str = datetime.utcfromtimestamp(occ['end_time']).strftime('%H:%M')
        print(f"   Train {occ['train']} from {start_str} to {end_str}")
    if i >= 3:  # limit preview
        break

#good to here

segments_df["From"] = segments_df["From"].astype(str).str.strip()
segments_df["To"] = segments_df["To"].astype(str).str.strip()


def simulate_test_train(train_id: str, test_trains_df: pd.DataFrame, 
                        segments_df: pd.DataFrame, segment_occupancy: dict) -> pd.DataFrame:
    segments_df["From"] = segments_df["From"].astype(str).str.strip()
    segments_df["To"] = segments_df["To"].astype(str).str.strip()

    df = test_trains_df[test_trains_df["ID"] == train_id].copy()

    df["Depart time"] = pd.to_datetime(df["Depart time"], format="%H:%M:%S").dt.time
    df["Departs_s"] = df["Depart time"].apply(lambda t: t.hour * 3600 + t.minute * 60 + t.second)
    df["Minimum Dwell"] = pd.to_numeric(df["Minimum Dwell"], errors="coerce").fillna(60).astype(int)
    df["arrivals"] = df["Departs_s"] - df["Minimum Dwell"]

    results = []
    current_time = df.iloc[0]["Departs_s"]
    last_passing_row = None
    last_depart_time = current_time

    for idx, row in df.iterrows():
        station = row["Station"].strip()
        dwell = row["Minimum Dwell"]
        passing = str(row["Passing"]).strip().upper()

        arrival_time = current_time
        depart_time = current_time + dwell

        results.append({
            "Station": station,
            "Arrival_s": arrival_time,
            "Arrival_str": datetime.fromtimestamp(arrival_time, tz=timezone.utc).strftime('%H:%M:%S'),
            "Depart_s": depart_time,
            "Depart_str": datetime.fromtimestamp(depart_time, tz=timezone.utc).strftime('%H:%M:%S'),
            "Delay": 0
        })

        if passing in ["Y", "END"]:
            if last_passing_row is not None:
                from_station = last_passing_row["Station"].strip()
                to_station = station

                from_idx = df.index.get_loc(last_passing_row.name)
                to_idx = df.index.get_loc(row.name)
                intermediate_rows = df.iloc[from_idx+1:to_idx+1]

                path_valid = True
                segment_sequence = []

                for i in range(len(intermediate_rows)):
                    s1 = df.iloc[from_idx + i]["Station"].strip()
                    s2 = df.iloc[from_idx + i + 1]["Station"].strip()
                    seg = segments_df[
                        ((segments_df["From"] == s1) & (segments_df["To"] == s2)) |
                        ((segments_df["From"] == s2) & (segments_df["To"] == s1))
                    ]
                    if seg.empty:
                        print(f"⚠️  Segment missing in segments_df: '{s1}' <-> '{s2}'")
                        path_valid = False
                        break
                    segment_sequence.append((s1, s2))

                if not path_valid:
                    midpoint = "Carrick loop"
                    seg1 = segments_df[
                        ((segments_df["From"] == from_station) & (segments_df["To"] == midpoint)) |
                        ((segments_df["From"] == midpoint) & (segments_df["To"] == from_station))
                    ]
                    seg2 = segments_df[
                        ((segments_df["From"] == midpoint) & (segments_df["To"] == to_station)) |
                        ((segments_df["From"] == to_station) & (segments_df["To"] == midpoint))
                    ]
                    if not seg1.empty and not seg2.empty:
                        segment_sequence = [(from_station, to_station)]
                    else:
                        print(f"⚠️  Segment path via midpoint '{midpoint}' missing: '{from_station}' <-> '{to_station}'")
                        raise ValueError(f"Segment path not found: {from_station} -> {to_station} via test train path")

                proposed_arrival = row["arrivals"]

                segment_key = (from_station, to_station) if (from_station, to_station) in segment_occupancy else (to_station, from_station)
                occs = segment_occupancy.get(segment_key, [])
                for occ in occs:
                    if not (proposed_arrival <= occ["start_time"] or last_depart_time >= occ["end_time"]):
                        print(f"🚧 Conflict on segment {segment_key} with train {occ['train']} from {occ['start_time']} to {occ['end_time']}")
                        last_depart_time = occ["end_time"]
                        proposed_arrival = last_depart_time

                current_time = proposed_arrival
                depart_time = current_time + dwell

                results[-1]["Arrival_s"] = current_time
                results[-1]["Depart_s"] = depart_time
                results[-1]["Arrival_str"] = datetime.fromtimestamp(current_time, tz=timezone.utc).strftime('%H:%M:%S')
                results[-1]["Depart_str"] = datetime.fromtimestamp(depart_time, tz=timezone.utc).strftime('%H:%M:%S')
                results[-1]["Delay"] = max(0, current_time - row["arrivals"])

                last_depart_time = depart_time
            last_passing_row = row
        else:
            current_time = depart_time

    return pd.DataFrame(results)

def print_simulation_result(df: pd.DataFrame, train_id: str):
    if df.empty:
        print(f"No simulation results for Test Train {train_id}")
        return

    print(f"\nSimulation Result for Test Train {train_id}")
    print("-" * 70)
    print(f"{'Station':25} | {'Arrival':8} | {'Depart':8} | {'Delay (s)':>10}")
    print("-" * 70)

    for _, row in df.iterrows():
        delay = row.get("Delay", 0)
        print(f"{row['Station']:25} | {row['Arrival_str']} | {row['Depart_str']} | {int(delay):>10}")

    total_time = df.iloc[-1]['Depart_s'] - df.iloc[0]['Arrival_s']
    print("-" * 70)
    print(f"Total Journey Time: {int(total_time // 60)} min {int(total_time % 60)} sec\n")


# Run the test
t1_result = simulate_test_train("T1", test_trains_df, segments_df, segment_occupancy)
print_simulation_result(t1_result, "T1")























# def simulate_test_train(train_id: str, test_trains_df: pd.DataFrame, 
#                         segments_df: pd.DataFrame, segment_occupancy: dict) -> pd.DataFrame:
#     segments_df["From"] = segments_df["From"].astype(str).str.strip()
#     segments_df["To"] = segments_df["To"].astype(str).str.strip()

#     df = test_trains_df[test_trains_df["ID"] == train_id].copy()

#     df["Depart time"] = pd.to_datetime(df["Depart time"], format="%H:%M:%S").dt.time
#     df["Departs_s"] = df["Depart time"].apply(lambda t: t.hour * 3600 + t.minute * 60 + t.second)
#     df["Minimum Dwell"] = pd.to_numeric(df["Minimum Dwell"], errors="coerce").fillna(60).astype(int)
#     df["arrivals"] = df["Departs_s"] - df["Minimum Dwell"]

#     results = []
#     current_time = df.iloc[0]["Departs_s"]
#     last_passing_row = None
#     last_depart_time = current_time

#     for idx, row in df.iterrows():
#         station = row["Station"].strip()
#         dwell = row["Minimum Dwell"]
#         passing = str(row["Passing"]).strip().upper()

#         arrival_time = current_time
#         depart_time = current_time + dwell

#         results.append({
#             "Station": station,
#             "Arrival_s": arrival_time,
#             "Arrival_str": datetime.fromtimestamp(arrival_time, tz=timezone.utc).strftime('%H:%M:%S'),
#             "Depart_s": depart_time,
#             "Depart_str": datetime.fromtimestamp(depart_time, tz=timezone.utc).strftime('%H:%M:%S'),
#             "Delay": 0
#         })

#         if passing in ["Y", "END"]:
#             if last_passing_row is not None:
#                 from_station = last_passing_row["Station"].strip()
#                 to_station = station

#                 from_idx = df.index.get_loc(last_passing_row.name)
#                 to_idx = df.index.get_loc(row.name)
#                 intermediate_rows = df.iloc[from_idx+1:to_idx+1]

#                 path_valid = True
#                 segment_sequence = []

#                 for i in range(len(intermediate_rows)):
#                     s1 = df.iloc[from_idx + i]["Station"].strip()
#                     s2 = df.iloc[from_idx + i + 1]["Station"].strip()
#                     seg = segments_df[
#                         ((segments_df["From"] == s1) & (segments_df["To"] == s2)) |
#                         ((segments_df["From"] == s2) & (segments_df["To"] == s1))
#                     ]
#                     if seg.empty:
#                         print(f"⚠️  Segment missing in segments_df: '{s1}' <-> '{s2}'")
#                         path_valid = False
#                         break
#                     segment_sequence.append((s1, s2))

#                 if not path_valid:
#                     midpoint = "Carrick loop"
#                     seg1 = segments_df[
#                         ((segments_df["From"] == from_station) & (segments_df["To"] == midpoint)) |
#                         ((segments_df["From"] == midpoint) & (segments_df["To"] == from_station))
#                     ]
#                     seg2 = segments_df[
#                         ((segments_df["From"] == midpoint) & (segments_df["To"] == to_station)) |
#                         ((segments_df["From"] == to_station) & (segments_df["To"] == midpoint))
#                     ]
#                     if not seg1.empty and not seg2.empty:
#                        # this is where we edited unsure
#                         segment_sequence = [(from_station, to_station)]  # for conflict checking

#                     else:
#                         print(f"⚠️  Segment path via midpoint '{midpoint}' missing: '{from_station}' <-> '{to_station}'")
#                         raise ValueError(f"Segment path not found: {from_station} -> {to_station} via test train path")

#                 proposed_arrival = row["arrivals"]

#                 segment_key = (from_station, to_station) if (from_station, to_station) in segment_occupancy else (to_station, from_station)
#                 occs = segment_occupancy.get(segment_key, [])
#                 for occ in occs:
#                     if not (proposed_arrival <= occ["start_time"] or last_depart_time >= occ["end_time"]):
#                         last_depart_time = occ["end_time"]
#                         proposed_arrival = last_depart_time

#                 current_time = proposed_arrival
#                 depart_time = current_time + dwell

#                 results[-1]["Arrival_s"] = current_time
#                 results[-1]["Depart_s"] = depart_time
#                 results[-1]["Arrival_str"] = datetime.fromtimestamp(current_time, tz=timezone.utc).strftime('%H:%M:%S')
#                 results[-1]["Depart_str"] = datetime.fromtimestamp(depart_time, tz=timezone.utc).strftime('%H:%M:%S')
#                 results[-1]["Delay"] = max(0, current_time - row["arrivals"])

#                 last_depart_time = depart_time
#             last_passing_row = row
#         else:
#             current_time = depart_time

#     return pd.DataFrame(results)

# def print_simulation_result(df: pd.DataFrame, train_id: str):
#     if df.empty:
#         print(f"No simulation results for Test Train {train_id}")
#         return

#     print(f"\nSimulation Result for Test Train {train_id}")
#     print("-" * 70)
#     print(f"{'Station':25} | {'Arrival':8} | {'Depart':8} | {'Delay (s)':>10}")
#     print("-" * 70)

#     for _, row in df.iterrows():
#         delay = row.get("Delay", 0)
#         print(f"{row['Station']:25} | {row['Arrival_str']} | {row['Depart_str']} | {int(delay):>10}")

#     total_time = df.iloc[-1]['Depart_s'] - df.iloc[0]['Arrival_s']
#     print("-" * 70)
#     print(f"Total Journey Time: {int(total_time // 60)} min {int(total_time % 60)} sec\n")


# # Run the test
# t1_result = simulate_test_train("T1", test_trains_df, segments_df, segment_occupancy)
# print_simulation_result(t1_result, "T1")






















































# def simulate_test_train(train_id: str, test_trains_df: pd.DataFrame, 
#                         segments_df: pd.DataFrame, segment_occupancy: dict) -> pd.DataFrame:
#     df = test_trains_df[test_trains_df["ID"] == train_id].copy()

#     df["Depart time"] = pd.to_datetime(df["Depart time"], format="%H:%M:%S").dt.time
#     df["Departs_s"] = df["Depart time"].apply(lambda t: t.hour * 3600 + t.minute * 60 + t.second)
#     df["Minimum Dwell"] = pd.to_numeric(df["Minimum Dwell"], errors="coerce").fillna(60).astype(int)
#     df["arrivals"] = df["Departs_s"] - df["Minimum Dwell"]

#     results = []
#     current_time = df.iloc[0]["Departs_s"]
#     last_passing_row = None
#     last_depart_time = current_time

#     for idx, row in df.iterrows():
#         station = row["Station"].strip()
#         dwell = row["Minimum Dwell"]
#         passing = str(row["Passing"]).strip().upper()

#         # Store result row for all stations (regardless of passing status)
#         arrival_time = current_time
#         depart_time = current_time + dwell

#         results.append({
#             "Station": station,
#             "Arrival_s": arrival_time,
#             "Arrival_str": datetime.fromtimestamp(arrival_time, tz=timezone.utc).strftime('%H:%M:%S'),
#             "Depart_s": depart_time,
#             "Depart_str": datetime.fromtimestamp(depart_time, tz=timezone.utc).strftime('%H:%M:%S'),
#             "Delay": 0
#         })

#         if passing in ["Y", "END"]:
#             if last_passing_row is not None:
#                 from_station = last_passing_row["Station"].strip()
#                 to_station = station

#                 segment_info = segments_df[
#                     ((segments_df["From"].str.strip() == from_station) &
#                      (segments_df["To"].str.strip() == to_station)) |
#                     ((segments_df["From"].str.strip() == to_station) &
#                      (segments_df["To"].str.strip() == from_station))
#                 ]

#                 if segment_info.empty:
#                     # Try a 2-hop lookup via known virtual midpoints
#                     virtual_midpoints = {"Carrick loop", "Ballymote", "Collooney"}
#                     found = False
#                     for mid in virtual_midpoints:
#                         part1 = segments_df[
#                             ((segments_df["From"].str.strip() == from_station) & (segments_df["To"].str.strip() == mid)) |
#                             ((segments_df["From"].str.strip() == mid) & (segments_df["To"].str.strip() == from_station))
#                         ]
#                         part2 = segments_df[
#                             ((segments_df["From"].str.strip() == mid) & (segments_df["To"].str.strip() == to_station)) |
#                             ((segments_df["From"].str.strip() == to_station) & (segments_df["To"].str.strip() == mid))
#                         ]
#                         if not part1.empty and not part2.empty:
#                             tight_run = int(part1.iloc[0]["Tight Run"]) + int(part2.iloc[0]["Tight Run"])
#                             is_single = (
#                                 part1.iloc[0]["Track"].strip().lower() == "single" or
#                                 part2.iloc[0]["Track"].strip().lower() == "single"
#                             )
#                             found = True
#                             break
#                     if not found:
#                         raise ValueError(f"Segment not found: {from_station} <-> {to_station}")
#                 else:
#                     tight_run = int(segment_info.iloc[0]["Tight Run"])
#                     is_single = segment_info.iloc[0]["Track"].strip().lower() == "single"

#                 proposed_arrival = last_depart_time + tight_run

#                 if is_single:
#                     segment_key = (
#                         (from_station, to_station)
#                         if (from_station, to_station) in segment_occupancy
#                         else (to_station, from_station)
#                     )
#                     occs = segment_occupancy.get(segment_key, [])
#                     for occ in occs:
#                         if not (proposed_arrival <= occ["start_time"] or last_depart_time >= occ["end_time"]):
#                             last_depart_time = occ["end_time"]
#                             proposed_arrival = last_depart_time + tight_run

#                 # Update times for this passing station
#                 current_time = proposed_arrival
#                 depart_time = current_time + dwell
#                 results[-1]["Arrival_s"] = current_time
#                 results[-1]["Depart_s"] = depart_time
#                 results[-1]["Arrival_str"] = datetime.fromtimestamp(current_time, tz=timezone.utc).strftime('%H:%M:%S')
#                 results[-1]["Depart_str"] = datetime.fromtimestamp(depart_time, tz=timezone.utc).strftime('%H:%M:%S')
#                 results[-1]["Delay"] = max(0, current_time - row["arrivals"])

#                 last_depart_time = depart_time
#             last_passing_row = row
#         else:
#             # For non-passing stations, just advance time by dwell
#             current_time = depart_time

#     return pd.DataFrame(results)

# def print_simulation_result(df: pd.DataFrame, train_id: str):
#     if df.empty:
#         print(f"No simulation results for Test Train {train_id}")
#         return

#     print(f"\nSimulation Result for Test Train {train_id}")
#     print("-" * 70)
#     print(f"{'Station':25} | {'Arrival':8} | {'Depart':8} | {'Delay (s)':>10}")
#     print("-" * 70)

#     for _, row in df.iterrows():
#         delay = row.get("Delay", 0)
#         print(f"{row['Station']:25} | {row['Arrival_str']} | {row['Depart_str']} | {int(delay):>10}")

#     total_time = df.iloc[-1]['Depart_s'] - df.iloc[0]['Arrival_s']
#     print("-" * 70)
#     print(f"Total Journey Time: {int(total_time // 60)} min {int(total_time % 60)} sec\n")


# # Run the test
# t1_result = simulate_test_train("T1", test_trains_df, segments_df, segment_occupancy)
# print_simulation_result(t1_result, "T1")


# Example usage
# t1_result = simulate_test_train("T1", test_trains_df, segments_df, segment_occupancy)
# print(t1_result)













# def save_segment_occupancy_csv(segment_occupancy: dict, base_dir: str):
#     output_dir = os.path.join(base_dir, "printouts")
#     os.makedirs(output_dir, exist_ok=True)

#     output_file = os.path.join(output_dir, "segment_occupancy.csv")

#     with open(output_file, mode="w", newline="") as f:
#         writer = csv.writer(f)
#         writer.writerow(["From", "To", "Train", "Start Time", "End Time"])

#         for (from_station, to_station), occ_list in segment_occupancy.items():
#             for occ in occ_list:
#                 writer.writerow([
#                     from_station,
#                     to_station,
#                     occ["train"],
#                     datetime.utcfromtimestamp(occ["start_time"]).strftime('%H:%M'),
#                     datetime.utcfromtimestamp(occ["end_time"]).strftime('%H:%M')
#                 ])

#     print(f"Segment occupancy saved to: {output_file}")

# save_segment_occupancy_csv(segment_occupancy, BASE_DIR)






























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


#speed slider
ax_slider = plt.axes([0.4, 0.05, 0.2, 0.07])
speed_slider = Slider(ax_slider, 'Speed Factor', valmin=1, valmax=10, valinit=3, valstep=1)

clock_text = ax.text(0.5, 0.95, "Time: 00:00", transform=ax.transAxes, ha='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.7))


# Global pause state
is_paused = False

# Add pause button
ax_pause = plt.axes([0.15, 0.02, 0.1, 0.05])  # [left, bottom, width, height]
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

#print("First 10 entries for Train", train_id)
#for i in range(min(10, len(times))):
 #   print(f"Station {i+1}: Depart={times[i]}, Arrival={arrivals[i]}, Dwell={actual_dwell[i]}")
start_time = timetable_df["Departs_s"].min() - 1000
end_time = timetable_df["Departs_s"].max() + 1000

FRAME_STEP = 5
FRAMES = 8000
INTERVAL = 15

sim_time = start_time




# # Time slider axes
# ax_time_slider = plt.axes([0.65, 0.05, 0.25, 0.03])  # [left, bottom, width, height]
# time_slider = Slider(ax_time_slider, 'Time', valmin=start_time, valmax=end_time, valinit=start_time, valstep=60)

# def seconds_to_timestr(seconds):
#     return datetime.utcfromtimestamp(seconds).strftime('%H:%M')

# # Update sim_time from slider
# def on_time_slider(val):
#     global sim_time, is_paused
#     is_paused = True
#     pause_button.label.set_text('Resume')
#     sim_time = val
#     fig.canvas.draw_idle()

# time_slider.on_changed(on_time_slider)





def animate_trains(frame):
   

    global sim_time, last_artists
    if is_paused:
        return last_artists  # return previous frame's state without change

    # rest of your animate_trains logic remains the same

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

        first_depart = departs[0]
        if first_depart - 300 <= t < first_depart:
            pos_pt = track.interpolate(positions[0])
            dot.set_data([pos_pt.x], [pos_pt.y])
            updated_artists.append(dot)
            continue

        # Final hold at last station
        final_depart = departs[-1]
        final_pos = positions[-1]
        if final_depart <= t <= final_depart + 300:
            pos_pt = track.interpolate(final_pos)
            dot.set_data([pos_pt.x], [pos_pt.y])
            updated_artists.append(dot)
            continue

        if t > final_depart + 300:
            dot.set_data([], [])
            continue

        idx = np.searchsorted(departs, t, side='right') - 1
        if idx < 0 or idx >= len(positions) - 1:
            dot.set_data([], [])
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
            continue

        pos_pt = track.interpolate(dist)
        dot.set_data([pos_pt.x], [pos_pt.y])
        updated_artists.append(dot)

    clock_text.set_text(f"Time: {datetime.utcfromtimestamp(t).strftime('%H:%M')}")
    updated_artists.append(clock_text)

    if frame == 0:
        ax.legend(loc='upper right')

    last_artists = updated_artists  # store for paused state
    return updated_artists


ani = FuncAnimation(fig, animate_trains, frames=FRAMES, interval=INTERVAL, blit=True)
plt.legend()
plt.show()
















