# segment_occupancy.py
import os
import pandas as pd
import numpy as np
from datetime import datetime

# Define base directory dynamically from the script location
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

# Define relative paths
RAW_DATA_DIR = os.path.join(BASE_DIR, "rawdata")
SHAPEFILE_DIR = os.path.join(BASE_DIR, "shapefiles")

# Paths to data files
SEGMENTS_FILE = os.path.join(RAW_DATA_DIR, "Segments2.xlsx")
TIMETABLE_FILE = os.path.join(RAW_DATA_DIR, "Timetable2.xlsx")

# --- Helper: Time parsing ---
def time_to_seconds(tstr):
    for fmt in ("%H:%M", "%H:%M:%S"):
        try:
            t = datetime.strptime(str(tstr).strip(), fmt)
            return t.hour * 3600 + t.minute * 60 + t.second
        except:
            continue
    return np.nan

# --- Load data ---
segments_df = pd.read_excel(SEGMENTS_FILE)
timetable_df = pd.read_excel(TIMETABLE_FILE)

# --- Clean and preprocess timetable ---
timetable_df["Passing"] = timetable_df["Passing"].astype(str).str.strip().str.upper()
timetable_df["Stop"] = timetable_df["Stop"].astype(str).str.strip().str.upper()
timetable_df["IsStop"] = timetable_df["Stop"].isin(["START", "Y", "END"])
timetable_df["Departs_s"] = timetable_df["Depart time"].apply(time_to_seconds)
timetable_df["Station"] = timetable_df["Station"].str.strip()
timetable_df = timetable_df.dropna(subset=["Station", "Departs_s"])

# --- Build segment occupancy records ---
records = []

for train_id in timetable_df['ID'].unique():
    train_data = timetable_df[timetable_df['ID'] == train_id].sort_values("Departs_s")
    stops = train_data["Station"].tolist()
    dep_times = train_data["Departs_s"].tolist()
    stop_flags = train_data["IsStop"].tolist()
    pass_flags = train_data["Passing"].tolist()

    for i in range(len(stops) - 1):
        from_station = stops[i]
        to_station = stops[i + 1]

        segment_row = segments_df[
            ((segments_df["From"] == from_station) & (segments_df["To"] == to_station)) |
            ((segments_df["From"] == to_station) & (segments_df["To"] == from_station))
        ]

        if segment_row.empty:
            continue

        segment_row = segment_row.iloc[0]
        track_type = segment_row["Track"].strip()
        tight_run = segment_row["Tight Run"]

        from_stop_flag = stop_flags[i]
        to_stop_flag = stop_flags[i + 1]

        dep_time = dep_times[i]
        arr_time = dep_times[i + 1]
        duration = arr_time - dep_time

        allow_passing = (
            track_type.lower() == "double" or
            pass_flags[i] == "Y" or
            pass_flags[i + 1] == "Y"
        )

        is_variable = (train_data.iloc[i]["Stop"] == "V" or train_data.iloc[i + 1]["Stop"] == "V")

        records.append({
            "TrainID": train_id,
            "SegmentFrom": from_station,
            "SegmentTo": to_station,
            "TrackType": track_type,
            "AllowPassing": allow_passing,
            "DepartTime": dep_time,
            "ArrivalTime": arr_time,
            "SegmentDuration": duration,
            "TightRun": tight_run,
            "IsStopAtFrom": from_stop_flag,
            "IsStopAtTo": to_stop_flag,
            "IsVariablePassThrough": is_variable
        })

segment_occupancy_df = pd.DataFrame(records)

def print_full_dataframe(df, rows=20):
    import pandas as pd
    with pd.option_context(
        'display.max_rows', rows,
        'display.max_columns', None,
        'display.width', 2000,
        'display.max_colwidth', None
    ):
        print(df.head(rows))

print_full_dataframe(segment_occupancy_df, rows=20)


