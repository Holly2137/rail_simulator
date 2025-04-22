import streamlit as st
import os
import pandas as pd
from PIL import Image

# Set page title and layout
st.set_page_config(page_title="Rail Simulation Viewer", layout="centered")
st.title("🚆 Rail Simulation Player")

# Paths to folders
ANIMATIONS_DIR = "animations"
FIGURES_DIR = "appfigs"
COMMENTARY_DIR = "commentary"
TABLES_DIR = "apptabs"

# List available .mp4 videos
video_files = [f for f in os.listdir(ANIMATIONS_DIR) if f.endswith(".mp4")]

if not video_files:
    st.error("No animations found in the 'animations' folder.")
else:
    # Map filenames to display names
    video_map = {f"Animation {i+1}": f for i, f in enumerate(video_files)}
    selected_label = st.selectbox("Select a simulation to play:", list(video_map.keys()))
    selected_video = video_map[selected_label]
    index = list(video_map.keys()).index(selected_label) + 1

    st.markdown(f"### Now Playing: {selected_label}")
    video_path = os.path.join(ANIMATIONS_DIR, selected_video)
    with open(video_path, 'rb') as f:
        st.video(f.read())

    # Columns for additional info
    col1, col2 = st.columns([1, 2])

    # Show graph if exists
    fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
    if os.path.exists(fig_path):
        col1.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_column_width=True)
    else:
        col1.warning("No graph available.")

    # Show commentary if exists
    com_path = os.path.join(COMMENTARY_DIR, f"com_{index}.txt")
    if os.path.exists(com_path):
        with open(com_path, 'r') as f:
            commentary = f.read()
        col2.markdown("**Commentary:**")
        col2.info(commentary)
    else:
        col2.warning("No commentary available.")

    # Show table if exists
    table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
    if os.path.exists(table_path):
        df = pd.read_csv(table_path)
        st.markdown("### Summary Table")
        st.dataframe(df)
    else:
        st.warning("No table data available.")



