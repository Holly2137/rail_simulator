
import streamlit as st
import os

# Set page title and layout
st.set_page_config(page_title="Rail Simulation Viewer", layout="centered")
st.title("🚆 Rail Simulation Player")

# Path to animations folder
ANIMATIONS_DIR = "animations"

# List available .mp4 videos
video_files = [f for f in os.listdir(ANIMATIONS_DIR) if f.endswith(".mp4")]

if not video_files:
    st.error("No animations found in the 'animations' folder.")
else:
    # Map filenames to display names
    video_map = {f"Animation {i+1}": f for i, f in enumerate(video_files)}
    selected_label = st.selectbox("Select a simulation to play:", list(video_map.keys()))
    selected_video = video_map[selected_label]

    st.markdown(f"### Now Playing: {selected_label}")
    video_path = os.path.join(ANIMATIONS_DIR, selected_video)
    with open(video_path, 'rb') as f:
        st.video(f.read())

    # Optional description box
    st.markdown("---")
    st.info("This is a pre-rendered simulation showing train movements across the Sligo-Dublin line.")



























