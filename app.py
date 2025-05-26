import streamlit as st
import os
import pandas as pd
from PIL import Image

# Set page config
st.set_page_config(
    page_title="Rail Simulation Viewer",
    layout="wide",
    page_icon="🚆"
)

st.markdown("""
<style>
    [data-testid="stAppViewContainer"] {
        background-color: #000000;
    }
    .stMarkdown, .stDataFrame, .stVideo, .stSelectbox label, .stSelectbox select, .stTextInput label, .stTextInput input {
        color: #ffffff !important;
    }
    h1, h2, h3, h4, h5, h6 {
        color: #ffffff !important;
    }
</style>
""", unsafe_allow_html=True)

# Title
st.title("🚆 Sligo-Dublin Rail Line Simulator 🚆")

# Paths
ANIMATIONS_DIR = "animations"
FIGURES_DIR = "appfigs"
COMMENTARY_DIR = "commentary"
TABLES_DIR = "apptables"
CHARTS_DIR = "appchart"

# Video selection
video_files = [f for f in os.listdir(ANIMATIONS_DIR) if f.endswith(".mp4")]
if not video_files:
    st.error("No animations found in the 'animations' folder.")
else:
    video_map = {f"Animation {i+1}": f for i, f in enumerate(video_files)}
    col1, _ = st.columns([1, 3])
    with col1:
        selected_label = st.selectbox(
            "🎬 Select a simulation to play:", 
            list(video_map.keys())
        )
    selected_video = video_map[selected_label]
    index = list(video_map.keys()).index(selected_label) + 1

    # 1. Video (left) + com1a (right)
    row1 = st.columns([3, 1], gap="large")
    with row1[0]:
        st.subheader("Playback")
        video_path = os.path.join(ANIMATIONS_DIR, selected_video)
        with open(video_path, 'rb') as f:
            st.video(f.read())
    with row1[1]:
        st.subheader("Description")
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}a.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())

    # 2. Stringline (left) + com1b (right)
    row2 = st.columns([3, 1], gap="large")
    with row2[0]:
        st.subheader("Stringline Diagram")
        fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
        if os.path.exists(fig_path):
            st.image(Image.open(fig_path), use_container_width=True)
    with row2[1]:
        st.subheader("")
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}b.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())

    # 3. Table (left) + com1c (right)
    row3 = st.columns([3, 1], gap="large")
    with row3[0]:
        st.subheader("Data Summary")
        table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
        if os.path.exists(table_path):
            df = pd.read_csv(table_path)
            st.dataframe(df, use_container_width=True)
    with row3[1]:
        st.subheader("")
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}c.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())

    # 4. Chart (left) + com1d (right)
    row4 = st.columns([3, 1], gap="large")
    with row4[0]:
        st.subheader("Analysis Chart")
        chart_path = os.path.join(CHARTS_DIR, f"chart_{index}.PNG")
        if os.path.exists(chart_path):
            st.image(Image.open(chart_path), use_container_width=True)
    with row4[1]:
        st.subheader("")
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}d.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())

    # Footer
    st.markdown("---")
    st.caption("🚉 Created by: Holly Briere-Edney,Shapefile © Heritage Studies Research Group, ATU Galway City (2023) CC BY-NC-SA 4.0, Built with Streamlit")
    
    
    