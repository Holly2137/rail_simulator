

# import streamlit as st
# import os
# import pandas as pd
# from PIL import Image

# # Set page title and layout
# st.set_page_config(page_title="Rail Simulation Viewer", layout="wide")
# st.title("🚆 Rail Simulation Player")

# # Custom CSS styling for better aesthetics
# st.markdown("""
#     <style>
#     .block-container {
#         background-color: #0b1e24;
#         color: white;
#     }
#     .box-wrap {
#         border: 3px solid #6fb7af;
#         border-radius: 15px;
#         padding: 20px;
#         box-shadow: 3px 3px 6px rgba(0,0,0,0.4);
#         margin-bottom: 25px;
#     }
#     h4 {
#         margin-top: 0;
#         margin-bottom: 10px;
#     }
#     </style>
# """, unsafe_allow_html=True)

# # Paths to folders
# ANIMATIONS_DIR = "animations"
# FIGURES_DIR = "appfigs"
# COMMENTARY_DIR = "commentary"
# TABLES_DIR = "apptables"

# # List available .mp4 videos
# video_files = [f for f in os.listdir(ANIMATIONS_DIR) if f.endswith(".mp4")]

# if not video_files:
#     st.error("No animations found in the 'animations' folder.")
# else:
#     # Map filenames to display names
#     video_map = {f"Animation {i+1}": f for i, f in enumerate(video_files)}
#     selected_label = st.selectbox("Select a simulation to play:", list(video_map.keys()))
#     selected_video = video_map[selected_label]
#     index = list(video_map.keys()).index(selected_label) + 1

#     # Top row: video and commentary (side by side)
#     top1, top2 = st.columns([3, 3])
#     with top1:
#         with st.container():
#             st.markdown("<div class='box-wrap'>", unsafe_allow_html=True)
#             st.subheader("Playback")
#             video_path = os.path.join(ANIMATIONS_DIR, selected_video)
#             with open(video_path, 'rb') as f:
#                 st.video(f.read())
#             st.markdown("</div>", unsafe_allow_html=True)

#     with top2:
#         with st.container():
#             st.markdown("<div class='box-wrap'>", unsafe_allow_html=True)
#             st.subheader("Insight")
#             com_path = os.path.join(COMMENTARY_DIR, f"com_{index}.txt")
#             if os.path.exists(com_path):
#                 with open(com_path, 'r') as f:
#                     commentary = f.read()
#                 st.markdown(commentary)
#             else:
#                 st.warning("No commentary available.")
#             st.markdown("</div>", unsafe_allow_html=True)

#     # Bottom row: graph and table
#     bottom1, bottom2 = st.columns([3, 3])
#     with bottom1:
#         with st.container():
#             st.markdown("<div class='box-wrap'>", unsafe_allow_html=True)
#             st.subheader("Route 7  Stringline Diagram")
#             fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
#             if os.path.exists(fig_path):
#                 st.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_container_width=True)
#             else:
#                 st.warning("No graph available.")
#             st.markdown("</div>", unsafe_allow_html=True)

#     with bottom2:
#         with st.container():
#             st.markdown("<div class='box-wrap'>", unsafe_allow_html=True)
#             st.subheader("Data Summary")
#             table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
#             if os.path.exists(table_path):
#                 df = pd.read_csv(table_path)
#                 st.dataframe(df, use_container_width=True)
#             else:
#                 st.warning("No table data available.")
#             st.markdown("</div>", unsafe_allow_html=True)


# import streamlit as st
# import os
# import pandas as pd
# from PIL import Image

# # Set page title and layout
# st.set_page_config(page_title="Rail Simulation Viewer", layout="wide")
# st.title("🚆 Rail Simulation Player")

# # Paths to folders
# ANIMATIONS_DIR = "animations"
# FIGURES_DIR = "appfigs"
# COMMENTARY_DIR = "commentary"
# TABLES_DIR = "apptables"

# # List available .mp4 videos
# video_files = [f for f in os.listdir(ANIMATIONS_DIR) if f.endswith(".mp4")]

# if not video_files:
#     st.error("No animations found in the 'animations' folder.")
# else:
#     # Map filenames to display names
#     video_map = {f"Animation {i+1}": f for i, f in enumerate(video_files)}
#     selected_label = st.selectbox("Select a simulation to play:", list(video_map.keys()))
#     selected_video = video_map[selected_label]
#     index = list(video_map.keys()).index(selected_label) + 1

#     # Top row: video and commentary (side by side)
#     top1, top2 = st.columns([2, 1],gap="large")
#     with top1:
#         st.subheader("Playback")
#         video_path = os.path.join(ANIMATIONS_DIR, selected_video)
#         with open(video_path, 'rb') as f:
#             st.video(f.read())

#     with top2:
#         st.subheader("Insight")
#         com_path = os.path.join(COMMENTARY_DIR, f"com_{index}.txt")
#         if os.path.exists(com_path):
#             with open(com_path, 'r') as f:
#                 commentary = f.read()
#             st.markdown(commentary)
#         else:
#             st.warning("No commentary available.")

#     # Bottom row: graph and table
#     bottom1, bottom2 = st.columns([2, 1],gap="large")
#     with bottom1:
#         st.subheader("Route 7 Stringline Diagram")
#         fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
#         if os.path.exists(fig_path):
#             st.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_container_width=True)
#         else:
#             st.warning("No graph available.")

#     with bottom2:
#         st.subheader("Data Summary")
#         table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
#         if os.path.exists(table_path):
#             df = pd.read_csv(table_path)
#             st.dataframe(df, use_container_width=True)
#         else:
#             st.warning("No table data available.")


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

# Title
st.title("🚆 Rail Simulation Player")

# Paths
ANIMATIONS_DIR = "animations"
FIGURES_DIR = "appfigs"
COMMENTARY_DIR = "commentary"
TABLES_DIR = "apptables"

# Video selection
video_files = [f for f in os.listdir(ANIMATIONS_DIR) if f.endswith(".mp4")]
if not video_files:
    st.error("No animations found in the 'animations' folder.")
else:
    video_map = {f"Animation {i+1}": f for i, f in enumerate(video_files)}
    selected_label = st.selectbox(
        "🎬 Select a simulation to play:", 
        list(video_map.keys())
    )
    selected_video = video_map[selected_label]
    index = list(video_map.keys()).index(selected_label) + 1

    # Section 1: Video and first commentary
    top1, top2 = st.columns([2, 1], gap="large")
    with top1:
        st.subheader("Playback")
        video_path = os.path.join(ANIMATIONS_DIR, selected_video)
        with open(video_path, 'rb') as f:
            st.video(f.read())

    with top2:
        st.subheader("Discription")
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}a.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())
        else:
            st.warning("No commentary available.")

    # Section 2: Table and second commentary
    mid1, mid2 = st.columns([2, 1], gap="large")
    with mid1:
        st.subheader("Data Summary")
        table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
        if os.path.exists(table_path):
            df = pd.read_csv(table_path)
            st.dataframe(df, use_container_width=True)
        else:
            st.warning("No table data available.")

    with mid2:
        
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}b.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())
        else:
            st.warning("No commentary available.")

    # Section 3: Stringline and third commentary
    bot1, bot2 = st.columns([2, 1], gap="large")
    with bot1:
        st.subheader("Route 7 Stringline Diagram")
        fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
        if os.path.exists(fig_path):
            st.image(Image.open(fig_path), use_container_width=True)
        else:
            st.warning("No graph available.")

    with bot2:
        
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}c.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                st.markdown(f.read())
        else:
            st.warning("No commentary available.")

    # Footer
    st.markdown("---")
    st.caption("🚉 Built with Streamlit")