# import streamlit as st
# import os
# import pandas as pd
# from PIL import Image

# # Set page title and layout
# st.set_page_config(page_title="Rail Simulation Viewer", layout="centered")
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

#     st.markdown(f"### Now Playing: {selected_label}")
#     video_path = os.path.join(ANIMATIONS_DIR, selected_video)
#     with open(video_path, 'rb') as f:
#         st.video(f.read())

#     # Columns for additional info
#     col1, col2 = st.columns([1, 2])

#     # Show graph if exists
#     fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
#     if os.path.exists(fig_path):
#         col1.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_column_width=True)
#     else:
#         col1.warning("No graph available.")

#     # Show commentary if exists
#     com_path = os.path.join(COMMENTARY_DIR, f"com_{index}.txt")
#     if os.path.exists(com_path):
#         with open(com_path, 'r') as f:
#             commentary = f.read()
#         col2.markdown("**Commentary:**")
#         col2.info(commentary)
#     else:
#         col2.warning("No commentary available.")

#     # Show table if exists
#     table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
#     if os.path.exists(table_path):
#         df = pd.read_csv(table_path)
#         st.markdown("### Summary Table")
#         st.dataframe(df)
#     else:
#         st.warning("No table data available.")



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

#     # Top row: video and commentary
#     st.markdown("<div style='display: flex; gap: 2rem;'>", unsafe_allow_html=True)

#     # Video block
#     st.markdown("""
#     <div style='flex: 2; border: 2px solid green; padding: 10px;'>
#         <h4>Animation</h4>
#     """, unsafe_allow_html=True)
#     video_path = os.path.join(ANIMATIONS_DIR, selected_video)
#     with open(video_path, 'rb') as f:
#         st.video(f.read())
#     st.markdown("</div>", unsafe_allow_html=True)

#     # Commentary block
#     st.markdown("""
#     <div style='flex: 1; border: 2px solid green; padding: 10px;'>
#         <h4>Commentary</h4>
#     """, unsafe_allow_html=True)
#     com_path = os.path.join(COMMENTARY_DIR, f"com_{index}.txt")
#     if os.path.exists(com_path):
#         with open(com_path, 'r') as f:
#             commentary = f.read()
#         st.markdown(commentary)
#     else:
#         st.warning("No commentary available.")
#     st.markdown("</div></div>", unsafe_allow_html=True)

#     # Bottom row: graph and table
#     bottom1, bottom2 = st.columns([3, 2])
#     with bottom1:
#         fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
#         if os.path.exists(fig_path):
#             st.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_container_width=True)
#             st.markdown("<div style='border: 2px solid green; padding: 10px;'>Graph Displayed Above</div>", unsafe_allow_html=True)
#         else:
#             st.warning("No graph available.")

#     with bottom2:
#         table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
#         if os.path.exists(table_path):
#             df = pd.read_csv(table_path)
#             st.markdown("### Summary Table")
#             st.dataframe(df, use_container_width=True)
#             st.markdown("<div style='border: 2px solid green; padding: 10px;'>Table Displayed Above</div>", unsafe_allow_html=True)
#         else:
#             st.warning("No table data available.")

# RULE: No code should be changed without explicit user approval (via 'x').

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
#     .box {
#         border: 3px solid #146d6d;
#         border-radius: 10px;
#         padding: 15px;
#         box-shadow: 2px 2px 5px rgba(0,0,0,0.3);
#         margin-bottom: 20px;
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
#             st.markdown("<div class='box'>", unsafe_allow_html=True)
#             st.markdown("#### Animation")
#             video_path = os.path.join(ANIMATIONS_DIR, selected_video)
#             with open(video_path, 'rb') as f:
#                 st.video(f.read())
#             st.markdown("</div>", unsafe_allow_html=True)

#     with top2:
#         with st.container():
#             st.markdown("<div class='box'>", unsafe_allow_html=True)
#             st.markdown("#### Commentary")
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
#             st.markdown("<div class='box'>", unsafe_allow_html=True)
#             fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
#             if os.path.exists(fig_path):
#                 st.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_container_width=True)
#             else:
#                 st.warning("No graph available.")
#             st.markdown("</div>", unsafe_allow_html=True)

#     with bottom2:
#         with st.container():
#             st.markdown("<div class='box'>", unsafe_allow_html=True)
#             table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
#             if os.path.exists(table_path):
#                 df = pd.read_csv(table_path)
#                 st.markdown("### Summary Table")
#                 st.dataframe(df, use_container_width=True)
#             else:
#                 st.warning("No table data available.")
#             st.markdown("</div>", unsafe_allow_html=True)
            

# RULE: No code should be changed without explicit user approval (via 'x').

import streamlit as st
import os
import pandas as pd
from PIL import Image

# Set page title and layout
st.set_page_config(page_title="Rail Simulation Viewer", layout="wide")
st.title("🚆 Rail Simulation Player")

# Custom CSS styling for better aesthetics
st.markdown("""
    <style>
    .block-container {
        background-color: #0b1e24;
        color: white;
    }
    .box-wrap {
        display: flex;
        flex-direction: column;
        border: 3px solid #6fb7af;
        border-radius: 15px;
        padding: 20px;
        box-shadow: 3px 3px 6px rgba(0,0,0,0.4);
        margin-bottom: 25px;
    }
    h4 {
        margin-top: 0;
        margin-bottom: 10px;
    }
    </style>
""", unsafe_allow_html=True)

# Paths to folders
ANIMATIONS_DIR = "animations"
FIGURES_DIR = "appfigs"
COMMENTARY_DIR = "commentary"
TABLES_DIR = "apptables"

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

    # Top row: video and commentary (side by side)
    top1, top2 = st.columns([3, 3])
    with top1:
        st.markdown(f"""
            <div class='box-wrap'>
                <h4>Playback</h4>
        """, unsafe_allow_html=True)
        video_path = os.path.join(ANIMATIONS_DIR, selected_video)
        with open(video_path, 'rb') as f:
            st.video(f.read())
        st.markdown("</div>", unsafe_allow_html=True)

    with top2:
        st.markdown(f"""
            <div class='box-wrap'>
                <h4>Insight</h4>
        """, unsafe_allow_html=True)
        com_path = os.path.join(COMMENTARY_DIR, f"com_{index}.txt")
        if os.path.exists(com_path):
            with open(com_path, 'r') as f:
                commentary = f.read()
            st.markdown(commentary)
        else:
            st.warning("No commentary available.")
        st.markdown("</div>", unsafe_allow_html=True)

    # Bottom row: graph and table
    bottom1, bottom2 = st.columns([3, 3])
    with bottom1:
        st.markdown(f"""
            <div class='box-wrap'>
                <h4>Route 7  Stringline Diagram</h4>
        """, unsafe_allow_html=True)
        fig_path = os.path.join(FIGURES_DIR, f"appfig_{index}.png")
        if os.path.exists(fig_path):
            st.image(Image.open(fig_path), caption=f"Graph for {selected_label}", use_container_width=True)
        else:
            st.warning("No graph available.")
        st.markdown("</div>", unsafe_allow_html=True)

    with bottom2:
        st.markdown(f"""
            <div class='box-wrap'>
                <h4>Data Summary</h4>
        """, unsafe_allow_html=True)
        table_path = os.path.join(TABLES_DIR, f"tab_{index}.csv")
        if os.path.exists(table_path):
            df = pd.read_csv(table_path)
            st.dataframe(df, use_container_width=True)
        else:
            st.warning("No table data available.")
        st.markdown("</div>", unsafe_allow_html=True)


