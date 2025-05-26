import os
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx

# --- 1. Setup directories ---
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
RAW_DATA_DIR = os.path.join(BASE_DIR, "rawdata")
SHAPEFILE_DIR = os.path.join(BASE_DIR, "shapefiles")

# --- 2. Define shapefile paths ---
LINE_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Line.shp")
STOPS_SHP = os.path.join(SHAPEFILE_DIR, "Sligo_Dublin_Stops.shp")

# --- 3. Load shapefiles ---
railline = gpd.read_file(LINE_SHP)
stations = gpd.read_file(STOPS_SHP)

# --- 4. Reproject to Web Mercator for basemap compatibility ---
railline = railline.to_crs(epsg=3857)
stations = stations.to_crs(epsg=3857)

# --- 5. Clean station names ---
stations["StopName"] = stations["StopName"].str.replace(" Train Station", "", regex=False).str.strip()

# --- 6. Plot the map ---
fig, ax = plt.subplots(figsize=(12, 12))
railline.plot(ax=ax, color='black', linewidth=4, label="Rail Line")
stations.plot(ax=ax, color='red', markersize=50, zorder=3, label="Stations")

# Add labels to stations
for _, row in stations.iterrows():
    ax.text(row.geometry.x + 2000, row.geometry.y + 4000, row["StopName"], fontsize=12, ha='left', color='darkslategray')

ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)

ax.set_title("Route 7: Sligo to Dublin (Connolly)", fontsize=16, pad=20)
ax.axis('off')
ax.legend(loc="lower left")

plt.tight_layout()
plt.show()
