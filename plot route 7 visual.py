
import os
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import matplotlib.patches as patches

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
    name = row["StopName"]

    # Rename Sligo
    if name == "Sligo":
        name = "Sligo (MacDiarmada)"

    # Offset Connolly differently to avoid overlap
    if name == "Connolly":
        x_offset = 2000
        y_offset = 1000
    else:
        x_offset = 1000
        y_offset = 3000

    ax.text(row.geometry.x + x_offset, row.geometry.y + y_offset, name,
            fontsize=12, weight='bold', ha='left', color='darkgoldenrod')

# Expand plot bounds to fit all labels
xlims = ax.get_xlim()
ylims = ax.get_ylim()
ax.set_xlim(xlims[0] - 2000, xlims[1] + 18000)
ax.set_ylim(ylims[0] - 10000, ylims[1] + 15000)

# Add basemap and border
#ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron)
ctx.add_basemap(ax, source=ctx.providers.CartoDB.Voyager)
#ctx.add_basemap(ax, source=ctx.providers.Esri.WorldTopoMap)
#ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.HOT)
#ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery)

# Add thick black border around the basemap extent
border_rect = patches.Rectangle(
    (ax.get_xlim()[0], ax.get_ylim()[0]),
    ax.get_xlim()[1] - ax.get_xlim()[0],
    ax.get_ylim()[1] - ax.get_ylim()[0],
    linewidth=5, edgecolor='black', facecolor='none', zorder=10
)
ax.add_patch(border_rect)

# Add credit box in bottom left (moved up slightly)
credit_text = ("Image by: Holly Briere-Edney\n"
                "Shapefile © Heritage Studies Research Group, ATU Galway City (2023), CC BY-NC-SA 4.0")
ax.text(0.01, 0.05, credit_text, transform=ax.transAxes, fontsize=8, ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4'))

# Add title inside plot area, center top
ax.text(0.5, 0.98, "Route 7: Sligo (MacDiarmada) to Dublin (Connolly)", transform=ax.transAxes,
        fontsize=16, weight='bold', ha='center', va='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

# Move legend higher and increase font size by 1.5x
legend = ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.09))
for text in legend.get_texts():
    text.set_fontsize(13.5)  # default is 9, so 9 * 1.5 = 13.5

ax.axis('off')
plt.tight_layout()
plt.show()

