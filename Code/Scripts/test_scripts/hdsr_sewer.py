import geopandas as gpd
import numpy as np
from rasterstats import zonal_stats

ahn_path = r"D:\Work\Project\P1414\GIS\AHN\AHN4_WSS_filled.tif"
knopen_path = r"D:\Work\Project\HL-23006\GIS\Riool\knopen_waterloop.shp"
bg_path = r"D:\Work\Project\HL-23006\GIS\Riool\Bemalingsgebieden.shp"
out_path = r"D:\Work\Project\HL-23006\GIS\Riool\comb.shp"

affected_columns = [
    "verhard_m2",
    "berging_m3",
    "aantal_inw",
    "dwa_inw",
    "dwa_bedrij",
    "geinst_cap",
    "poc",
]

knopen_gdf = gpd.read_file(knopen_path)
bg_gdf = gpd.read_file(bg_path)
stats = zonal_stats(
    bg_gdf.geometry,
    ahn_path,
    stats=["mean"],
    # add_stats={"nanmean": np.nanmean},
    nodata=np.nan,
)
# print(stats)
bg_gdf["ahn_mean"] = [stat["mean"] for stat in stats]

bg_gdf = bg_gdf.loc[
    :,
    [
        "id",
        "naam",
        "stelsel",
        "gemeente",
        "verhard_m2",
        "berging_m3",
        "berging_mm",
        "aantal_inw",
        "dwa_inw",
        "dwa_bedrij",
        "geinst_cap",
        "poc",
        "idrwzi",
        "ahn_mean",
        "geometry",
    ],
]
knopen_mask = (
    (knopen_gdf["functie"] == "Externe overstort")
    | (knopen_gdf["functie"] == "Regenwateroverstort")
    | (knopen_gdf["functie"] == "Regenwateruitlaat")
    | (knopen_gdf["functie"] == "uitlaat/lozingswerk")
)
knopen_gdf = knopen_gdf.loc[knopen_mask, :]
knopen_gdf = knopen_gdf.loc[
    :,
    [
        "objectid_1",
        "knoop_id",
        "type",
        "subtype",
        "functie",
        "afvoertype",
        "drempelbre",
        "drempelniv",
        "waterloop",
        "geometry",
    ],
]


comb_gdf = gpd.sjoin(left_df=knopen_gdf, right_df=bg_gdf, how="left", predicate="within")
g_gdf = comb_gdf.groupby("id")
for name, group in g_gdf:
    n = group.shape[0]
    if n > 1:
        mask = comb_gdf["id"] == name
        comb_gdf.loc[mask, affected_columns] = (
            comb_gdf.loc[mask, affected_columns].astype(float).divide(n)
        )

comb_gdf.to_file(out_path)
