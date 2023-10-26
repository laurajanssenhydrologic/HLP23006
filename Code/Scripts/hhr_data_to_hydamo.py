import geopandas as gpd
import pandas as pd


def combine_wl_and_profile(wl_gdf: gpd.GeoDataFrame, wl_data_df: pd.DataFrame) -> gpd.GeoDataFrame:
    wl_data_df["CODE_WTRG"] = wl_data_df["CODE_PROF"].apply(
        lambda code: code.removeprefix("PROF_")
    )
    wl_gdf = gpd.GeoDataFrame(pd.merge(left=wl_gdf, right=wl_data_df, on="CODE_WTRG"))

    wl_gdf = wl_gdf.rename(
        columns=dict([("BREEDTE", "bodembreedte"), ("BODEMHOOGTE", "bodemhoogte benedenstrooms")])
    )
    wl_gdf = wl_gdf.drop(columns="CODE_PROF")
    wl_gdf["bodemhoogte bovenstrooms"] = wl_gdf["bodemhoogte benedenstrooms"]
    wl_gdf["hoogte insteek linkerzijde"] = wl_gdf["bodemhoogte benedenstrooms"] + 1
    wl_gdf["hoogte insteek rechterzijde"] = wl_gdf["bodemhoogte benedenstrooms"] + 1
    return wl_gdf


if __name__ == "__main__":
    kunstwerk_shp_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\kunstwerk.shp"
    pomp_data_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\pompen.csv"
    profiel_shp_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\profiel.shp"
    profiel_data_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\profielen.csv"
    retentie_shp_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\retentie.shp"
    retentie_data_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\retentie.csv"
    stuw_data_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\stuwen.csv"
    watergang_shp_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\watergang.shp"

    out_path = r"D:\Work\Project\HL-23029\GIS\polderwatergangen\polderwatergang.gpkg"

    kw_gdf = gpd.read_file(kunstwerk_shp_path)

    pomp_data_df = pd.read_csv(pomp_data_path)
    pomp_codes = pomp_data_df["CODE_KNSTW"].values[:]
    pomp_gdf = kw_gdf.loc[kw_gdf["CODE_KNSTW"].isin(pomp_codes), :]
    pomp_gdf = gpd.GeoDataFrame(pd.merge(left=pomp_gdf, right=pomp_data_df, on="CODE_KNSTW"))
    print(pomp_gdf.head())
    pomp_gdf.to_file(out_path, layer="gemaal")

    stuw_data_df = pd.read_csv(stuw_data_path)
    stuw_codes = stuw_data_df["CODE_KNSTW"].values[:]
    stuw_gdf = kw_gdf.loc[kw_gdf["CODE_KNSTW"].isin(stuw_codes), :]
    stuw_gdf = gpd.GeoDataFrame(pd.merge(left=stuw_gdf, right=stuw_data_df, on="CODE_KNSTW"))
    print(stuw_gdf.head())
    stuw_gdf.to_file(out_path, layer="stuw")

    ret_gdf = gpd.read_file(retentie_shp_path)
    ret_data_df = pd.read_csv(retentie_data_path)
    ret_gdf = gpd.GeoDataFrame(pd.merge(left=ret_gdf, right=ret_data_df, on="CODE_RET"))
    print(ret_gdf.head())
    ret_gdf.to_file(out_path, layer="retentie")

    wl_gdf = gpd.read_file(watergang_shp_path)
    wl_data_df = pd.read_csv(profiel_data_path)
    wl_gdf = combine_wl_and_profile(wl_gdf=wl_gdf, wl_data_df=wl_data_df)
    print(wl_gdf.head())
    wl_gdf.to_file(out_path, layer="waterloop")
