from copy import copy

import geopandas as gpd
import numpy as np
import pandas as pd


def _clip_structures_by_branches(dm, features: list, buffer: float = 1, min_overlap: float = 0.95):

    buffered_branches = (
        gpd.GeoDataFrame(
            dm.waterloop.buffer(buffer, cap_style=2),
            columns=["geometry"],
            crs=dm.waterloop.crs,
            geometry="geometry",
        )
        .dissolve(by=None)
        .explode(ignore_index=True)
    )

    for feature in features:
        if (
            (feature == "waterloop")
            or (feature == "profielpunt")
            or (feature == "keringen")
            or (feature == "peilgebieden")
            # or (feature == "rivier_profielen_data")
            # or (feature == "rivier_profielen_ruwheid")
        ):
            continue

        ddm_feature = getattr(dm, feature)
        ddm_feature = gpd.GeoDataFrame(
            ddm_feature[ddm_feature["geometry"].notna()],
            crs=ddm_feature.crs,
            geometry="geometry",
        )

        if feature == "profiellijn":
            _buffered_branches = gpd.GeoDataFrame(
                dm.waterloop.dissolve(by=None).buffer(20).explode(ignore_index=True),
                columns=["geometry"],
                crs=ddm_feature.crs,
                geometry="geometry",
            )
            # clipped_gdf = gpd.overlay(
            #     ddm_feature, _buffered_branches, how="intersection", keep_geom_type=True
            # )
            print(feature)
            print(ddm_feature.shape[0])
            columns = ddm_feature.columns
            clipped_gdf = gpd.sjoin(
                left_df=ddm_feature,
                right_df=buffered_branches,
                how="left",
                predicate="intersects",
            )
            clipped_gdf = clipped_gdf.loc[clipped_gdf["index_right"].notnull(), :]
            clipped_gdf = clipped_gdf[columns]
            print(clipped_gdf.shape[0])
            setattr(dm, feature, clipped_gdf)

            if dm.profielpunt is not None:
                profiel_punt_gdf = dm.profielpunt

                profiel_punt_gdf_out = None

                for ix, (name, line) in enumerate(clipped_gdf.iterrows()):
                    points = profiel_punt_gdf.loc[
                        profiel_punt_gdf["profiellijnid"] == line["globalid"], :
                    ]
                    if profiel_punt_gdf_out is None:
                        profiel_punt_gdf_out = pd.DataFrame(data=points)
                    else:
                        profiel_punt_gdf_out = pd.concat([profiel_punt_gdf_out, points])

                dm.profielpunt = gpd.GeoDataFrame(
                    profiel_punt_gdf_out, geometry="geometry", crs=profiel_punt_gdf.crs
                )

        elif ddm_feature.shape[0] > 0:
            print(feature)
            print(ddm_feature.shape[0])
            # clipped_gdf = gpd.overlay(
            #     ddm_feature, buffered_branches, how="intersection", keep_geom_type=True
            # )
            columns = ddm_feature.columns
            clipped_gdf = gpd.sjoin(
                left_df=ddm_feature, right_df=buffered_branches, how="left", predicate="within"
            )
            clipped_gdf = clipped_gdf.loc[clipped_gdf["index_right"].notnull(), :]
            clipped_gdf = clipped_gdf[columns]
            print(clipped_gdf.shape[0])

            # mls_struct_bool = clipped_gdf.geometry.type == "MultiLineString"
            # if np.sum(mls_struct_bool) > 0:
            #     _clipped_gdf = copy(clipped_gdf)

            #     for ix, struct in clipped_gdf.iterrows():
            #         new_length = np.sum(clipped_gdf.loc[struct.name].geometry.length)
            #         old_length = np.sum(ddm_feature.loc[struct.name].geometry.length)

            #         if mls_struct_bool[ix]:
            #             if (new_length / old_length) >= min_overlap:
            #                 _clipped_gdf.loc[struct.name] = ddm_feature.loc[struct.name]
            #             else:
            #                 _clipped_gdf.drop(index=struct.name, inplace=True)

            #         elif not mls_struct_bool[ix]:
            #             if (new_length / old_length) < min_overlap:
            #                 _clipped_gdf.drop(index=struct.name, inplace=True)
            #     clipped_gdf = _clipped_gdf

            setattr(dm, feature, clipped_gdf)
    return dm
