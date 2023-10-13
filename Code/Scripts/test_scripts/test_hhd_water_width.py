from multiprocessing import Pool

import geopandas as gpd


def find_branch_width(branch_geom, buffer_list, name, watervlak_geometry):
    for jx, buffer in enumerate(buffer_list):
        buffered_branches = branch_geom.buffer(buffer, cap_style=2)
        polygon = buffered_branches.intersection(watervlak_geometry)

        if jx > 0:
            overlap_area = polygon.area / buffered_branches.area
            if overlap_area < 0.9:
                width = buffer_list[jx - 1] * 2
                break
            else:
                width = buffer * 2

    return dict([("id", name), ("overlap", overlap_area), ("width", width)])


def main():
    p_folder = r"D:\Work\Project\P1414\GIS"
    branches_path = p_folder + r"\Uitgesneden watergangen\HHD_v7_test.shp"
    watervlak_path = (
        p_folder
        + r"\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Watervoerend deel.shp"
    )

    branches_gdf = gpd.read_file(branches_path, geometry="geometry")
    watervlak_gdf = gpd.read_file(watervlak_path, geometry="geometry")
    watervlak_geometry = watervlak_gdf.dissolve(by=None).geometry.values[0]

    buffer_list = [1.25, 2.5, 5, 10, 15, 20, 25, 50, 100]

    cpus = 8
    pool = Pool(processes=cpus)
    results = []
    for ix, (name, branch) in enumerate(branches_gdf.iterrows()):

        branch_geom = branch.geometry

        # result = find_branch_width(
        #     branch_geom=branch_geom,
        #     buffer_list=buffer_list,
        #     name=name,
        #     watervlak_geometry=watervlak_geometry,
        # )
        # results.append(result)
        kwds = dict(
            [
                ("branch_geom", branch_geom),
                ("buffer_list", buffer_list),
                ("name", name),
                ("watervlak_geometry", watervlak_geometry),
            ]
        )
        results.append(pool.apply_async(find_branch_width, kwds=kwds))

    for r in results:
        res = r.get()
        print(res)

    pool.close()
    pool.join()


if __name__ == "__main__":
    main()
