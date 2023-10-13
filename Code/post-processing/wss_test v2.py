import numpy as np
import rasterio
from rasterio.transform import from_bounds
from rasterio.windows import Window

raster_path = r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\fig_wl\480.tiff"
tile_path = r"D:\Work\Project\P1414\Models\Combined\V21_WBD_HD\run_model_changed_coords\dflowfm\output\fig_wl\tiles"

max_area = 500e6
max_x = max_y = np.sqrt(max_area)

csv = """ahn_version,3
scenario_type,3
scenario_calc_type,max
scenario_damage_table,dt.cfg
event_name,waterlevel,floodtime,repairtime_roads,repairtime_buildings,floodmonth,repetition_time
"""
csv_file_name = "index.csv"


with rasterio.open(raster_path) as src:
    bounds = src.bounds
    shape = src.shape
    transform = src.transform

    xorg = transform[2]
    yorg = transform[5]
    print(xorg, yorg)

    xres = transform[0]
    yres = transform[4]

    x_pixels = np.floor(max_x / xres)
    y_pixels = np.floor(max_y / yres)

    print(x_pixels, y_pixels)

    nx = np.ceil(shape[1] / x_pixels).astype(int)
    ny = np.ceil(shape[0] / y_pixels).astype(int)

    print(nx, ny)

    for x_slice in range(nx):
        for y_slice in range(ny):
            col_offset = x_slice * x_pixels
            row_offset = y_slice * y_pixels
            window = Window(
                col_off=col_offset, row_off=row_offset, width=x_pixels, height=y_pixels
            )
            wimage = src.read(1, window=window)

            transform = from_bounds(
                west=col_offset * xres + xorg,
                east=(col_offset + x_pixels) * xres + xorg,
                south=(row_offset + y_pixels) * yres + yorg,
                north=row_offset * yres + yorg,
                width=x_pixels,
                height=y_pixels,
            )
            name = r"{}_{}".format(x_slice, y_slice)
            file_name = r"{}_{}.tiff".format(x_slice, y_slice)
            entry = name + "," + file_name + ",24,10,10,9\n"
            csv = csv + entry

            filepath = tile_path + r"\\" + file_name
            with rasterio.open(
                filepath,
                "w",
                count=1,
                driver="GTiff",
                dtype=wimage.dtype,
                transform=transform,
                height=wimage.shape[0],
                width=wimage.shape[1],
            ) as dst:
                dst.write(wimage, 1)

with open(tile_path + r"\\" + csv_file_name, "w") as dst:
    for line in csv:
        dst.write(line)
