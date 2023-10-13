import re
from scipy.spatial import KDTree
import numpy as np
import os
import shutil
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import json 
from shapely.geometry import Point, LineString
import geopandas as gpd
from ipyleaflet import Map, Marker, projections, basemaps
import ipyleaflet as ipl
import pyproj
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points

class ModifyDambreakSettingsDHYDRO(WidgetStyling):
    def __init__(self, path):
        self.path = os.path.join(path, 'dflowfm\structures.ini')
        self.structure_textfile, self.dambreak_settings = self.load_dambreak_settings()
        self.strucutre_template_no_dambreak = self.structure_textfile
        self.wrote_backup = False
        
        self.set_default_layout_and_styling()

        self.settings_to_modify = ['crestLevelIni', 't0', 'timeToBreachToMaximumDepth',
                                  'crestLevelMin', 'breachWidthIni', 'f1', 'f2', 'uCrit']

        self.widgets = {}
        for setting in self.settings_to_modify:
            self.widgets[setting] = ipy.FloatText(
                value = float(self.dambreak_settings[setting]),
                description=f'{setting}:',
                disabled=False
                )
            self.widgets[setting].layout = self.item_layout
            self.widgets[setting].style = self.item_style
        
        self.widgets_to_display = [widget for widget in self.widgets.values()]
        
        
    def load_dambreak_settings(self):
        # find occurences of dambreak structure type
        pattern = r"type\s*=\s*dambreak"
        matches = {}
        with open(self.path, "r") as f:
            structures_textfile = f.readlines()  
        for i, line in enumerate(structures_textfile):
            if re.search(pattern, line):
                matches[i] = line

        # if there is only one occurence, start going up to look for [Structure]
        if len(matches) != 1:
            raise Exception(f"Error, {len(matches)} dambreaks found in structures.ini")
        found = False
        start_index = list(matches.keys())[0]
        pattern = r"\[Structure\]"
        look_up_dist = 0
        while not found and look_up_dist < 10:
            look_up_dist = look_up_dist + 1
            i = start_index - look_up_dist
            if re.search(pattern, structures_textfile[i]):
                found = True
                dambreak_settings = structures_textfile[i:]
                structures_textfile = structures_textfile[:i]
        if not found:
            raise Exception(f"Error, did not find dambreak structure in structures.ini")
            
        # find end of dambreak strucutre
        for i in range(len(dambreak_settings)):
            if dambreak_settings[i] == '\n':
                dambreak_settings = dambreak_settings[:i]
                structures_textfile = structures_textfile + dambreak_settings[i:]
                break

        # read paramters from structure file
        pattern = r"(\w+)\s*=\s*([^\n#]+)(?:\s*#\s*(.*))?"
        dambreak_settings_dict = {}
        # Find all key-value pairs in the text data
        for line in dambreak_settings:
            match = re.search(pattern, line)
            if match:
                key = match[1].strip()
                value = match[2].strip()
                dambreak_settings_dict[key] = value
        return structures_textfile, dambreak_settings_dict
    
    def modify_dikebreak_location(self, db_gdf, fw_gdf, max_dist = 100, z_default = 0):
        fw_point_list = []
        for ix, branch in fw_gdf.iterrows():
            coords = branch.geometry.coords[:]

            for x, y, z in coords:
                fw_point_list.append([x, y, z])

        fw_kdtree = KDTree(data=fw_point_list)

        for name, db in db_gdf.iterrows():
            geom = db.geometry
            if not isinstance(geom, LineString):
                raise ValueError("I need a LineString...")

            coords = np.array(geom.coords[:])
            centroid = geom.interpolate(0.5, normalized=True)
            if len(coords[0]) == 2:
                mid_point = [centroid.x, centroid.y, z_default]
            elif len(coords[0]) == 3:
                mid_point = [centroid.x, centroid.y, centroid.z]
            else:
                raise ValueError("can't cope with {} coordinates".format(len(coords[0])))

            _, ix_point_in_kering = fw_kdtree.query(mid_point, distance_upper_bound=max_dist)
            point_in_kering = fw_point_list[ix_point_in_kering]
            self.dambreak_settings['crestLevelIni'] = round(point_in_kering[2] - 4, 2) # - db.crestlevelini
            self.dambreak_settings['t0'] = 7200
            
            self.dambreak_settings['numCoordinates'] = coords.shape[0]
            self.dambreak_settings['startLocationX'] = point_in_kering[0]
            self.dambreak_settings['startLocationY'] = point_in_kering[1]
            self.dambreak_settings['xCoordinates'] = ' '.join(str(x) for x in coords[:, 0].tolist())
            self.dambreak_settings['yCoordinates'] = ' '.join(str(x) for x in coords[:, 1].tolist())

#             These are removed from the sas service, since it is not required. This way D-HYDRO will determine this.
#             self.dambreak_settings['waterLevelDownstreamLocationX'] = coords[-1, 0]
#             self.dambreak_settings['waterLevelDownstreamLocationY'] = coords[-1, 1]
#             self.dambreak_settings['waterLevelUpstreamLocationX'] = coords[0, 0]
#             self.dambreak_settings['waterLevelUpstreamLocationY'] = coords[0, 1]
            keys_to_remove = ['waterLevelDownstreamLocationX', 'waterLevelDownstreamLocationY',
                             'waterLevelUpstreamLocationX', 'waterLevelUpstreamLocationY']
            for key in keys_to_remove:
                if key in  self.dambreak_settings.keys():
                    self.dambreak_settings.pop(key)
     
    def write_to_structures(self, write_output = True, backup_original = True):
        structures = self.strucutre_template_no_dambreak[:]
        structures.append('\n')
        structures.append('[Structure]\n')
        
        for key, value in self.dambreak_settings.items():
            structures.append(f"    {key.ljust(35)}= {value}\n")
        
        if write_output:
            if backup_original and self.wrote_backup == False:
                backup_file_loc = os.path.join(os.path.dirname(self.path), 'stuctures_backup.ini')
                shutil.copyfile(self.path, backup_file_loc)
                self.wrote_backup == True
                
            with open(self.path, 'w') as file:
                for line in structures:
                    file.write(line)
        return structures

    def update_widget_values(self):
        for setting in self.settings_to_modify:
            self.widgets[setting].value = self.dambreak_settings[setting]
        
    def display_widgets(self):
        self.update_widget_values()
        items = ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style)
        display(items)
        
        button = ipy.Button(description="Update settings", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.update_settings_widget)

    def update_settings_widget(self, b):
        for setting in self.settings_to_modify:
            self.dambreak_settings[setting] = self.widgets[setting].value
        
        lines = self.write_to_structures(write_output = True)
        _, dambreak_settings_from_file = self.load_dambreak_settings() # read from structures.ini, so you are sure you have correct values displayed
        
        clear_output(wait=True)
        self.display_widgets()
        display("Dike breach settings are:")
        print_settings_dict = {k: dambreak_settings_from_file[k] for k in self.settings_to_modify}
        print(json.dumps(print_settings_dict, indent=4))
    


class DambreakWidget(WidgetStyling):
    def __init__(self):
        # Define the two coordinate systems
        self.crs_rd =28992
        self.crs_map = 4326
        crs_rd = pyproj.CRS.from_epsg(self.crs_rd)
        crs_map = pyproj.CRS.from_epsg(self.crs_map)

        # Create a transformer object to convert between the two coordinate systems
        self.map_to_rd = pyproj.Transformer.from_crs(crs_map, crs_rd)
        self.rd_to_map = pyproj.Transformer.from_crs(crs_rd, crs_map)
        
        self.path_keringen = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data\Combined_test_v14_WBD.gpkg')
        self.get_keringen()
        self.dambreak_layer = None

        self.set_default_layout_and_styling()
        
    def get_keringen(self):
        keringen = gpd.read_file(self.path_keringen,
                    layer = 'keringen',
                    crs=self.crs_rd)

        strings_to_match = ['hhsk', 'hhd', 'hdsr', 'hhr', 'wagv', 'ark']
        self.keringen = keringen[keringen['code'].str.startswith(tuple(strings_to_match))]
        self.keringen_plot = self.keringen.to_crs(4326)
        self.geojson_keringen = self.gdf_to_geojson(self.keringen_plot, layer_name = "keringen")
    
    def draw_map(self, center = [51.970682, 4.64013599]):
        self.m = Map(center=center, zoom=12)
        self.marker = Marker(location=center, draggable=True)
        self.m.add_layer(self.marker);
        self.m.add_layer(self.geojson_keringen)
        display(self.m)
        
        button_confirm = ipy.Button(description="Confirm dike breach location", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button_confirm, output)
        button_confirm.on_click(self.snap_marker)
        
        button_reset = ipy.Button(description="Reset dike breach location", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button_reset, output)
        button_reset.on_click(self.reset_map)
    
    def reset_map(self, b):
        clear_output(wait = True)
        if "self.dambreak" in locals():
            del self.dambreak
        
        if "self.marker.location" in locals():
            center = self.marker.location
            self.draw_map(center = center)
        else:
            self.draw_map()
    
    def gdf_to_geojson(self, input_gdf, layer_name):
        json_data = json.loads(input_gdf.to_json())
        geojson = ipl.GeoJSON(data=json_data, name = layer_name)
        return geojson
    
    def snap_marker(self, b):
        snap_rd, snap_map = self.marker_to_kering()
        self.marker.location = (snap_map.x, snap_map.y)
    
    def marker_to_kering(self):
        lat, lon = self.marker.location
        x, y = self.map_to_rd.transform(lat, lon)

        # find the nearest point on the line to the point
        point = Point(x, y)
        lines = self.keringen

        # build a spatial index on the line GeoDataFrame
        closest_line = lines.iloc[lines.distance(point).argmin()]
        nearest = nearest_points(point, closest_line['geometry'])

        snap_rd = nearest[1]
        lat, lon = self.rd_to_map.transform(snap_rd.x, snap_rd.y)
        snap_map = Point(lat, lon)
        self.draw_perpendicular_line(closest_line, snap_rd)
        return [snap_rd, snap_map]

    def draw_perpendicular_line(self, closest_line, snap_point):
        best_dist = 99999999
        best_index = None
        for i, coord in enumerate(closest_line.geometry.coords):
            dist = ((snap_point.x - coord[0])**2 + (snap_point.y - coord[1])**2)**0.5
            if dist < best_dist:
                best_dist = dist
                best_index = i
        points_of_interest = closest_line.geometry.coords[max(best_index - 1, 0): min(len(closest_line.geometry.coords), best_index+1)]
        
        points_of_interest = [[p[0], p[1]] for p in points_of_interest]

        point_on_line = [snap_point.x, snap_point.y]

        perp_line = self.calc_perpendicular_line(points_of_interest, point_on_line, 500)
        
        line_geom = LineString(perp_line)
        gdf = gpd.GeoDataFrame(geometry=[line_geom], crs = self.crs_rd)
        self.dambreak = gdf
        gdf = gdf.to_crs(self.crs_map)
        
        if self.dambreak_layer is None:
            self.geojson_dambreak = self.gdf_to_geojson(gdf, layer_name = "dambreak")
            self.dambreak_layer = self.m.add_layer(self.geojson_dambreak)
        else:
            self.geojson_dambreak = self.gdf_to_geojson(gdf, layer_name = "dambreak")
            dambreak_layer = self.m.add_layer(self.geojson_dambreak)

    def calc_perpendicular_line(self, points_of_interest, point_on_line, len_perpendicular):
        x1, y1 = points_of_interest[0]
        x2, y2 = points_of_interest[1]
        x_point, y_point = point_on_line
        if x2 - x1 == 0:
            slope = 999999999
            perp_slope = 0
        elif y2-y1 == 0:
            slope = 0
            perp_slope = 999999999
        else:
            slope = (y2 - y1) / (x2 - x1)
            perp_slope = -1/slope

        # Normalize the vector representing the line
        dx = 1 / ((1 + perp_slope ** 2)**0.5)
        dy = perp_slope * dx

        # Calculate the endpoints of the perpendicular line segment
        x_start = x_point - dx/2 * len_perpendicular
        y_start = y_point - dy/2* len_perpendicular
        x_end = x_point + dx/2* len_perpendicular
        y_end = y_point + dy/2* len_perpendicular
        return [[x_start, y_start], [x_end, y_end]]
    


# 1.3
# from notebooks.background_scripting.v1.modify_dambreak import DambreakWidget

# dambreak_input = DambreakWidget()   
# dambreak_input.draw_map()

# 1.4
# from notebooks.background_scripting.v1.modify_dambreak import ModifyDambreakSettingsDHYDRO

# dambreak_settings_dhydro = ModifyDambreakSettingsDHYDRO(model_path)
# if 'dambreak_input' in globals():
#     if hasattr(dambreak_input, 'dambreak'):
#         dambreak_settings_dhydro.modify_dikebreak_location(dambreak_input.dambreak, dambreak_input.keringen)
# dambreak_settings_dhydro.display_widgets()