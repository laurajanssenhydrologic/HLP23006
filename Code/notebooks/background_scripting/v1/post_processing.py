from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import os
import ipywidgets as ipy
from IPython.display import display, clear_output

from flowmeshreader import load_meta_data, load_map_data, mesh_to_tiff, load_mesh
from plotting import raster_plot_with_context

from IPython.display import display, clear_output
import netCDF4 as nc
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def read_output_folder(model_folder):
    files = os.listdir(os.path.join(model_folder, 'dflowfm'))
    mdu_file = [file for file in files if file.endswith('.mdu')][0]
    mdu_file = os.path.join(model_folder, 'dflowfm', mdu_file)
    model_name =  os.path.basename(mdu_file).split('.')[0]
    fm_folder = os.path.dirname(mdu_file)

    # find the outputdir according to the mdu file
    output_dir_lines = []
    with open(mdu_file, 'r') as f:
        for line in f.readlines():
            if line.lower().strip().startswith('outputdir'):
                output_dir_lines.append(line)

    if len(output_dir_lines) > 1:
        raise Exception("Found multiple output dir in mdu file!")

    if len(output_dir_lines) == 0:
        output_dir = os.path.join(fm_folder, "DFM_OUTPUT_{}".format(model_name))
    elif len(output_dir_lines) == 1:
        output_dir_line = output_dir_lines[0]
        if '#' in output_dir_line:
            output_dir_line = output_dir_line.split('#')[0]
        output_dir_name = output_dir_line.split('=')[-1].strip()
        if output_dir_name == '.':
            output_dir = fm_folder
        else:
            output_dir = os.path.join(fm_folder, output_dir_name)
    return output_dir


class PlotSettingsBreach(WidgetStyling):
    """
    Class that contains all the widgetes and settings for plotting the dike breach information.
    """
    def __init__(self, model_folder):
        output_folder = read_output_folder(model_folder)
        self.settings = {}
        his_files = [x for x in os.listdir(output_folder) if x.endswith('his.nc')]
        if len(his_files) != 1:
            raise Exception(f"Found {len(his_files)} his files in {output_folder}, should be 1. Are you sure the model has finished computing? If so, check if you can find the output in the model folder, which should be {model_folder}")
        self.settings['his_path'] = os.path.join(output_folder, his_files[0])
        self.settings['output_file_path'] = os.path.join(output_folder,  "post_processing")

        plot_variables = self.load_plot_variables()
        self.settings['variables'] = plot_variables

        default_variables = [var for var in plot_variables if var.endswith('discharge')]
        self.settings['plot_variables'] = default_variables

        self.settings_to_print = ['his_path', 'output_file_path', 'plot_variables']

        # list of widgets
        self.set_default_layout_and_styling()

        self.output_folder_images =  ipy.Text(
            value= self.settings['output_file_path'],
            placeholder= self.settings['output_file_path'],
            description='Path to save images',
            disabled=False,
        )
        
        self.variables_to_plot = ipy.SelectMultiple(
            options=self.settings['variables'],
            value=self.settings['plot_variables'],
            description='Dike breach variables to plot:',
            disabled=False, 
            rows = len(self.settings['variables']),
            )

        self.widgets_to_display = [self.output_folder_images, self.variables_to_plot]

        for widget in self.widgets_to_display:
            widget.layout = self.item_layout
            widget.style = self.item_style

    def load_plot_variables(self):
        plot_variables = []
        nc_file = nc.Dataset(self.settings['his_path'])
        variables = nc_file.variables.keys()
        t = nc_file.variables['time'][:]
        for var in variables:
            if var.startswith('dambreak'):
                data = nc_file.variables[var][:]
                try:
                    if len(data) == len(t):
                        plot_variables.append(var)
                except:
                    pass
        nc_file.close()
        return plot_variables

    def display_widgets(self):
        items_map = ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style)
        display(items_map)

        button = ipy.Button(description="Plot breach information", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.plot_map)

    def update_settings_dict(self):
        self.settings['output_file_path'] = self.output_folder_images.value
        self.settings['plot_variables'] = self.variables_to_plot.value
    
    def plot_map(self, b):
        self.update_settings_dict()
        clear_output(wait=True)
        self.display_widgets()
        display("Plot settings are:")
        print_settings_dict = {k: v for k, v in self.settings.items() if k in self.settings_to_print}
        print(json.dumps(print_settings_dict, sort_keys=True, indent=4))
        BreachPlotter(self.settings)

class PlotSettingsMap(WidgetStyling):
    """
    Class that contains all the widgetes and settings for plotting model results on the 2D map.
    """
    def __init__(self, model_folder):
        output_folder = read_output_folder(model_folder)
        self.settings = {}
        map_files = [x for x in os.listdir(output_folder) if x.endswith('map.nc')]
        if len(map_files) != 1:
            raise Exception(f"Found {len(map_files)} map files in {output_folder}, should be 1. Are you sure the model has finished computing? If so, check if you can find the output in the model folder, which should be {model_folder}")
        self.settings['map_path'] = os.path.join(output_folder, map_files[0])
        self.settings['output_file_path'] = os.path.join(output_folder,  "post_processing")

        self.settings['mesh_resolution'] = self.find_minimum_resolution_from_mesh()
        self.settings['timesteps_in_hours'] = self.get_timesteps()
        
        self.settings['timestep'] = self.settings['timesteps_in_hours'][-1]
        
        self.settings['variables'] = load_meta_data(self.settings['map_path'])
        self.variables_map = [x for x in self.settings['variables'] if 'flowelem' not in x.lower()]
        self.settings['plot_variable'] = 'Mesh2d_waterdepth' if 'Mesh2d_waterdepth' in self.settings['variables'] else self.settings['variables'][0]
        
        self.settings['aggregation_type_options'] =  ['maximum', 'minimum', 'timestep']
        self.settings['aggregation_type'] = 'timestep'

        self.settings['min_value_legend'] = 0
        self.settings['max_value_legend'] = 1
        self.settings['color_map_options'] = ['viridis', 'hsv', 'rainbow', 'plasma', 'inferno', 'magma', 'Blues', 'Reds']
        self.settings['color_map'] = self.settings['color_map_options'][0]

        self.settings_to_print = ['output_file_path', 'mesh_resolution', 'timestep', 'plot_variable', 'aggregation_type', 'min_value_legend', 'max_value_legend', 'color_map']

        # list of widgets
        self.set_default_layout_and_styling()

        self.output_folder_images =  ipy.Text(
            value= self.settings['output_file_path'],
            placeholder= self.settings['output_file_path'],
            description='Path to save images',
            disabled=False,
        )
        
        self.nc_variable_widget = ipy.Dropdown(
            options=self.settings['variables'],
            value=self.settings['plot_variable'],
            description='Variable to plot:',
            disabled=False,
        )

        self.output_resolution = ipy.FloatText(
            value = self.settings['mesh_resolution'],
            description='Output resolution (m):',
            disabled=False,
        )

        self.aggregation_type = ipy.Dropdown(
            options= self.settings['aggregation_type_options'],
            value=self.settings['aggregation_type'] ,
            description='Min, max or timestep:',
            disabled=False,
        )

        self.timestep_selector = ipy.FloatSlider(
            value= self.settings['timestep'],
            min=self.settings['timesteps_in_hours'][0],
            max=self.settings['timesteps_in_hours'][-1],
            step=self.settings['timesteps_in_hours'][1] - self.settings['timesteps_in_hours'][0],
            description='Timestep (hours): ',
            disabled= False,
        )

        self.min_value_legend = ipy.FloatText(
            value= self.settings['min_value_legend'],
            description='Minimum value legend:',
            disabled= False,
        )

        self.max_value_legend = ipy.FloatText(
            value= self.settings['max_value_legend'],
            description='Maximum value legend:',
            disabled= False,
        )

        self.colormap_picker = ipy.Dropdown(
            options=self.settings['color_map_options'],
            value=self.settings['color_map'],
            description='Colormap:',
            disabled=False,
        )

        self.widgets_to_display_map = [self.output_folder_images, self.nc_variable_widget, self.output_resolution, self.aggregation_type, self.timestep_selector,
                                       self.min_value_legend, self.max_value_legend, self.colormap_picker]

        for widget in self.widgets_to_display_map:
            widget.layout = self.item_layout
            widget.style = self.item_style

        
    def display_widgets(self):
        items_map = ipy.VBox(children=self.widgets_to_display_map, layout=self.box_layout, style = self.box_style)
        display(items_map)

        button_map = ipy.Button(description="Plot map", layout = self.button_layout, style = self.button_style)
        output = ipy.Output()
        display(button_map, output)
        button_map.on_click(self.plot_map)

    def find_minimum_resolution_from_mesh(self, nr_points:int = 500):
        """
        Find the distance between cells in the mesh.
        Finding distance for all cells would take a short while, so only use the first x nr_points.
        You will probablty find the minimum at these first points, or at least a good estimate 
        This is used as a default value for the mesh resolution.
        """
        x, y = load_mesh(self.settings['map_path'])
        points = np.array(tuple(zip(x,y)))[:nr_points]
        # find the distances between all points. This is a vectorized method to do this:
        distances = np.sqrt(np.sum((points[:,np.newaxis,:] - points[np.newaxis,:,:])**2, axis=-1))
        distances = np.array(distances).flatten()
        distances = distances[distances != 0]
        return np.min(distances)

    def get_timesteps(self):
        """
        Get the timesteps that are present in the netcdf file
        """
        nc_file = nc.Dataset(self.settings['map_path'])
        timesteps_in_seconds = nc_file.variables['time']
        timesteps_in_hours = np.divide(timesteps_in_seconds, 3600)
        return timesteps_in_hours

    def update_settings_dict(self):
        """
        Update the settings dictionary based on the values of the widget
        """
        self.settings['output_file_path'] = self.output_folder_images.value
        self.settings['mesh_resolution'] = self.output_resolution.value
        self.settings['timestep'] = round(self.timestep_selector.value,2)
        self.settings['plot_variable'] = self.nc_variable_widget.value
        self.settings['aggregation_type'] = self.aggregation_type.value
        self.settings['min_value_legend'] = self.min_value_legend.value
        self.settings['max_value_legend'] = self.max_value_legend.value
        self.settings['color_map'] = self.colormap_picker.value
    
    def plot_map(self, b):
        """
        Funciton activated on button press. Calls the map plotter. 
        """
        self.update_settings_dict()
        clear_output(wait=True)
        self.display_widgets()
        display("Plot settings are:")
        print_settings_dict = {k: v for k, v in self.settings.items() if k in self.settings_to_print}
        if self.settings['aggregation_type'] != 'timestep':
            print_settings_dict['timestep'] = "not applicable" 
        print(json.dumps(print_settings_dict, sort_keys=True, indent=4))
        MapPlotter(self.settings)
      

class BreachPlotter():
    """
    Class that is responsible for generating the plot of the dikebreach
    """
    def __init__(self, plot_settings):
        self.settings = plot_settings
        self.plot_output()   

    def plot_output(self):
        """
        Plot the dikebreach characteristics that the user selected
        """
        nc_file = nc.Dataset(self.settings['his_path'])
        t = nc_file.variables['time'][:]
        
        fig, axes = plt.subplots(
            figsize = (12, 5 + 4 * len(self.settings['plot_variables'])),
            nrows = len(self.settings['plot_variables']),
            sharex = True)

        if len(self.settings['plot_variables']) == 1:
            axes = [axes]
        elif len(self.settings['plot_variables']) == 0:
            return 0

        for i, var in enumerate(self.settings['plot_variables']):
            data = nc_file.variables[var][:]
            axes[i].plot(t, data, color = '#3587A4')
            axes[i].set_title(var)
            axes[i].set_xlim([min(t), max(t)])
            if min(data) == 0.0:
                axes[i].set_ylim(bottom = 0)
            axes[i].set_title(var)

        if os.path.exists(self.settings['output_file_path']) == False:
            os.makedirs(self.settings['output_file_path'])
        plt.savefig(os.path.join(self.settings['output_file_path'], 'breach.png'))
        nc_file.close()

class MapPlotter():
    """
    Class responsible for generating plots based on plot settings
    """
    def __init__(self, plot_settings):
        self.settings = plot_settings
        self.distance_tol = self.settings['mesh_resolution'] * 4
        self.interpolation = r"nearest"
        
        self.plot_output()       
        
    def plot_output(self):
        """
        Plot the output that is wanted by the user.
        """
        # load mesh coordinates and data from netCDF 
        map_data = load_map_data(self.settings['map_path'], self.settings['plot_variable'])
        map_data[map_data < 0.01] = np.nan
        out = ipy.HTML("<p>Generating plot... (first plot can take a while)</p>")
        display(out)

        if self.settings['aggregation_type'] == 'minimum':
            data = np.nanmin(map_data, axis = 0)
            start_depth = map_data[1,:]
            data[np.where(start_depth > 0.01)] = np.nan
            title = f"{self.settings['aggregation_type']} {self.settings['plot_variable']}"
        elif self.settings['aggregation_type'] == 'maximum':
            data = np.nanmax(map_data, axis = 0)
            start_depth = map_data[1,:]
            data[np.where(start_depth > 0.01)] = np.nan
            title = f"{self.settings['aggregation_type']} {self.settings['plot_variable']}"
        elif self.settings['aggregation_type'] == 'timestep':
            t_index = (np.abs(self.settings['timesteps_in_hours'] - self.settings['timestep'])).argmin()
            data = map_data[t_index, :]
            title = f"{self.settings['plot_variable']} after {round(self.settings['timestep'],2)} hours"
        else:
            raise Exception("Error in selection of plot settings")

        raster_path = os.path.join(self.settings['output_file_path'], f"{self.settings['plot_variable']}.tiff")
        if os.path.exists(self.settings['output_file_path']) == False:
            os.makedirs(self.settings['output_file_path'])

        # convert to raster and save as tiff
        _, _, grid_data = mesh_to_tiff(
            data,
            self.settings['map_path'],
            raster_path,
            self.settings['mesh_resolution'],
            self.distance_tol,
            interpolation=self.interpolation,
        )
        fig, ax = raster_plot_with_context(
            raster_path = raster_path,
            epsg = 28992, 
            clabel = f"{self.settings['plot_variable']}", 
            cmap = self.settings['color_map'], 
            vmin = self.settings['min_value_legend'],
            vmax = self.settings['max_value_legend']
            )
        plt.savefig(os.path.join(self.settings['output_file_path'], f"{self.settings['plot_variable']}.png"))
        out.value = ""