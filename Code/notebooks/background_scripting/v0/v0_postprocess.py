import os
import sys
from pathlib import Path

from flowmeshreader import load_meta_data, load_map_data, mesh_to_tiff, load_mesh
from plotting import raster_plot_with_context

import ipywidgets as widgets
from IPython.display import display, clear_output
from netCDF4 import Dataset
import json
import numpy as np


class PlotSettings():
    """
    Class that contains all the widgetes and settings for plotting model results.
    Also contians the functions that communicate with the widgets from the jupyter notebook. 
    """
    def __init__(self, model_container):
        self.settings = {}
        self.settings['input_file_path'] = model_container.output_folder + r"\dflowfm\DFM_OUTPUT_test\test_map.nc"
        self.settings['output_file_path'] = model_container.output_folder + r"\post_processing/output.tiff"
        self.settings['mesh_resolution'] = self.find_minimum_resolution_from_mesh()
        self.settings['timesteps_in_hours'] = self.get_timesteps()
        self.settings['timstep'] = self.settings['timesteps_in_hours'][-1]
        
        self.settings['variables'] = load_meta_data(self.settings['input_file_path'])
        self.variables_map = [x for x in self.settings['variables'] if 'flowelem' not in x.lower()]
        self.settings['plot_variable'] = 'Mesh2d_waterdepth' if 'Mesh2d_waterdepth' in self.settings['variables'] else self.settings['variables'][0]
        
        self.settings['aggregation_type_options'] =  ['maximum', 'minimum', 'timestep']
        self.settings['aggregation_type'] = 'timestep'
        self.settings_to_print = ['output_file_path', 'mesh_resolution', 'timestep', 'plot_variable', 'aggregation_type']
        
        # list of widgets
        layout = widgets.Layout(width='50%')
        style = {'description_width': '50%'}

        self.output_folder_images =  widgets.Text(
            value= self.settings['output_file_path'],
            placeholder= self.settings['output_file_path'],
            description='Path to GIS folder',
            disabled=False,
            layout = layout,
            style = style
        )
        
        self.nc_variable_widget = widgets.Dropdown(
            options=self.settings['variables'],
            value=self.settings['plot_variable'],
            description='Variable to plot:',
            disabled=False,
            layout = layout,
            style = style
        )

        self.output_resolution = widgets.FloatText(
            value = self.settings['mesh_resolution'],
            description='Output resolution (m):',
            disabled=False,
            layout = layout,
            style = style
        )

        self.aggregation_type = widgets.Dropdown(
            options= self.settings['aggregation_type_options'],
            value=self.settings['aggregation_type'] ,
            description='Min, max or timestep:',
            disabled=False,
            layout = layout,
            style = style
        )

        self.timestep_selector = widgets.FloatSlider(
            value= self.settings['timstep'],
            min=self.settings['timesteps_in_hours'][0],
            max=self.settings['timesteps_in_hours'][-1],
            step=self.settings['timesteps_in_hours'][1] - self.settings['timesteps_in_hours'][0],
            description='Timestep (hours): ',
            disabled= False,
            layout = layout,
            style = style
        )
        
    def display_widgets(self):
        display(self.nc_variable_widget)
        display(self.output_resolution)
        display(self.aggregation_type)
        display(self.timestep_selector)
        button = widgets.Button(description="Plot results")
        output = widgets.Output()
        display(button, output)
        button.on_click(self.update_settings_dict)

    def find_minimum_resolution_from_mesh(self, nr_points:int = 500):
        """
        Find the distance between cells in the mesh.
        Finding distance for all cells would take a short while, so only use the first x nr_points.
        You will probablty find the minimum at these first points, or at least a good estimate 
        This is used as a default value for the mesh resolution.
        """
        x, y = load_mesh(self.settings['input_file_path'])
        points = np.array(tuple(zip(x,y)))[:nr_points]
        # find the distances between all points. This is a vectorized method to do this:
        distances = np.sqrt(np.sum((points[:,np.newaxis,:] - points[np.newaxis,:,:])**2, axis=-1))
        distances = np.array(distances).flatten()
        distances = distances[distances != 0]
        return np.min(distances)

    def get_timesteps(self):
        nc = Dataset(self.settings['input_file_path'])
        timesteps_in_seconds = nc.variables['time']
        timesteps_in_hours = np.divide(timesteps_in_seconds, 3600)
        return timesteps_in_hours

    def update_settings_dict(self, b):
        self.settings['output_file_path'] = self.output_folder_images.value
        self.settings['mesh_resolution'] = self.output_resolution.value
        self.settings['timestep'] = round(self.timestep_selector.value,2)
        self.settings['plot_variable'] = self.nc_variable_widget.value
        self.settings['aggregation_type'] = self.aggregation_type.value
        
        clear_output(wait=True)
        self.display_widgets()
        display("Plot settings are:")
        print_settings_dict = {k: v for k, v in self.settings.items() if k in self.settings_to_print}
        if self.settings['aggregation_type'] != 'timestep':
            print_settings_dict['timestep'] = "not applicable" 
        print(json.dumps(print_settings_dict, sort_keys=True, indent=4))
        ModelPlotter(self.settings)
        
class ModelPlotter():
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
        map_data = load_map_data(self.settings['input_file_path'], self.settings['plot_variable'])
        out = widgets.HTML("<h3>Generating plot...</h3>")
        display(out)

        if self.settings['aggregation_type'] == 'minimum':
            data = np.min(map_data, axis = 0)
            title = f"{self.settings['aggregation_type']} {self.settings['plot_variable']}"
        elif self.settings['aggregation_type'] == 'maximum':
            data = np.max(map_data, axis = 0)
            title = f"{self.settings['aggregation_type']} {self.settings['plot_variable']}"
        elif self.settings['aggregation_type'] == 'timestep':
            t_index = (np.abs(self.settings['timesteps_in_hours'] - self.settings['timestep'])).argmin()
            data = map_data[t_index, :]
            title = f"{self.settings['plot_variable']} after {round(self.settings['timestep'],2)} hours"
        else:
            raise Exception("Error in selection of plot settings")

        # convert to raster and save as tiff
        _, _, grid_data = mesh_to_tiff(
            data,
            self.settings['input_file_path'],
            self.settings['output_file_path'],
            self.settings['mesh_resolution'],
            self.distance_tol,
            interpolation=self.interpolation,
        )
        fig, ax = raster_plot_with_context(
            raster_path = self.settings['output_file_path'], 
            epsg = 28992, 
            clabel = f"{self.settings['plot_variable']}", 
            cmap = "Blues", 
            title = title,
            )
        out.value = ""

  