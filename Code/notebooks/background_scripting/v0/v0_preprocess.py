import ipywidgets as widgets
import os
import json
from IPython.display import display, clear_output
from data_structures.dhydro_data import DHydroData
import subprocess
import warnings
warnings.filterwarnings('ignore')


class ModelSettings:
    """
    Class that contains all the settings that need to be used to create the model.
    Also contians the functions that communicate with the widgets from the jupyter notebook. 
    """
    def __init__(self):
        self.options_dict = {'HHSK':'HHSK_clipped_wnp', 'HDSR':'HDSR_clipped_test', 'HHD':'HHD_clipped_wnp', 'HHR':'HHR_clipped_wnp',
               "WAGV":'WAGV_clipped', 'ARKNZK':'ARKNZK', 'Rijntakken': 'Rijntakken', 'RMM':'RMM'}
        self.config_dict = {'HHSK':'hhsk', 'HDSR':'hdsr', 'HHD':'hhd', 'HHR':'hhr',
                       "WAGV":'wagv', 'ARKNZK':'ark_nzk', 'Rijntakken': 'rijntakken', 'RMM':'rijnmaasmonding'}
        

        # initiate default values for settings
        self.settings = {}
        self.settings['folder'] = r"D:\Work\Project\P1414"
        self.settings['output_folder'] = r"D:\Work\Project\P1414"
        self.settings['build_database'] = False
        self.settings['load_gpkgs'] = False
        self.settings['build_model'] = False
        self.settings['gpkgs_to_include'] = list(self.options_dict.keys())
        self.settings_to_print = ['folder', 'build_database', 'load_gpkgs', 'build_model', 'gpkgs_to_include', "output_folder"]
        self.update_gpkgs()

        # list of widgets
        layout = widgets.Layout(width='50%')
        style = {'description_width': '50%'}
        
        self.folder_widget = widgets.Text(
            value= self.settings['folder'],
            placeholder= self.settings['folder'],
            description='Path to GIS folder:',
            disabled=False,
            layout = layout,
            style = style
        )

        self.output_folder_widget = widgets.Text(
            value= self.settings['output_folder'],
            placeholder= self.settings['output_folder'],
            description='Export folder for model:',
            disabled=False,
            layout = layout,
            style = style
        )

        self.loc_choice = widgets.SelectMultiple(options=list(self.options_dict.keys()),
            value=self.settings['gpkgs_to_include'],
            description='gpkgs to include:',
            disabled=False, 
            rows = len(self.settings['gpkgs_to_include']),
            layout = layout,
            style = style
            )

        self.build_database_widget = widgets.Checkbox(
            value=self.settings['build_database'],
            description='Build database',
            disabled=False,
            layout = layout,
            style = style
        )

        self.load_gpkgs_widget = widgets.Checkbox(
            value=self.settings['load_gpkgs'],
            description='Load gpkgs',
            disabled=False,
            layout = layout,
            style = style
        )

        self.build_model_widget = widgets.Checkbox(
            value=self.settings['build_model'],
            description='Build model',
            disabled=False,
            layout = layout,
            style = style
        )
        
    def update_gpkgs(self):
        self.settings['gpkgs_list']= [os.path.join(self.settings['folder'], 'GIS\HYDAMO', f"{self.options_dict[loc]}.gpkg") for loc in self.settings['gpkgs_to_include']]
        self.settings['config_list'] = [f"{self.config_dict[loc]}_config" for loc in self.settings['gpkgs_to_include']]
        

    def update_settings_dict(self, b):
        self.settings['folder'] = self.folder_widget.value
        self.settings['output_folder'] = self.output_folder_widget.value
        self.settings['build_database'] = self.build_database_widget.value
        self.settings['load_gpkgs'] = self.load_gpkgs_widget.value
        self.settings['build_model'] = self.build_model_widget.value
        self.settings['gpkgs_to_include'] = self.loc_choice.value
        self.update_gpkgs()
        
        clear_output(wait=True)
        self.display_widgets()
        display("Model settings are:")
        print_settings_dict = {k: self.settings[k] for k in self.settings_to_print}
        print(json.dumps(print_settings_dict, indent=4))
        
    def display_widgets(self):
        display(self.folder_widget)
        display(self.output_folder_widget)
        display(self.loc_choice)
        display(self.build_database_widget)
        display(self.load_gpkgs_widget)
        display(self.build_model_widget)
        button = widgets.Button(description="Update settings")
        output = widgets.Output()
        display(button, output)
        button.on_click(self.update_settings_dict)


class ModelFromNotebook:
    """
    Class that can build the database, load gpkgs, and build the model. 
    This class immediately starts performing these activities based on the model settings.
    This class also can run the model on the users pc
    """
    def __init__(self, model_settings):
        self.settings = model_settings
        self.gpkg_file = self.settings['folder'] + r"\GIS\HYDAMO\Combined_test_v7.6.gpkg"
        self.output_folder =  self.settings['output_folder'] + r"\Model"
        self.config_dhydro = r"combined_config"
        self.snap_dist_list = [0, 0, 10, 10, 50, 10, 10, 100]
        self.defaults = r"defaults"

    def build_model(self):

        if self.settings['build_database']:
            out = widgets.HTML()
            progress = widgets.IntProgress(min=0, max=len(self.settings['config_list']), description = "Progress:")
            display(out)
            display(progress)

            dhd = DHydroData()
            for ix, config in enumerate(self.settings['config_list']):
                out.value = f"<h3>Building database for config: {config} ({ix+1} out of {len(self.settings['config_list'])})</h3>"
                progress.value = ix
                print(config)

                dhd.hydamo_from_raw_data(
                    defaults=self.defaults, config=config, branch_snap_dist=self.snap_dist_list[ix]
                )
                try:
                    dhd.fixed_weirs_from_raw_data(config=config, defaults=self.defaults)
                except AttributeError:
                    pass
            progress.value += 1
            dhd.clip_structures_by_branches()
            dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=self.defaults, min_length=500)
            dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=self.defaults, min_length=500)
            dhd.hydamo_to_gpkg(output_gpkg=self.gpkg_file)

        if self.settings['load_gpkgs']:
            progress = widgets.IntProgress(min=0, max=len(self.settings['gpkgs_list']), description = "Progress:")
            out = widgets.HTML()
            display(progress)
            display(out)

            dhd = DHydroData()
            for ix, gpkg in enumerate(self.settings['gpkgs_list']):
                out.value = f"<h3>Loading gpkgs for: {gpkg} ({ix+1} out of {len(self.settings['gpkgs_list'])})</h3>"
                progress.value = ix

                # 2. load data
                dhd.hydamo_from_gpkg(gpkg, branch_snap_dist=self.snap_dist_list[ix])
            progress.value += 1
            dhd.fixed_weirs_from_raw_data(config="wegen_config", defaults=self.defaults, min_length=500)
            dhd.fixed_weirs_from_raw_data(config="relief_config", defaults=self.defaults, min_length=500)
            dhd.hydamo_to_gpkg(output_gpkg=self.gpkg_file)

        if  self.settings['build_model']:
            out = widgets.HTML("<h3>Building model in D-HYDRO...</h3>")
            display(out)
            # 1. initialize an instance of DHydamoData
            dhd = DHydroData()

            # 2. load data
            dhd.hydamo_from_gpkg(self.gpkg_file)

            # remove brug as it needs a cs
            del dhd.ddm.brug
            dhd.features.remove("brug")

            # 3. save as dhydro model
            dhd.to_dhydro(config=self.config_dhydro , output_folder=self.output_folder)

    def run_model(self):
        nr_log_lines = 10
        full_logs = []
        out_title = widgets.HTML(f"<p>Last {nr_log_lines} logging lines: </p>")
        out = widgets.HTML("")
        display(out_title)
        display(out)

        filepath = os.path.join(self.output_folder)
        with subprocess.Popen(os.path.join(filepath, 'run.bat'), cwd = filepath ,  stdout=subprocess.PIPE, bufsize=1, universal_newlines=True) as p:
            for b in p.stdout:
                full_logs.append(b)
                new_logging = ""
                for line in full_logs[-nr_log_lines:]:
                    new_logging += f"<p>{line}</p>"
                out.value = new_logging
        
        if p.returncode != 0:
            raise subprocess.CalledProcessError(p.returncode, p.args)
        