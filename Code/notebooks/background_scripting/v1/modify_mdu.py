import os
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import json
from datetime import datetime
import shutil

def copy_model(path_model:str, scenario_name:str):
    """
    Copy a model from the database to the model_run folder, such that it can be used and modified

    Args:
        path_model (str): path to the model
        scenario_name (str): name of the scenario

    Returns:
        str: path to the new folder of the model
    """
    now = datetime.now()
    formatted_datetime = now.strftime("%Y-%m-%dT%H-%M-%S")
      
    new_name = os.path.basename(path_model) + f"_{formatted_datetime}_{scenario_name}"
    model_run_path =  os.path.join(os.path.dirname(os.path.dirname(path_model)), 'Model_runs')
    if os.path.exists(model_run_path) == False:
        os.mkdir(model_run_path)
    new_path =os.path.join(model_run_path, new_name)
    destination = shutil.copytree(path_model, new_path, copy_function = shutil.copy) 
    return new_path

class ModifyMDU(WidgetStyling):
    """
    Class for modifying the MDU settings of the model
    """
    def __init__(self, model_folder):
        files = os.listdir(os.path.join(model_folder, 'dflowfm'))
        mdu_file = [file for file in files if file.endswith('.mdu')][0]
        self.mdu_path = os.path.join(model_folder, 'dflowfm', mdu_file)

        self.run_bat_file = os.path.join(model_folder, 'run.bat')
        backup_run_bat = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), r'data\\run.bat')
        if os.path.isfile(self.run_bat_file) == False:
            shutil.copy(backup_run_bat, self.run_bat_file)

        self.mdu_lines = self.read_mdu()
        
        self.modify_parameter_mdu(parameter = "statsinterval", new_value = 1) # always do this

        self.settings_to_modify = ['tStart', 'tStop', 'mapInterval']
        self.settings_names = {
            'tStart': 'Model start time (tStart in hours)',
            'tStop': "Model end time (tStop in hours)", 
            'mapInterval': "Interval to write map file (mapInterval in minutes)"
            }
        self.settings_in_hours = ['tStart', 'tStop']
        self.settings_in_minutes = ['mapInterval']
        self.settings = {}
        self.settings['refDate'] = datetime.strptime(self.read_parameter_mdu('refDate'), '%Y%m%d')
        
        self.set_default_layout_and_styling()

        self.widgets = {}
        for setting in self.settings_to_modify:
            self.settings[setting] =  float(self.read_parameter_mdu(setting))
            val = self.convert_to_sas(setting, self.settings[setting])
            self.widgets[setting] = ipy.IntText(
                value = val,
                description=f'{self.settings_names[setting]}:',
                disabled=False
                )
            
            self.widgets[setting].layout = self.item_layout
            self.widgets[setting].style = self.item_style
        
        self.settings['DHYDRO location'] = self.read_run_bat()
        self.widgets['DHYDRO location'] = ipy.Text(
            value = self.settings['DHYDRO location'],
            description=f'Location of dhydro run_dimr:',
            disabled=False,
            style = self.item_style,
            layout = self.item_layout
            )

        self.widgets_to_display = [self.widgets[setting] for setting in self.widgets.keys()]
    
    def convert_to_sas(self, key, value):
        """
        Convert unit from seconds to hours or minutes
        """
        if key in self.settings_in_hours:
            return value / 60 / 60
        if key in self.settings_in_minutes:
            return value / 60
        return value
    
    def convert_to_mdu(self, key, value):
        """
        Convert unit from hours or minutes to seconds
        """
        if key in self.settings_in_hours:
            return value * 60 * 60
        if key in self.settings_in_minutes:
            return value * 60
        return value

    def read_run_bat(self):
        """
        read run.bat line that calls dhydro
        """
        with open(self.run_bat_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("call "):
                path =  line.split('"')[1]
        return path

    def modify_run_bat(self):
        """
        modify run.bat line that calls dhydro
        """
        with open(self.run_bat_file, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("call "):
                lines[i] = 'call \"{}\"\n'.format(self.settings['DHYDRO location'])
        with open(self.run_bat_file, 'w') as f:
            for line in lines:
                f.write(line)

    def search_string_in_file(self, string:str, lines:list):
        """
        Function to find the lines that start with a certain string

        Args:
            string (str): string to find
            lines (list): lines to find string in 

        Returns:
            indices and string that is found
        """
        line_numbers, strings = [], []
        for i, line in enumerate(lines):
            if line.strip(' ').lower().startswith(string.lower()):
                line_numbers.append(i)
                strings.append(line)
        if len(strings) == 1:
            return line_numbers[0], strings[0]
        else:
            return None
        
    def read_mdu(self):
        """
        read all lines of the mdu file
        """
        with open(os.path.join(self.mdu_path), 'r') as f:
            lines = f.readlines()
        return lines
    
    def save_mdu(self):
        """
        save all lines of the mdu file
        """
        with open(os.path.join(self.mdu_path), 'w') as f:
            for line in self.mdu_lines:
                f.write(line)
    
    def modify_parameter_mdu(self, parameter:str, new_value:str):
        """
        Modify a parameter in the mdu dictionary

        Args:
            parameter (str): parameter name
            new_value (str): new value
        """
        search_string = f"{parameter} " # add spaces to prevent accidentaly having wrong value
        index, line = self.search_string_in_file(search_string, self.mdu_lines)
        key = line.split("=")[0]
        value = f"= {new_value}"
        self.mdu_lines[index] = f"{key}{value}\n"
    
    def read_parameter_mdu(self, parameter):
        search_string = f"{parameter} " # add spaces to prevent accidentaly having wrong value
        index, line = self.search_string_in_file(search_string, self.mdu_lines)
        value = line.split("=")[-1].strip(' ')
        if '\n' in value:
            value = value.replace('\n', '')
        return value

    def display_widgets(self):
        items = ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style)
        display(items)
        
        button = ipy.Button(description="Update settings", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.update_settings_widget)
    
    def update_settings_widget(self, b):
        for setting in self.settings_to_modify:
            widget_value = self.widgets[setting].value 
            self.settings[setting] = self.convert_to_mdu(setting, widget_value)
            self.modify_parameter_mdu(parameter = setting, new_value = self.settings[setting])

        self.settings['DHYDRO location'] = self.widgets['DHYDRO location'].value.strip('"')
        self.widgets['DHYDRO location'].value = self.widgets['DHYDRO location'].value.strip('"')
        self.modify_run_bat()    
        
        self.save_mdu()
        clear_output(wait=True)
        self.display_widgets()
        display("MDU settings are:")
        print_settings_dict = {k: self.settings[k] for k in self.widgets.keys()}
        for key in print_settings_dict.keys():
            if key in self.settings_to_modify:
                print_settings_dict[key] = self.convert_to_sas(key, self.settings[key])
        print(json.dumps(print_settings_dict, indent=4))
