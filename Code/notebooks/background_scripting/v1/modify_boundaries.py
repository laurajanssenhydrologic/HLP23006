import os
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import ipywidgets as ipy
from IPython.display import display, clear_output
import json
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt

class ModifyBoundaries(WidgetStyling):
    """
    Class for modifying the boundary conditions of the model
    """
    def __init__(self, input_path, mdu_settings, dambreak_settings):
        self.modify_rhine_discharge = ModifyRhineDischarge(input_path, mdu_settings, dambreak_settings)

        self.model_folder = input_path
        files = os.listdir(os.path.join(self.model_folder, 'dflowfm'))
        for file in files:
            if file.lower().endswith('.ini') and 'initial' in file.lower() and 'water' in file.lower() and 'level' in file.lower():
                ic_file = file
            if 'bnd.ext' in file.lower():
                bnd_file = file
        self.ic_folder = os.path.join(self.model_folder, 'dflowfm', ic_file)
        self.bnd_folder = os.path.join(self.model_folder, 'dflowfm', bnd_file)

        self.line_indices = {}

        # init noordzee
        start_string_noord = 'branchId=noor_'
        self.line_indices['Waterlevel North Sea'] = self.get_inital_conditions_lines(start_string_noord)

        # init markermeer
        start_string_mark = 'branchId=mark_'
        self.line_indices['Waterlevel Markermeer'] = self.get_inital_conditions_lines(start_string_mark)

        self.set_default_layout_and_styling()

        self.read_functions = {
            'Waterlevel North Sea' : self.read_wl_initial,
            'Waterlevel Markermeer' : self.read_wl_initial
        }

        self.write_functions = {
            'Waterlevel North Sea' : self.write_wl_initial,
            'Waterlevel Markermeer' : self.write_wl_initial
        }
        
        self.settings = {}
        self.widgets = {}

        for setting in self.read_functions.keys():
            self.settings[setting] = self.read_functions[setting](self.line_indices[setting])

            self.widgets[setting] = ipy.FloatText(
                value =  self.settings[setting],
                description= setting,
                disabled=False
                )
            self.widgets[setting].layout = self.item_layout
            self.widgets[setting].style = self.item_style

        self.widgets_to_display = [widget for widget in self.widgets.values()]

    def find_start_string_in_lines(self, lines: list, start_string: str, stop_on_first_match:bool = False):
        """
        Find the occurence of a certain string in a list of strings

        Args:
            lines (list): lines to check
            start_string (str): wich line should the string start with
            stop_on_first_match (bool, optional): Wether to stop search when there is one match Defaults to False.

        Returns:
            indices of the lines in which the string is found
        """
        lines_found = []
        for i, line in enumerate(lines):
            if line.replace(' ', '').startswith(start_string):
                lines_found.append(i)
                if stop_on_first_match:
                    return i
        return lines_found

    def write_lines_to_file(self, file_path:str, lines:list):
        """Write lines to a text file

        Args:
            file_path (str): Path of the file to write
            lines (list): lines to write to the file

        Raises:
            Exception: in case the lines you provide are not the same length as the existin lines in the file
        """
        with open(file_path, 'r') as f:
            old_lines = f.readlines()

        if len(old_lines) == len(lines): #small check
            with open(file_path, 'w') as f:
                    for line in lines:
                        f.write(line)
        else:
            raise Exception("error")

    def get_inital_conditions_lines(self, start_string:str):
        """Find the line indices of the initial conditions that you want

        Args:
            start_string (str): which string the inital conditon should start with

        Raises:
            Exception: if the string was not found

        Returns:
            _type_: indices at which the inital condition is located
        """
        with open(self.ic_folder) as f:
            lines = f.readlines()

        lines_with_noordzee_branch = self.find_start_string_in_lines(
            lines,
            start_string = start_string)
        if len(lines_with_noordzee_branch) == 0:
            raise Exception(f"{len(lines_with_noordzee_branch)} occurences found with start string: {start_string}")

        indices_noordzee = []
        for line_index in lines_with_noordzee_branch:
            # search for discharge, start searching from the point where the branch name was found
            start_string = 'values='
            max_lookahead = 4
            search_in_lines = lines[line_index:line_index+max_lookahead]
            index_wl = self.find_start_string_in_lines(  
                search_in_lines,
                start_string = start_string,
                stop_on_first_match=True)
            indices_noordzee.append(index_wl + line_index)

        return indices_noordzee

    def read_wl_initial(self, line_indices:list):
        """

        Args:
            line_indices (list): indices to read

        Raises:
            Exception: small error check to see if all initial water levels are the same

        Returns:
            waterlevel
        """
        with open(self.ic_folder, 'r') as f:
            lines = f.readlines()
        wls = []
        for index in line_indices:
            wl_line = lines[index]
            if '#' in wl_line:
                wl_line = wl_line.split('#')[0]
            values = wl_line.split('=')[-1].strip().split(' ')
            for value in values:
                wls.append(float(value))

        if [wls[0] == wl for wl in wls]:
            return wls[0]
        else:
            raise Exception("Wls are not the same")

    def write_wl_initial(self, line_indices, new_value):
        with open(self.ic_folder) as f:
            lines = f.readlines()

        for index in line_indices:
            wl_line = lines[index]
            key = wl_line.split('=')[0]
            value = f"= {new_value} {new_value}"
            if '#' in wl_line:
                comment = wl_line.split('#')[-1]
                lines[index] = f"{key}{value}    #{comment}"
            else:
                lines[index] = f"{key}{value}\n"
        self.write_lines_to_file(self.ic_folder, lines)

    def display_widgets(self):
        """
        Display widgets and the update button
        """
        all_widgets = self.widgets_to_display + self.modify_rhine_discharge.widgets_to_display
        items = ipy.VBox(children=all_widgets, layout=self.box_layout, style = self.box_style)
        display(items)
        
        button = ipy.Button(description="Update settings", style = self.button_style, layout = self.button_layout)
        output = ipy.Output()
        display(button, output)
        button.on_click(self.update_settings_widget)
    
    def update_settings_widget(self, b):
        """
        update the settings based on the widget values, and modify the d-hydro values accordingly

        Args:
            b (_type_): only required to maket his funciton callable by button widget
        """
        print_settings = {}
        for setting in self.read_functions.keys():
            self.settings[setting] = self.widgets[setting].value
            self.write_functions[setting](self.line_indices[setting], self.settings[setting])
            print_settings[setting] = self.read_functions[setting](self.line_indices[setting])
        
        for setting in self.modify_rhine_discharge.widgets.keys():
            self.modify_rhine_discharge.settings[setting] = self.modify_rhine_discharge.widgets[setting].value
            print_settings[setting] = self.modify_rhine_discharge.settings[setting]
            
        self.modify_rhine_discharge.update()
        
        clear_output(wait=True)
        self.display_widgets()
        display("Boundary conditions are:")
        print_settings_dict = print_settings
        print(json.dumps(print_settings_dict, indent=4))



class ModifyRhineDischarge(ModifyBoundaries, WidgetStyling):
    """
    Class for modifying the discharge of the Rhine
    """
    def __init__(self, model_path, mdu_settings, dambreak_settings):
        self.bnd_folder = os.path.join(model_path, r'dflowfm\bnd.ext')
        self.rhine_bc_loc = os.path.join(model_path, r'dflowfm\rhine.bc')
        self.branch_name = 'rijn_wl_DuitseRijn'
        self.id = 'LateralSource_1D_1'

        self.model_start_time = mdu_settings['refDate'] + timedelta(seconds = mdu_settings['tStart'])
        self.tStart = int(mdu_settings['tStart'])
        self.tStop = int(mdu_settings['tStop'])
        self.tBreach = int(dambreak_settings['t0'])
        self.model_duration = self.tStop - self.tStart

        self.settings = {}
        self.settings['Rhine basic discharge'] = 5000
        self.settings['Rhine peak discharge'] = 15000
        self.settings['Rhine event start (hour)'] = 0
        self.settings['Rhine event duration (hours)'] = 24

        if os.path.exists(self.rhine_bc_loc):
            pass
        else:
            self.modify_rhine_discharge_to_timeseries()
        
        self.update()

        self.set_default_layout_and_styling()
        self.widgets = {}
        for setting in self.settings:
            self.widgets[setting] = ipy.FloatText(
                value = float(self.settings[setting]),
                description=f'{setting}:',
                disabled=False
                )
            self.widgets[setting].layout = self.item_layout
            self.widgets[setting].style = self.item_style
        self.widgets_to_display = [widget for widget in self.widgets.values()]
    
    
    def modify_rhine_discharge_to_timeseries(self):
        """
        if currently the dischrage is fixed, replace it by a file reference.
        """
        with open(self.bnd_folder) as f:
            lines = f.readlines()
        start_string = f'branchId={self.branch_name}'
        lines_found = self.find_start_string_in_lines(lines, start_string = start_string)
        
        if len(lines_found) != 1:
            raise Exception(f"{len(lines_found)} occurences of startstring {start_string}")
        else:
            start_search_index = lines_found[0]

        # search for discharge, start searching from the point where the branch name was found
        max_lookahead = 4
        search_in_lines = lines[start_search_index: start_search_index+max_lookahead]
        start_string = 'discharge='
        index_discharge = self.find_start_string_in_lines(
            search_in_lines,
            start_string = start_string,
            stop_on_first_match=True)
        
        backward_search_distance= 4
        index_id = self.find_start_string_in_lines(
            lines[start_search_index-backward_search_distance:start_search_index],
            'id=',
            stop_on_first_match=True
        )
        self.id = lines[start_search_index + index_id - backward_search_distance].split('=')[-1].strip(' ')
        discharge_line = lines[index_discharge + start_search_index]
        discharge_default = discharge_line.split('=')[-1].strip(' ')
        if discharge_default == 'rhine.bc':
            raise Exception("Restart the notebook")
        key = discharge_line.split('=')[0]
        value = f"= {os.path.basename(self.rhine_bc_loc)}"
        lines[index_discharge + start_search_index] = f"{key}{value}\n"
        self.write_lines_to_file(self.bnd_folder, lines)


    def write_discharge_rhine(self):
        """
        write the rhine discharge to a file
        """
        header = [
            "[General]",
            "    fileVersion           = 1.01"      ,          
            "    fileType              = boundConds",
        ]
        
        content_rhine = [
            "[forcing]",
            f"    name                  = {self.id}",            
            "    function              = timeseries",        
            "    timeInterpolation     = linear",           
            "    quantity              = time",      
            f"    unit                  = minutes since {self.model_start_time.strftime('%Y-%m-%d %H:%M:%S')}",
            "    quantity              = lateral_discharge",   
            "    unit                  = m³/s"
        ]

        timeseries = self.timeseries_lines
        all_lines = header + content_rhine + timeseries

        with open(self.rhine_bc_loc, 'w') as f:
            for line in all_lines:
                f.write(line)
                f.write('\n')

    def generate_discharge_timeseries(self):
        """
        Generate discharge timeseries based on the user input paramters
        """
        if self.tStop % 300 != 0 or self.tStart % 300 != 0:
            raise Exception("tstart or tstop is not valid")
        self.T = []
        self.Q = []
        duration = self.settings['Rhine event duration (hours)'] * 60 * 60
        offset = self.settings['Rhine event start (hour)'] * 60 * 60
        for t in range(self.tStart, self.tStop + 300, 300):
            self.T.append(t)
            if t < offset + duration and t > offset:
                Q_calc =  self.settings['Rhine basic discharge'] + (self.settings['Rhine peak discharge']-self.settings['Rhine basic discharge']) * (1 + np.cos(np.pi + 2*np.pi*((t-offset)/duration))) * 1/2
            else:
                Q_calc = self.settings['Rhine basic discharge']
            self.Q.append(Q_calc)
        
        self.timeseries_lines = []
        for i in range(len(self.T)):
            self.timeseries_lines.append(f"    {self.T[i]} {self.Q[i]}")
        
    def plot_timeseries(self):
        """
        Generate a plot of the discharge timeseries, together with tstart, tstop and tbreach.
        """
        fig, ax = plt.subplots(figsize = (13, 5))
        T_plot = [x/3600 for x in self.T]
        ax.plot(T_plot, self.Q, color = '#3587A4', lw = 3, label = 'Discharge Rhine')
        range = max(T_plot) - min(T_plot) 
        ax.set_xlim([min(T_plot) - range * 0.01, max(T_plot) + range * 0.01])
        ax.set_ylim(ax.get_ylim())
        ax.set_title("Discharge Rhine")
        ax.set_xlabel("T (hours)")
        ax.set_ylabel("Discharge (m3/s)")


        ax.plot([self.tStart/3600, self.tStart/3600], [-10000, 100000], color = 'orange', label = "Start of simulation")
        ax.plot([self.tStop/3600, self.tStop/3600], [-10000, 100000], color = 'orange', label = "End of simulation")
        ax.plot([self.tBreach/3600, self.tBreach/3600], [-10000, 100000], color = 'red', label = "Timestep dikebreach")
        ax.legend()
        ax.set_ylim(bottom = 0)


    def update(self):
        """
        function to update the timeseries, write it to file, and plot it
        """
        self.generate_discharge_timeseries()
        self.write_discharge_rhine()
        self.plot_timeseries()