import ipywidgets as ipy
from IPython.display import display, clear_output
import json
from notebooks.background_scripting.v1.widget_styling import WidgetStyling
import os

class ModelSettings(WidgetStyling):
    """
    Class that contains all the settings that need to be used to create the model.
    Also contians the functions that communicate with the widgets from the jupyter notebook. 
    """
    def __init__(self, git_dir):
        # initiate default values for settings
        self.settings = {}
        self.settings['folder'] = os.path.join(git_dir, 'Model_database')
        try: 
            self.settings['model_options'] = self.select_models(self.settings['folder'])
            self.settings['model'] = self.settings['model_options'][0]
        except:
            self.settings['model_options'] = ['The folder you selected is not valid!', 'Please update the folder above and hit "Update settings".']
            self.settings['model'] = self.settings['model_options'][0]

        self.settings['scenario_name'] = 'run'
        self.settings_to_print = ['folder', 'model', 'scenario_name']
      
        self.folder_widget = ipy.Text(
            value= self.settings['folder'],
            placeholder= self.settings['folder'],
            description='Path to folder with models:',
            disabled=False
        )

        self.select_model_widget = ipy.Dropdown(
            options=self.settings['model_options'],
            value=self.settings['model'],
            description='Model:',
            disabled=False
        )

        self.scenario_name_widget = ipy.Text(
            value=self.settings['scenario_name'],
            placeholder=self.settings['scenario_name'],
            description='Scenario name:',
            disabled=False   
        )

        self.set_default_layout_and_styling()
        

        self.widgets_to_display = [self.folder_widget, self.select_model_widget, self.scenario_name_widget]
        for widget in self.widgets_to_display:
            widget.layout = self.item_layout
            widget.style = self.item_style

    def select_models(self, folder):
        """
        Function that checks for each folder if the folder contains a model (contains a .xml file) or not
        """
        valid_folders = []
        for potential_folder in os.listdir(folder):
            for f in os.listdir(os.path.join(folder, potential_folder)):
                if f.endswith('.xml'):
                    valid_folders.append(potential_folder)
                    break
        return valid_folders



    def update_settings_dict(self, b):
        """
        Function to update settings based on the widget values after the update button is pressed.
        """
        if self.settings['folder'] != self.folder_widget.value.strip(' '):
            self.settings['folder'] = self.folder_widget.value

            try:
                self.settings['model_options'] = self.select_models(self.settings['folder'])
                self.select_model_widget.options = self.settings['model_options']
                self.settings['model'] = self.settings['model_options'][0]
                self.select_model_widget.value = self.settings['model']
            except:
                self.settings['model_options'] =  ['The folder you selected is not valid!', 'Please update the folder above and hit "Update settings".']
                self.select_model_widget.options = self.settings['model_options']
                self.settings['model'] = self.settings['model_options'][0]
                self.select_model_widget.value = self.settings['model']
        else:
            self.settings['model'] = self.select_model_widget.value
            self.settings['scenario_name'] = self.scenario_name_widget.value

        clear_output(wait=True)
        self.display_widgets()
        display("Model settings are:")
        print_settings_dict = {k: self.settings[k] for k in self.settings_to_print}
        print(json.dumps(print_settings_dict, indent=4))
        
    def display_widgets(self):
        """
        Function to display all widgets, and the update button.
        """
        items = ipy.VBox(children=self.widgets_to_display, layout=self.box_layout, style = self.box_style)
        display(items)
        button = ipy.Button(
            description="Update settings", 
            style = self.button_style, 
            layout = self.button_layout,
            )
        output = ipy.Output()
        display(button, output)
        button.on_click(self.update_settings_dict)
