import ipywidgets as ipy

class WidgetStyling():
    """
    Class that contains the default widget styling. This class is inherited such that the default layout and styling can be applied. 
    """
    def set_default_layout_and_styling(self):
        self.box_layout = ipy.Layout(
            border='2px solid #3587A4',
            width='100%',
            height='',
        )

        self.box_style = {}

        self.item_layout = ipy.Layout(
            width = '95%',
            margin = '5px'
            )

        self.item_style = {
            'description_width': '40%',
            }
    
        self.item_style_wide_description = {
            'description_width': '65%',
            }

        self.button_style = {
            'button_color' :'#3587A4',
            'text_color' : 'white'
        }

        self.button_layout =ipy.Layout(
            width = '99%',
            height = '35px'
        )

        self.progress_style = {
            'description_width':'initial'
        }

        self.progress_layout = ipy.Layout(
            width = '99%'
        )

    

