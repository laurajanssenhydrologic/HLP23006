from notebooks.background_scripting.v1.widget_styling import WidgetStyling

class TestConnection(WidgetStyling):
    """
    Class that does some small checks to see if things are installed correctly. 
    Tries to import some packages, and a class. 
    """
    def __init__(self):
        from flowmeshreader import load_meta_data, load_map_data, mesh_to_tiff, load_mesh
        from plotting import raster_plot_with_context
        self.set_default_layout_and_styling()
        print("Imports are succesfull!")

