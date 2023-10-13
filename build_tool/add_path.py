
with open("../jupyter_book/conf.py", "a") as myfile:
    myfile.write("import os\n")
    myfile.write("import sys\n")
    myfile.write("sys.path.insert(0, os.path.abspath('..\HydroLogic_Inundation_toolbox\Readers'))\n")
    myfile.write("sys.path.insert(0, os.path.abspath('..\HydroLogic_Inundation_toolbox'))\n")
    myfile.write("sys.path.insert(0, os.path.abspath('..\Code'))\n")
    myfile.write("sys.path.insert(0, os.path.abspath('..\Code\data_structures'))\n")
    myfile.write("napoleon_custom_sections = [('Returns', 'params_style')]\n")