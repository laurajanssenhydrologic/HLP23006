# Installation Guide
The easiest way to install D-HydroLogic is by setting up a Conda environment and by downloading D-HydroLogic directly from GitHub.

## 1. Creating a Python environment
In order to succesfully install the D-HydroLogic package, a Python environment is required. 
An easy way to create such an environment is by creating a Conda environment using the .yml file that can be found on the [GitHub repository](https://github.com/HydroLogicBV/P1414/blob/main/environment_sas.yml).
This can be done using the code snippets below.
````
conda env create -f environment_sas.yml
````

````
conda activate inundation_toolbox_sas
````

:::{note}
The Conda manual provides more information on how to create a new Python environment from a .yml file. 
This manual can be found [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).
:::

## 2. Installing from GitHub
In order to install D-HydroLogic from GitHub, first clone the repository, then install it using pip from the downloaded directory; per below.

````
git clone https://github.com/HydroLogicBV/P1414.git
````

````
pip install ./hydrolib_adapted
````

## Troubleshooting

Incase you get an error that includes the following: "Error: Could not install packages due to an OSError: \[WinError 5\] Access is denied ...".
Try the following command to finish installing all required packages:
````
pip install ugrid imod ipympl isort tqdm matplotlib_scalebar
````