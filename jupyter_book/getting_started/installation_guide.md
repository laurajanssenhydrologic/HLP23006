# Installation Guide
The easiest way to install D-HydroLogic is through downloading it directly from GitHub.

:::{important}
In order to install the D-HydroLogic package, a Conda environment is required.
The most recent environment specification can be found [here](https://github.com/HydroLogicBV/P1414/blob/main/environment.yml).
:::

## Installing from GitHub
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