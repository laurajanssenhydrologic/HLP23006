# ROI Toolbox
| ![Schematic overview](./Software_diagram_v2.png) |
|:--:|
| Figure 1:  schematic overview of the workflows in [DhydroData](dhydro_data_api.rst). |

The toolchain in the ROI Toolbox to build D-HYDRO models is shown in the figure below.
The following API-references are provided:

::::{grid}
:gutter: 4

:::{grid-item-card} 
:columns: 4
:link: hydamo_helpers_api
:link-type: doc
:text-align: center

**HYDAMO Helper Functions**
^^^
Helper functions to convert raw-data to a [ROI Data Model](roi_data_api.rst) compliant format.
:::

:::{grid-item-card} 
:columns: 4
:link: dhydro_data_api
:link-type: doc
:text-align: center

**DhydroData**
^^^
Main class to read and write HYDAMO geopackages, convert raw data into HYDAMO format, and build D-HYDRO models.
:::

:::{grid-item-card} 
:columns: 4
:link: to_dhydro_helpers_api
:link-type: doc
:text-align: center

**To Dhydro Helper Functions**
^^^
Helper functions that can convert a [ROI Data Model](roi_data_api.rst) into a D-HYDRO model.
:::

:::{grid-item-card} 
:columns: 4
:link: hydamo_data_api
:link-type: doc
:text-align: center

**HYDAMO Data Model**
^^^
HYDAMO data model.
:::

:::{grid-item-card} 
:columns: 4
:link: roi_data_api
:link-type: doc
:text-align: center

**ROI Data Model**
^^^
Extension of the [HYDAMO Data Model](hydamo_data_api.rst).
:::

:::{grid-item-card} 
:columns: 4
:link: fixed_weirs_helpers_api
:link-type: doc
:text-align: center

**Fixed Weirs Helper Functions**
^^^
Helper functions to convert fixed weirs to a [ROI Data Model](roi_data_api.rst) compliant format.
:::

::::{grid}