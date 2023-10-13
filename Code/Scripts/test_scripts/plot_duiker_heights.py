import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interpn
from scipy.stats import linregress


def density_scatter(x, y, sort=True, bins=20, **kwargs):
    """
    Scatter plot colored by 2d histogram
    """

    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn(
        (0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])),
        data,
        np.vstack([x, y]).T,
        method="splinef2d",
        bounds_error=False,
    )

    # To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    return x, y, z


duikers_path = r"D:\Work\Project\HL-23006\GIS\Legger\Kokers_Lijnen_filled.shp"
duikers_gdf = gpd.read_file(duikers_path).to_crs("EPSG:28992")

plt.figure()
x = duikers_gdf["AHN1"].values
y = duikers_gdf["HOOGTEBOK1"].values.astype(float)
select_bool = (np.isnan(x)) | (np.isnan(y))

x = x[~select_bool]
y = y[~select_bool]
x, y, z = density_scatter(x=x, y=y, sort=True, bins=50)
alpha = (z - np.min(z)) / (np.max(z) - np.min(z))
plt.scatter(x=x, y=y, c=z, alpha=alpha, label="Duiker")

# p = np.polyfit(x=x, y=y, deg=1)


slope, intercept, r_value, p_value, std_err = linregress(x, y)
x_p = [np.min(x), np.max(x)]
y_p = [np.min(x) * slope + intercept, np.max(x) * slope + intercept]
print(slope, intercept, r_value, p_value)

plt.plot(
    x_p, y_p, label="linear fit: {:.2f}x + {:.2f}\n$R^2$: {:.2f}".format(slope, intercept, r_value)
)

plt.colorbar(label="Density (-)")

plt.legend()
plt.xlabel("AHN (m+NAP)")
plt.ylabel("HOOGTEBOK (m+NAP)")

plt.figure()
x = duikers_gdf["AHN2"].values
y = duikers_gdf["HOOGTEBOK2"].values.astype(float)
select_bool = (np.isnan(x)) | (np.isnan(y))

x = x[~select_bool]
y = y[~select_bool]
x, y, z = density_scatter(x=x, y=y, sort=True, bins=50)
alpha = (z - np.min(z)) / (np.max(z) - np.min(z))
plt.scatter(x=x, y=y, c=z, alpha=alpha, label="Duiker")

slope, intercept, r_value, p_value, std_err = linregress(x, y)
x_p = [np.min(x), np.max(x)]
y_p = [np.min(x) * slope + intercept, np.max(x) * slope + intercept]
print(slope, intercept, r_value, p_value)


plt.plot(
    x_p, y_p, label="linear fit: {:.2f}x + {:.2f}\n$R^2$: {:.2f}".format(slope, intercept, r_value)
)

plt.colorbar(label="Density (-)")

plt.legend()
plt.xlabel("AHN (m+NAP)")
plt.ylabel("HOOGTEBOK (m+NAP)")
plt.show()
