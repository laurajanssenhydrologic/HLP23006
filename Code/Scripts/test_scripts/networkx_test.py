import math
import warnings

import geopandas as gpd
import libpysal
import matplotlib.pyplot as plt
import momepy
import networkx as nx
import numpy as np

momepy.gdf_to_nx


def _angle(a, b, c):
    """
    Measure angle between a-b, b-c. In degrees.
    Helper for gdf_to_nx.
    Adapted from cityseer's implementation.
    """
    a1 = math.degrees(math.atan2(b[1] - a[1], b[0] - a[0]))
    a2 = math.degrees(math.atan2(c[1] - b[1], c[0] - b[0]))
    return abs((a2 - a1 + 180) % 360 - 180)


def _generate_primal(G, gdf_network, fields, multigraph):
    """
    Generate primal graph.
    Helper for gdf_to_nx.
    """
    G.graph["approach"] = "primal"

    msg = "%s. This can lead to unexpected behaviour. The intended usage of the conversion function is with networks made of LineStrings only."

    if not "LineString" in gdf_network.geom_type.unique():
        warnings.warn(
            message=msg % "The given network does not contain any LineString.",
            category=RuntimeWarning,
        )

    if len(gdf_network.geom_type.unique()) > 1:
        warnings.warn(
            message=msg % "The given network consists of multiple geometry types.",
            category=RuntimeWarning,
        )

    key = 0
    for row in gdf_network.itertuples():
        first = row.geometry.coords[0]
        last = row.geometry.coords[-1]

        # Round location
        # first = (np.around(first[0], -1), np.around(first[1], -1))
        # last = (np.around(last[0], -1), np.around(last[1], -1))
        first = (np.around(first[0], 5), np.around(first[1], 5))
        last = (np.around(last[0], 5), np.around(last[1], 5))

        data = [r for r in row][1:]
        attributes = dict(zip(fields, data))
        if multigraph:
            G.add_edge(first, last, key=key, **attributes)
            key += 1
        else:
            G.add_edge(first, last, **attributes)


def _generate_dual(G, gdf_network, fields, angles, multigraph, angle):
    """
    Generate dual graph
    Helper for gdf_to_nx.
    """
    G.graph["approach"] = "dual"
    key = 0

    sw = libpysal.weights.Queen.from_dataframe(gdf_network, silence_warnings=True)
    cent = gdf_network.geometry.centroid
    gdf_network["temp_x_coords"] = cent.x
    gdf_network["temp_y_coords"] = cent.y

    for i, row in enumerate(gdf_network.itertuples()):
        centroid = (row.temp_x_coords, row.temp_y_coords)
        data = [f for f in row][1:-2]
        attributes = dict(zip(fields, data))
        G.add_node(centroid, **attributes)

        if sw.cardinalities[i] > 0:
            for n in sw.neighbors[i]:
                start = centroid
                end = (
                    gdf_network["temp_x_coords"].iloc[n],
                    gdf_network["temp_y_coords"].iloc[n],
                )
                p0 = row.geometry.coords[0]
                p1 = row.geometry.coords[-1]
                geom = gdf_network.geometry.iloc[n]
                p2 = geom.coords[0]
                p3 = geom.coords[-1]
                points = [p0, p1, p2, p3]
                shared = [x for x in points if points.count(x) > 1]
                if shared:  # fix for non-planar graph
                    remaining = [e for e in points if e not in [shared[0]]]
                    if len(remaining) == 2:
                        if angles:
                            angle_value = _angle(remaining[0], shared[0], remaining[1])
                            if multigraph:
                                G.add_edge(start, end, key=0, **{angle: angle_value})
                                key += 1
                            else:
                                G.add_edge(start, end, **{angle: angle_value})
                        else:
                            if multigraph:
                                G.add_edge(start, end, key=0)
                                key += 1
                            else:
                                G.add_edge(start, end)


def gdf_to_nx(
    gdf_network,
    approach="primal",
    length="mm_len",
    multigraph=True,
    directed=False,
    angles=True,
    angle="angle",
):
    """
    Convert LineString GeoDataFrame to networkx.MultiGraph or other Graph as per
    specification.

    Preserves columns as edge or node attributes (depending on the ``approach``).
    Index is not preserved.

    See the User Guide page :doc:`../../user_guide/graph/convert` for details.

    Parameters
    ----------
    gdf_network : GeoDataFrame
        GeoDataFrame containing objects to convert
    approach : str, default 'primal'
        Allowed options are ``'primal'`` or ``'dual'``. Primal graph represents
        endpoints as nodes and LineStrings as edges, dual graph represents
        LineStrings as nodes and their topological relation as edges. In such a
        case, it can encode an angle between LineStrings as an edge attribute.
    length : str, default 'mm_len'
        name of attribute of segment length (geographical) which will be saved to graph
    multigraph : bool, default True
        create ``MultiGraph`` of ``Graph`` (potentially directed). ``MutliGraph``
        allows multiple
        edges between any pair of nodes, which is a common case in street networks.
    directed : bool, default False
        create directed graph (``DiGraph`` or ``MultiDiGraph``). Directionality follows
        the order of LineString coordinates.
    angles : bool, default True
        capture angles between LineStrings as an attribute of a dual graph. Ignored if
        ``approach="primal"``.
    angle : str, default 'angle'
        name of attribute of angle between LineStrings which will be saved to graph.
        Ignored if ``approach="primal"``.

    Returns
    -------
    networkx.Graph,
    networkx.MultiGraph,
    networkx.DiGraph,
    networkx.MultiDiGraph
        Graph as per specification

    See also
    --------
    nx_to_gdf

    Examples
    --------
    >>> import geopandas as gpd
    >>> df = gpd.read_file(momepy.datasets.get_path('bubenec'), layer='streets')
    >>> df.head(5)
                                                geometry
    0  LINESTRING (1603585.640 6464428.774, 1603413.2...
    1  LINESTRING (1603268.502 6464060.781, 1603296.8...
    2  LINESTRING (1603607.303 6464181.853, 1603592.8...
    3  LINESTRING (1603678.970 6464477.215, 1603675.6...
    4  LINESTRING (1603537.194 6464558.112, 1603557.6...

    Primal graph:

    >>> G = momepy.gdf_to_nx(df)
    >>> G
    <networkx.classes.multigraph.MultiGraph object at 0x7f8cf90fad50>

    >>> G_directed = momepy.gdf_to_nx(df, directed=True)
    >>> G_directed
    <networkx.classes.multidigraph.MultiDiGraph object at 0x7f8cf90f56d0>

    >>> G_digraph = momepy.gdf_to_nx(df, multigraph=False, directed=True)
    >>> G_digraph
    <networkx.classes.digraph.DiGraph object at 0x7f8cf9150c10>

    >>> G_graph = momepy.gdf_to_nx(df, multigraph=False, directed=False)
    >>> G_graph
    <networkx.classes.graph.Graph object at 0x7f8cf90facd0>

    Dual graph:

    >>> G_dual = momepy.gdf_to_nx(df, approach="dual")
    >>> G_dual
    <networkx.classes.multigraph.MultiGraph object at 0x7f8cf9150fd0>


    """
    gdf_network = gdf_network.copy()
    if "key" in gdf_network.columns:
        gdf_network.rename(columns={"key": "__key"}, inplace=True)

    if multigraph and directed:
        net = nx.MultiDiGraph()
    elif multigraph and not directed:
        net = nx.MultiGraph()
    elif not multigraph and directed:
        net = nx.DiGraph()
    else:
        net = nx.Graph()

    net.graph["crs"] = gdf_network.crs
    gdf_network[length] = gdf_network.geometry.length
    fields = list(gdf_network.columns)

    if approach == "primal":
        _generate_primal(net, gdf_network, fields, multigraph)

    elif approach == "dual":
        if directed:
            raise ValueError("Directed graphs are not supported in dual approach.")

        _generate_dual(net, gdf_network, fields, angles=angles, multigraph=multigraph, angle=angle)

    else:
        raise ValueError(f"Approach {approach} is not supported. Use 'primal' or 'dual'.")

    return net


# gdf = gpd.read_file(
#     r"D:\Work\Project\P1414\GIS\HDSR\Legger\Hydro_Objecten(2)\HydroObject.shp"
# ).explode()
# gdf = gdf.loc[gdf["IWS_W_WATB"] > 5, :]
# gdf = gpd.read_file(
#     r"D:\Work\Project\P1414\GIS\HHDelfland\Legger_Delfland_shp\Oppervlaktewaterlichamen\Primair water.shp"
# ).explode()
# gdf = gpd.read_file(
#     r"D:\Work\Project\P1414\GIS\HHRijnland\Legger\Watergang\Watergang_as.shp"
# ).explode(ignore_index=True)
# gdf = gdf.loc[gdf["BODEMBREED"] > 5, :]
# gdf = gpd.read_file(r"D:\Work\Project\P1414\GIS\HHSK\Legger\Hoofdwatergang.shp").explode()
# gdf = gdf.loc[gdf["WATERBREED"] > 5, :]
gdf = gpd.read_file(r"D:\Work\Project\P1414\GIS\WAGV\hydrovak\hydrovak.shp").explode()
gdf = gdf.loc[gdf["IWS_W_WATB"] > 5, :]

print("rounding coords")
G = gdf_to_nx(gdf)
# G = nx.Graph(G)
component_list = sorted(nx.connected_components(G), key=len, reverse=True)

for n in range(3):
    sG = G.subgraph(component_list[n])

    new_lengths = old_lengths = 0
    for x, y, data in sG.edges(data=True):
        # new_lengths += edge
        print(data)
        break

# print(G.edges[:, :, :]["OBJECTID"])
# mering three most impotant components to network

# print(H.edges)
# print(
#     H.get_edge_data(
#         u=(5.294029083000055, 51.97453857700003), v=(5.302407304000042, 51.97215656900005)
#     )
# )
# print(nx.degree_histogram(G))
print(nx.degree_histogram(H))

# trim non-conected nodes
# to_be_removed = [x for x in G.nodes() if G.degree(x) <= 1]
# for x in to_be_removed:
#     G.remove_node(x)

# to_be_removed = [x for x in G.nodes() if G.degree(x) < 1]
# for x in to_be_removed:
#     G.remove_node(x)


# print(H.nodes)
# Plot
positions = {n: [n[0], n[1]] for n in list(H.nodes)}
f, ax = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
gdf.plot(color="k", ax=ax[0])
for i, facet in enumerate(ax):
    facet.set_title(("Rivers", "Graph")[i])
    facet.axis("off")
nx.draw(H, positions, ax=ax[1], node_size=5)
plt.show()
