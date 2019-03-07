import glob
import os
import pathlib
import re

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np


# include the path here to a folder full of rasters that QGIS can load
path_to_rasters = "."

# folder for output sections as PNGs
path_for_outputs = "."


def parse_raster_layer_depths(rfilename):
    '''Parse raster layer depths from filename.

    Args:
        rfilename (str): the filename of the raster.

    Returns: a numpy array [top, bottom] of the depths that the
    raster layer represents.

    e.g. the file named Con006_doi_gm_013.2-017m.ers
    represents the conductivity at depths between 13.2 and 17 m,
    so this function returns np.array([13.2, 17.0]).

    '''
    m = re.search(r"gm_(\d\d\d\.\d)-(\d\d\d\.\d)m", rfilename)
    assert m
    return np.array([float(m.group(1)), float(m.group(2))])


def get_relevant_depth_slice(depth, raster_layers):
    for rlayer in raster_layers:
        if rlayer["top"] <= depth and rlayer["bottom"] >= depth:
            return rlayer


def imrot(arr):
    return np.flipud(np.rot90(arr))


def run(
    path,
    output_path,
    dl,
    dz,
    bottom=None,
    name_field="Id",
    figsize=None,
    cmaps=("magma"),
):
    '''Run the script.

    Args:
        path (str): path to the raster files
        output_path (str): where to store cross section pngs
        dl (float): section x-axis spacing (arbitrary, smaller is slower)
        dz (float): section z-axis spacing (arbitrary, smaller is slower)
        name_field (str): name of an attribute from the vector polyline layer
            which define the section line. defaults to the polyline "Id".
        figsize (tuple): for matplotlib.pyplot.figure()
        cmaps (list): for matplotlib.imshow - you can generate plots for
            multiple colormaps, or just the one.

    This function does:

    1. reads in all rasters from path
    2. parses a top and bottom z value for each raster
    3. iterates over all selected polylines in the QGIS map canvas - these
    are the section lines
    4. extracts values from the rasters along each section line
    5. plots a section for each section line.

    '''
    rlayers = []
    for rfilename in glob.glob(os.path.join(path, "*.ers")):
        name = os.path.split(rfilename)[1]
        try:
            depths = parse_raster_layer_depths(rfilename)
        except AssertionError:
            pass
        else:
            print("loading {}...".format(name), end="")
            rlayer = {
                "layer": QgsRasterLayer(rfilename),
                "filename": rfilename,
                "name": name,
                "top": depths[0],
                "bottom": depths[1],
            }
            print(" {} to {} m".format(rlayer["top"], rlayer["bottom"]))
            rlayers.append(rlayer)

    get_rlayer = lambda z: get_relevant_depth_slice(z, rlayers)

    tops = np.array([r["top"] for r in rlayers])
    bottoms = np.array([r["bottom"] for r in rlayers])
    if bottom is None:
        bottom = max([r["bottom"] for r in rlayers])

    slines = []
    for sline_feature in iface.activeLayer().selectedFeatures():
        geom = sline_feature.geometry()
        if geom.type() == QgsWkbTypes.LineGeometry:
            lines = geom.asMultiPolyline()
            for p0, p1 in lines:
                line = p1 - p0
                if np.abs(line.x()) > np.abs(line.y()):
                    sline = {"sect-x-func": lambda pt: pt.x(), "sect-x-name": "Easting"}
                else:
                    sline = {
                        "sect-x-func": lambda pt: pt.y(),
                        "sect-x-name": "Northing",
                    }
                sline.update(
                    {
                        "attrs": dict(
                            zip(
                                sline_feature.fields().names(),
                                sline_feature.attributes(),
                            )
                        ),
                        "length": line.length(),
                        "feature": sline_feature,
                        "vector": line,
                        "origin": p0,
                        "end": p1,
                    }
                )

                print("found section_line: {}".format(sline))
                slines.append(sline)

    profiles = []
    for sline in slines:
        nl = int(np.ceil(sline["length"] / dl))
        nz = int(np.ceil(bottom / dz))
        arr = np.zeros((nl, nz))
        profile = {"sect": sline, "name": sline["attrs"][name_field]}
        print("extracting section {}={}...".format(name_field, profile["name"]), end="")
        point_features = []
        for il in range(nl):
            l = il * dl
            v = sline["vector"].normalized() * l
            p = sline["origin"] + v
            point = QgsFeature()
            point.setGeometry(QgsGeometry.fromPointXY(p))
            point_features.append(point)
            for iz in range(nz):
                z = iz * dz
                rlayer = get_rlayer(z)
                if rlayer:
                    val, res = rlayer["layer"].dataProvider().sample(p, 1)
                    val = float(val)
                else:
                    val = np.nan
                arr[il, iz] = val
        profile.update({"data": arr})
        profiles.append(profile)
        print(" done")

        # Plot the arr points on the active layer (temporary/scratch QGIS layer)
        # Totally optional

        # (res, outFeats) = iface.activeLayer().dataProvider().addFeatures(point_features)

    # Set colour scale from data

    vmin = min([np.nanmin(pr["data"]) for pr in profiles])
    vmax = max([np.nanmax(pr["data"]) for pr in profiles])
    
    # or manually
    vmin = 10 ** 0
    vmax = 10 ** 3

    for cmap in cmaps:
        for pr in profiles:
            fig = plt.figure(figsize=(7, 4))
            ax = fig.add_subplot(111)
            dfunc = pr["sect"]["sect-x-func"]
            left = dfunc(pr["sect"]["origin"])
            right = dfunc(pr["sect"]["end"])
            extent = [left, right, bottom, 0]
            im = ax.imshow(
                imrot(pr["data"]),
                aspect="auto",
                extent=extent,
                norm=colors.LogNorm(),
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                interpolation="nearest",
            )
            cbar = plt.colorbar(im)
            cbar.set_label("Conductivity (mS/m)")
            x0, x1 = ax.get_xlim()
            ax.set_xlim(min((x0, x1)), max((x0, x1)))
            ax.set_title("{}: {}".format(name_field, pr["name"]))
            ax.set_xlabel(pr["sect"]["sect-x-name"])
            ax.set_ylabel("Depth (m)")
            fig.tight_layout()
            fig.savefig(
                str(
                    pathlib.Path(output_path)
                    / "vertsec_{}_{}_{}.png".format(cmap, name_field, pr["name"])
                ),
                dpi=120,
            )


run(
    path_to_rasters,
    path_for_outputs,
    dl=50,
    dz=5,
    bottom=150,
    name_field="x-sect-id",
    cmaps=("magma", ),
)

# optionally...

plt.show()

# QGIS keeps figures open. Use this to clear it up:

# plt.close("all")
