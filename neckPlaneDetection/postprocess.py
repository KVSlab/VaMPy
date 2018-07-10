#!/usr/bin/env python

from remove_aneurysm import *
from common import *
from os import path, listdir
from extractaneurysmneckplanesection import extractaneurysmneckplanesection_run
from computeaneurysmsaccenterline import computeaneurysmsaccenterline_run
from computesurfacesandvolumes import computesurfacesandvolumes_run
from extractparentvesselparameters import extractparentvesselparameters_run


def get_data(dirpath):
    # Naming convention
    domain_name = "Case" + dirpath.split(path.sep)[-1]
    domain_path = dirpath

    # Check if case is already analyzed
    if path.exists(path.join(domain_path, "csvfile", domain_name + "_angles.csv")):
        return

    # Models
    model_path = path.join(domain_path, domain_name + ".vtp")
    centerline_recon_path = path.join(domain_path, domain_name + "_reconstructedmodel_cl.vtp")
    voronoi_path = path.join(domain_path, domain_name + "_reconstructedmodel_vor.vtp")

    # Parameters
    anu_num = 0  # TODO: Read from info.txt
    smooth = True
    smooth_factor = 0.25
    bif = False #True
    lower = True #False
    addPoint = True
    cylinder_factor = 7.0
    aneurysm_type = "terminal"

    s = ""
    s += "" if not bif else "_bif"
    s += "" if not smooth else "_smooth"
    s += "" if not addPoint else "_extraPoint"
    s += "" if not lower else "_lower"
    s += "" if cylinder_factor == 7.0 else "_cyl%s" % cylinder_factor
    model_new_surface = path.join(dirpath, domain_name + "_anu"+s+".vtp")

    if not path.exists(model_new_surface):
        new_surface = remove(dirpath, smooth, smooth_factor, bif, addPoint, lower,
                             cylinder_factor, aneurysm_type, anu_num,
                             domain_name)
    else:
        new_surface = ReadPolyData(model_new_surface)

    if new_surface.GetNumberOfPoints() == 0:
        print("ERROR: No points in the reconstructed surface")
        return ##raise ValueError

    # Compute new centerline
    inlet, outlets = get_centers(None, dirpath)
    compute_centerlines(inlet, outlets, centerline_recon_path, new_surface,
                        resampling=0.1)

    makeVoronoiDiagram(new_surface, voronoi_path)

    # Double check that the new voronoi diagram does not contain extreme points
    # this is unlikely, but have been observed
    vor = ReadPolyData(voronoi_path)
    vor = remove_extreme_points(vor, ReadPolyData(model_path.replace(".vtp",
                                                                     "_voronoi_anu.vtp")))
    WritePolyData(vor, voronoi_path)

    if not path.exists(path.join(domain_path, "csvfile", domain_name + "_parameters.csv")):
        print("extractaneurysmneckplanesection_run")
        extractaneurysmneckplanesection_run(domain_path, domain_name, "terminal")
    if not path.exists(path.join(domain_path, "csvfile", domain_name + "_centerline.csv")):
        print("computeaneurysmsaccenterline_run")
        computeaneurysmsaccenterline_run(domain_path, domain_name)
    if not path.exists(path.join(domain_path, "csvfile", domain_name + "_surfacesandvolumes.csv")):
        print("computesurfacesandvolumes_run")
        computesurfacesandvolumes_run(domain_path, domain_name)
    if not path.exists(path.join(domain_path, "csvfile", domain_name + "_angles.csv")):
        print("extractparentvesselparameters_run")
        extractparentvesselparameters_run(domain_path, domain_name, "terminal")


if __name__ == "__main__":
    # input, path to folder
    # TODO: Add argumentparser, and make the other parameters availeble as well
    get_data("Case73")
