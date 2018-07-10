from common import *
from argparse import ArgumentParser
from os import path, listdir
from time import time
from copy import deepcopy

import sys
import math

# Local import
from patchandinterpolatecenterlines import *
#from vmtkClassifyBifucation import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
import ToolRepairSTL


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()
    # TODO: Change to --feature --no-feature for bool agruments
    # Ex. 
    # feature_parser = parser.add_mutually_exclusive_group(required=False)
    # feature_parser.add_argument('--feature', dest='feature', action='store_true')
    # feature_parser.add_argument('--no-feature', dest='feature', action='store_false')
    # parser.set_defaults(feature=True)

    # Also add choise when that is approporiate

    # Could also make bool be int and add choise 0, 1 to be true false...

    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('--case', type=str, default=None, help="Choose case")
    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                        "smoothed before it is manipulated", metavar="smooth")
    parser.add_argument('--smooth_factor', type=float, default=0.25,
                         help="If smooth option is true then each voronoi point" + \
                         " that has a radius less then MISR*(1-smooth_factor) at" + \
                         " the closest centerline point is removes",
                         metavar="smoothening_factor")
    parser.add_argument("--bif", type=bool, default=False,
                        help="interpolate bif as well")
    parser.add_argument("--addPoint", type=bool, default=False,
                        help="Add additional point for integration")
    parser.add_argument("--lower", type=bool, default=False,
                        help="Make a fourth line to interpolate along that" + \
                             " is lower than the other bif line.")
    parser.add_argument("--cylinder_factor", type=float, default=7.0,
                        help="Factor for choosing the smaller cylinder")
    parser.add_argument("--type", type=str, default="terminal",
                        help="Which aneurysm type, either 'lateral' or" + \
                        " 'terminal'")
    parser.add_argument("--anu_num", type=int, default=0,
                        help="If multiple aneurysms, choise one")

    args = parser.parse_args()

    if args.type not in ["lateral", "terminal"]:
        RuntimeError(("Aneurysm type: %s is not a valid choise, chose either" + \
                     "'lateral' or 'terminal'.") % args.type)

    return args.s, args.smooth_factor, args.bif, args.addPoint, args.d,\
            args.case, args.lower, args.cylinder_factor, args.type, args.anu_num


def get_points(data, key, bif=False):
    if bif:
        div_points = np.asarray([data["dau1"][key], data["dau2"][key]])
    else:
        div_points = np.asarray([data["par"][key], data["dau1"][key], data["dau2"][key]])

    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)

    return points, div_points


def get_startpoint(centerline):
    line = ExtractSingleLine(centerline, 0)
    return line.GetPoints().GetPoint(0)


def merge_cl(centerline, end_point, div_point):
    merge = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    N_lines = centerline.GetNumberOfLines()

    arrays = []
    N_, names = get_number_of_arrays(centerline)
    for i in range(N_):
        tmp = centerline.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        array = get_vtk_array(names[i], tmp_comp, centerline.GetNumberOfPoints())
        arrays.append(array)

    # Find lines to merge
    lines = [ExtractSingleLine(centerline, i) for i in range(N_lines)]
    locators = [get_locator(lines[i]) for i in range(N_lines)]
    div_ID = [locators[i].FindClosestPoint(div_point[0]) for i in range(N_lines)]
    end_ID = [locators[i].FindClosestPoint(end_point[0]) for i in range(N_lines)]
    dist = [np.sum(lines[i].GetPoint(end_ID[i]) - end_point[0]) for i in range(N_lines)]
    #change = [j for j in range(N_lines) if dist[j] != 0 or dist[j] > 5]

    # Find the direction of each line
    map_other = {0: 1, 1: 0}
    ID0 = locators[0].FindClosestPoint(end_point[1])
    ID1 = locators[1].FindClosestPoint(end_point[1])
    dist0 = math.sqrt(np.sum((np.asarray(lines[0].GetPoint(ID0)) - end_point[1])**2))
    dist1 = math.sqrt(np.sum((np.asarray(lines[1].GetPoint(ID1)) - end_point[1])**2))
    end1 = 0 if dist0 < dist1 else 1
    end2 = int(not end1)
    for i in range(2, N_lines):
        ID1 = locators[i].FindClosestPoint(end_point[1])
        ID2 = locators[i].FindClosestPoint(end_point[2])
        dist1 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID1)) - end_point[1])**2))
        dist2 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID2)) - end_point[2])**2))
        map_other[i] = end1 if dist1 > dist2 else end2

    counter = 0
    for i in range(centerline.GetNumberOfLines()):
        line = lines[i]

        # Check if it should be merged
        loc = get_locator(line)
        clipp_id = loc.FindClosestPoint(end_point[0])
        div_id = loc.FindClosestPoint(div_point[0])
        clipp_dist = distance(line.GetPoint(clipp_id), end_point[0])
        div_dist = distance(line.GetPoint(div_id), div_point[0])
        tol = getTolerance(line)*3
        merge_bool = True
        if clipp_dist > tol or div_dist > tol:
            merge_bool = False

        # Get the other line
        other = lines[map_other[i]]
        N = line.GetNumberOfPoints()
        cellArray.InsertNextCell(N)

        for j in range(N):
            # Add point
            if div_ID[i] < j < end_ID[i] and merge_bool: # and i in change:
                new = (np.asarray(other.GetPoint(j)) +
                    np.asarray(line.GetPoint(j))) / 2.
                points.InsertNextPoint(new)
            else:
                points.InsertNextPoint(line.GetPoint(j))

            cellArray.InsertCellPoint(counter)

            # Add array
            for k in range(N_):
                num = arrays[k].GetNumberOfComponents()
                if num == 1:
                    tmp = line.GetPointData().GetArray(names[k]).GetTuple1(j)
                    arrays[k].SetTuple1(counter, tmp)
                elif num == 3:
                    tmp = line.GetPointData().GetArray(names[k]).GetTuple3(j)
                    arrays[k].SetTuple3(counter, tmp[0], tmp[1], tmp[2])
                else:
                    print("Add more options")
                    sys.exit(0)

            counter += 1

    # Insert points, lines and arrays
    merge.SetPoints(points)
    merge.SetLines(cellArray)
    for i in range(N_):
        merge.GetPointData().AddArray(arrays[i])

    return merge


def remove(dirpath, smooth, smooth_factor, bif, addPoint, lower,
         cylinder_factor, aneurysm_type, anu_num, name):
    # Input filenames
    model_path = path.join(dirpath, name + ".vtp")

    # Output names
    # Surface
    model_smoothed_path = path.join(dirpath, name + "_smooth.vtp")

    # Centerliens
    centerline_par_path = path.join(dirpath, name + "_centerline_par.vtp")
    centerline_par_ordered_path = path.join(dirpath, name + "_centerline_orderd_par.vtp")
    centerline_aneurysm_path = path.join(dirpath, name + "_centerline_aneurysm.vtp")
    centerline_dau1_path = path.join(dirpath, name + "_centerline_dau1.vtp")
    centerline_dau2_path = path.join(dirpath, name + "_centerline_dau2.vtp")
    centerline_bif_path = path.join(dirpath, name + "_centerline_bif.vtp")
    centerline_complete_path = path.join(dirpath, name + "_centerline_complete.vtp")
    centerline_clipped_path = path.join(dirpath, name + "_centerline_clipped_anu.vtp")
    centerline_clipped_bif_path = path.join(dirpath, name + "_centerline_clipped_bif_anu.vtp")
    centerline_new_path = path.join(dirpath, name + "_centerline_interpolated_anu.vtp")
    centerline_new_bif_path = path.join(dirpath, name + "_centerline_interpolated_bif_anu.vtp")
    centerline_new_bif_lower_path = path.join(dirpath, name + "_centerline_interpolated_bif_lower_anu.vtp")
    centerline_relevant_outlets_path = path.join(dirpath, name + "_centerline_relevant_outlets.vtp")

    # Voronoi diagrams
    voronoi_path = path.join(dirpath, name +"_voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, name +"_voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, name +"_voronoi_clipped_anu.vtp")
    voronoi_anu_path = path.join(dirpath, name + "_voronoi_anu.vtp")

    # Points
    points_clipp_path = path.join(dirpath, name + "_clippingpoints.vtp")
    points_div_path = path.join(dirpath, name +"_divergingpoints.vtp")

    # Naming based on different options
    s = ""
    s += "" if not bif else "_bif"
    s += "" if not smooth else "_smooth"
    s += "" if not addPoint else "_extraPoint"
    s += "" if not lower else "_lower"
    s += "" if cylinder_factor == 7.0 else "_cyl%s" % cylinder_factor
    model_new_surface = path.join(dirpath, name + "_anu"+s+".vtp")

    # Read and check model
    if not path.exists(model_path):
        RuntimeError("The given directory: %s did not contain the file: model.vtp" % dirpath)

    # Get aneurysm type
    parameters = get_parameters(dirpath)
    if "aneurysm_type" in parameters.keys():
        aneurysm_type = parameters["aneurysm_type"]
        print("Aneurysm type read from info.txt file: %s" % aneurysm_type)

    # Clean surface
    surface = ReadPolyData(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    # Check the mesh if there is redundant nodes or NaN triangles
    if not "check_surface" in parameters.keys():
        ToolRepairSTL.surfaceOverview(surface)
        ToolRepairSTL.foundAndDeleteNaNTriangles(surface)
        surface = ToolRepairSTL.cleanTheSurface(surface)
        foundNaN = ToolRepairSTL.foundAndDeleteNaNTriangles(surface)
        if foundNaN == True:
            raise RuntimeError('There is an issue with the surface. Nan' + \
                                'coordinates or some other shenanigans.')
        else:
            parameters["check_surface"] = True
            write_parameters(parameters, dirpath)

        # Check connectivity and only choose the surface with the largest area
        connectivity = getConnectivity(surface)
        region_array = get_array("RegionId", connectivity)
        if region_array.max() > 0:
            print("WARNING: The surface is not connected, therefore only the" + \
                    " connected surface with the larges area will be keept.")
            surfaces = []
            for i in range(region_array):
                surfaces.append(threshold(connectivity, "RegionId", lower=(i-0.1),
                                        upper=(i+0.1), type="between", source=0))
                area.append(compute_area(surfaces[-1]))

            surface = surfaces[area.index(max(area))]

    # Get a capped and uncapped version of the surface
    open_surface = surface
    capped_surface = ReadPolyData(model_path.replace(".vtp", "_closed.vtp"))

    # Get aneurysm "end point"
    aneurysm = get_aneurysm_dome(capped_surface, dirpath, anu_num)
    outlet1, outlet2 = get_relevant_outlets(capped_surface, dirpath)

    # Get inlet and outlets
    inlet, outlets = get_centers(open_surface, dirpath)

    # Sort outlets
    tmp_outlets = np.array(outlets).reshape(len(outlets)//3, 3)
    outlet1_index = np.argsort(np.sum((tmp_outlets - outlet1)**2, axis=1))[0]
    outlet2_index = np.argsort(np.sum((tmp_outlets - outlet2)**2, axis=1))[0]
    tmp_outlets = tmp_outlets.tolist()
    if max(outlet1_index, outlet2_index) == outlet1_index:
        outlet1 = tmp_outlets.pop(outlet1_index)
        outlet2 = tmp_outlets.pop(outlet2_index)
    else:
        outlet2 = tmp_outlets.pop(outlet2_index)
        outlet1 = tmp_outlets.pop(outlet1_index)
    outlet_rest = (np.array(tmp_outlets).flatten()).tolist()
    outlets = outlet1 + outlet2 + outlet_rest
    data = {}
    for i in range(len(outlets)//3):
        data["outlet"+str(i)] = outlets[3*i:3*(i+1)]
    write_parameters(data, dirpath)

    # There is a wide range of different centerlines that is computed and
    # it is probably to many. However, to keep the current algorithm
    # and robustness of the method it is necessary..

    # Compute parent artery and aneurysm centerline
    centerline_aneurysm = compute_centerlines(inlet, aneurysm,
                                              centerline_aneurysm_path,
                                              capped_surface, resampling=0.1)
    centerline_par = compute_centerlines(inlet, outlets,
                                          centerline_par_path,
                                          capped_surface, resampling=0.1)
    centerlines_complete = compute_centerlines(inlet, outlets + aneurysm,
                                               centerline_complete_path,
                                               capped_surface, resampling=0.1)

    # Compute centerlines for landmaring aneurysm
    centerline_dau1 = compute_centerlines(outlet2, aneurysm + outlet1,
                                          centerline_dau1_path, capped_surface,
                                         resampling=0.1)
    centerline_dau2 = compute_centerlines(outlet1, aneurysm + outlet2,
                                          centerline_dau2_path, capped_surface,
                                         resampling=0.1)

    # Additional centerline for bifurcation
    if aneurysm_type == "terminal":
        centerline_relevant_outlets = compute_centerlines(inlet, outlet1 + outlet2,
                                                          centerline_relevant_outlets_path,
                                                          capped_surface,
                                                          resampling=0.1)
        centerline_bif = compute_centerlines(outlet1, outlet2,
                                             centerline_bif_path,
                                             capped_surface, resampling=0.1)
    else:
        centerline_relevat_outlets = None
        centerline_bif = None

    # Create a tolerance for diverging
    tolerance = getTolerance(centerline_par)

    # Get data from centerlines
    data = getData(centerline_relevant_outlets, centerline_dau1, centerline_dau2, tolerance, aneurysm_type)
    write_parameters(data, dirpath)

    # Compute and smooth voornoi diagram (not aneurysm)
    print("Compute voronoi diagram")
    voronoi = makeVoronoiDiagram(surface, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        parameters = get_parameters(dirpath)
        number_of_aneurysms = len([a for a in parameters.keys() if "aneurysm_" in a])
        if number_of_aneurysms == 1:
            voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerline_par, 0.25)
        else:
            aneu_centerline = ExtractSingleLine(centerline_complete,
                                                centerline_complete.GetNumberOfCells() - 1)
            div_aneu_id = []
            for i in range(centerline_complete.GetNumberOfCells()-1):
                div_aneu_id.append(centerline_div(aneu_centerline,
                                                  ExtractSingleLine(centerline_complete, i)))
            div_aneu_id = max(div_aneu_id)
            aneu_centerline = ExtractSingleLine(aneu_centerline, start=div_aneu_id)
            voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi,
                                                          centerline_par, 0.25,
                                                          no_smooth=aneu_centerline)

        voronoi_smoothed = remove_extreme_points(voronoi_smoothed, voronoi)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)

        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)

    # divpoints and endpoints, for bif or lower
    key = "div_point"
    div_points = get_points(data, key, bif=aneurysm_type=="lateral")
    if aneurysm_type == "terminal":
        div_points_bif = get_points(data, key, bif=True)

    key = "end_point"
    end_points = get_points(data, key, bif=aneurysm_type=="lateral")
    if aneurysm_type == "terminal":
        end_points_bif = get_points(data, key, bif=True)

    center = ((3/9.)*div_points[1][0] + (3/9.)*div_points[1][1] + \
                            (3/9.)*div_points[1][2]).tolist()

    write_points(div_points[0], points_div_path)
    write_points(end_points[0], points_clipp_path)

    print("Clipping centerlines and voronoi diagram.")
    patch_cl = CreateParentArteryPatches(centerline_par, end_points[0])
    WritePolyData(patch_cl, centerline_clipped_path)

    if aneurysm_type == "terminal":
        patch_bif_cl = CreateParentArteryPatches(centerline_bif, end_points_bif[0])
        WritePolyData(patch_bif_cl, centerline_clipped_bif_path)

    # Clipp the voronoi diagram
    masked_voronoi = MaskVoronoiDiagram(voronoi, patch_cl)
    voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    WritePolyData(voronoi_clipped, voronoi_clipped_path)

    # Interpolate the centerline
    print("Interpolate centerlines and voronoi diagram.")
    interpolated_cl = InterpolatePatchCenterlines(patch_cl, centerline_par,
                                                  div_points[0].GetPoint(0),
                                                  None, False)
    WritePolyData(interpolated_cl, centerline_new_path.replace(".vtp", "1.vtp"))

    if aneurysm_type == "terminal":
        if bif:
            print("Start interpolate bif")
            # TODO: Fix interpolate bif
            interpolated_bif = InterpolatePatchCenterlines(patch_bif_cl, centerline_bif,
                                                            None, "bif", True)
            WritePolyData(interpolated_bif, centerline_new_bif_path)

        if lower:
            print("Start interpolate lower")
            div_points_bif[0].SetPoint(0, center[0], center[1], center[2])
            interpolated_bif_lower = InterpolatePatchCenterlines(patch_bif_cl, centerline_bif,
                                                                div_points_bif[0].GetPoint(0),
                                                                "lower", True)
            WritePolyData(interpolated_bif_lower, centerline_new_bif_lower_path)

    print("Start merge")
    interpolated_cl = merge_cl(interpolated_cl, div_points[1],
                               end_points[1])
    WritePolyData(interpolated_cl, centerline_new_path)

    # Interpolate voronoi diagram
    bif_ = []
    if aneurysm_type == "terminal":
        if lower and bif:
            bif_ = [interpolated_bif, interpolated_bif_lower, patch_bif_cl]
        elif bif:
            bif_ = [interpolated_bif, patch_bif_cl]
        elif lower:
            bif_ = [interpolated_bif_lower, patch_bif_cl]

    print("Start interpolate voronoi diagram")
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, patch_cl,
                                                       voronoi_clipped,
                                                       end_points[0],
                                                       bif_, cylinder_factor)
    WritePolyData(interpolated_voronoi, voronoi_anu_path.replace(".vtp", "_remove.vtp"))
    interpolated_voronoi = remove_extreme_points(interpolated_voronoi, voronoi_clipped)
    WritePolyData(interpolated_voronoi, voronoi_anu_path)

    # Write a new surface from the new voronoi diagram
    print("Create new surface")
    new_surface = create_new_surface(interpolated_voronoi)
    WritePolyData(new_surface, model_new_surface)

    return new_surface


if  __name__ == "__main__":
    smooth, smooth_factor, bif, addPoint, basedir, case, lower, \
    cylinder_factor, aneurysm_type, anu_num = read_command_line()
    #main(path.join(path.dirname(path.abspath(__file__)), "C0037"), smooth,
    #      smooth_factor, bif, addPoint, lower, cylinder_factor, aneurysm_type)
    folders = listdir(basedir) if case is None else [case]
    for folder in folders:
        if folder[:2] in ["P0", "C0"]:
            #if not path.exists(path.join(basedir, folder, "model_anu_bif_smooth_extraPoint_lower.vtp")):
            print("Working on case", folder)
            main(path.join(basedir, folder), smooth, smooth_factor,
                    bif, addPoint, lower, cylinder_factor, aneurysm_type, anu_num)
