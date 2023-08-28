## Adhoc landmarking script of the left atrium
## method is a generalized version of the method found in:
## Tobon-Gomez, Catalina, et al. "Benchmark for algorithms segmenting
##   the left atrium from 3D CT and MRI datasets." IEEE transactions on medical
##   imaging 34.7 (2015): 1460-1473.

## Test data can be aquired from the Left Atrium Segmentation Challenge
## http://www.cardiacatlas.org/challenges/left-atrium-segmentation-challenge/

## Writen by Aslak W. Bergersen, 2019
## Modified by Henrik A. Kjeldsberg, 2022

import argparse
import time
from os import getcwd, makedirs, path
from pathlib import Path

try:
    from morphman.common import *
    from vmtk import vmtkscripts
except:
    raise ImportError("The scipt is dependent on morphMan, for install instructions see" + \
                      " https://morphman.readthedocs.io/en/latest/installation.html")

try:
    from vmtkpointselector import *
except ImportError:
    pass


def get_surface_closest_to_point(clipped, point):
    """Check the connectivty of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        clipped (vtkPolyData): The clipped segments of the surface.
        point (list): The point of interest. Keep region closest to this point

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """
    connectivity = vtk_compute_connectivity(clipped, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return clipped

    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0))
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(point))
        distances.append(get_distance(region_point, point))

    # Remove the region with the closest distance
    region_of_interest = regions[distances.index(min(distances))]

    return region_of_interest


def provide_region_points(surface, provided_points, dir_path=None):
    """
    Get relevant region points from user selected points on a input surface.

    Args:
        provided_points (ndarray): Point(s) representing area to refine
        surface (vtkPolyData): Surface model.
        dir_path (str): Location of info.json file

    Returns:
        points (list): List of relevant outlet IDs
    """
    # Fix surface
    cleaned_surface = vtk_clean_polydata(surface)
    triangulated_surface = vtk_triangulate_surface(cleaned_surface)

    if provided_points is None:
        # Select seeds
        print("--- Please select regions to refine in rendered window")
        SeedSelector = vmtkPickPointSeedSelector()
        SeedSelector.SetSurface(triangulated_surface)
        SeedSelector.Execute()

        regionSeedIds = SeedSelector.GetTargetSeedIds()
        get_point = triangulated_surface.GetPoints().GetPoint
        points = [list(get_point(regionSeedIds.GetId(i))) for i in range(regionSeedIds.GetNumberOfIds())]
    else:
        surface_locator = get_vtk_point_locator(surface)
        provided_points = [[provided_points[3 * i], provided_points[3 * i + 1], provided_points[3 * i + 2]]
                           for i in range(len(provided_points) // 3)]

        points = [list(surface.GetPoint(surface_locator.FindClosestPoint(p))) for p in provided_points]

    if dir_path is not None:
        info = {"number_of_regions": len(points)}

        for i in range(len(points)):
            info["region_%d" % i] = points[i]
            write_parameters(info, dir_path)

    return points


def extract_LA_and_LAA(input_path, input_path_2):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along centerlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    base_path = get_path_names(input_path_2)
    model_name = base_path.split("/")[-1]
    if "_" in model_name:
        model_name = model_name.split("_")[0]
        base_path = '/'.join(base_path.split("/")[:-1] + [model_name])

    centerline_path = base_path + "_centerline.vtp"
    # Change the path creating a new directory for the LA + LAA (without PVs and MV)
    """common_path_vtp = path.join(base_path, 'LA_and_LAA_without_PVs_and_MV')
    if not path.exists(common_path_vtp):
        makedirs(common_path_vtp)
    original_vtp_path = common_path_vtp"""
    #base_path = Path(common_path_vtp)
    la_and_laa_path = base_path + "_la_and_laa.vtp"

    # Open the surface file.
    print("--- Load model file\n")
    surface = read_polydata(input_path)

    # Open the volume file
    print("--- Load 3D geometry file\n")
    volume = read_polydata(input_path_2)

    if is_surface_capped(surface)[0]:
        capped_surface = surface
        surface = get_uncapped_surface(surface, gradients_limit=0.01, area_limit=20, circleness_limit=5)
    else:
        capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface, base_path)
    p_outlet = np.array(inlet)

    # Make centerlines
    # Check if voronoi and pole_ids exists
    centerlines, _, _ = compute_centerlines(inlet, outlets, centerline_path, capped_surface,
                                            resampling=0.1, smooth=False,
                                            base_path=base_path)

    # Clip PVs
    print("--- Clipping PVs")
    for i in range((len(outlets) // 3)):
        line_tmp = extract_single_line(centerlines, i)
        line = extract_single_line(line_tmp, 0, start_id=20, end_id=line_tmp.GetNumberOfPoints() - 20)
        line = compute_splined_centerline(line, nknots=10, isline=True)

        # Resample line
        l = get_curvilinear_coordinate(line)
        step = 5 * np.mean(l[1:] - l[:-1])
        line = vmtk_resample_centerline(line, step)

        area, sections = vmtk_compute_centerline_sections(capped_surface, line)

        # Get arrays
        a = get_point_data_array("CenterlineSectionArea", area)
        n = get_point_data_array("FrenetTangent", area, k=3)
        l = get_curvilinear_coordinate(area)

        # Compute 'derivative' of the area
        dAdX = (a[1:, 0] - a[:-1, 0]) / (l[1:] - l[:-1])

        # Check only from "middle" of lumen and towards PV
        half_dAdX = int(len(dAdX) / 2)
        dAdX = dAdX[half_dAdX:]

        # tolerance for clipping the PVs and MV!!! --> Vary the number after dAdX < - "Number2vary"
        # if it is C_167 ->>> stop_id = np.nonzero(dAdX < -80)[0][-1] + half_dAdX + 3
        # else ->>> stop_id = np.nonzero(dAdX < -100)[0][-1] + half_dAdX + 3
        stop_id = np.nonzero(dAdX < -70)[0][-1] + half_dAdX + 3

        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip PVs of the geometry
        plane = vtk_plane(center, normal)
        volume, clipped = vtk_clip_polydata(volume, plane)

        # Find part to keep
        volume = vtk_clean_polydata(volume)#, input_path_2)
        clipped = vtk_clean_polydata(clipped)#, input_path_2)
        p_boundary = np.array(line_tmp.GetPoint(line_tmp.GetNumberOfPoints() - 1))

        surf_loc = get_vtk_point_locator(volume)
        clip_loc = get_vtk_point_locator(clipped)
        id_surf = surf_loc.FindClosestPoint(p_boundary)
        id_clip = clip_loc.FindClosestPoint(p_boundary)
        p_surface = np.array(volume.GetPoint(id_surf))
        p_clipped = np.array(clipped.GetPoint(id_clip))
        dist_surface = np.linalg.norm(p_surface - p_boundary)
        dist_clipped = np.linalg.norm(p_clipped - p_boundary)
        if dist_surface < dist_clipped:
            volume, clipped = clipped, volume

        volume = attach_clipped_regions_to_surface(volume, clipped, center)

    # Clip MV
    print("--- Clipping MV")
    line = extract_single_line(centerlines, 0)
    line = extract_single_line(line, 0, start_id=100, end_id=line.GetNumberOfPoints() - 40)
    line = compute_splined_centerline(line, nknots=10, isline=True)

    l = get_curvilinear_coordinate(line)
    step = 5 * np.mean(l[1:] - l[:-1])
    line = vmtk_resample_centerline(line, step)

    area, sections = vmtk_compute_centerline_sections(capped_surface, line)

    # Get arrays
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)
    l = get_curvilinear_coordinate(area)

    # Compute 'derivative' of the area
    dAdX = (a[1:, 0] - a[:-1, 0]) / (l[1:] - l[:-1])
    stop_id = np.nonzero(dAdX > 20)[0][0]

    normal = -n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    volume, clipped = vtk_clip_polydata(volume, plane)#, input_path_2)

    # Find part to keep
    volume = vtk_clean_polydata(volume)
    clipped = vtk_clean_polydata(clipped)
    p_boundary = p_outlet

    surf_loc = get_vtk_point_locator(volume)
    clip_loc = get_vtk_point_locator(clipped)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(volume.GetPoint(id_surf))
    p_clipped = np.array(clipped.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    if dist_surface < dist_clipped:
        volume, clipped = clipped, volume

    volume = attach_clipped_regions_to_surface(volume, clipped, center)

    print("--- Saving LA and LAA to: {}".format(la_and_laa_path))
    write_polydata(volume, la_and_laa_path)


def extract_LAA_and_LA_body(input_path, laa_point):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along centerlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    base_path = get_path_names(input_path)
    model_name = base_path.split("/")[-1]
    if "_" in model_name:
        model_name = model_name.split("_")[0]
        base_path = '/'.join(base_path.split("/")[:-1] + [model_name])

    # Change the path creating a new directory for the LA + LAA (without PVs and MV)
    """common_path_vtp = path.join(base_path, 'LAbody_and_LAA_isolated')
    if not path.exists(common_path_vtp):
        makedirs(common_path_vtp)
    original_vtp_path = base_path
    base_path = common_path_vtp"""


    laa_model_path = base_path + "_laa.vtp"
    laa_centerline_path = base_path + "_laa_centerline.vtp"

    # Open the surface file.
    print("--- Load model file\n")
    input_path = base_path + "_la_and_laa.vtp"
    surface = read_polydata(input_path)

    if is_surface_capped(surface)[0]:
        capped_surface = surface
        surface = get_uncapped_surface(surface, gradients_limit=0.01, area_limit=20, circleness_limit=5)
    else:
        capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface, base_path)
    p_outlet = inlet

    # Get area and corresponding centers
    parameters = get_parameters(base_path)

    # Check if LAA exists in parameters
    if "region_0" in parameters.keys():
        appendage_point = parameters["region_0"]
    else:
        p_laa = provide_region_points(capped_surface, laa_point, None)
        appendage_point = p_laa[0]

    print("--- LAA defined at point: {:.6f} {:.6f} {:.6f}"
          .format(appendage_point[0], appendage_point[1], appendage_point[2]))

    # Compute centerline to orifice from MV to get tangent
    laa_centerlines, _, _ = compute_centerlines(p_outlet, appendage_point, laa_centerline_path, capped_surface,
                                                resampling=0.1, smooth=False, base_path=base_path)

    id_start = int(laa_centerlines.GetNumberOfPoints() * 0.2)
    id_stop = int(laa_centerlines.GetNumberOfPoints() * 0.8)
    line = extract_single_line(laa_centerlines, 0, start_id=id_start, end_id=id_stop)
    write_polydata(line, base_path + "_cl_to_check.vtp")
    laa_l = get_curvilinear_coordinate(line)
    step = 10 * np.mean(laa_l[1:] - laa_l[:-1])
    line = vmtk_resample_centerline(line, step)
    line = compute_splined_centerline(line, nknots=10, isline=True)
    area, sections = vmtk_compute_centerline_sections(surface, line)

    # Get arrays
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)
    l = get_curvilinear_coordinate(area)

    # Compute 'derivative' of the area
    dAdX = (a[1:, 0] - a[:-1, 0]) / (l[1:] - l[:-1])

    # Check only from "middle" of lumen and towards PV
    half_dAdX = int(len(dAdX) / 4)
    dAdX = dAdX[half_dAdX:]

    # Stopping criteria
    tolerance_for_clipping_laa = 3  # FIXME: Input parameter? --> Vary the number after dAdX < - "Number2vary"
    # if it is C_167: stop_id = np.nonzero(dAdX < -600)[0][-1] + half_dAdX + tolerance_for_clipping_laa # WORKED
    stop_id = np.nonzero(dAdX < -200)[0][-1] + half_dAdX + tolerance_for_clipping_laa
    normal = n[stop_id]
    center = area.GetPoint(stop_id)

    total_surface = surface # I need to save this surface to use in the LA body
    # Clip the LAA model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)

    # Find part to keep
    surface = vtk_clean_polydata(surface)
    clipped = vtk_clean_polydata(clipped)
    p_boundary = appendage_point

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    if dist_surface > dist_clipped:
        surface, clipped = clipped, surface

    surface = get_surface_closest_to_point(surface, center)

    print("--- Saving LAA to: {}".format(laa_model_path))
    # Fix LAA surface
    cleaned_LAA_surface = vtk_clean_polydata(surface)
    triangulated_LAA_surface = vtk_triangulate_surface(cleaned_LAA_surface)
    write_polydata(triangulated_LAA_surface, laa_model_path)


    la_body_model_path = base_path + "_la_body.vtp"
    print("--- Saving LA body to: {}".format(la_body_model_path))
    # Fix LA body surface
    cleaned_LA_body_surface = vtk_clean_polydata(clipped)
    triangulated_LA_body_surface = vtk_triangulate_surface(cleaned_LA_body_surface)
    write_polydata(triangulated_LA_body_surface, la_body_model_path)





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--case_uncapped_geometry2compute_centerlines")
    parser.add_argument("--case_geometry2be_clipped")
    parser.add_argument("--laa", default=None, type=float, nargs="+")
    parser.add_argument("--includes_laa_and_la_body", default=0)
    args = parser.parse_args()
    case_path   = args.case_uncapped_geometry2compute_centerlines
    case_path_2 = args.case_geometry2be_clipped
    laa_point = args.laa

    print("--- Extracting geometry from: {} to compute the centerlines".format(case_path))
    print("--- Extracting geometry from: {} to clip its geometry".format(case_path_2))

    scale = 1  # Get seconds
    t0 = time.time()
    HI = ["TAWSS", "OSI", "RRT", "ECAP", "TWSSG"]
    original_vtp_path = case_path_2
    for hi in HI:
        case_path_2 = path.join(original_vtp_path, '{}'.format(hi))
        if not path.exists(case_path_2):
            makedirs(case_path_2)
        case_path_2 = Path(case_path_2)
        case_path_2 = case_path_2 / "{}.vtp".format(hi)
        extract_LA_and_LAA(str(case_path), str(case_path_2))
        t1 = time.time()
        print("--- LA Extraction complete")
        print("--- Time spent extracting LA & LAA: {:.3f} s".format((t1 - t0) / scale))

        if args.includes_laa_and_la_body != 0:
            extract_LAA_and_LA_body(str(case_path_2), laa_point)
            t2 = time.time()
            print("--- LAA and LA body Extraction complete")
            print("--- Time spent extracting LAA and LA body: {:.3f} s".format((t2 - t1) / scale))