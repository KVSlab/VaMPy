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
        get_point = surface.GetPoints().GetPoint
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


def extract_LA_and_LAA(input_path):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    base_path = get_path_names(input_path)

    centerline_path = base_path + "_centerline.vtp"
    la_and_laa_path = base_path + "_la_and_laa.vtp"

    # Open the surface file.
    print("--- Load model file\n")
    surface = read_polydata(input_path)

    is_capped, _ = is_surface_capped(surface)
    if not is_capped:
        capped_surface = vmtk_cap_polydata(surface)
    else:
        capped_surface = surface

    # Get area and corresponding centers
    parameters = get_parameters(base_path)
    p_outlet = parameters['outlet']  # Get point at MV outlet

    # Get inlets
    inlets = []
    areas = []
    k = 0
    while True:
        try:
            inlet_k = parameters['inlet{}'.format(k)]
            area_k = parameters['inlet{}_area'.format(k)]
            inlets.append(inlet_k)
            areas.append(area_k)
        except:
            print("--- Found {} inlets".format(k))
            break
        k += 1

    # Sort inlets and get four largest
    sorted_inlets = [x for _, x in sorted(zip(areas, inlets))][-4:]
    inlets_aslist = np.array(sorted_inlets).flatten().tolist()

    # Make centerlines
    # Check if voronoi and pole_ids exists
    centerlines, _, _ = compute_centerlines(p_outlet, inlets_aslist,
                                            centerline_path, capped_surface,
                                            resampling=0.1, smooth=False,
                                            base_path=base_path)

    # Clip PVs
    for i in range(4):
        line_tmp = extract_single_line(centerlines, i)
        line = extract_single_line(line_tmp, 0, start_id=40, end_id=line_tmp.GetNumberOfPoints() - 100)
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

        stop_id = np.nonzero(dAdX < -100)[0][-1] + half_dAdX + 3

        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        surface, clipped = vtk_clip_polydata(surface, plane)

        # Find part to keep
        surface = vtk_clean_polydata(surface)
        clipped = vtk_clean_polydata(clipped)
        p_boundary = np.array(line_tmp.GetPoint(line_tmp.GetNumberOfPoints() - 1))

        surf_loc = get_vtk_point_locator(surface)
        clip_loc = get_vtk_point_locator(clipped)
        id_surf = surf_loc.FindClosestPoint(p_boundary)
        id_clip = clip_loc.FindClosestPoint(p_boundary)
        p_surface = np.array(surface.GetPoint(id_surf))
        p_clipped = np.array(surface.GetPoint(id_clip))
        dist_surface = np.linalg.norm(p_surface - p_boundary)
        dist_clipped = np.linalg.norm(p_clipped - p_boundary)
        if dist_surface < dist_clipped:
            surface, clipped = clipped, surface

        surface = attach_clipped_regions_to_surface(surface, clipped, center)

    # Clip MV
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
    # Stopping criteria
    stop_id = np.nonzero(dAdX > 50)[0][0] + 3
    normal = -n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)
    if surface.GetNumberOfCells() < clipped.GetNumberOfCells():
        surface, clipped = clipped, surface

    surface = attach_clipped_regions_to_surface(surface, clipped, center)

    print("--- Saving LAA to: {}".format(la_and_laa_path))
    write_polydata(surface, la_and_laa_path)


def extract_LAA(input_path, laa_point):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    base_path = get_path_names(input_path)

    laa_model_path = base_path + "_laa.vtp"
    laa_centerline_path = base_path + "_laa_centerline.vtp"

    # Open the surface file.
    print("--- Load model file\n")
    surface = read_polydata(input_path)

    is_capped, _ = is_surface_capped(surface)
    if not is_capped:
        capped_surface = vmtk_cap_polydata(surface)
    else:
        capped_surface = surface

    # Get area and corresponding centers
    parameters = get_parameters(base_path)

    p_outlet = parameters['outlet']  # Get point at MV outlet

    p_laa = provide_region_points(surface, laa_point, None)
    appendage_point = p_laa[0]
    print("--- LAA defined at point: {:.6f} {:.6f} {:.6f}"
          .format(appendage_point[0], appendage_point[1], appendage_point[2]))

    # Compute centerline to orifice from MV to get tangent
    laa_centerlines, _, _ = compute_centerlines(p_outlet, appendage_point, laa_centerline_path, capped_surface,
                                                resampling=0.1, smooth=False, base_path=base_path)

    line = extract_single_line(laa_centerlines, 0, start_id=50, end_id=laa_centerlines.GetNumberOfPoints() - 50)
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

    # Stopping criteria
    stop_id = np.nonzero(dAdX < -500)[0][0] + 6
    normal = n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)
    if surface.GetNumberOfCells() > clipped.GetNumberOfCells():
        surface, clipped = clipped, surface

    surface = get_surface_closest_to_point(surface, center)

    print("--- Saving LAA to: {}".format(laa_model_path))
    write_polydata(surface, laa_model_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--case")
    parser.add_argument("--laa", default=None, type=float, nargs="+")
    args = parser.parse_args()
    case_path = args.case
    laa_point = args.laa

    t0 = time.time()
    # extract_LA_and_LAA(case_path)
    t1 = time.time()
    extract_LAA(case_path, laa_point)
    t2 = time.time()
    print("--- Extraction complete")
    scale = 1  # Get seconds
    print("Time spent extracting LA & LAA: {:.3f} s".format((t1 - t0) / scale))
    print("Time spent extracting LAA: {:.3f} s".format((t2 - t1) / scale))
