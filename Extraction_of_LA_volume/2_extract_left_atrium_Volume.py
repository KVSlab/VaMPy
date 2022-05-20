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
from os import path, makedirs, getcwd
from mpi4py import MPI

try:
    from morphman.common import *
    from vmtk import vmtkscripts
    from vmtk import vtkvmtk
    import vmtk
    import vtk
except:
    raise ImportError("The scipt is dependent on morphMan, for install instructions see" + \
                      " https://morphman.readthedocs.io/en/latest/installation.html")

try:
    from vmtkpointselector import *
except ImportError:
    pass


def vtk_merge_UnstructuredGridData(inputs):
    """
    Appends one or more UnstructuredGrid
    datates together into a single
    UnstructuredGrid dataset.

    Args:
        inputs (list): List of vtkUnstructuredGridData objects.

    Returns:
        merged_data (vtkUnstructuredGridData): Single UnstructuredGrid dataset.
    """
    append_filter = vtk.vtkAppendFilter()
    for input_ in inputs:
        append_filter.AddInputData(input_)
    append_filter.Update()
    merged_data = append_filter.GetOutput()

    return merged_data


def vtk_compute_connectivity_VOLUME(volume, mode="All", closest_point=None, show_color_regions=True,
                             mark_visited_points=False):
    """Wrapper of vtkPolyDataConnectivityFilter. Compute connectivity.

    Args:
        show_color_regions (bool): Turn on/off the coloring of connected regions.
        mark_visited_points (bool): Specify whether to record input point ids that appear in the output.
        volume (vtkUnstructuredGridData): Input volume data.
        mode (str): Type of connectivity filter.
        closest_point (list): Point to be used for mode='Closest'
    """
    connectivity = vtk.vtkConnectivityFilter()
    connectivity.SetInputData(volume)

    # Mark each region with "RegionId"
    if mode == "All":
        connectivity.SetExtractionModeToAllRegions()
    elif mode == "Largest":
        connectivity.SetExtractionModeToLargestRegion()
    elif mode == "Closest":
        if closest_point is None:
            print("ERROR: point not set for extracting closest region")
            sys.exit(0)
        connectivity.SetExtractionModeToClosestPointRegion()
        connectivity.SetClosestPoint(closest_point)

    if show_color_regions:
        connectivity.ColorRegionsOn()

    if mark_visited_points:
        connectivity.MarkVisitedPointIdsOn()

    connectivity.Update()
    output = connectivity.GetOutput()

    return output

def attach_clipped_regions_to_volume(volume, clipped, center):
    """Check the connectivty of a clipped volume, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        volume  (vtkUnstructuredGridData):
        clipped (vtkUnstructuredGridData): The clipped segments of the volume.
        center (list): The center of the clipping point

    Returns:
        volume (vtkUnstructuredGridData): The volume where only one segment has been removed.
    """

    connectivity = vtk_compute_connectivity_VOLUME(clipped, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return volume
    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0))
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(get_distance(region_point, center))

    # Remove the region with the closest distance
    print("distances.index", distances.index(min(distances)))
    regions.pop(distances.index(min(distances)))

    # Add the other regions back to the surface
    volume = vtk_merge_UnstructuredGridData(regions + [volume])
    #volume = vtk_clean_polydata(volume)
    #volume = vtk_triangulate_surface(volume)

    return volume

def vtk_clip_UnstructuredGridData(volume, cutter=None, value=0.0, get_inside_out=False, generate_clip_scalars=False):
    """Clip the input vtkUnstructuredGridData object with a cutter function (plane, box, etc)

    Args:
        generate_clip_scalars (bool): If True, output scalar values will be interpolated from implicit function values.
        get_inside_out (bool): Get inside out, default is False
        volume (vtkUnstructuredGridData): Input vtkUnstructuredGridData for clipping
        cutter (vtkBox, vtkPlane): Function for cutting the polydata (default None).
        value (float): Distance to the ImplicteFunction or scalar value to clip.

    Returns:
        clipper (vtkUnstructuredGridData): The clipped volume
    """
    clipper = vtk.vtkClipDataSet()
    clipper.SetInputData(volume)
    if cutter is None:
        clipper.GenerateClipScalarsOff()
    else:
        clipper.SetClipFunction(cutter)
    if get_inside_out:
        clipper.InsideOutOn()
    if generate_clip_scalars and cutter is not None:
        clipper.GenerateClipScalarsOn()
    clipper.SetValue(value)
    clipper.Update()

    return clipper.GetOutput()

#def get_surface_closest_to_point(clipped, point):
    """Check the connectivty of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        clipped (vtkPolyData): The clipped segments of the surface.
        point (list): The point of interest. Keep region closest to this point

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """
    """connectivity = vtk_compute_connectivity(clipped, mode="All")
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

    return region_of_interest"""

def get_surface_closest_to_point(volume, clipped, center):
    """Check the connectivty of a clipped volume, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        volume (vtkUnstructuredGridData):
        clipped (vtkUnstructuredGridData): The clipped segments of the volume.
        center (list): The center of the clipping point

    Returns:
        volume (vtkUnstructuredGridData): The volume where only one segment has been removed.
    """

    connectivity = vtk_compute_connectivity_VOLUME(volume, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return volume
    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0))
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(get_distance(region_point, center))
        #write_polydata(regions[i], la_and_laa_path+str(i)+".vtp")

    

    # Add the other regions back to the surface
    volume = regions[distances.index(min(distances))] 

    # Remove the region with the closest distance
    regions.pop(distances.index(min(distances)))

    clipped = vtk_merge_UnstructuredGridData(regions + [clipped])
    #clipped = vtk_clean_polydata(clipped)
    #clipped = vtk_triangulate_surface(clipped)

    return volume, clipped



def get_surface_closest_to_point_SURFACE(surface, clipped, center):
    """Check the connectivty of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        surface (vtkPolyData):
        clipped (vtkPolyData): The clipped segments of the surface.
        center (list): The center of the clipping point

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """

    connectivity = vtk_compute_connectivity(surface, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return surface
    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0))
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(get_distance(region_point, center))
        regions[i] = vtk_convert_unstructured_grid_to_polydata(regions[i])
        #write_polydata(regions[i], la_and_laa_path+str(i)+".vtp")

    

    # Add the other regions back to the surface
    surface = vtk_convert_unstructured_grid_to_polydata( regions[distances.index(min(distances))] )
    surface = vtk_clean_polydata(surface)
    surface = vtk_triangulate_surface(surface)

    # Remove the region with the closest distance
    print("distances_surface.index", distances.index(min(distances)))
    regions.pop(distances.index(min(distances)))

    clipped = vtk_convert_unstructured_grid_to_polydata(clipped)
    clipped = vtk_merge_polydata(regions + [clipped])
    clipped = vtk_clean_polydata(clipped)
    clipped = vtk_triangulate_surface(clipped)

    return surface, clipped


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
    base_path = get_path_names(input_path)
    model_name = base_path.split("/")[-1]

    if "_" in model_name:
        model_name = model_name.split("_")[0]
        tstep=get_path_names(input_path_2).split("/")[-1].split("_")[-1]
        print("tstep=", tstep)
        base_path = '/'.join(base_path.split("/")[:-1] + [model_name])

    centerline_path = base_path + "_centerline.vtp"
    la_and_laa_path = base_path + "_la_and_laa_{}.vtu".format(str(tstep))
    la_and_laa_path_sur = base_path + "_la_and_laa.vtp"

    #la_and_laa_path_clipped = base_path + "_la_and_laa_clipped.vtu"
    #la_and_laa_path_sur_clipped = base_path + "_la_and_laa_clipped.vtp"


    # Open the surface file.
    print("--- Load model file\n")
    surface = read_polydata(input_path)

    # Open the volume file
    print("--- Load 3D geometry file\n")
    volume = read_polydata(input_path_2)
    volume_original = volume

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
#    for i in range(1):
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
        print("center", center)
        print("normal", normal)
        plane = vtk_plane(center, normal)
        plane_2 = vtk_plane(center, -1*normal)
        clipped= vtk_clip_UnstructuredGridData(volume, plane_2)  # new not fancy solution
        volume = vtk_clip_UnstructuredGridData(volume,   plane)             # new not fancy solution
        surface, clipped_sur = vtk_clip_polydata(surface,plane)

        # Find part to keep
        #volume = vtk_clean_polydata(volume, input_path_2)
        #clipped = vtk_clean_polydata(clipped, input_path_2)
        surface = vtk_clean_polydata(surface)
        clipped_sur = vtk_clean_polydata(clipped_sur)
        p_boundary = np.array(line_tmp.GetPoint(line_tmp.GetNumberOfPoints() - 1))

        surf_loc = get_vtk_point_locator(surface)
        clip_loc = get_vtk_point_locator(clipped_sur)
        id_surf = surf_loc.FindClosestPoint(p_boundary)
        id_clip = clip_loc.FindClosestPoint(p_boundary)
        p_surface = np.array(surface.GetPoint(id_surf))
        p_clipped = np.array(clipped_sur.GetPoint(id_clip))
        dist_surface = np.linalg.norm(p_surface - p_boundary)
        dist_clipped = np.linalg.norm(p_clipped - p_boundary)
        if dist_surface < dist_clipped:
            volume, clipped = clipped, volume
            surface, clipped_sur = clipped_sur, surface

        volume = attach_clipped_regions_to_volume(volume, clipped, center)
        surface = vtk_clean_polydata(surface)
        clipped_sur = vtk_clean_polydata(clipped_sur)
        surface = attach_clipped_regions_to_surface(surface, clipped_sur, center)

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
    print("center", center)
    print("normal", normal)
    plane = vtk_plane(center, normal)
    plane_2 = vtk_plane(center, -1*normal)
    #volume, clipped = vtk_clip_polydata(volume,  plane)
    clipped = vtk_clip_UnstructuredGridData(volume, plane_2)  # new not fancy solution
    volume  = vtk_clip_UnstructuredGridData(volume,   plane)  # new not fancy solution
    surface, clipped_sur = vtk_clip_polydata(surface, plane)

    # Find part to keep
    #volume = vtk_clean_polydata(volume)
    #clipped = vtk_clean_polydata(clipped)
    surface = vtk_clean_polydata(surface)
    clipped_sur = vtk_clean_polydata(clipped_sur)
    p_boundary = p_outlet

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped_sur)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped_sur.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    if dist_surface < dist_clipped:
        volume, clipped = clipped, volume
        surface, clipped_sur = clipped_sur, surface

    volume = attach_clipped_regions_to_volume(volume, clipped, center)
    surface = vtk_clean_polydata(surface)
    clipped_sur = vtk_clean_polydata(clipped_sur)
    surface = attach_clipped_regions_to_surface(surface, clipped_sur, center)

    print("--- Saving LA and LAA to: {}".format(la_and_laa_path))
    write_polydata(volume, la_and_laa_path)

    #print("--- Saving LA and LAA to: {}".format(la_and_laa_path_clipped))
    #write_polydata(clipped, la_and_laa_path_clipped)

    print("--- Saving LA and LAA to: {}".format(la_and_laa_path_sur))
    write_polydata(surface, la_and_laa_path_sur)

    #print("--- Saving LA and LAA to: {}".format(la_and_laa_path_sur_clipped))
    #write_polydata(clipped_sur, la_and_laa_path_sur_clipped)


def extract_LAA_and_LA_body(input_path, input_path_2, laa_point):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along centerlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    #---------------------------------------------------------------------------------------
    # File paths --> vtp
    input_path_original = input_path
    base_path = get_path_names(input_path)
    model_name = base_path.split("/")[-1]
    if "_" in model_name:
        model_name = model_name.split("_")[0]
        base_path = '/'.join(base_path.split("/")[:-1] + [model_name])


    laa_model_path = base_path + "_laa.vtp"
    laa_centerline_path = base_path + "_laa_centerline.vtp"

    # Open the surface file.
    print("--- Load model file\n")
    input_path = base_path + "_la_and_laa.vtp"
    surface = read_polydata(input_path)
    #---------------------------------------------------------------------------------------
    # Create a folder for laa --> vtu
    #common_path_laa = path.join(get_path_names(input_path_original.split("/")[:-1]), "LAA")    
    #if not path.exists(common_path_laa):
    #    makedirs(common_path_laa)

    # Define tstep from the vtu file
    tstep=get_path_names(input_path_2).split("/")[-1].split("_")[-1]
    print("tstep=", tstep)

    laa_model_path_2 = base_path + "_laa_{}.vtu".format(str(tstep))
    #laa_model_path_2 = base_path + "_laa_{}.vtu".format(str(tstep))
    # Open the surface file.
    print("--- Load model file\n")
    input_path_2 = base_path + "_la_and_laa_{}.vtu".format(str(tstep))
    volume = read_polydata(input_path_2)
    #---------------------------------------------------------------------------------------
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
    plane_2 = vtk_plane(center, -1*normal)
    clipped = vtk_clip_UnstructuredGridData(volume, plane_2)  # new not fancy solution
    volume  = vtk_clip_UnstructuredGridData(volume,   plane)  # new not fancy solution
    surface, clipped_sur = vtk_clip_polydata(surface, plane)

    # Find part to keep
    surface = vtk_clean_polydata(surface)
    clipped_sur = vtk_clean_polydata(clipped_sur)
    p_boundary = appendage_point

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped_sur)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped_sur.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    if dist_surface > dist_clipped:
        surface, clipped_sur = clipped_sur, surface
        volume, clipped = clipped, volume

    surface     = vtk_convert_unstructured_grid_to_polydata(surface)
    clipped_sur = vtk_convert_unstructured_grid_to_polydata(clipped_sur)
    surface, clipped_sur = get_surface_closest_to_point_SURFACE(surface, clipped_sur, center)
    volume, clipped = get_surface_closest_to_point(volume, clipped, center)

    #-----------------------------------------------------------------------------
    # Save vtp files
    print("--- Saving LAA to: {}".format(laa_model_path))
    # Fix LAA surface
    cleaned_LAA_surface = vtk_clean_polydata(vtk_convert_unstructured_grid_to_polydata(surface))
    triangulated_LAA_surface = vtk_triangulate_surface(cleaned_LAA_surface)
    write_polydata(triangulated_LAA_surface, laa_model_path)


    la_body_model_path = base_path + "_la_body.vtp"
    print("--- Saving LA body to: {}".format(la_body_model_path))
    
    # Fix LA body surface
    cleaned_LA_body_surface = vtk_clean_polydata(clipped_sur)
    triangulated_LA_body_surface = vtk_triangulate_surface(cleaned_LA_body_surface)
    write_polydata(triangulated_LA_body_surface, la_body_model_path)

    #-------------------------------------------------------------------------------
    # Save vtu files
    print("--- Saving LAA to: {}".format(laa_model_path_2))
    # Fix LAA surface
    write_polydata(volume, laa_model_path_2)


    la_body_model_path_2 = base_path + "_la_body_{}.vtu".format(str(tstep))
    print("--- Saving LA body to: {}".format(la_body_model_path_2))
    # Fix LA body surface
    write_polydata(clipped, la_body_model_path_2)

    # #-------------------------------------------------------------------------------
    # Compute and save quartiles from vtu files
    """print("--- Computing quartiles from LAA to: {}".format(laa_model_path_2))
    # Fix LAA surface
    computeQuartiles_LAA = vtk.vtkComputeQuartiles()
    computeQuartiles_LAA.SetInputData(volume)

    computeQuartiles_LAA.GetTable()
    q_LAA=computeQuartiles_LAA.GetOutput()
    print("q_LAA", q_LAA)"""

    





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    """parser.add_argument("--case_uncapped_geometry2compute_centerlines")
    parser.add_argument("--case_geometry2be_clipped")"""
    parser.add_argument("--laa", default=None, type=float, nargs="+")
    parser.add_argument("--includes_laa_and_la_body", default=0)
    args = parser.parse_args()
    """case_path   = args.case_uncapped_geometry2compute_centerlines
    case_path_2 = args.case_geometry2be_clipped"""
    laa_point = args.laa

    N_tsteps = 200
    t_init=0
    case_path   =  "/home/sergionio/Downloads/test_nu/200_tsteps/C_167_remeshed_surface.vtp"


    for i in range(N_tsteps):
        tstep = t_init + i
        case_path_2 =  "/home/sergionio/Downloads/test_nu/200_tsteps/nu_relative_last_tstep_{}.vtu".format(str(tstep))

        print("--- Extracting geometry from: {} to compute the centerlines".format(case_path))
        print("--- Extracting geometry from: {} to clip its geometry".format(case_path_2))

        scale = 1  # Get seconds
        t0 = time.time()
        extract_LA_and_LAA(case_path, case_path_2)
        t1 = time.time()
        print("--- LA Extraction complete")
        print("--- Time spent extracting LA & LAA: {:.3f} s".format((t1 - t0) / scale))

        if args.includes_laa_and_la_body != 0:
            extract_LAA_and_LA_body(case_path, case_path_2, laa_point)
            t2 = time.time()
            print("--- LAA and LA body Extraction complete")
            print("--- Time spent extracting LAA and LA body: {:.3f} s".format((t2 - t1) / scale))

    """print("--- Extracting geometry from: {} to compute the centerlines".format(case_path))
    print("--- Extracting geometry from: {} to clip its geometry".format(case_path_2))

    scale = 1  # Get seconds
    t0 = time.time()
    extract_LA_and_LAA(case_path, case_path_2)
    t1 = time.time()
    print("--- LA Extraction complete")
    print("--- Time spent extracting LA & LAA: {:.3f} s".format((t1 - t0) / scale))

    if args.includes_laa_and_la_body != 0:
        extract_LAA_and_LA_body(case_path, case_path_2, laa_point)
        t2 = time.time()
        print("--- LAA and LA body Extraction complete")
        print("--- Time spent extracting LAA and LA body: {:.3f} s".format((t2 - t1) / scale))"""