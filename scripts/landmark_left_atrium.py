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


def extract_LA_and_LAA(case, condition, is_local, clip_volume=False):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """

    # File paths
    print("--- Load model file\n")
    if is_local:
        save_path = f"/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics"
        save_path_vtu = f"/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics"
        model_path = f'/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_{condition.lower()}/{case}'
    else:
        save_path = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics"
        save_path_vtu = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/FlowMetrics"
        model_path = f'/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/'

    input_path = path.join(save_path, f"hemodynamics.vtp")
    clipped_model = input_path.replace(".vtp", "_la_and_laa.vtp")

    input_path_brt = path.join(save_path_vtu, f"blood_residence_time.vtu")
    input_path_energy = path.join(save_path_vtu, f"energy.vtu")
    volume_brt = read_polydata(input_path_brt)
    volume_energy = read_polydata(input_path_energy)
    clipped_model_brt = input_path_brt.replace(".vtu", "_la_and_laa.vtu")
    clipped_model_energy = input_path_energy.replace(".vtu", "_la_and_laa.vtu")

    new_cl = input_path.replace(".vtp", "_centerline.vtp")
    surface = read_polydata(input_path)
    open_surface_path = path.join(model_path, 'model_remeshed_surface.vtp')
    centerline_path = path.join(model_path, 'model_centerlines.vtp')
    open_surface = read_polydata(open_surface_path)
    capped_surface = surface

    # Centers
    inlet, outlets = compute_centers(open_surface)
    p_outlet = np.array(inlet)

    # Centerline
    capped_surface = vtk_clean_polydata(capped_surface)
    la_centerlines, _, _ = compute_centerlines(inlet, outlets, new_cl, capped_surface, resampling=0.1, smooth=True)
    # Read original
    la_centerline = read_polydata(centerline_path)

    # Clip PVs
    N_PVs = 4
    for i in range(N_PVs):
        print(f"--- Clipping PV ({i + 1})")
        if case == "0008":
            la_centerline_i = extract_single_line(la_centerline, i)
        else:
            la_centerline_i = extract_single_line(la_centerlines, i)
        start = int(la_centerline_i.GetNumberOfPoints() * 0.5)
        stop = int(la_centerline_i.GetNumberOfPoints() * 0.95)

        line = extract_single_line(la_centerline_i, 0, start_id=start, end_id=stop)
        la_l = get_curvilinear_coordinate(line)
        step = 5 * np.mean(la_l[1:] - la_l[:-1])
        line = vmtk_resample_centerline(line, step)
        line = compute_splined_centerline(line, nknots=10, isline=True)
        area, sections = vmtk_compute_centerline_sections(capped_surface, line)

        # Get arrays
        a = get_point_data_array("CenterlineSectionArea", area)
        n = get_point_data_array("FrenetTangent", area, k=3)

        # Compute 'derivative' of the area
        dAdX = np.gradient(a.T[0], step)

        # Find the largest change in cross-sectional area
        tol = 5
        lim = -50
        ID = -1
        if case == "2022" and condition == "sr":
            tol = 10

        stop_id = np.nonzero(dAdX < lim)[0][ID] + tol
        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        if clip_volume:
            plane_inv = vtk_plane(center, -normal)
            # BRT
            clipped_volume_brt, _ = vtk_clip_polydata(volume_brt, plane_inv, clip_volume=True)
            volume_brt, _ = vtk_clip_polydata(volume_brt, plane, clip_volume=True)

            # Energy
            clipped_volume_energy, _ = vtk_clip_polydata(volume_energy, plane_inv, clip_volume=True)
            volume_energy, _ = vtk_clip_polydata(volume_energy, plane, clip_volume=True)

        surface, clipped = vtk_clip_polydata(surface, plane, clip_volume=False)

        # Find part to keep
        surface = vtk_clean_polydata(surface)
        clipped = vtk_clean_polydata(clipped)
        p_boundary = np.array(la_centerline_i.GetPoint(la_centerline_i.GetNumberOfPoints() - 1))

        surf_loc = get_vtk_point_locator(surface)
        clip_loc = get_vtk_point_locator(clipped)
        id_surf = surf_loc.FindClosestPoint(p_boundary)
        id_clip = clip_loc.FindClosestPoint(p_boundary)
        p_surface = np.array(surface.GetPoint(id_surf))
        p_clipped = np.array(clipped.GetPoint(id_clip))
        dist_surface = np.linalg.norm(p_surface - p_boundary)
        dist_clipped = np.linalg.norm(p_clipped - p_boundary)
        if dist_surface < dist_clipped:
            surface, clipped = clipped, surface
            if clip_volume:
                volume_brt, clipped_volume_brt = clipped_volume_brt, volume_brt
                volume_energy, clipped_volume_energy = clipped_volume_energy, volume_energy
        surface = attach_clipped_regions_to_surface(surface, clipped, center)
        if clip_volume:
            volume_brt = attach_clipped_regions_to_surface(volume_brt, clipped_volume_brt, center, clip_volume=True)
            volume_energy = attach_clipped_regions_to_surface(volume_energy, clipped_volume_energy, center,
                                                              clip_volume=True)

    # Clip MV
    print("--- Clipping MV")
    # Load Initial centerline
    la_centerline_v0 = extract_single_line(la_centerline, 0)
    p0 = np.array(la_centerline_v0.GetPoint(0))

    la_centerline_v1 = extract_single_line(la_centerlines, 0)
    p1 = np.array(la_centerline_v1.GetPoint(0))

    # Define normal vector
    n = p0 - p1
    n /= np.linalg.norm(n)

    normal = -n
    center = p0

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)
    if clip_volume:
        plane_inv = vtk_plane(center, -normal)
        clipped_volume_brt, _ = vtk_clip_polydata(volume_brt, plane_inv, clip_volume=True)
        volume_brt, _ = vtk_clip_polydata(volume_brt, plane, clip_volume=True)

        clipped_volume_energy, _ = vtk_clip_polydata(volume_energy, plane_inv, clip_volume=True)
        volume_energy, _ = vtk_clip_polydata(volume_energy, plane, clip_volume=True)
    # Find part to keep
    surface = vtk_clean_polydata(surface)
    clipped = vtk_clean_polydata(clipped)
    p_boundary = p_outlet

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    if dist_surface < dist_clipped:
        surface, clipped = clipped, surface
        if clip_volume:
            volume_brt, clipped_volume_brt = clipped_volume_brt, volume_brt
            volume_energy, clipped_volume_energy = clipped_volume_energy, volume_energy

    print("--- Saving LA and LAA to: {}".format(clipped_model))
    surface = attach_clipped_regions_to_surface(surface, clipped, center)
    # write_polydata(surface, clipped_model)

    if clip_volume:
        print("--- Saving BRT LA and LAA volume to: {}".format(clipped_model_brt))
        volume_brt = attach_clipped_regions_to_surface(volume_brt, clipped_volume_brt, center, clip_volume=True)
        # write_polydata(volume_brt, clipped_model_brt)

        print("--- Saving Energy LA and LAA volume to: {}".format(clipped_model_energy))
        volume_energy = attach_clipped_regions_to_surface(volume_energy, clipped_volume_energy, center,
                                                          clip_volume=True)
        # write_polydata(volume_energy, clipped_model_energy)
    # def separate_LA_and_LAA(case, condition, cycle, laa_point, is_local, clip_volume=False):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # Case 1: Postprocessing
    # folder = PATH/Hemodynamics/CASE or PATH/FlowMetrics/CASE (CASE_INDEX_cycle_N.vtp)

    # Case 2: Clipping model (CASE_clipped.vtp or similar)
    # folder = PATH/CASE

    # File paths
    print("--- Load model file\n")
    if is_local:
        save_path = f"/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics"
        save_path_vtu = f"/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics"
        model_path = f'/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_{condition.lower()}/{case}'
    else:
        save_path = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics"
        save_path_vtu = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/FlowMetrics"
        model_path = f'/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/'

    input_path = path.join(save_path, f"hemodynamics.vtp")
    # clipped_model = input_path.replace(".vtp", "_la_and_laa.vtp")
    laa_model_path = input_path.replace('.vtp', '_laa.vtp')
    la_model_path = input_path.replace('.vtp', '_la.vtp')

    # surface = read_polydata(clipped_model)
    if clip_volume:
        input_path_brt = path.join(save_path_vtu, f"blood_residence_time.vtu")
        input_path_energy = path.join(save_path_vtu, f"energy.vtu")

        clipped_model_brt = input_path_brt.replace('.vtu', '_la_and_laa.vtu')
        laa_model_path_brt = input_path_brt.replace('.vtu', '_laa.vtu')
        la_model_path_brt = input_path_brt.replace('.vtu', '_la.vtu')

        clipped_model_energy = input_path_energy.replace('.vtu', '_la_and_laa.vtu')
        laa_model_path_energy = input_path_energy.replace('.vtu', '_laa.vtu')
        la_model_path_energy = input_path_energy.replace('.vtu', '_la.vtu')

        # volume_brt = read_polydata(clipped_model_brt)
        # volume_energy = read_polydata(clipped_model_energy)
    print("-- Reading centerlines")
    capped_surface = vmtk_cap_polydata(surface)

    laa_centerlines_path = path.join(model_path, "model_region_centerline_0.vtp")
    laa_centerlines = read_polydata(laa_centerlines_path)
    id_start = int(laa_centerlines.GetNumberOfPoints() * 0.1)
    id_stop = int(laa_centerlines.GetNumberOfPoints() * 0.8)


    print("-- Resample, Spline and Section computation" )
    line = extract_single_line(laa_centerlines, 0, start_id=id_start, end_id=id_stop)
    laa_l = get_curvilinear_coordinate(line)
    step = 2 * np.mean(laa_l[1:] - laa_l[:-1])
    line = vmtk_resample_centerline(line, step)
    line = compute_splined_centerline(line, nknots=20, isline=True)
    area, sections = vmtk_compute_centerline_sections(capped_surface, line)

    # Get arrays
    print("-- Extract point arrays" )
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)

    # Stopping criteria
    dAdX = np.gradient(a.T[0], step)
    lim = -200
    tol = 5
    ID = -1
    if case == "0008":
        tol = 20
    if case == "0081":
        tol = 1
    if condition == "sr":
        if case == "0007":
            lim = -100
        if case == "0026":
            ID = 0
            tol = 20
    if condition == "af":
        if case == "0024":
            tol = 0
        if case == "0028":
            ID = -3
            tol = 0
        if case == "0080":
            tol = 6
        if case == "1029":
            ID = -2
        if case == "1039":
            lim = -100
            tol = 7
        if case == "0026":
            tol = 20
        if case == "0004":
            tol = 10
        if case == "0003":
            tol = 0
        if case == "0035":
            lim = -150

    stop_id = np.nonzero(dAdX < lim)[0][ID] + tol
    normal = n[stop_id]
    center = area.GetPoint(stop_id)
    if case == "0081":
        center = np.array([22.36, 100.057, -162.79])
        normal = np.array([0.355, -0.234, 0.904])
    if case == "0022" and condition == "af":
        center = np.array([15.99, 115.17, -169.2])
        normal = np.array([0.3368, 0.5127, 0.789])
    if case == "0078" and condition == "af":
        center = np.array([35.33, 101.25, -172.5])
        normal = np.array([0.647, -0.1069, 0.7549])

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)
    if clip_volume:
        plane_inv = vtk_plane(center, -normal)

        # BRT
        clipped_volume_brt, _ = vtk_clip_polydata(volume_brt, plane_inv, clip_volume=True)
        volume_brt, _ = vtk_clip_polydata(volume_brt, plane, clip_volume=True)

        # Energy
        clipped_volume_energy, _ = vtk_clip_polydata(volume_energy, plane_inv, clip_volume=True)
        volume_energy, _ = vtk_clip_polydata(volume_energy, plane, clip_volume=True)

    # Find part to keep
    surface = vtk_clean_polydata(surface)
    clipped = vtk_clean_polydata(clipped)
    p_boundary = laa_point

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    # Extract LAA:
    if dist_surface > dist_clipped:
        surface, clipped = clipped, surface
        if clip_volume:
            volume_brt, clipped_volume_brt = clipped_volume_brt, volume_brt
            volume_energy, clipped_volume_energy = clipped_volume_energy, volume_energy

    surface_whole = attach_clipped_regions_to_surface(surface, clipped, center)
    laa_surface = get_surface_closest_to_point(surface_whole, center)
    if clip_volume:
        # BRT
        laa_volume_brt = attach_clipped_regions_to_surface(volume_brt, clipped_volume_brt, center, clip_volume=True)
        # laa_volume_brt = merge_dataset(laa_volume_brt)
        laa_volume_brt = get_surface_closest_to_point(laa_volume_brt, center, is_volume=True)
        # Energy
        laa_volume_energy = attach_clipped_regions_to_surface(volume_energy, clipped_volume_energy, center,
                                                              clip_volume=True)
        # laa_volume_energy = merge_dataset(laa_volume_energy)
        laa_volume_energy = get_surface_closest_to_point(laa_volume_energy, center, is_volume=True)

    # Extract LA:
    if dist_surface < dist_clipped:
        surface, clipped = clipped, surface
        if clip_volume:
            volume_brt, clipped_volume_brt = clipped_volume_brt, volume_brt
            volume_energy, clipped_volume_energy = clipped_volume_energy, volume_energy

    surface_whole = attach_clipped_regions_to_surface(surface, clipped, center)
    la_surface = get_surface_closest_to_point(surface_whole, center)
    if clip_volume:
        la_volume_brt = attach_clipped_regions_to_surface(volume_brt, clipped_volume_brt, center, clip_volume=True)
        # la_volume_brt = get_surface_closest_to_point(la_volume_brt_tmp, center, is_volume=True)

        la_volume_energy = attach_clipped_regions_to_surface(volume_energy, clipped_volume_energy, center,
                                                             clip_volume=True)
        # la_volume_energy = get_surface_closest_to_point(la_volume_energy, center, is_volume=True)
    print("--- Saving LAA to: {}".format(laa_model_path))
    write_polydata(laa_surface, laa_model_path)

    print("--- Saving LA to: {}".format(la_model_path))
    write_polydata(la_surface, la_model_path)

    if clip_volume:
        print("--- Saving BRT LAA (vtu) to: {}".format(laa_model_path_brt))
        write_polydata(laa_volume_brt, laa_model_path_brt)

        print("--- Saving BRT LA (vtu) to: {}".format(la_model_path_brt))
        write_polydata(la_volume_brt, la_model_path_brt)

        print("--- Saving Energy LAA (vtu) to: {}".format(laa_model_path_energy))
        write_polydata(laa_volume_energy, laa_model_path_energy)

        print("--- Saving Energy LA (vtu) to: {}".format(la_model_path_energy))
        write_polydata(la_volume_energy, la_model_path_energy)


def vtk_clip_polydata(surface, cutter=None, value=0, get_inside_out=False, generate_clip_scalars=False,
                      clip_volume=False):
    """Clip the input vtkPolyData object with a cutter function (plane, box, etc)

    Args:
        generate_clip_scalars (bool): If True, output scalar values will be interpolated from implicit function values.
        get_inside_out (bool): Get inside out, default is False
        surface (vtkPolyData): Input vtkPolyData for clipping
        cutter (vtkBox, vtkPlane): Function for cutting the polydata (default None).
        value (float): Distance to the ImplicteFunction or scalar value to clip.
        clip_volume (bool): Clips volumetric surface if true.

    Returns:
        clipper (vtkPolyData): The clipped surface
    """
    clipper = vtk.vtkClipDataSet() if clip_volume else vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    if cutter is None:
        clipper.GenerateClipScalarsOff()
    else:
        clipper.SetClipFunction(cutter)
    if get_inside_out:
        clipper.InsideOutOn()
    if generate_clip_scalars and cutter is not None:
        clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(value)
    clipper.Update()

    return clipper.GetOutput(), clipper.GetClippedOutput()


def attach_clipped_regions_to_surface(surface, clipped, center, clip_volume=False):
    """Check the connectivity of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        surface (vtkPolyData):
        clipped (vtkPolyData): The clipped segments of the surface.
        center (list): The center of the clipping point
        clip_volume (bool): Clips volumetric surface if True

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """
    connectivity = vtk_compute_connectivity(clipped, mode="All", is_volume=clip_volume)
    if connectivity.GetNumberOfPoints() == 0:
        return surface
    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(
            vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0, volume=clip_volume))
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(get_distance(region_point, center))

    # Remove the region with the closest distance
    regions.pop(distances.index(min(distances)))

    # Add the other regions back to the surface
    surface = vtk_merge_polydata(regions + [surface], is_volume=clip_volume)
    if not clip_volume:
        surface = vtk_clean_polydata(surface)
        surface = vtk_triangulate_surface(surface)

    return surface


def vtk_merge_polydata(inputs, is_volume=False):
    """
    Appends one or more polygonal
    datasets together into a single
    polygonal dataset.

    Args:
        inputs (list): List of vtkPolyData objects.

    Returns:
        merged_data (vtkPolyData): Single polygonal dataset.
    """
    append_filter = vtk.vtkAppendFilter() if is_volume else vtk.vtkAppendPolyData()
    for input_ in inputs:
        append_filter.AddInputData(input_)
    append_filter.Update()
    merged_data = append_filter.GetOutput()

    return merged_data


def vtk_compute_threshold(surface, name, lower=0, upper=1, threshold_type="between", source=1, volume=False):
    """Wrapper for vtkThreshold. Extract a section of a surface given a criteria.

    Args:
        surface (vtkPolyData): The input data to be extracted.
        name (str): Name of scalar array.
        lower (float): Lower bound.
        upper (float): Upper bound.
        threshold_type (str): Type of threshold (lower, upper, between)
        source (int): PointData or CellData.

    Returns:
        surface (vtkPolyData): The extracted surface based on the lower and upper limit.
    """
    # source = 1 uses cell data as input
    # source = 0 uses point data as input

    # Apply threshold
    vtk_threshold = vtk.vtkThreshold()
    vtk_threshold.SetInputData(surface)
    if threshold_type == "between":
        vtk_threshold.SetLowerThreshold(lower)
        vtk_threshold.SetUpperThreshold(upper)
    elif threshold_type == "lower":
        vtk_threshold.SetLowerThreshold(lower)
    elif threshold_type == "upper":
        vtk_threshold.SetUpperThreshold(upper)
    else:
        print((("%s is not a threshold type. Pleace chose from: upper, lower" +
                ", or between") % threshold_type))
        sys.exit(0)

    vtk_threshold.SetInputArrayToProcess(0, 0, 0, source, name)
    vtk_threshold.Update()
    surface = vtk_threshold.GetOutput()

    # Convert to polydata
    if not volume:
        surface = vtk_convert_unstructured_grid_to_polydata(surface)

    return surface


def vtk_compute_connectivity(surface, mode="All", closest_point=None, show_color_regions=True,
                             mark_visited_points=False, is_volume=False):
    """Wrapper of vtkPolyDataConnectivityFilter. Compute connectivity.

    Args:
        show_color_regions (bool): Turn on/off the coloring of connected regions.
        mark_visited_points (bool): Specify whether to record input point ids that appear in the output.
        surface (vtkPolyData): Input surface data.
        mode (str): Type of connectivity filter.
        closest_point (list): Point to be used for mode='Closest'
    """
    connectivity = vtk.vtkConnectivityFilter() if is_volume else vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)

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


def get_surface_closest_to_point(clipped, point, is_volume=False):
    """Check the connectivty of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        clipped (vtkPolyData): The clipped segments of the surface.
        point (list): The point of interest. Keep region closest to this point

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """
    connectivity = vtk_compute_connectivity(clipped, mode="All", is_volume=is_volume)
    if connectivity.GetNumberOfPoints() == 0:
        return clipped

    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(
            vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1,
                                  source=0, volume=is_volume))
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


if __name__ == "__main__":
    failed = []
    clip_volume = True

    parser = argparse.ArgumentParser()
    parser.add_argument('--case', help='Description for foo argument', required=True)
    parser.add_argument('--condition', help='Description for bar argument', required=True)
    args = parser.parse_args()

    conditions = [args.condition.lower()]
    cases = [args.case]
    local = False
    for condition in conditions:
        for case in cases:
            print("--- Extracting case: {}".format(case))
            # Get LAA point
            print(f"--- Loading LAA points for {case}, condition {condition}")
            if local:
                # Local
                laa_apex_point_path = f"/Users/henriakj/PhD/Code/VaMPy/models/models_inria/laa_apex_points_{condition}.json"
            else:
                # Oracle
                laa_apex_point_path = f"/home/opc/Simulation40/laa_apex_points_{condition}.json"

            with open(laa_apex_point_path) as f:
                info = json.load(f)
            laa_point = info[case][0]
            t0 = time.time()

            t1 = time.time()
            try:
                extract_LA_and_LAA(case, condition, local, clip_volume)
                # separate_LA_and_LAA(case, condition, cycle, laa_point, local, clip_volume)
            except Exception as e:
                print(f"--- FAILED for case {case}, condition {condition}, Error: {e}")
                failed.append(f"{condition}:{case}")

            t2 = time.time()
            print(f"--- Time spent extracting LAA: {t2 - t1:.3f} s")
