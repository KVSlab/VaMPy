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
from os import mkdir
from os.path import isdir

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


def extract_LA_and_LAA(folder, index, cycle, clip_volume=False):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    case = folder.split("/")[-1]
    cyclename = "_cycle_{:02d}".format(int(cycle)) if cycle is not None else ""
    filename = "{}_{}{}".format(case,index, cyclename)
    if clip_volume:
        filetype = ".vtu"
        input_path = path.join(folder,"..", "VTU", filename + filetype)
    else:
        filetype = "vtp"
        input_path = path.join(folder, "VTP", filename + filetype)
    save_path = path.join(folder, "CLIPPED")
    #filename = folder.split("/")[-1].split("_")[0]
    #save_path = folder.rsplit("/", 1)[0]
    #input_path = path.join(folder)  # , "VTU", filename + filetype)

    if not isdir(save_path):
        mkdir(save_path)

    centerline_path = path.join(save_path, filename + "_centerline.vtp")
    la_and_laa_path = path.join(save_path, filename + "_la_and_laa.vtp")
    la_and_laa_path_vtu = path.join(save_path, filename + "_la_and_laa.vtu")

    # Open the surface file.
    print("--- Load model file\n")
    if clip_volume:
        # TODO: Add as input parameter
        surface = read_polydata(
            "/Users/henriakj/PhD/Code/VaMPy/models/models_for_convergence_study_upf_af_n_1/LA_20CYCLE_3M/LA_remeshed_surface.vtp")
        volume = read_polydata(input_path)
    else:
        surface = read_polydata(input_path)

    if is_surface_capped(surface)[0]:
        capped_surface = surface
        surface = get_uncapped_surface(capped_surface, gradients_limit=0.035, area_limit=15, circleness_limit=4)
    else:
        capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface, save_path)
    p_outlet = np.array(inlet)

    # Make centerlines
    # Check if voronoi and pole_ids exists
    centerlines, _, _ = compute_centerlines(inlet, outlets, centerline_path, capped_surface,
                                            resampling=0.1, smooth=False, base_path=save_path)

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

        stop_id = np.nonzero(dAdX < -100)[0][-1] + half_dAdX + 3  # 10

        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        if clip_volume:
            plane_inv = vtk_plane(center, -normal)
            clipped_volume, _ = vtk_clip_polydata(volume, plane_inv, clip_volume=True)
            volume, _ = vtk_clip_polydata(volume, plane, clip_volume=True)

        surface, clipped = vtk_clip_polydata(surface, plane, clip_volume=False)

        # Find part to keep
        surface = vtk_clean_polydata(surface)
        clipped = vtk_clean_polydata(clipped)
        p_boundary = np.array(line_tmp.GetPoint(line_tmp.GetNumberOfPoints() - 1))

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
                volume, clipped_volume = clipped_volume, volume
        surface = attach_clipped_regions_to_surface(surface, clipped, center)
        if clip_volume:
            volume = attach_clipped_regions_to_surface(volume, clipped_volume, center, clip_volume=True)

    # LA006, LA023
    # Clip MV
    print("--- Clipping MV")
    line = extract_single_line(centerlines, 0)
    line = extract_single_line(line, 0, start_id=0, end_id=line.GetNumberOfPoints() - 40)
    line = compute_splined_centerline(line, nknots=10, isline=True)

    l = get_curvilinear_coordinate(line)
    dx = np.mean(l[1:] - l[:-1])
    step = 5 * dx
    line = vmtk_resample_centerline(line, step)

    area, sections = vmtk_compute_centerline_sections(capped_surface, line)

    # Get arrays
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)
    l = get_curvilinear_coordinate(area)

    # Compute 'derivative' of the area
    dAdX = (a[1:, 0] - a[:-1, 0]) / (l[1:] - l[:-1])
    stop_id = 40  # np.nonzero(dAdX > 50)[0][0]

    normal = n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)
    if clip_volume:
        plane_inv = vtk_plane(center, -normal)
        clipped_volume, _ = vtk_clip_polydata(volume, plane_inv, clip_volume=True)
        volume, _ = vtk_clip_polydata(volume, plane, clip_volume=True)

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
            volumee, clipped_volume = clipped_volume, volume

    print("--- Saving LA and LAA to: {}".format(la_and_laa_path))
    surface = attach_clipped_regions_to_surface(surface, clipped, center)
    write_polydata(surface, la_and_laa_path)

    if clip_volume:
        volume = attach_clipped_regions_to_surface(volume, clipped_volume, center, clip_volume=True)
        write_polydata(volume, la_and_laa_path_vtu)


def extract_LA_or_LAA(folder, laa_point, index, cycle, clip_volume=False):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    case = folder.split("/")[-1]
    cyclename = "_cycle_{:02d}".format(int(cycle)) if cycle is not None else ""
    filename = "{}_{}{}".format(case,index, cyclename)

    if clip_volume:
        filetype = ".vtu"
    else:
        filetype = ".vtp"
    save_path = path.join(folder, "CLIPPED")

    # filename = folder.split("/")[-1].split("_")[0]
    # save_path = folder.rsplit("/", 1)[0]
    # input_path = path.join(folder)  # , "VTU", filename + filetype)

    if not isdir(save_path):
        mkdir(save_path)

    laa_centerline_path = path.join(save_path, filename + "_laa_centerline.vtp")
    clipped_model = path.join(save_path, filename + "_la_and_laa.vtp")
    clipped_model_vtu = path.join(save_path, filename + "_la_and_laa.vtu")
    la_model_path = path.join(save_path, filename + "_la" + filetype)
    laa_model_path = path.join(save_path, filename + "_laa" + filetype)

    # Open the surface file.
    print("--- Load model file\n")
    surface = read_polydata(clipped_model)
    if clip_volume:
        volume = read_polydata(clipped_model_vtu)

    if is_surface_capped(surface)[0]:
        capped_surface = surface
        surface = get_uncapped_surface(surface, gradients_limit=0.01, area_limit=20, circleness_limit=5)
    else:
        capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface, save_path)
    p_outlet = inlet

    # Get area and corresponding centers
    parameters = get_parameters(save_path)

    # Check if LAA exists in parameters
    if "region_0" in parameters.keys():
        appendage_point = parameters["region_0"]
    else:
        if laa_point is None:
            model = save_path.split("/")[-1]
            print("Loading LAA points for model {}".format(model))
            laa_points_path = save_path.rsplit("/", 1)[0] + "/LA20_LAA_POINTS.json"
            with open(laa_points_path, "r") as f:
                laa_points = json.load(f)
            appendage_point = laa_points[model]
        else:
            p_laa = provide_region_points(capped_surface, laa_point, None)
            appendage_point = p_laa[0]

    print("--- LAA defined at point: {:.6f} {:.6f} {:.6f}"
          .format(appendage_point[0], appendage_point[1], appendage_point[2]))

    # Compute centerline to orifice from MV to get tangent
    laa_centerlines, _, _ = compute_centerlines(p_outlet, appendage_point, laa_centerline_path, capped_surface,
                                                resampling=0.1, smooth=False, base_path=save_path)
    id_start = int(laa_centerlines.GetNumberOfPoints() * 0.25)
    id_stop = int(laa_centerlines.GetNumberOfPoints() * 0.9)
    # for la013 id_stop = int(laa_centerlines.GetNumberOfPoints() * 0.6)
    line = extract_single_line(laa_centerlines, 0, start_id=id_start, end_id=id_stop)
    laa_l = get_curvilinear_coordinate(line)
    step = 2.5 * np.mean(laa_l[1:] - laa_l[:-1])
    line = vmtk_resample_centerline(line, step)
    line = compute_splined_centerline(line, nknots=10, isline=True)
    area, sections = vmtk_compute_centerline_sections(surface, line)

    # Get arrays
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)
    l = get_curvilinear_coordinate(area)

    # Stopping criteria
    # tolerance = int(0.025 * len(laa_l))  # FIXME: Input parameter?

    dAdX = np.gradient(a.T[0], np.mean(l[1:] - l[:-1]))
    stop_id = np.nonzero(dAdX < -150)[0][-1] + 10  # + tolerance
    normal = n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)
    if clip_volume:
        plane_inv = vtk_plane(center, -normal)
        clipped_volume, _ = vtk_clip_polydata(volume, plane_inv, clip_volume=True)
        volume, _ = vtk_clip_polydata(volume, plane, clip_volume=True)

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

    # Extract LAA:
    if dist_surface > dist_clipped:
        surface, clipped = clipped, surface
        if clip_volume:
            volume, clipped_volume = clipped_volume, volume

    laa_surface = get_surface_closest_to_point(surface, center)
    if clip_volume:
        laa_volume = get_surface_closest_to_point(volume, center, is_volume=True)

    # Extract LA:
    if dist_surface < dist_clipped:
        surface, clipped = clipped, surface
        if clip_volume:
            volume, clipped_volume = clipped_volume, volume

    surface_whole = attach_clipped_regions_to_surface(surface, clipped, center)
    la_surface = get_surface_closest_to_point(surface_whole, center)
    if clip_volume:
        la_volume = attach_clipped_regions_to_surface(volume, clipped_volume, center, clip_volume=True)

    print("--- Saving LAA to: {}".format(laa_model_path))
    if clip_volume:
        write_polydata(laa_volume, laa_model_path)
    else:
        write_polydata(laa_surface, laa_model_path)

    print("--- Saving LA to: {}".format(la_model_path))
    if clip_volume:
        write_polydata(la_volume, la_model_path)
    else:
        write_polydata(la_surface, la_model_path)


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
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder")
    parser.add_argument("--index")
    parser.add_argument("--cycle", default=None)
    parser.add_argument("--laa", default=None, type=float, nargs="+")
    parser.add_argument("--volume", default=0, type=int)
    args = parser.parse_args()
    folder = args.folder
    laa_point = args.laa
    index = args.index
    cycle = args.cycle
    clip_volume = True if int(args.volume) == 1 else False

    print("--- Extracting from: {}".format(folder))

    scale = 1  # Get seconds
    t0 = time.time()
    extract_LA_and_LAA(folder, index, cycle, clip_volume)
    t1 = time.time()
    print("--- LA Extraction complete")
    print("--- Time spent extracting LA & LAA: {:.3f} s".format((t1 - t0) / scale))

    extract_LA_or_LAA(folder, laa_point, index, cycle, clip_volume)
    t2 = time.time()
    print("--- LAA Extraction complete")
    print("--- Time spent extracting LAA: {:.3f} s".format((t2 - t1) / scale))
