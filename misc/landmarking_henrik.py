## Adhoc landmarking script of the left atrium
## method is a generalized version of the method found in:
## Tobon-Gomez, Catalina, et al. "Benchmark for algorithms segmenting
##   the left atrium from 3D CT and MRI datasets." IEEE transactions on medical
##   imaging 34.7 (2015): 1460-1473.

## Test data can be aquired from the Left Atrium Segmentation Challenge
## http://www.cardiacatlas.org/challenges/left-atrium-segmentation-challenge/

## Writen py Aslak W. Bergersen, 2019

import argparse
import time

from IPython import embed

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
        return surface

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


def group_pvs(outlets):
    """Group the pulmonary veins into two groups 'clusters' (left/right).

    Args:
        outlets (list): Flat list of coordinates of the outlets

    Returns:
        inlets0 (list): Flat list of the coordinates of the outlets in group 1.
        inlets1 (list): Flat list of the coordinates of the outlets in group 2.
    """
    # Group the PVs
    inlets = [outlets[i * 3:(i + 1) * 3] for i in range(len(outlets) // 3)][0]
    inlets = inlets.tolist()
    inlet0 = inlets[0]
    inlet1 = inlets[np.argmax(np.sqrt(np.sum((np.array(inlets) - np.array(inlet0)) ** 2, axis=1)))]
    inlets0 = []
    inlets1 = []
    for inlet in inlets:
        dist1 = get_distance(inlet0, inlet)
        dist2 = get_distance(inlet1, inlet)
        if dist1 > dist2:
            inlets1 += inlet
        else:
            inlets0 += inlet

    return inlets0, inlets1


# Overwrite the smooth function from morphMan
def smooth_voronoi_diagram(voronoi, centerlines, remove, smoothing_factor):
    """Implementation of the smoothing algorithm presnted in Ford et al. 2009 Br J Radiol.
    Breifly, the algorithm finds the closest centerline point for each point in the
    Voronoi diagram. Only the Voronoi points where the radius is > (1-smoothing_factor)
    * the maximal inscribed sphere radius of the centerline point. In this implementation
    there is no smoothing at the inlet and outlet, and there is also a threshold on the
    distance to the centerline as well. To gain more control, there is also possible to
    provide an additional centerline (remove) for which all Voronoi points cloest to this
    will be removed.

    Args:
        voronoi (vtkPolyData): Voronoi diagram to be smoothed.
        centerlines (vtkPolyData): Centerline of the geometry.
        remove (vtkPolyData): Centerline to remove the Voronoi from.
        smoothing_factor (float): Smooththing factor.

    Returns:
        smoothed_diagram (vtkPolyData): The smoothed Voronoi diagram.
    """
    number_of_points = voronoi.GetNumberOfPoints()
    thresholds = get_point_data_array(radiusArrayName, centerlines) * (1 - smoothing_factor)

    # Do not smooth inlet and outlets, set threshold to -1
    start = 0
    end = 0
    for i in range(centerlines.GetNumberOfLines()):
        line = extract_single_line(centerlines, i)
        length = get_curvilinear_coordinate(line)
        end_ = line.GetNumberOfPoints() - 1
        end += end_

        # Point buffer start
        end_id = end_ - np.argmin(np.abs(-(length - length.max()) - thresholds[end]))
        start_id = np.argmin(np.abs(length - thresholds[start]))

        diff_start = start_id // 2
        thresholds[start:start + diff_start] = -1
        thresholds[end - end_id // 2:end] = -1
        start += end_ + 1
        end += 1

    # Get locators
    locator = get_vtk_point_locator(centerlines)
    locator_remove = get_vtk_point_locator(remove)

    # Holders for the new data
    smoothed_diagram = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array_numpy = np.zeros(number_of_points)

    # Loop over each Voronoi point
    count = 0
    for i in range(number_of_points):
        point = voronoi.GetPoint(i)
        radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        id_ = locator.FindClosestPoint(point)
        id_r = locator_remove.FindClosestPoint(point)
        r_point = remove.GetPoint(id_r)
        cl_point = centerlines.GetPoint(id_)

        if ((radius >= thresholds[id_] and get_distance(point, cl_point) < thresholds[id_]) \
            or thresholds[id_] < 0) and get_distance(point, cl_point) < get_distance(point, r_point):
            points.InsertNextPoint(point)
            cell_array.InsertNextCell(1)
            cell_array.InsertCellPoint(count)
            radius_array_numpy[count] = radius
            count += 1

    # Add new data to the Voronoi diagram
    radius_array = create_vtk_array(radius_array_numpy, radiusArrayName, k=1)
    smoothed_diagram.SetPoints(points)
    smoothed_diagram.SetVerts(cell_array)
    smoothed_diagram.GetPointData().AddArray(radius_array)

    return smoothed_diagram


def clipp_surface(surface, name, value=0.5):
    """ Wrapper around vtkClipPolyData

    Args:
        surface (vtkPolyData): Surface to clip
        name (str): Name of the array in surface to clip along.
        value (float): Value of the array to clip.

    Returns:
        surface (vtkPolyData): Clipped surface.
        clipped (vtkPolyData): The section that has been removed.
    """
    surface.GetPointData().SetActiveScalars(name)

    surfaceClipper = vtk.vtkClipPolyData()
    surfaceClipper.SetInputData(surface)
    surfaceClipper.SetValue(value)
    surfaceClipper.InsideOutOff()
    surfaceClipper.SetGenerateClippedOutput(1)
    surfaceClipper.Update()

    surface = surfaceClipper.GetOutput()
    clipped = surfaceClipper.GetClippedOutput()

    return surface, clipped


def landmark_atrium(input_path, input_laa_path, output_path):
    """Algorithm for landmarking the left atrium based on the cross-sectional area along
    centerlines, and a smoothed Voronoi diagram similar to what is introduced in Ford et
    al. (2009). The algorithm assumes that there are only four PVs and if there are
    multiple bifurcations included, that the four largest outlets belong to each PV.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """
    # File paths
    base_path = get_path_names(input_path)
    if input_path.endswith(".vtk"):
        surf = read_polydata(input_path)
        if isinstance(surf, vtk.vtkUnstructuredGrid):
            surf = geometry_filter(surf)
        write_polydata(surf, input_path.replace(".vtk", ".vtp"))
        input_path = input_path.replace(".vtk", ".vtp")

    centerline_path = base_path + "_centerline.vtp"
    laa_centerline_path = base_path + "_laa_centerline.vtp"
    smooth_voronoi_path = base_path + "_smooth_voronoi.vtp"
    smoothed_surface_path = base_path + "_smooth_surface.vtp"
    full_centerline_path = base_path + "_full_centerline.vtp"
    only_four_pvs_path = base_path + "_only_four_pvs.vtp"

    # Get an open and closed version of the surface
    surface, capped_surface = prepare_surface(base_path, input_path)
    surface_with_laa, _ = prepare_surface(base_path, input_path)

    # Create two set of centerlines one from PV's to mitral plane, and another from left
    # to right PVs  and Voronoi diagram
    inlet, outlets = compute_centers(surface, base_path)

    # Get area and corresponding centers
    parameters = get_parameters(base_path)
    counter = 0
    area = []
    main_outlets = []
    while "outlet" + str(counter) in parameters.keys():
        area.append(parameters["outlet" + str(counter) + "_area"])
        main_outlets.append(parameters["outlet" + str(counter)])
        counter += 1

    # Get centers of the four largest outlets
    sorted_outlets = np.array(main_outlets)[np.argsort(area)]
    area = np.array(main_outlets)[np.argsort(area)]
    main_outlets = sorted_outlets[-4:]

    # Make centerlines
    # Check if voronoi and pole_ids exists
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, main_outlets.flatten().tolist(),
                                                         centerline_path, capped_surface,
                                                         resampling=0.1, smooth=False,
                                                         base_path=base_path)
    if len(sorted_outlets) > 4:
        full_centerlines, _, _ = compute_centerlines(inlet, outlets, full_centerline_path,
                                                     capped_surface, resampling=0.1,
                                                     smooth=False, pole_ids=pole_ids,
                                                     voronoi=voronoi, base_path=base_path)
        locator = get_vtk_point_locator(centerlines)
        remove = []
        for i in range(len(sorted_outlets)):
            line = extract_single_line(full_centerlines, i)
            for j in range(line.GetNumberOfPoints()):
                point = line.GetPoint(j)
                cl_point = centerlines.GetPoint(locator.FindClosestPoint(point))
                dist = get_distance(cl_point, point)
                if dist > 2:
                    break
            remove.append(extract_single_line(line, 0, start_id=j))
        remove = vtk_merge_polydata(remove)
    else:
        remove = None

    # Smooth Voronoi diagram along centerlines
    if not path.exists(smooth_voronoi_path):
        new_voronoi = smooth_voronoi_diagram(voronoi, centerlines, remove, 0.5)
        write_polydata(new_voronoi, smooth_voronoi_path)
    else:
        new_voronoi = read_polydata(smooth_voronoi_path)

    # Create a new surface
    if not path.exists(smoothed_surface_path):
        new_surface = create_new_surface(new_voronoi, poly_ball_size=[150, 150, 150])
        new_surface = vtk_compute_connectivity(new_surface, mode="Largest")
        write_polydata(new_surface, smoothed_surface_path)
    else:
        new_surface = read_polydata(smoothed_surface_path)

    # Clip outlets
    boundary = vtk_extract_feature_edges(surface)
    normals = []
    for i in range(5):
        if i == 4:
            point = inlet
        else:
            point = main_outlets[i]
        tmp_boundary = vtk_compute_connectivity(boundary, mode="Closest",
                                                closest_point=point)
        p1 = tmp_boundary.GetPoint(0)
        p2 = tmp_boundary.GetPoint(tmp_boundary.GetNumberOfPoints() // 3)
        p3 = tmp_boundary.GetPoint(tmp_boundary.GetNumberOfPoints() * 2 // 3)
        n = np.cross(np.array(p2) - np.array(p1), np.array(p3) - np.array(p1))

        plane = vtk_plane(point, n)
        new_surface, clipped = vtk_clip_polydata(new_surface, plane)
        if new_surface.GetNumberOfCells() < clipped.GetNumberOfCells():
            new_surface, clipped = clipped, new_surface
        new_surface = attach_clipped_regions_to_surface(new_surface, clipped, point)

    # Landmark PV's
    new_surface_capped = vmtk_cap_polydata(new_surface)
    write_polydata(new_surface_capped, only_four_pvs_path)
    for i in range(4):
        line = extract_single_line(centerlines, i)
        line = extract_single_line(line, 0, start_id=40, end_id=line.GetNumberOfPoints() - 100)
        line = compute_splined_centerline(line, nknots=10, isline=True)
        area, sections = vmtk_compute_centerline_sections(new_surface_capped, line)

        # Get arrays
        a = get_point_data_array("CenterlineSectionArea", area)
        n = get_point_data_array("FrenetTangent", area, k=3)
        l = get_curvilinear_coordinate(area)

        # Compute 'derivative' of the area
        dAdX = (a[1:, 0] - a[:-1, 0]) / (l[1:] - l[:-1])

        # Stopping criteria
        stop_id = np.nonzero(dAdX < -50)[0][-1] + 3
        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        surface, clipped = vtk_clip_polydata(surface, plane)
        if surface.GetNumberOfCells() < clipped.GetNumberOfCells():
            surface, clipped = clipped, surface

        surface = attach_clipped_regions_to_surface(surface, clipped, center)

    # Landmark MV
    line = extract_single_line(centerlines, 0)
    line = extract_single_line(line, 0, start_id=100, end_id=line.GetNumberOfPoints() - 40)
    line = compute_splined_centerline(line, nknots=10, isline=True)
    area, sections = vmtk_compute_centerline_sections(new_surface_capped, line)

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

    write_polydata(surface, base_path + "_no_PVs_or_MV.vtp")
    # Compute distances
    distanceFilter = vmtkscripts.vmtkSurfaceDistance()
    distanceFilter.ReferenceSurface = new_surface
    distanceFilter.DistanceArrayName = "DistanceLAA"
    distanceFilter.Surface = surface
    distanceFilter.Execute()
    distances = distanceFilter.Surface

    LAA_array = distances.GetPointData().GetArray("DistanceLAA")
    surface.GetPointData().AddArray(LAA_array)

    LAA_center = np.asarray([-2, 35, -10])  # laa_surface.GetCenter()
    clipped, surface = clipp_surface(surface, "DistanceLAA")
    # clipped, surface = clipp_surface(surface, "DistanceLAA")
    # TODO: Find LAA center without providing LAA
    # LAA = get_connectivity(clipped,
    write_polydata(clipped, base_path + "_clipped.vtp")
    write_polydata(surface, base_path + "_surface_clipped.vtp")
    surface = attach_clipped_regions_to_surface(surface, clipped, LAA_center)

    # Store LAA
    connectivity = vtk_compute_connectivity(clipped)
    region_id = get_point_data_array("RegionId", connectivity)
    LAA = [None]
    LAA_surface_area = 1E-16
    for i in range(int(region_id.max()) + 1):
        region = vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0)
        region_surface_area = vtk_compute_mass_properties(region)
        if region_surface_area > LAA_surface_area:
            LAA_surface_area = region_surface_area
            LAA[0] = region

    LAA_uncapped = LAA[0]
    LAA_capped = vmtk_cap_polydata(LAA[0])
    write_polydata(LAA_uncapped, input_path.replace(".vtp", "_LAA.vtp"))

    tol = 1E-12
    LAA_orifice = []
    for i in range(LAA_capped.GetNumberOfPoints()):
        p = LAA_capped.GetPoint(i)
        p_u = LAA_uncapped.GetPoint(i)
        p = np.array(p)
        p_u = np.array(p_u)
        d = np.linalg.norm(p - p_u)

        if d > tol:
            LAA_orifice.append(p)

    if len(LAA_orifice) > 1:
        LAA_orifice = np.array(LAA_orifice)
        mean_LAA_orifice = np.average(LAA_orifice, axis=0)
    else:
        mean_LAA_orifice = LAA_orifice[0]

    # Clip LAA
    test_automated_clipping = False
    if test_automated_clipping:

        # Compute centerline to orifice from MV to get tangent
        line = extract_single_line(centerlines, 0)
        p_mv = np.array(line.GetPoint(0))
        laa_centerlines, _, _ = compute_centerlines(p_mv.tolist(), mean_LAA_orifice.tolist(), laa_centerline_path,
                                                    capped_surface, resampling=0.1,
                                                    smooth=False, pole_ids=pole_ids,
                                                    voronoi=voronoi, base_path=base_path)
        laa_centerlines = vmtk_compute_geometric_features(laa_centerlines, False)

        n_laa = get_point_data_array("FrenetTangent", laa_centerlines, k=3)
        l_laa = np.linspace(0, 1, len(n_laa))
        n_x = n_laa[:, 0]
        diff = (n_x[1:] - n_x[:-1]) / (l_laa[1:] - l_laa[:-1])

        skip = int(len(diff) * 0.01)  # Ignore part close to LAA edge
        stop_id = np.nonzero(abs(diff[skip:]) > 50)[0][0] + 50
        normal = n_laa[stop_id]
        center = mean_LAA_orifice

        plane = vtk_plane(center, normal)
        new_surface, clipped = vtk_clip_polydata(surface_with_laa, plane)
        write_polydata(new_surface, "newsurf.vtp")
        write_polydata(clipped, "newsurf_clip.vtp")
        if new_surface.GetNumberOfCells() < clipped.GetNumberOfCells():
            new_surface, clipped = clipped, new_surface
        new_surface = attach_clipped_regions_to_surface(new_surface, clipped, center)

        write_polydata(new_surface, "newsurf_att.vtp")
        write_polydata(clipped, "newsurf_clip_att.vtp")

    else:
        laa_point = provide_region_points(surface_with_laa, None, None)
        la_point = laa_point[0]
        # Compute centerline to orifice from MV to get tangent
        line = extract_single_line(centerlines, 0)
        p_mv = np.array(line.GetPoint(0))
        laa_centerlines, _, _ = compute_centerlines(p_mv.tolist(), la_point, laa_centerline_path,
                                                    capped_surface, resampling=0.1,
                                                    smooth=False, pole_ids=pole_ids,
                                                    voronoi=voronoi, base_path=base_path)

        line = extract_single_line(laa_centerlines, 0, start_id=50, end_id=laa_centerlines.GetNumberOfPoints() - 50)
        laa_l = get_curvilinear_coordinate(area)
        step = 10 * np.mean(laa_l[1:] - laa_l[:-1])
        line = vmtk_resample_centerline(line, step)
        line = compute_splined_centerline(line, nknots=10, isline=True)
        area, sections = vmtk_compute_centerline_sections(surface_with_laa, line)
        write_polydata(sections, "tmp_area_LAA.vtp")

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

        np.savetxt("laa_a.txt", a)
        np.savetxt("laa_l.txt", l)
        np.savetxt("laa_dadx.txt", dAdX)

        # Clip the model
        plane = vtk_plane(center, normal)
        surface, clipped = vtk_clip_polydata(surface_with_laa, plane)
        if surface.GetNumberOfCells() > clipped.GetNumberOfCells():
            surface, clipped = clipped, surface

        surface = get_surface_closest_to_point(surface, center)

        write_polydata(surface, base_path + "_isolated_LAA.vtp")

        embed()

    # Compute distances
    old_surface = read_polydata(input_path)
    distanceFilter = vmtkscripts.vmtkSurfaceDistance()
    distanceFilter.ReferenceSurface = surface
    distanceFilter.DistanceArrayName = "Distance"
    distanceFilter.Surface = old_surface
    distanceFilter.Execute()
    distances = distanceFilter.Surface

    LA_array = distances.GetPointData().GetArray("Distance")

    LA_full_volume = []
    LA_volume = []
    LAA_volume = []
    PV_volume = []

    # Apply landmarking to all models
    for i in range(10):
        input_path = "landmark/0" + str(i) + "_registered.vtp"
        print(input_path)
        if path.exists(input_path):
            surface = read_polydata(input_path)
        else:
            surface = read_polydata(input_path.replace(".vtp", ".vtk"))
            if isinstance(surface, vtk.vtkUnstructuredGrid):
                surface = vmtk_compute_geometric_features(surface)
            write_polydata(surface, input_path)

        surface.GetPointData().AddArray(LA_array)
        surface = clean_surface(surface)
        surface = triangulate_surface(surface)

        # Store full LA as vtp
        LA_full_capped = vmtk_cap_surface(surface)
        clipped, surface = clipp_surface(surface, "Distance", value=1e-3)

        # Store LA
        LA_capped = capp_surface(surface)
        write_polydata(surface, input_path.replace(".vtp", "_only_LA.vtp"))

        # Store LAA
        connectivity = vtk_compute_connectivity(clipped)
        region_id = get_point_data_array("RegionId", connectivity)
        regions = []
        centers = []
        for i in range(int(region_id.max()) + 1):
            region = vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0)
            centers.append(region.GetCenter())
            regions.append(region)

        LAA_id = int(np.argmin(np.sqrt(np.sum((np.array(centers) - LAA_center) ** 2, axis=1))))
        LAA = regions.pop(LAA_id)
        LAA_capped = capp_surface(LAA)
        write_polydata(LAA, input_path.replace(".vtp", "_only_LAA.vtp"))

        # Store PVs
        PVs = merge_data(regions)
        PVs_capped = capp_surface(PVs)
        write_polydata(PVs, input_path.replace(".vtp", "_only_PVs.vtp"))

        LA_full_volume.append(compute_volume(LA_full_capped))
        LA_volume.append(compute_volume(LA_capped))
        LAA_volume.append(compute_volume(LAA_capped))
        PV_volume.append(compute_volume(PVs_capped))

    return LA_full_volume, LA_volume, LAA_volume, PV_volume


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

    for i in range(4):
        line = extract_single_line(centerlines, i)
        line = extract_single_line(line, 0, start_id=40, end_id=line.GetNumberOfPoints() - 100)
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

        stop_id = np.nonzero(dAdX < -50)[0][-1]  # - 3
        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        surface, clipped = vtk_clip_polydata(surface, plane)
        if surface.GetNumberOfCells() < clipped.GetNumberOfCells():
            surface, clipped = clipped, surface

        surface = attach_clipped_regions_to_surface(surface, clipped, center)

    # Landmark MV
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

    write_polydata(surface, la_and_laa_path)


def extract_LAA(input_path):
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

    p_laa = provide_region_points(surface, None, None)
    la_point = p_laa[0]

    # Compute centerline to orifice from MV to get tangent
    laa_centerlines, _, _ = compute_centerlines(p_outlet, la_point, laa_centerline_path, capped_surface,
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
    args = parser.parse_args()
    case_path = args.case

    t0 = time.time()
    extract_LA_and_LAA(case_path)
    t1 = time.time()
    extract_LAA(case_path)
    t2 = time.time()
    print("--- Extraction complete")
    scale = 1  # Get seconds
    print("Time spent extracting LA & LAA: {:.3f} s".format((t1 - t0) / scale))
    print("Time spent extracting LAA: {:.3f} s".format((t2 - t1) / scale))
