from morphman.common import *
from vtk.numpy_interface import dataset_adapter as dsa

import ImportData
from NetworkBoundaryConditions import FlowSplitting

try:
    from vmtkpointselector import *
except ImportError:
    pass
import numpy as np
from os import path

# Global array names
radiusArrayName = 'MaximumInscribedSphereRadius'
parallelTransportNormalsArrayName = 'ParallelTransportNormals'
groupIDsArrayName = "GroupIds"
abscissasArrayName = 'Abscissas'
clippingArrayName = 'ClippingArray'
branchClippingArrayName = 'BranchClippingArray'
distanceToTubeArrayName = 'DistanceToTubeFunction'
closedArrayName = 'ClosedSection'
eikonalSolutionArrayName = 'EikonalSolutionArray'
edgeArrayName = 'EdgeArray'
edgePCoordArrayName = 'EdgePCoordArray'
costFunctionArrayName = 'CostFunctionArray'
distanceToSpheresArrayName = "DistanceToSpheres"
cellEntityArrayName = "CellEntityIds"

# Options not available from commandline
divergingRatioToSpacingTolerance = 2.0
interpolationHalfSize = 3
voronoiCoreCutOffThreshold = 0.75
numberOfSplineAnalyzedPoints = 40
phiValues = [float(i) for i in range(2, 43, 2)]
thetaStep = 2.0

# Shortcuts
version = vtk.vtkVersion().GetVTKMajorVersion()


def get_regions_to_refine(surface, provided_points, dir_path):
    """
    Determines which regions to refine.

    Args:
        surface (vtkPoylData): Surface model
        provided_points (ndarray): Points defining regions to refine
        dir_path (str): Path to save directory

    Returns:
        region_points (list): List of points representing regions to refine
    """
    # Check if info exists
    if not path.isfile(path.join(dir_path, dir_path + ".txt")):
        provide_region_points(surface, provided_points, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    dome = []
    for key, value in parameters.items():
        if key.startswith("region_"):
            dome.append(value)

    if not dome:
        dome = provide_region_points(surface, provided_points, dir_path)

    # Flatten list
    return [item for sublist in dome for item in sublist]


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


def make_voronoi_diagram(surface, file_path):
    """
    Compute the voronoi diagram of surface model.

    Args:
        surface (polydata): Capped surface model to create a Voronoi diagram of.
        file_path (str): Absolute path to surface model path.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    if path.isfile(file_path):
        return read_polydata(file_path)

    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = 0
    voronoi.Execute()

    write_polydata(voronoi.VoronoiDiagram, file_path)

    return voronoi.VoronoiDiagram


def compute_centers_for_meshing(surface, is_atrium, case_path=None, test_capped=False):
    """
    Compute the center of all the openings in the surface. The inlet is chosen based on
    the largest area for arteries (or aneurysm). However, for atrium, the outlet is chosen based on
    the largest area (new).

    Args:
        test_capped (bool): Check if surface is capped.
        surface (vtkPolyData): centers of the openings.
        case_path (str): path to case directory.
        is_atrium (bool): Check if it is an atrium model.

    Returns:
        inlet_center (list): Inlet center.
        outlet_centers (list): A flattened list with all the outlet centers.
    """
    # Get cells which are open
    cells = vtk_extract_feature_edges(surface)

    if cells.GetNumberOfCells() == 0 and not test_capped:
        print("WARNING: The model is capped, so it is uncapped, but the method is experimental.")
        uncapped_surface = get_uncapped_surface(surface)
        compute_centers_for_meshing(uncapped_surface, is_atrium, case_path, test_capped)
    elif cells.GetNumberOfCells() == 0 and test_capped:
        return False, 0

    # Compute connectivity of the cells
    outputs = vtk_compute_connectivity(cells)

    # Get connectivity array
    region_array = get_point_data_array("RegionId", outputs)

    if test_capped:
        return region_array.max() >= 1, region_array.max()

    # Get points
    points = []
    get_point = outputs.GetPoints().GetPoint
    for i in range(region_array.shape[0]):
        points.append(get_point(i))
    points = np.asarray(points)

    # Get area and center
    area = []
    center = []
    for i in range(int(region_array.max()) + 1):
        # Compute area
        boundary = vtk_compute_threshold(outputs, "RegionId", lower=i - 0.1, upper=i + 0.1, threshold_type="between",
                                         source=0)

        delaunay_filter = vtk.vtkDelaunay2D()
        delaunay_filter.SetInputData(boundary)
        delaunay_filter.Update()
        area.append(vtk_compute_mass_properties(delaunay_filter.GetOutput()))

        # Get center
        center.append(np.mean(points[(region_array == i).nonzero()[0]], axis=0))

    # Assume multiple inlets for atrium, and multiple outlets for arteries
    boundary_name = "outlet" if is_atrium else "inlet"
    boundaries_name = "inlet%d" if is_atrium else "outlet%d"

    boundary_area_name = boundary_name + "_area"
    boundaries_area_name = boundaries_name + "_area"
    boundary_id = area.index(max(area))
    if case_path is not None:
        info = {boundary_name: center[boundary_id].tolist(), boundary_area_name: area[boundary_id]}
        p = 0
        for i in range(len(area)):
            if i == boundary_id:
                p = -1
                continue
            info[boundaries_name % (i + p)] = center[i].tolist()
            info[boundaries_area_name % (i + p)] = area[i]

        write_parameters(info, case_path)

    boundary_center = center[boundary_id].tolist()  # center of the outlet
    center.pop(boundary_id)

    center_ = [item for sublist in center for item in sublist]  # centers of the inlets

    return center_, boundary_center


def get_centers_for_meshing(surface, is_atrium, dir_path, use_flow_extensions=False):
    """
    Get the centers of the inlet and outlets.

    Args:
        surface (vtkPolyData): An open surface.
        is_atrium (bool): True if model represents atrium
        dir_path (str): Path to the case file.
        use_flow_extensions (bool): Turn on/off flow extension.

    Returns:
        inlet (list): A flatt list with the point of the inlet
        outlet (list): A flatt list with the points of all the outlets.
    """

    # Check if info exists
    if use_flow_extensions or not path.isfile(path.join(dir_path, dir_path + ".json")):
        compute_centers_for_meshing(surface, is_atrium, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    outlets = []
    inlets = []
    for key, value in parameters.items():
        if "area" not in key and "id" not in key and "relevant" not in key:
            if "outlet" in key:
                outlets += value
            elif "inlet" in key:
                inlets += value

    num_outlets = len(outlets) // 3
    num_inlets = len(inlets) // 3

    if is_atrium and num_inlets != 0:
        inlets = []
        for i in range(num_inlets):
            inlets += parameters["inlet%d" % i]
    if not is_atrium and num_outlets != 0:
        outlets = []
        for i in range(num_outlets):
            outlets += parameters["outlet%d" % i]

    # FIXIT: atrium case has several inlets (instead of inlet) and only one outlet (instead of outlets).
    if inlets == [] and outlets == []:
        inlets, outlets = compute_centers_for_meshing(surface, is_atrium, dir_path)

        print("The number of outlets =", len(outlets) // 3)
        print("The number of inlets =", len(inlets) // 3)
        print()

    return inlets, outlets


def dist_sphere_curvature(surface, centerlines, region_center, misr_max, save_path, factor):
    """
    Determines the target edge length for each cell on the surface, including
    potential refinement or coarsening of certain user specified areas.
    Level of refinement/coarseness is determined based on surface curvature.

    Args:
        surface (vtkPolyData): Input surface model
        centerlines (vtkPolyData): Centerlines of input model
        region_center (list): Point representing region to refine
        misr_max (list): Maximum inscribed sphere radius in region of refinement
        save_path (str): Location to store processed surface
        factor (float): Coarsening factor, determining the level of refinement (<1) or coarsening (>1)

    Returns:
        surface (vtkPolyData): Processed surface model with info on cell specific target edge length
    """

    # Get longest centerline
    length = []
    for i in range(centerlines.GetNumberOfLines()):
        line = extract_single_line(centerlines, i)
    length.append(get_curvilinear_coordinate(line)[-1])
    ind_longest = length.index(max(length))

    # Get all bifurcations along the longest centerline
    bif = []
    bif_id = []
    longest_line = extract_single_line(centerlines, ind_longest)
    tol = get_centerline_tolerance(centerlines)

    for i in range(centerlines.GetNumberOfLines()):
        if i == ind_longest: continue
        comp_line = extract_single_line(centerlines, i)
        for j in range(comp_line.GetNumberOfPoints()):
            pnt1 = longest_line.GetPoints().GetPoint(j)
            pnt2 = comp_line.GetPoints().GetPoint(j)
            if get_distance(pnt1, pnt2) > tol:
                bif.append(pnt1)
                bif_id.append(j)
                break

    # Remove bifurcations detected twice
    pop = []
    for i in range(len(bif)):
        for j in range(i + 1, len(bif)):
            dist = get_distance(bif[i], bif[j])
            if dist < tol * 6:
                pop.append(j)

    for i in sorted(set(pop))[::-1]:
        bif.pop(i)
        bif_id.pop(i)

    # Set siphon to be the location of the first branch, assumes ICA as
    # inlet, and that the ophthalmic is segmented.
    bif_id.sort()
    siphon_point = line.GetPoints().GetPoint(bif_id[0])
    distance_to_sphere = compute_distance_to_sphere(surface, siphon_point)

    # Add the center of the sac
    for i in range(len(region_center)):
        distance_to_sphere = compute_distance_to_sphere(distance_to_sphere, region_center[i],
                                                        distance_scale=0.2 / (misr_max[i] * 2.5))

    # Compute curvature
    curvatureFilter = vmtkscripts.vmtkSurfaceCurvature()
    curvatureFilter.AbsoluteCurvature = 1
    curvatureFilter.MedianFiltering = 1
    curvatureFilter.CurvatureType = "gaussian"
    curvatureFilter.Offset = 0.15
    curvatureFilter.BoundedReciprocal = 1
    curvatureFilter.Surface = distance_to_sphere
    curvatureFilter.Execute()

    # Multiple the surface
    curvatureSurface = curvatureFilter.Surface
    curvatureArray = get_point_data_array("Curvature", curvatureSurface)
    distance_to_sphere_array = get_point_data_array("DistanceToSpheres", distance_to_sphere)
    size_array = curvatureArray * distance_to_sphere_array * factor

    size_vtk_array = create_vtk_array(size_array, "Size")
    curvatureSurface.GetPointData().AddArray(size_vtk_array)

    write_polydata(curvatureSurface, save_path)

    return distance_to_sphere


def dist_sphere_constant(surface, centerlines, region_center, misr_max, save_path, edge_length):
    """
    Determines the target edge length for each cell on the surface, including
    potential refinement or coarsening of certain user specified areas.
    Level of refinement/coarseness is determined based on user selected region, otherwise a constant target edge length
    is selected.

    Args:
        surface (vtkPolyData): Input surface model
        centerlines (vtkPolyData): Centerlines of input model
        region_center (list): Point representing region to refine
        misr_max (list): Maximum inscribed sphere radius in region of refinement
        save_path (str): Location to store processed surface
        edge_length (float): Target edge length

    Returns:
        surface (vtkPolyData): Processed surface model with info on cell specific target edge length
    """
    # Constant meshing method with possible refined area.
    # --- Compute the distanceToCenterlines
    distToCenterlines = vmtkscripts.vmtkDistanceToCenterlines()
    distToCenterlines.Surface = surface
    distToCenterlines.Centerlines = centerlines
    distToCenterlines.Execute()
    distance_to_sphere = distToCenterlines.Surface

    # Reduce element size in region
    coarsen = False
    if coarsen:
        F = 3 / 4  # For 3xMISIR
    else:
        F = 1
    use_laa = False
    if use_laa:
        # Test LAA specific refinement
        laa_path = save_path.rsplit("/", 1)[0] + "/LA015_laa.vtp"
        laa = read_polydata(laa_path)
        for i in range(len(region_center)):
            distance_to_sphere = compute_distance_to_surface(distance_to_sphere, laa, min_distance=edge_length / 3,
                                                         max_distance=edge_length)
    else:
        for i in range(len(region_center)):
            distance_to_sphere = compute_distance_to_sphere(distance_to_sphere,
                                                        region_center[i],
                                                        min_distance=edge_length / 3,
                                                        max_distance=edge_length,
                                                        distance_scale=F * edge_length * 3 / 4 / (misr_max[i] * 2.))

    element_size = edge_length + np.zeros((surface.GetNumberOfPoints(), 1))
    if len(region_center) != 0:
        distance_to_spheres_array = get_point_data_array("DistanceToSpheres", distance_to_sphere)
        if coarsen:
            element_size = np.maximum(element_size, distance_to_spheres_array)
        else:
            element_size = np.minimum(element_size, distance_to_spheres_array)

    vtk_array = create_vtk_array(element_size, "Size")
    distance_to_sphere.GetPointData().AddArray(vtk_array)
    write_polydata(distance_to_sphere, save_path)

    return distance_to_sphere


def dist_sphere_diam(surface, centerlines, region_center, misr_max, save_path, factor):
    """
    Determines the target edge length for each cell on the surface, including
    potential refinement or coarsening of certain user specified areas.
    Level of refinement/coarseness is determined based on the distance to the centerline.

    Args:
        surface (vtkPolyData): Input surface model
        centerlines (vtkPolyData): Centerlines of input model
        region_center (list): Point representing region to refine
        misr_max (list): Maximum inscribed sphere radius in region of refinement
        save_path (str): Location to store processed surface
        factor (float): Coarsening factor, determining the level of refinement (<1) or coarsening (>1)

    Returns:
        surface (vtkPolyData): Processed surface model with info on cell specific target edge length
    """
    # Meshing method following Owais way.
    # --- Compute the distanceToCenterlines
    distToCenterlines = vmtkscripts.vmtkDistanceToCenterlines()
    distToCenterlines.Surface = surface
    distToCenterlines.Centerlines = centerlines
    distToCenterlines.Execute()
    distance_to_sphere = distToCenterlines.Surface

    # Compute element size based on diameter
    upper = 20
    lower = 6
    diameter_array = 2 * get_point_data_array("DistanceToCenterlines", distance_to_sphere)
    element_size = 13. / 35 * diameter_array ** 2 + lower
    element_size[element_size > upper] = upper
    element_size[element_size < lower] = lower
    elements_vtk = create_vtk_array(element_size, "Num elements")
    distance_to_sphere.GetPointData().AddArray(elements_vtk)
    element_size = diameter_array / element_size

    # Reduce element size in aneurysm
    for i in range(len(region_center)):
        distance_to_sphere = compute_distance_to_sphere(distance_to_sphere,
                                                        region_center[i],
                                                        max_distance=100,
                                                        distance_scale=0.2 / (misr_max[i] * 2.))
    if len(region_center) == 0:
        element_size *= factor
    else:
        distance_to_spheres_array = get_point_data_array("DistanceToSpheres", distance_to_sphere)
        element_size = np.minimum(element_size, distance_to_spheres_array) * factor

    vtk_array = create_vtk_array(element_size, "Size")
    distance_to_sphere.GetPointData().AddArray(vtk_array)
    write_polydata(distance_to_sphere, save_path)

    return distance_to_sphere


def mesh_alternative(surface):
    """
    Subdivides and smoothes the input surface model, preparing it for volumetric meshing.

    Args:
        surface (vtkPolyData): Input surface model to be meshed alternatively

    Returns:
        surface (vtkPolyData): Smoothed model
    """
    print("--- Meshing failed.")
    print("--- Proceeding with surface smooting and meshing.")
    surface = vmtk_smooth_surface(surface, "laplace", iterations=500)

    surfaceSubdivision = vmtkscripts.vmtkSurfaceSubdivision()
    surfaceSubdivision.Surface = surface
    surfaceSubdivision.Method = "butterfly"
    surfaceSubdivision.Execute()
    surface = surfaceSubdivision.Surface

    return vmtk_smooth_surface(surface, "laplace", iterations=500)


def compute_distance_to_sphere(surface, center_sphere, radius_sphere=0.0, distance_offset=0.0, distance_scale=0.01,
                               min_distance=0.2, max_distance=0.3):
    """
    Computes cell specific target edge length (distances) based on input criterion.

    Args:
        surface (vtkPolyData): Input surface model
        center_sphere (list): Point representing region to refine
        radius_sphere (list): Maximum inscribed sphere radius in region of refinement
        distance_offset (float): Offsets the relevant region by a offset
        distance_scale (float): Determines how the distance is scaled based on the distance from the relevant region
        min_distance (float): Minimum distance away from the relevant region before scaling starts
        max_distance (float): Maximum distance away from the relevant region before scaling stops

    Returns:
        surface (vtkPolyData): Modified surface model with distances
    """
    # Check if there allready exists a distance to spheres
    N = surface.GetNumberOfPoints()
    number, names = get_number_of_arrays(surface)
    add = False
    if distanceToSpheresArrayName not in names:
        add = True

    # Get array
    if add:
        dist_array = get_vtk_array(distanceToSpheresArrayName, 1, N)
        surface.GetPointData().AddArray(dist_array)
    else:
        dist_array = surface.GetPointData().GetArray("DistanceToSpheres")

    # Compute distance
    for i in range(N):
        distanceToSphere = dist_array.GetComponent(i, 0)

        # Get distance, but factor in size of sphere
        newDist = get_distance(center_sphere, surface.GetPoints().GetPoint(i)) - radius_sphere

        # Set offset and scale distance
        newDist = distance_offset + newDist * distance_scale

        # Capp to min distance
        if newDist < min_distance:
            newDist = min_distance

        # Capp to max distance
        if newDist > max_distance:
            newDist = max_distance

        # Keep smallest distance
        newDist = min(newDist, distanceToSphere) if not add else newDist

        # TODO: Edit if coarsening
        coarsen = False
        if coarsen:
            s = 4.9  # 2 and 3 works
            # Coarsening:
            dist_array.SetComponent(i, 0, s * 1.9 - (s - 1) * newDist)
        else:
            dist_array.SetComponent(i, 0, newDist)

    return surface


def compute_distance_to_surface(surface, subsurface, min_distance=0.2, max_distance=0.3):
    """
    Computes cell specific target edge length (distances) based on input criterion.

    Args:
        surface (vtkPolyData): Input surface model
        subsurface (vtkPolyData): Surface model to compare with
        min_distance (float): Minimum distance away from the relevant region before scaling starts
        max_distance (float): Maximum distance away from the relevant region before scaling stops

    Returns:
        surface (vtkPolyData): Modified surface model with distances
    """
    # Check if there allready exists a distance to spheres
    N = surface.GetNumberOfPoints()
    number, names = get_number_of_arrays(surface)
    add = False
    if distanceToSpheresArrayName not in names:
        add = True

    # Get array
    if add:
        dist_array = get_vtk_array(distanceToSpheresArrayName, 1, N)
        surface.GetPointData().AddArray(dist_array)
    else:
        dist_array = surface.GetPointData().GetArray("DistanceToSpheres")

    # Define locators
    subsurface_locator = get_vtk_point_locator(subsurface)

    # Compute distance
    distances = []
    for i in range(N):
        distanceToSphere = dist_array.GetComponent(i, 0)

        p_surface = np.array(surface.GetPoint(i))
        id_subsurface = subsurface_locator.FindClosestPoint(p_surface)
        p_subsurface = subsurface.GetPoint(id_subsurface)
        newDist = get_distance(p_surface, p_subsurface)

        tolDist = 1.5
        distances.append(newDist)

        if newDist < tolDist:
            newDist = min_distance
        else:
            dD = ((newDist - tolDist) / newDist) ** (3 / 2)
            newDist = max(min_distance, dD * max_distance)

        # Keep smallest distance
        newDist = min(newDist, distanceToSphere) if not add else newDist

        dist_array.SetComponent(i, 0, newDist)

    return surface


def generate_mesh(surface, add_boundary_layer):
    """
    Generates a mesh suitable for CFD from a input surface model.

    Args:
        surface (vtkPolyData): Surface model to be meshed.
        add_boundary_layer (bool): Adds boundary layer if true

    Returns:
        mesh (vtkUnstructuredGrid): Output mesh
        remeshedsurface (vtkPolyData): Remeshed version of the input model
    """
    # Compute the mesh.
    meshGenerator = vmtkscripts.vmtkMeshGenerator()
    meshGenerator.Surface = surface
    meshGenerator.ElementSizeMode = "edgelengtharray"  # Variable size mesh
    meshGenerator.TargetEdgeLengthArrayName = "Size"  # Variable size mesh#
    if add_boundary_layer:
        meshGenerator.BoundaryLayer = 1
        meshGenerator.BoundaryLayerOnCaps = 0
        # n_layers = [1,2,3,4,6,8,10]
        # meshGenerator.NumberOfSubLayers = n_layers[6]
        # normal, down down up up
        # Cases 1 2 3 4 5 6. nr 3 = default
        cases = [0.55, 0.7, 0.85, 1.0, 1.15, 1.3]
        layer_ratios = [0.25, 0.40, 0.55, 0.70, 0.85, 1]
        # meshGenerator.SubLayerRatio = layer_ratios[-3]
        # Do 0 - 2 - 4 - 10 layer study
        # Note: Defaults: 0.85, 4, 0.75
        meshGenerator.BoundaryLayerThicknessFactor = 0.85  # ID=2 = default
        meshGenerator.NumberOfSubLayers = 4
        meshGenerator.SubLayerRatio = 0.75

        meshGenerator.Tetrahedralize = 1
        meshGenerator.VolumeElementScaleFactor = 0.8
        meshGenerator.EndcapsEdgeLengthFactor = 1.0
    else:
        meshGenerator.BoundaryLayer = 0
        meshGenerator.BoundaryLayerOnCaps = 1

    # Mesh
    meshGenerator.Execute()

    # Remeshed surface, store for later
    remeshSurface = meshGenerator.RemeshedSurface

    # Full mesh
    mesh = meshGenerator.Mesh

    return mesh, remeshSurface


def find_boundaries(model_path, mean_inflow_rate, network, mesh, verbose_print, is_atrium):
    """
    Finds inlet and outlet boundary IDs after complete meshing, including
    mean flow ratio and area ratios between outlets or inlets (determined by type of model)

    Args:
        model_path (str): Path to model files
        mean_inflow_rate (float): Flow rate
        network (Network): Flow splitting network based on network boundary condition
        mesh (vtkUnstructuredGrid): Volumetric mesh
        verbose_print (bool): Prints additional info if True
        is_atrium (bool): Determines if model represents atrium or artery
    """
    # Extract the surface mesh of the wall
    wallMesh = vtk_compute_threshold(mesh, "CellEntityIds", lower=0.5, upper=1.5)
    boundaryReferenceSystems = vmtkscripts.vmtkBoundaryReferenceSystems()
    boundaryReferenceSystems.Surface = wallMesh
    boundaryReferenceSystems.Execute()
    refSystem = boundaryReferenceSystems.ReferenceSystems
    cellEntityIdsArray = get_vtk_array('CellEntityIds', 0, refSystem.GetNumberOfPoints())
    refSystem.GetPointData().AddArray(cellEntityIdsArray)

    # Extract the surface mesh of the end caps
    boundarySurface = vtk_compute_threshold(mesh, "CellEntityIds", upper=1.5, threshold_type="upper")
    pointCells = vtk.vtkIdList()
    surfaceCellEntityIdsArray = vtk.vtkIntArray()
    surfaceCellEntityIdsArray.DeepCopy(boundarySurface.GetCellData().GetArray('CellEntityIds'))

    # Find the corresponding couple (mesh outlet ID, network ID).
    ids = []
    for i in range(refSystem.GetNumberOfPoints()):
        distancePoints = 10000000
        pointId = boundarySurface.FindPoint(refSystem.GetPoint(i))
        boundarySurface.GetPointCells(pointId, pointCells)
        cellId = pointCells.GetId(0)
        cellEntityId = surfaceCellEntityIdsArray.GetValue(cellId)
        cellEntityIdsArray.SetValue(i, cellEntityId)

        meshPoint = refSystem.GetPoint(i)
        for element in network.elements:
            if element.IsAnOutlet():
                networkPoint = element.GetOutPointsx1()[0]
            if element.IsAnInlet():
                networkPoint = element.GetInPointsx0()[0]
            if vtk.vtkMath.Distance2BetweenPoints(meshPoint, networkPoint) < distancePoints:
                distancePoints = vtk.vtkMath.Distance2BetweenPoints(meshPoint, networkPoint)
                closest = element.GetId()
        if network.elements[closest].IsAnInlet():
            verbose_print('I am the inlet, Sup?')
            verbose_print(network.elements[closest].GetInPointsx0()[0])
            ids.insert(0, [cellEntityId, mean_inflow_rate])
        else:
            beta = network.elements[closest].GetBeta()
            ids.append([cellEntityId, beta])
            verbose_print(beta)
            verbose_print(network.elements[closest].GetOutPointsx1()[0])
        verbose_print('CellEntityId: %d\n' % cellEntityId)
        verbose_print('meshPoint: %f, %f, %f\n' % (meshPoint[0], meshPoint[1], meshPoint[2]))
        verbose_print(ids)

    # Store information for the solver.
    inlet_id = [ids[0][0] - 1]
    outlet_ids = []
    area_ratios = []
    for k in range(1, refSystem.GetNumberOfPoints()):
        outlet_ids.append(ids[k][0] - 1)
        area_ratios.append(ids[k][1])

    info = {
        "mean_flow_rate": mean_inflow_rate,
        "area_ratio": area_ratios
    }

    # Swap outlet with inlet if meshing atrium model
    if is_atrium:
        info['inlet_ids'] = outlet_ids
        info['outlet_id'] = inlet_id
    else:
        info['inlet_id'] = inlet_id
        info['outlet_ids'] = outlet_ids

    write_parameters(info, model_path)


def setup_model_network(centerlines, file_name_probe_points, region_center, verbose_print):
    """
    Sets up network used for network boundary condition model.

    Args:
        centerlines (vtkPolyData): Centerlines representing meshed model
        file_name_probe_points (str): Save path of probe points
        region_center (list): List of points representing region of refinement
        verbose_print (bool): Prints additional info if True

    Returns:
        network (Network): Network model
        probe_points (ndarray): Probe points where velocity and pressure is to be sampled
    """
    # Set the network object used in the scripts for
    # boundary conditions and probes.
    network = ImportData.Network()
    centerlinesBranches = ImportData.SetNetworkStructure(centerlines, network, verbose_print)

    if not path.isfile(file_name_probe_points):
        # Get the list of coordinates for the probe points along the network centerline.
        listProbePoints = ImportData.GetListProbePoints(centerlinesBranches, network, verbose_print)
        listProbePoints += region_center

        print("--- Saving probes points in: %s\n" % file_name_probe_points)
        probe_points = np.array(listProbePoints)
        probe_points.dump(file_name_probe_points)
    else:
        probe_points = np.load(file_name_probe_points, allow_pickle=True)

    # Set the flow split and inlet boundary condition
    # Compute the outlet boundary condition percentages.
    flowSplitting = FlowSplitting()
    flowSplitting.ComputeAlphas(network, verbose_print)
    flowSplitting.ComputeBetas(network, verbose_print)
    flowSplitting.CheckTotalFlowRate(network, verbose_print)

    return network, probe_points


def compute_flow_rate(is_atrium, inlet, parameters):
    """
    Computes mean flow rate used as boundary condition for pulsatile flow condition

    Args:
        is_atrium (bool): Determines if model is atrium or artery
        inlet (list): List of points representing midpoint of boundaries
        parameters (dict): Dictionary containing model parameters

    Returns:
        mean_inflow_rate (float): Mean inflow rate
    """
    # FIXME: Add plausible boundary conditions for atrial flow
    flow_rate_factor = 0.27
    if is_atrium:
        Total_inlet_area = 0
        num_inlets = len(inlet) // 3
        for i in range(num_inlets):
            Total_inlet_area += parameters["inlet%s_area" % i]
        mean_inflow_rate = flow_rate_factor * Total_inlet_area
    else:
        mean_inflow_rate = flow_rate_factor * parameters["inlet_area"]

    return mean_inflow_rate


def write_mesh(compress_mesh, file_name_surface_name, file_name_vtu_mesh, file_name_xml_mesh, mesh, remeshed_surface):
    """
    Writes the mesh to DOLFIN format, and compresses to .gz format

    Args:
        compress_mesh (bool): Compressed mesh to zipped format
        file_name_surface_name (str): Path to remeshed surface model
        file_name_vtu_mesh (str): Path to VTK mesh
        file_name_xml_mesh (str): Path to XML mesh
        mesh (vtuUnstructuredGrid): Meshed surface model
        remeshed_surface (vtkPolyData): Remeshed surface model
    """
    # Write mesh in VTU format
    write_polydata(remeshed_surface, file_name_surface_name)
    write_polydata(mesh, file_name_vtu_mesh)

    # Write mesh to FEniCS to format
    meshWriter = vmtkscripts.vmtkMeshWriter()
    meshWriter.CellEntityIdsArrayName = "CellEntityIds"
    meshWriter.Mesh = mesh
    meshWriter.Mode = "ascii"
    meshWriter.Compressed = compress_mesh
    meshWriter.OutputFileName = file_name_xml_mesh
    meshWriter.Execute()


def add_flow_extension(surface, centerlines, include_outlet, extension_length=2.0,
                       extension_mode="boundarynormal"):
    """
    Adds flow extensions to either all inlets or all outlets with specified extension length.

    Args:
        surface (vtkPolyData): Surface model to extend
        centerlines (vtkPolyData): Centerlines in model
        include_outlet (bool): Determines if outlet should be included or not
        extension_length (float): Determines length of flow extensions. Factor is multiplied with MISR at relevant boundary
        extension_mode (str): Determines how extensions are place, either normal to boundary or following centerline direction

    Returns:
        surface_extended (vtkPolyData): Extended surface model
    """
    # Mimick behaviour of vmtkflowextensionfilter
    boundaryExtractor = vtkvmtk.vtkvmtkPolyDataBoundaryExtractor()
    boundaryExtractor.SetInputData(surface)
    boundaryExtractor.Update()
    boundaries = boundaryExtractor.GetOutput()

    # Find outlet
    lengths = []
    for i in range(boundaries.GetNumberOfCells()):
        lengths.append(get_curvilinear_coordinate(boundaries.GetCell(i))[-1])
    outlet_id = lengths.index(max(lengths))

    # Exclude outlet or inlets
    boundaryIds = vtk.vtkIdList()
    for i in range(centerlines.GetNumberOfLines() + 1):
        if include_outlet and i == outlet_id:
            boundaryIds.InsertNextId(i)
        if not include_outlet and i != outlet_id:
            boundaryIds.InsertNextId(i)

    flowExtensionsFilter = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
    flowExtensionsFilter.SetInputData(surface)
    flowExtensionsFilter.SetCenterlines(centerlines)
    flowExtensionsFilter.SetAdaptiveExtensionLength(1)
    flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(1)
    flowExtensionsFilter.SetExtensionRatio(extension_length)
    flowExtensionsFilter.SetTransitionRatio(1.0)
    flowExtensionsFilter.SetCenterlineNormalEstimationDistanceRatio(1.0)
    if extension_mode == "centerlinedirection":
        flowExtensionsFilter.SetExtensionModeToUseCenterlineDirection()
    if extension_mode == "boundarynormal":
        flowExtensionsFilter.SetExtensionModeToUseNormalToBoundary()
    flowExtensionsFilter.SetInterpolationModeToThinPlateSpline()
    flowExtensionsFilter.SetBoundaryIds(boundaryIds)
    flowExtensionsFilter.SetSigma(1.0)
    flowExtensionsFilter.Update()

    surface_extended = flowExtensionsFilter.GetOutput()

    return surface_extended

def tmp():
    search_id = 0
    f2_factors = [0.95, 1.25, 1.42, 1.11]
    for i in range(centerlines.GetNumberOfLines()):
        if include_outlet and i == outlet_id:
            boundaryIds = vtk.vtkIdList()
            boundaryIds.InsertNextId(i)
        if not include_outlet:
            boundaryIds = vtk.vtkIdList()
            if i == 3:
                boundaryIds.InsertNextId(1)
            else:
                boundaryIds.InsertNextId(0)
            centerline = extract_single_line(centerlines, i)
            r = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(0 + search_id)
            search_id += centerline.GetNumberOfPoints()
            factor = extension_length / r * f2_factors[i]  # in [mm]

            flowExtensionsFilter = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
            flowExtensionsFilter.SetInputData(surface)
            flowExtensionsFilter.SetCenterlines(centerline)
            flowExtensionsFilter.SetAdaptiveExtensionLength(1)
            flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(1)
            flowExtensionsFilter.SetExtensionRatio(factor)
            flowExtensionsFilter.SetTransitionRatio(1.0)
            flowExtensionsFilter.SetCenterlineNormalEstimationDistanceRatio(1.0)
            if extension_mode == "centerlinedirection":
                flowExtensionsFilter.SetExtensionModeToUseCenterlineDirection()
            if extension_mode == "boundarynormal":
                flowExtensionsFilter.SetExtensionModeToUseNormalToBoundary()
            flowExtensionsFilter.SetInterpolationModeToThinPlateSpline()
            flowExtensionsFilter.SetBoundaryIds(boundaryIds)
            flowExtensionsFilter.SetSigma(1.0)
            flowExtensionsFilter.Update()

            surface = flowExtensionsFilter.GetOutput()

    if include_outlet:
        flowExtensionsFilter = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
        flowExtensionsFilter.SetInputData(surface)
        flowExtensionsFilter.SetCenterlines(centerlines)
        flowExtensionsFilter.SetAdaptiveExtensionLength(1)
        flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(1)
        flowExtensionsFilter.SetExtensionRatio(extension_length)
        flowExtensionsFilter.SetTransitionRatio(1.0)
        flowExtensionsFilter.SetCenterlineNormalEstimationDistanceRatio(1.0)
        if extension_mode == "centerlinedirection":
            flowExtensionsFilter.SetExtensionModeToUseCenterlineDirection()
        if extension_mode == "boundarynormal":
            flowExtensionsFilter.SetExtensionModeToUseNormalToBoundary()
        flowExtensionsFilter.SetInterpolationModeToThinPlateSpline()
        flowExtensionsFilter.SetBoundaryIds(boundaryIds)
        flowExtensionsFilter.SetSigma(1.0)
        flowExtensionsFilter.Update()

        surface = flowExtensionsFilter.GetOutput()

    return surface

def remesh_surface(surface, edge_length, exclude=None):
    surface = dsa.WrapDataObject(surface)
    if cellEntityArrayName not in surface.CellData.keys():
        surface.CellData.append(np.zeros(surface.VTKObject.GetNumberOfCells()) + 1, cellEntityArrayName)
    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
    remeshing.Surface = surface.VTKObject
    remeshing.CellEntityIdsArrayName = cellEntityArrayName
    remeshing.TargetEdgeLength = edge_length
    remeshing.MaxEdgeLength = 1e6
    remeshing.MinEdgeLength = 0.0
    remeshing.TargetEdgeLengthFactor = 1.0
    remeshing.TargetEdgeLengthArrayName = ""
    remeshing.TriangleSplitFactor = 5.0
    remeshing.ElementSizeMode = "edgelength"
    if exclude is not None:
        remeshing.ExcludeEntityIds = exclude

    remeshing.Execute()

    remeshed_surface = remeshing.Surface

    return remeshed_surface
