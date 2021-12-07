from morphman.common import *

import ImportData
from NetworkBoundaryConditions import FlowSplitting

try:
    from vmtkpointselector import *
except ImportError:
    pass
import numpy as np
from os import path, remove

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

    Args:
        surface:
        provided_points:
        dir_path:

    Returns:

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


def compute_centers_for_meshing(polyData, atrium_present, case_path=None, test_capped=False):
    """
    Compute the center of all the openings in the surface. The inlet is chosen based on
    the largest area for arteries (or aneurysm). However, for atrium, the outlet is chosen based on
    the largest area (new).

    Args:
        test_capped (bool): Check if surface is capped.
        polyData (vtkPolyData): centers of the openings.
        case_path (str): path to case directory.
        atrium_present (bool): Check if it is an atrium model.

    Returns:
        inlet_center (list): Inlet center.
        outlet_centers (list): A flattened list with all the outlet centers.
    """
    # Get cells which are open
    cells = vtk_extract_feature_edges(polyData)

    if cells.GetNumberOfCells() == 0 and not test_capped:
        print("WARNING: The model is capped, so it is uncapped, but the method is experimental.")
        uncapped_surface = get_uncapped_surface(polyData)
        compute_centers_for_meshing(uncapped_surface, atrium_present, case_path, test_capped)
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
    if atrium_present:
        # Store the center and area
        boundary_name = "outlet"
        boundaries_name = "inlet%d"
    else:
        boundary_name = "inlet"
        boundaries_name = "outlet%d"

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


def get_centers_for_meshing(surface, atrium_present, dir_path, flowext=False):
    """
    Get the centers of the inlet and outlets.

    Args:
        surface (vtkPolyData): An open surface.
        dir_path (str): Path to the case file.
        flowext (bool): Turn on/off flow extension.

    Returns:
        inlet (list): A flatt list with the point of the inlet
        outlet (list): A flatt list with the points of all the outlets.
    """

    # Check if info exists
    if flowext or not path.isfile(path.join(dir_path, dir_path + ".json")):
        compute_centers_for_meshing(surface, atrium_present, dir_path)

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

    if atrium_present and num_inlets != 0:
        inlets = []
        for i in range(num_inlets):
            inlets += parameters["inlet%d" % i]
    if not atrium_present and num_outlets != 0:
        outlets = []
        for i in range(num_outlets):
            outlets += parameters["outlet%d" % i]

    # FIXIT: atrium case has several inlets (instead of inlet) and only one outlet (instead of outlets).
    if inlets == [] and outlets == []:
        inlets, outlets = compute_centers_for_meshing(surface, atrium_present, dir_path)

        print("The number of outlets =", len(outlets) // 3)
        print("The number of inlets =", len(inlets) // 3)
        print()

    return inlets, outlets


def dist_sphere_curv(surface, centerlines, sac_center, misr_max, fileName, factor):
    """

    Args:
        surface:
        centerlines:
        sac_center:
        misr_max:
        fileName:
        factor:

    Returns:

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
    for i in range(len(sac_center)):
        distance_to_sphere = compute_distance_to_sphere(distance_to_sphere, sac_center[i],
                                                        distanceScale=0.2 / (misr_max[i] * 2.5))

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

    write_polydata(curvatureSurface, fileName)

    return distance_to_sphere


def dist_sphere_constant(surface, centerlines, sac_center, misr_max, fileName, edge_length):
    """

    Args:
        surface:
        centerlines:
        sac_center:
        misr_max:
        fileName:
        edge_length:

    Returns:

    """
    # Constant meshing method with possible refined area.
    # --- Compute the distanceToCenterlines
    distTocenterlines = vmtkscripts.vmtkDistanceToCenterlines()
    distTocenterlines.Surface = surface
    distTocenterlines.Centerlines = centerlines
    distTocenterlines.Execute()
    distance_to_sphere = distTocenterlines.Surface

    # Reduce element size in region
    for i in range(len(sac_center)):
        distance_to_sphere = compute_distance_to_sphere(distance_to_sphere,
                                                        sac_center[i],
                                                        minDistance=edge_length / 3,
                                                        maxDistance=edge_length,
                                                        distanceScale=edge_length * 3 / 4 / (misr_max[i] * 2.))

    element_size = edge_length + np.zeros((surface.GetNumberOfPoints(), 1))
    if len(sac_center) != 0:
        distance_to_spheres_array = get_point_data_array("DistanceToSpheres", distance_to_sphere)
        element_size = np.minimum(element_size, distance_to_spheres_array)

    vtk_array = create_vtk_array(element_size, "Size")
    distance_to_sphere.GetPointData().AddArray(vtk_array)
    write_polydata(distance_to_sphere, fileName)

    return distance_to_sphere


def dist_sphere_diam(surface, centerlines, sac_center, misr_max, fileName, factor):
    """

    Args:
        surface:
        centerlines:
        sac_center:
        misr_max:
        fileName:
        factor:

    Returns:

    """
    # Meshing method following Owais way.
    # --- Compute the distanceToCenterlines
    distTocenterlines = vmtkscripts.vmtkDistanceToCenterlines()
    distTocenterlines.Surface = surface
    distTocenterlines.Centerlines = centerlines
    distTocenterlines.Execute()
    distance_to_sphere = distTocenterlines.Surface

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
    # element_size[element_size < 0.12] = 0.12

    # Reduce element size in aneurysm
    for i in range(len(sac_center)):
        distance_to_sphere = compute_distance_to_sphere(distance_to_sphere,
                                                        sac_center[i],
                                                        maxDistance=100,
                                                        distanceScale=0.2 / (misr_max[i] * 2.))
    if len(sac_center) == 0:
        element_size *= factor
    else:
        distance_to_spheres_array = get_point_data_array("DistanceToSpheres", distance_to_sphere)
        element_size = np.minimum(element_size, distance_to_spheres_array) * factor

    vtk_array = create_vtk_array(element_size, "Size")
    distance_to_sphere.GetPointData().AddArray(vtk_array)
    write_polydata(distance_to_sphere, fileName)

    return distance_to_sphere


def mesh_alternative(surface):
    """

    Args:
        surface:

    Returns:

    """
    print("--- Meshing failed.")
    print("--- Proceeding with surface smooting and meshing.")
    surface = vmtk_smooth_surface(surface, "laplace", iterations=500)

    subdiv = vmtkscripts.vmtkSurfaceSubdivision()
    subdiv.Surface = surface
    subdiv.Method = "butterfly"
    subdiv.Execute()
    surface = subdiv.Surface

    return vmtk_smooth_surface(surface, "laplace", iterations=500)


def compute_distance_to_sphere(surface, centerSphere, radiusSphere=0.0,
                               distanceOffset=0.0, distanceScale=0.01,
                               minDistance=0.2, maxDistance=0.3,
                               distanceToSpheresArrayName="DistanceToSpheres"):
    """

    Args:
        surface:
        centerSphere:
        radiusSphere:
        distanceOffset:
        distanceScale:
        minDistance:
        maxDistance:
        distanceToSpheresArrayName:

    Returns:

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
        newDist = get_distance(centerSphere, surface.GetPoints().GetPoint(i)) - radiusSphere

        # Set offset and scale distance
        newDist = distanceOffset + newDist * distanceScale

        # Capp to min distance
        if newDist < minDistance:
            newDist = minDistance

        # Capp to max distance
        if newDist > maxDistance:
            newDist = maxDistance

        # Keep smallest distance
        newDist = min(newDist, distanceToSphere) if not add else newDist

        dist_array.SetComponent(i, 0, newDist)

    return surface


def generate_mesh(surface):
    """

    Args:
        surface:

    Returns:

    """
    # Compute the mesh.
    meshGenerator = vmtkscripts.vmtkMeshGenerator()
    meshGenerator.Surface = surface
    meshGenerator.ElementSizeMode = "edgelengtharray"  # Variable size mesh
    meshGenerator.TargetEdgeLengthArrayName = "Size"  # Variable size mesh
    meshGenerator.BoundaryLayer = 1
    meshGenerator.NumberOfSubLayers = 4
    meshGenerator.BoundaryLayerOnCaps = 0
    meshGenerator.BoundaryLayerThicknessFactor = 0.85
    meshGenerator.SubLayerRatio = 0.75
    meshGenerator.Tetrahedralize = 1
    meshGenerator.VolumeElementScaleFactor = 0.8
    meshGenerator.EndcapsEdgeLengthFactor = 1.0

    # Mesh
    meshGenerator.Execute()

    # Remeshed surface, store for later
    remeshSurface = meshGenerator.RemeshedSurface

    # Full mesh
    mesh = meshGenerator.Mesh

    return mesh, remeshSurface


def find_boundaries(model_path, mean_inflow_rate, network, polyDataVolMesh, verbose_print, is_atrium):
    """

    Args:
        model_path:
        mean_inflow_rate:
        network:
        polyDataVolMesh:
        verbose_print:
        is_atrium:

    Returns:

    """
    # Extract the surface mesh of the wall
    wallMesh = vtk_compute_threshold(polyDataVolMesh, "CellEntityIds", lower=0.5, upper=1.5)
    boundaryReferenceSystems = vmtkscripts.vmtkBoundaryReferenceSystems()
    boundaryReferenceSystems.Surface = wallMesh
    boundaryReferenceSystems.Execute()
    refSystem = boundaryReferenceSystems.ReferenceSystems
    cellEntityIdsArray = get_vtk_array('CellEntityIds', 0, refSystem.GetNumberOfPoints())
    refSystem.GetPointData().AddArray(cellEntityIdsArray)

    # Extract the surface mesh of the end caps
    boundarySurface = vtk_compute_threshold(polyDataVolMesh, "CellEntityIds", upper=1.5, threshold_type="upper")
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

    Args:
        centerlines:
        file_name_probe_points:
        region_center:
        verbose_print:

    Returns:

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


def compute_flow_rate(atrium_present, inlet, parameters):
    """

    Args:
        atrium_present:
        inlet:
        parameters:

    Returns:

    """
    # FIXME: Add plausible boundary conditions for atrial flow
    flow_rate_factor = 0.27
    if atrium_present:
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

    Args:
        compress_mesh:
        file_name_surface_name:
        file_name_vtu_mesh:
        file_name_xml_mesh:
        mesh:
        remeshed_surface:

    Returns:

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
    polyDataVolMesh = mesh
    return polyDataVolMesh


def add_flow_extension(surface, centerlines, include_outlet, extension_length=2.0,
                       extension_mode="boundarynormal"):
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
    flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(0)
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
