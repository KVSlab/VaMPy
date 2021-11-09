import argparse
import os
import shutil

from morphman.common import *
from vtk.numpy_interface import dataset_adapter as dsa

from common import get_centers_for_meshing, find_boundaries, setup_model_network, compute_flow_rate
from visualize import visualize

cell_id_name = "CellEntityIds"


def main(case_path, move_surface, add_extensions, edge_length, patient_specific, recompute_mesh, recompute_all,
         flow_extension_length):
    """
    Automatically generate movement and  mesh of surface model in .vtu and .xml format.
    Assumes the user either has a set of displaced models or models with mapped displacement fields,
    located in moved_path or mapped_path, respectively.

    Can add patient-specific or arbitrary movement to model.

    Args:
        case_path (str): Path to case
        move_surface (bool): To move surface or not
        add_extensions (bool): To add flow extensions or not
        edge_length (float): Mesh resolution, characteristic edge length
        patient_specific (bool): If case has patient-specific movement or not
        recompute_mesh (bool): Computes mesh if true
        recompute_all (bool): Computes everything if true
    """
    # Find model_path
    if "vtp" in case_path:
        model_path = case_path.replace(".vtp", "")
    elif "stl" in case_path:
        model_path = case_path.replace(".stl", "")

    cl_path = model_path + "_cl.vtp"
    case = model_path.split("/")[-1]
    mapped_path = model_path + "_mapped"
    moved_path = model_path + "_moved"
    mesh_path = model_path + ".vtu"
    mesh_xml_path = mesh_path.replace(".vtu", ".xml")

    if recompute_all:
        remove_preprocessing_files(case_path, model_path, moved_path)

    if not path.exists(moved_path):
        os.mkdir(moved_path)

    # Compute centerlines and get center of mitral valve as new origin
    surface = read_polydata(case_path)
    surface_to_mesh = None

    # Cap surface with flow extensions
    capped_surface = vmtk_cap_polydata(surface)
    outlets, inlet = get_centers_for_meshing(surface, True, model_path)
    centerlines, _, _ = compute_centerlines(inlet, outlets, cl_path, capped_surface, resampling=0.01)
    centerline = extract_single_line(centerlines, 0)
    origin = centerline.GetPoint(0)

    # Get movement
    if move_surface:
        print("-- Moving surface --")
        if patient_specific:
            move_atrium_real(case_path, mapped_path, moved_path, case)
        else:
            # Use constructed movement
            move_atrium(case_path, origin, moved_path, case)

    # Add flow extensions
    if add_extensions and path.exists(moved_path):
        surface_to_mesh = add_flow_extensions(surface, model_path, moved_path, centerlines, flow_extension_length,
                                              edge_length)

    if (not path.exists(mesh_path) or recompute_mesh) and surface_to_mesh is not None:
        mesh = generate_mesh(surface_to_mesh, mesh_path, mesh_xml_path, edge_length)

    if path.exists(mesh_path):
        mesh = read_polydata(mesh_path)

    if path.exists(mesh_path):
        find_ids(surface_to_mesh, model_path, mesh)


def remove_preprocessing_files(case_path, model_path, moved_path):
    folder = model_path.rsplit("/", 1)[0]
    folder_content = os.listdir(folder)
    files_to_remove = [f for f in folder_content if not case_path.endswith(f) and not moved_path.endswith(f)]
    for f in files_to_remove:
        directory = path.join(folder, f)
        if path.isdir(directory):
            shutil.rmtree(directory)
        else:
            os.remove(directory)


def find_ids(surface_to_mesh, model_path, mesh):
    file_name_flow_centerlines = model_path + "_cl_ext.vtp"
    probe_path = model_path + "_probes"

    print("-- Finding boundary IDs --")
    inlet, outlets = get_centers_for_meshing(surface_to_mesh, True, model_path, flowext=True)
    capped_surface = vmtk_cap_polydata(surface_to_mesh)
    centerlines, _, _ = compute_centerlines(outlets, inlet, file_name_flow_centerlines, capped_surface, resampling=0.1)
    network, probe_points = setup_model_network(centerlines, probe_path, [], lambda *a: None)

    parameters = get_parameters(model_path)

    mean_inflow_rate = compute_flow_rate(True, inlet, parameters)

    find_boundaries(model_path, mean_inflow_rate, network, mesh, lambda *a: None)

    # Display the flow split at the outlets, inlet flow rate, and probes.
    visualize(network.elements, probe_points, surface_to_mesh, mean_inflow_rate)


def IdealVolume(t):
    LA_volume = [36858.89622880263, 42041.397558417586, 47203.72790128924, 51709.730141809414, 56494.613640032476,
                 53466.224048278644, 46739.80937044214, 45723.76234837754, 46107.69142568748, 34075.82037837897]

    time = np.linspace(0, 1, len(LA_volume))
    LA_smooth = splrep(time, LA_volume, s=1e6, per=True)
    vmin = 37184.998997815936
    vmax = 19490.21405487303
    volume = (splev(t, LA_smooth) - vmin) / vmax
    return volume


def move_atrium_real(case_path, mapped_path, moved_path, case):
    surface = read_polydata(case_path)
    surface = dsa.WrapDataObject(surface)
    mapped_surfaces = sorted(os.listdir(mapped_path))
    n_frames = len(mapped_surfaces)
    for frame in range(1, n_frames):
        displaced_surface = read_polydata(path.join(mapped_path, "%s_%02d.vtp" % (case, frame)))
        displaced_surface = dsa.WrapDataObject(displaced_surface)
        displacement = displaced_surface.PointData["displacement"]
        displaced_surface.Points += displacement
        surface.Points += displacement
        write_polydata(surface.VTKObject, path.join(moved_path, "%s_%02d.vtp" % (case, frame)))
        surface.Points -= displacement


def move_atrium(case_path, origin, moved_path, case, cycle=1.0, n_frames=20):
    # Params
    A = 25 / 3 * 2
    t_array = np.linspace(0, cycle, n_frames)
    volumes = []

    surface = read_polydata(case_path)
    write_polydata(surface, path.join(moved_path, "%s_00.vtp" % case))

    def SineVolume(t):
        return -np.sin(np.pi * t)

    for i, t in enumerate(t_array):
        surface = read_polydata(case_path)
        surface = dsa.WrapDataObject(surface)
        points = surface.Points

        for j in range(len(points)):
            p = points[j]
            # displacement = IdealVolume(t)
            displacement = SineVolume(t)

            # Axial movement
            x_o = origin[0]
            y_o = origin[1]
            x_0 = p[0]
            y_0 = p[1]

            scaling_x = (x_0 - x_o)
            scaling_y = (y_0 - y_o)
            x_new = A / 100 * scaling_x * displacement
            y_new = A / 100 * scaling_y * displacement
            z_new = A / 100 * scaling_y * displacement
            #z_new = 0 # A * displacement

            # Longitudinal movement
            p_new = np.array([x_new, y_new, z_new])
            points[j] += p_new

        surface.SetPoints(points)

        write_polydata(surface.VTKObject, path.join(moved_path, "%s_%02d.vtp" % (case, i + 1)))
        capped_surface = vmtk_cap_polydata(surface.VTKObject)

        volume = vtk_compute_mass_properties(capped_surface, compute_volume=True)
        volumes.append(volume)


def capp_surface(remeshed_extended, offset=1):
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = remeshed_extended
    capper.Interactive = 0
    capper.Method = "centerpoint"
    capper.TriangleOutput = 0
    capper.CellEntityIdOffset = offset
    capper.Execute()
    surface = capper.Surface

    return surface


def add_flow_extensions(surface, model_path, moved_path, centerline, flow_extension_length, resolution=1.9):
    # Create result paths
    points_path = model_path + "_points.np"
    remeshed_path = model_path + "_remeshed.vtp"
    extended_path = model_path + "_extended"
    remeshed_extended_path = model_path + "_remeshed_extended.vtp"
    if not path.exists(extended_path):
        os.mkdir(extended_path)

    # remeshed = original
    if path.exists(remeshed_path):
        print("-- Remeshing --")
        remeshed = read_polydata(remeshed_path)
    else:
        remeshed = remesh_surface(surface, resolution)
        remeshed = vtk_clean_polydata(remeshed)
        write_polydata(remeshed, remeshed_path)

    # Create surface extensions on the original surface
    if flow_extension_length is not None:
        print("-- Adding flow extensions --")
        remeshed_extended = add_flow_extension(remeshed, centerline, include_outlet=False,
                                               extension_length=flow_extension_length)
        remeshed_extended = add_flow_extension(remeshed_extended, centerline, include_outlet=True,
                                               extension_length=flow_extension_length)

        # Smooth at edges
        remeshed_extended = vmtk_smooth_surface(remeshed_extended, "laplace", iterations=50)
    else:
        print("Skipping flow extensions")
        remeshed_extended = vmtk_smooth_surface(remeshed, "laplace", iterations=25)

    write_polydata(remeshed_extended, remeshed_extended_path)

    # Get a point mapper
    distance, point_map = get_point_map(remeshed, remeshed_extended)

    # Add extents to all surfaces
    extended_surfaces = sorted([f for f in os.listdir(moved_path) if f[:2] == "LA"])
    n_surfaces = len(extended_surfaces)

    print("-- Projecting surfaces --")
    points = np.zeros((remeshed_extended.GetNumberOfPoints(), 3, n_surfaces))
    for i in range(n_surfaces):
        model_path = path.join(moved_path, extended_surfaces[i])
        if i == 0:
            points[:, :, i] = dsa.WrapDataObject(remeshed_extended).Points
            continue

        tmp_surface = read_polydata(model_path)
        new_path = path.join(extended_path, model_path.split("/")[-1])
        if not path.exists(new_path):
            move_surface_model(tmp_surface, surface, remeshed, remeshed_extended, distance, point_map, new_path, i,
                               points)

    # Resample points and write to file
    N = 200
    points[:, :, -1] = points[:, :, 0]
    time = np.linspace(0, 1, points.shape[2])
    N2 = N + N // (time.shape[0] - 1)
    move = np.zeros((points.shape[0], points.shape[1], N + 1))
    move[:, 0, :] = resample(points[:, 0, :], N2, time, axis=1)[0][:, :N - N2 + 1]
    move[:, 1, :] = resample(points[:, 1, :], N2, time, axis=1)[0][:, :N - N2 + 1]
    move[:, 2, :] = resample(points[:, 2, :], N2, time, axis=1)[0][:, :N - N2 + 1]

    points = move
    points.dump(points_path)

    return remeshed_extended


def generate_mesh(surface, mesh_path, mesh_xml_path, resolution):
    # Cap mitral valve
    print("-- Meshing surface --")
    surface = dsa.WrapDataObject(surface)
    surface.CellData.append(np.zeros(surface.VTKObject.GetNumberOfCells()) + 1,
                            cell_id_name)
    remeshed_all_capped = capp_surface(surface.VTKObject)
    remeshed_all_capped = remesh_surface(remeshed_all_capped, resolution, exclude=[1])

    # Mesh volumetric
    sizingFunction = vtkvmtk.vtkvmtkPolyDataSizingFunction()
    sizingFunction.SetInputData(remeshed_all_capped)
    sizingFunction.SetSizingFunctionArrayName("Volume")
    sizingFunction.SetScaleFactor(0.8)
    sizingFunction.Update()

    surfaceToMesh = vmtkscripts.vmtkSurfaceToMesh()
    surfaceToMesh.Surface = sizingFunction.GetOutput()
    surfaceToMesh.Execute()

    tetgen = vmtkscripts.vmtkTetGen()
    tetgen.Mesh = surfaceToMesh.Mesh
    tetgen.GenerateCaps = 0
    tetgen.UseSizingFunction = 1
    tetgen.SizingFunctionArrayName = "Volume"
    tetgen.CellEntityIdsArrayName = cell_id_name
    tetgen.Order = 1
    tetgen.Quality = 1
    tetgen.PLC = 1
    tetgen.NoBoundarySplit = 1
    tetgen.RemoveSliver = 1
    tetgen.OutputSurfaceElements = 1
    tetgen.OutputVolumeElements = 1
    tetgen.Execute()

    mesh = tetgen.Mesh
    write_polydata(mesh, mesh_path)

    meshWriter = vmtkscripts.vmtkMeshWriter()
    meshWriter.CellEntityIdsArrayName = "CellEntityIds"
    meshWriter.Mesh = mesh
    meshWriter.Mode = "ascii"
    meshWriter.Compressed = True
    meshWriter.OutputFileName = mesh_xml_path
    meshWriter.Execute()

    return mesh


def move_surface_model(surface, original, remeshed, remeshed_extended, distance, point_map, file_path, i, points):
    surface = dsa.WrapDataObject(surface)
    original = dsa.WrapDataObject(original)
    remeshed = dsa.WrapDataObject(remeshed)
    remeshed_extended = dsa.WrapDataObject(remeshed_extended)

    if "displacement" in original.PointData.keys():
        original.VTKObject.GetPointData().RemoveArray("displacement")

    if "displacement" in remeshed_extended.PointData.keys():
        remeshed_extended.VTKObject.GetPointData().RemoveArray("displacement")

    # Get displacement field
    original.PointData.append(surface.Points - original.Points, "displacement")

    # Get
    projector = vmtkscripts.vmtkSurfaceProjection()
    projector.Surface = remeshed_extended.VTKObject
    projector.ReferenceSurface = original.VTKObject
    projector.Execute()

    # New surface
    new_surface = projector.Surface
    new_surface = dsa.WrapDataObject(new_surface)

    # Manipulate displacement in the extensions
    displacement = new_surface.PointData["displacement"]
    displacement[remeshed.Points.shape[0]:] = distance * displacement[point_map]

    # Move the mesh points
    new_surface.Points += displacement
    write_polydata(new_surface.VTKObject, file_path)
    points[:, :, i] = new_surface.Points.copy()
    new_surface.Points -= displacement


def get_point_map(remeshed, remeshed_extended):
    remeshed = dsa.WrapDataObject(remeshed)
    remeshed_extended = dsa.WrapDataObject(remeshed_extended)

    # Get lengths
    num_re = remeshed.Points.shape[0]
    num_ext = remeshed_extended.Points.shape[0] - remeshed.Points.shape[0]

    # Get locators
    inner_feature = vtk_compute_connectivity(vtk_extract_feature_edges(remeshed.VTKObject))
    outer_feature = vtk_compute_connectivity(vtk_extract_feature_edges(remeshed_extended.VTKObject))
    locator_remeshed = get_vtk_point_locator(remeshed.VTKObject)

    n_features = outer_feature.GetPointData().GetArray("RegionId").GetValue(outer_feature.GetNumberOfPoints() - 1)
    inner_features = np.array(
        [vtk_compute_threshold(inner_feature, "RegionId", i, i + 0.1, source=0) for i in range(n_features + 1)])
    outer_features = np.array(
        [vtk_compute_threshold(outer_feature, "RegionId", i, i + 0.1, source=0) for i in range(n_features + 1)])

    inner_regions = [dsa.WrapDataObject(feature) for feature in inner_features]
    inner_locators = [get_vtk_point_locator(feature) for feature in inner_features]
    inner_points = [feature.GetNumberOfPoints() for feature in inner_features]

    outer_regions = [dsa.WrapDataObject(feature) for feature in outer_features]
    outer_locators = [get_vtk_point_locator(feature) for feature in outer_features]
    outer_points = [feature.GetNumberOfPoints() for feature in outer_features]
    boundary_map = {i: j for i, j in zip(np.argsort(inner_points), np.argsort(outer_points))}

    # Get distance and point map
    distances = np.zeros(num_ext)
    point_map = np.zeros(num_ext)

    for i in range(num_ext):
        point = remeshed_extended.Points[num_re + i]
        tmp_dist = 1E16
        tmp_id = -1
        # Some hacks to find the correct corresponding points
        for region_id in range(len(outer_features)):
            region_id_out = boundary_map[region_id]
            id_i = inner_locators[region_id].FindClosestPoint(point)
            id_o = outer_locators[region_id_out].FindClosestPoint(point)

            p_i = inner_features[region_id].GetPoint(id_i)
            p_o = outer_features[region_id_out].GetPoint(id_o)

            dist_o = get_distance(point, p_i)
            dist_i = get_distance(point, p_o)
            dist_total = dist_i + dist_o
            if dist_total < tmp_dist:
                tmp_dist = dist_total
                tmp_id = region_id

        regionId = tmp_id
        regionId_out = boundary_map[regionId]
        inner_id = inner_locators[regionId].FindClosestPoint(point)
        outer_id = outer_locators[regionId_out].FindClosestPoint(point)

        dist_to_boundary = get_distance(point, outer_regions[regionId_out].Points[outer_id])
        dist_between_boundaries = get_distance(inner_regions[regionId].Points[inner_id],
                                               outer_regions[regionId_out].Points[outer_id])
        distances[i] = dist_to_boundary / dist_between_boundaries
        point_map[i] = locator_remeshed.FindClosestPoint(inner_regions[regionId].Points[inner_id])

    # Let the points corresponding to the caps have distance 0
    point_map = point_map.astype(int)
    return distances, point_map


def add_flow_extension(surface, centerlines, include_outlet, extension_length=2.0):
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
        else:
            if i != outlet_id:
                boundaryIds.InsertNextId(i)

    flowExtensionsFilter = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
    flowExtensionsFilter.SetInputData(surface)
    flowExtensionsFilter.SetCenterlines(centerlines)
    flowExtensionsFilter.SetAdaptiveExtensionLength(1)
    flowExtensionsFilter.SetAdaptiveExtensionRadius(1)
    flowExtensionsFilter.SetAdaptiveExtensionLength(1)
    flowExtensionsFilter.SetAdaptiveNumberOfBoundaryPoints(1)
    flowExtensionsFilter.SetExtensionRatio(extension_length)
    flowExtensionsFilter.SetExtensionRadius(1.0)
    flowExtensionsFilter.SetTransitionRatio(1.0)
    flowExtensionsFilter.SetCenterlineNormalEstimationDistanceRatio(1.0)
    flowExtensionsFilter.SetExtensionModeToUseCenterlineDirection()
    flowExtensionsFilter.SetInterpolationModeToThinPlateSpline()
    flowExtensionsFilter.SetBoundaryIds(boundaryIds)
    flowExtensionsFilter.SetSigma(1.0)
    flowExtensionsFilter.Update()

    surface = flowExtensionsFilter.GetOutput()

    return surface


def remesh_surface(surface, edge_length, exclude=None):
    surface = dsa.WrapDataObject(surface)
    if cell_id_name not in surface.CellData.keys():
        surface.CellData.append(np.zeros(surface.VTKObject.GetNumberOfCells()) + 1, cell_id_name)
    remeshing = vmtkscripts.vmtkSurfaceRemeshing()
    remeshing.Surface = surface.VTKObject
    remeshing.CellEntityIdsArrayName = cell_id_name
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


def read_command_line():
    """
    Read arguments from commandline and return all values in a dictionary.
    """
    '''Command-line arguments.'''
    parser = argparse.ArgumentParser(
        description="Add rigid flow extensions to Atrial models.")

    parser.add_argument('-i', '--inputModel', type=str, required=True, dest='fileNameModel',
                        help="Input file containing the 3D model.")
    parser.add_argument('-m', '--movement', type=str2bool, required=False, default=False, dest="moveSurface",
                        help="Add movement to input surface.")
    parser.add_argument('-e', '--extension', type=str2bool, required=False, default=False, dest="addExtensions",
                        help="Add extension to moved surface.")
    parser.add_argument('-el', '--edge-length', type=float, required=False, default=1.9, dest="edgeLength",
                        help="Edge length resolution for meshing.")
    parser.add_argument('-p', '--patient-specific', type=str2bool, required=False, default=False,
                        dest="patientSpecific", help="Use patient specific data or constructed movement.")
    parser.add_argument('-r', '--recompute-mesh', type=str2bool, required=False, default=False,
                        dest="recomputeMesh", help="Recomputes mesh if true.")
    parser.add_argument('-ra', '--recompute-all', type=str2bool, required=False, default=False,
                        dest="recomputeAll", help="Recomputes everything if true.")
    parser.add_argument('-fl', '--flowextlen', dest="flowExtLen", default=None, type=float,
                        help="Adds flow extensions of given length.")

    args, _ = parser.parse_known_args()

    return dict(case_path=args.fileNameModel, move_surface=args.moveSurface, add_extensions=args.addExtensions,
                edge_length=args.edgeLength, patient_specific=args.patientSpecific, recompute_mesh=args.recomputeMesh,
                recompute_all=args.recomputeAll, flow_extension_length=args.flowExtLen)


if __name__ == '__main__':
    main(**read_command_line())
