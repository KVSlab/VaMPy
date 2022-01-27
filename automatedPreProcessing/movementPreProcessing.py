import argparse
import os
import shutil

from morphman.common import *
from vtk.numpy_interface import dataset_adapter as dsa

from common import get_centers_for_meshing, find_boundaries, setup_model_network, compute_flow_rate, add_flow_extension
from moving_common import create_funnel, get_point_map, move_surface_model
from visualize import visualize

cell_id_name = "CellEntityIds"


def main(case_path, move_surface, add_extensions, edge_length, patient_specific, recompute_mesh, recompute_all,
         flow_extension_length, clamp_boundaries):
    """
    Automatically generate movement and  mesh of surface model in .vtu and .xml format.
    Assumes the user either has a set of displaced models or models with mapped displacement fields,
    located in moved_path or mapped_path, respectively.

    Can add patient-specific or arbitrary movement to model.

    Args:
        flow_extension_length (float): Length of flow extensions, factor of MISR
        clamp_boundaries (bool): Clamps inlet(s) and outlet if true
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
    remeshed_path = model_path + "_remeshed.vtp"
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

    # # remeshed = original
    # if path.exists(remeshed_path):
    #     print("-- Remeshing --")
    #     surface = read_polydata(remeshed_path)
    # else:
    #     surface = remesh_surface(surface, edge_length)
    #     surface = vtk_clean_polydata(surface)
    #     write_polydata(surface, remeshed_path)

    # Get movement
    if move_surface:
        print("-- Moving surface --")
        if patient_specific:
            move_atrium_real(case_path, mapped_path, moved_path, case)
        else:
            # Use constructed movement
            move_atrium_generic(case_path, origin, moved_path, case, centerlines)
    exit()
    # Add flow extensions
    if add_extensions and path.exists(moved_path):
        surface_to_mesh = add_flow_extensions(surface, model_path, moved_path, centerlines, flow_extension_length,
                                              clamp_boundaries, edge_length)

    if (not path.exists(mesh_path) or recompute_mesh) and surface_to_mesh is not None:
        mesh = generate_mesh_without_layers(surface_to_mesh, mesh_path, mesh_xml_path, edge_length)

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


def move_atrium_generic(remeshed_path, origin, moved_path, case, centerlines, cycle=1.0, n_frames=20):
    area_path = remeshed_path.replace(".stl", "") + "_area.txt"
    volume_path = remeshed_path.replace(".stl", "") + "_volume.txt"
    inlet_points = []
    for i in range(centerlines.GetNumberOfLines()):
        centerline = extract_single_line(centerlines, i)
        ip = centerline.GetPoint(centerline.GetNumberOfPoints() - 1)
        inlet_points.append(np.array(ip))
    inlet_points = np.array(inlet_points)

    x_i = np.mean(inlet_points[:, 0])
    y_i = np.mean(inlet_points[:, 1])
    z_i = np.mean(inlet_points[:, 2])
    p0 = np.array([x_i, y_i, z_i])  # At mean between PVs

    x_o = origin[0]
    y_o = origin[1]
    z_o = origin[2]
    p1 = origin  # At MV

    # Get normal vector
    n = n_z = (p1 - p0) / la.norm(p1 - p0)
    n_x = np.array([0, -n[2], n[1]])
    n_y = np.cross(n_z, n_x)

    # Params
    A = 1 / 4 * 1 / 3.5
    t_array = np.linspace(0, cycle, n_frames)
    volumes = []
    areas = []
    surface = read_polydata(remeshed_path)
    write_polydata(surface, path.join(moved_path, "%s_00.vtp" % case))

    def gaussian(t, sd, mu):
        return 1 / np.sqrt(2 * np.pi * sd ** 2) * np.exp((-(t - mu) ** 2 / (2 * sd ** 2)))

    def GaussianVolume(t, with_af=True):
        scale = 2
        t_systole = 0.75 / scale
        t_a_wave = 1.5 / scale

        sd_s = 0.3 / scale
        sd_a = 0.2 / scale

        factor = 1 / 2
        if with_af:
            return gaussian(t, sd_s, t_systole)
        else:
            return gaussian(t, sd_s, t_systole) + factor * gaussian(t, sd_a, t_a_wave)

    def SineVolume(t):
        return np.sin(np.pi * t)

    def get_scale(p):
        return 0.5 * 0.4

    is_af = True
    for i, t in enumerate(t_array):
        surface = read_polydata(remeshed_path)
        surface = dsa.WrapDataObject(surface)
        points = surface.Points

        for j in range(len(points)):
            p = points[j]

            # Get displacement profile
            # displacement = IdealVolume(t)
            d = GaussianVolume(t, is_af)

            d_z = d * get_scale(p)

            x_mark = (p - p0).dot(n_x)
            y_mark = (p - p0).dot(n_y)
            z_mark = (p - p0).dot(n_z)

            x_new = x_mark * d * A
            y_new = y_mark * d * A
            z_new = z_mark * d_z * A

            X_mark = np.array([x_new, y_new, z_new])
            R_inv = np.array([n_x, n_y, n_z]).T

            p_new = R_inv.dot(X_mark)
            points[j] += p_new

        surface.SetPoints(points)

        write_polydata(surface.VTKObject, path.join(moved_path, "%s_%02d.vtp" % (case, i + 1)))
        capped_surface = vmtk_cap_polydata(surface.VTKObject)

        surf_area = vtk_compute_mass_properties(capped_surface, compute_volume=False, compute_surface_area=True)
        volume = vtk_compute_mass_properties(capped_surface, compute_volume=True, compute_surface_area=False)
        volumes.append(volume)
        areas.append(surf_area)
    np.savetxt(area_path, areas)
    np.savetxt(volume_path, volumes)


def move_atrium(case_path, origin, moved_path, case, centerlines, cycle=1.0, n_frames=20):
    inlet_points = []
    for i in range(centerlines.GetNumberOfLines()):
        centerline = extract_single_line(centerlines, i)
        ip = centerline.GetPoint(centerline.GetNumberOfPoints() - 1)
        inlet_points.append(np.array(ip))
    inlet_points = np.array(inlet_points)

    x_i = np.mean(inlet_points[:, 0])
    y_i = np.mean(inlet_points[:, 1])
    z_i = np.mean(inlet_points[:, 1])

    x_o = origin[0]
    y_o = origin[1]
    z_o = origin[2]

    # Params
    A = 25 / 3 * 2
    t_array = np.linspace(0, cycle, n_frames)
    volumes = []

    surface = read_polydata(case_path)
    write_polydata(surface, path.join(moved_path, "%s_00.vtp" % case))

    def SineVolume(t):
        return np.sin(np.pi * t)

    for i, t in enumerate(t_array):
        surface = read_polydata(case_path)
        surface = dsa.WrapDataObject(surface)
        points = surface.Points

        for j in range(len(points)):
            p = points[j]
            # displacement = IdealVolume(t)
            d = SineVolume(t)

            # Axial movement
            x_0 = p[0]
            y_0 = p[1]

            scaling_x = (x_0 - x_o)
            scaling_y = (y_0 - y_o)
            x_new = A / 100 * scaling_x * d
            y_new = A / 100 * scaling_y * d
            z_new = A / 100 * scaling_y * d
            # z_new = 0 # A * displacement

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


def add_flow_extensions(surface, model_path, moved_path, centerline, flow_extension_length, clamp_boundaries,
                        resolution=1.9):
    # Create result paths
    points_path = model_path + "_points.np"
    extended_path = model_path + "_extended"
    remeshed_path = model_path + "_remeshed.vtp"
    remeshed_extended_path = model_path + "_remeshed_extended.vtp"
    extended_cl_path = model_path + "_cl_ext.vtp"
    funnel_path = model_path + "_funnel.vtp"
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
                                               extension_length=flow_extension_length * 3,
                                               extension_mode="boundarynormal")

        # Smooth at edges
        remeshed_extended = vmtk_smooth_surface(remeshed_extended, "laplace", iterations=50)
    else:
        print("-- Skipping flow extensions --")
        remeshed_extended = vmtk_smooth_surface(remeshed, "laplace", iterations=50)

    write_polydata(remeshed_extended, remeshed_extended_path)

    capped_surface = vmtk_cap_polydata(remeshed_extended)
    outlets, inlet = get_centers_for_meshing(remeshed_extended, True, model_path)
    centerlines_extended, _, _ = compute_centerlines(inlet, outlets, extended_cl_path, capped_surface, resampling=0.1)
    cl0 = extract_single_line(centerline, 0)
    cl0_ext = extract_single_line(centerlines_extended, 1)

    remeshed_extended = create_funnel(remeshed_extended, cl0, cl0_ext, 1)
    write_polydata(remeshed_extended, funnel_path)

    # Get a point mapper
    distance, point_map = get_point_map(remeshed, remeshed_extended)

    # Add extents to all surfaces
    extended_surfaces = sorted([f for f in os.listdir(moved_path) if f[:2] == "LA" or f[:2] == "Fu"])
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
                               points, clamp_boundaries)
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


def generate_mesh_without_layers(surface, mesh_path, mesh_xml_path, resolution):
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
    parser.add_argument('-c', '--clamp-boundaries', dest="clamp", default=False, type=str2bool,
                        help="Clamps boundaries at inlet(s) and outlet if true.")

    args, _ = parser.parse_known_args()

    return dict(case_path=args.fileNameModel, move_surface=args.moveSurface, add_extensions=args.addExtensions,
                edge_length=args.edgeLength, patient_specific=args.patientSpecific, recompute_mesh=args.recomputeMesh,
                recompute_all=args.recomputeAll, flow_extension_length=args.flowExtLen, clamp_boundaries=args.clamp)


if __name__ == '__main__':
    main(**read_command_line())
