# Python 2 / 3 support
from __future__ import print_function

import argparse
import json

import ImportData
import ToolRepairSTL
from NetworkBoundaryConditions import FlowSplitting
# Local imports
from common import *
from simulate import run_simulation
from visualize import visualize


def run_pre_processing(filename_model, verbose_print, smoothing_method, smoothing_factor, smooth_aneurysm,
                       meshing_method, aneurysm_present, create_flow_extensions, viz, config_path, number_of_sac_points,
                       coarsening_factor, compress_mesh=True):
    """
    Automatically generate mesh of surface model in .vtu and .xml format, including prescribed
    flow rates at inlet and outlet based on flow network model.

    Runs simulation of meshed case on a remote ssh server if server configuration is provided.

    Args:
        filename_model (str): Name of case
        verbose_print (bool): Toggles verbose mode
        smoothing_method (str): Method for surface smoothing
        smoothing_factor (float): Smoothing parameter
        smooth_aneurysm (bool): Toggles smoothing of aneurysm (if present)
        meshing_method (str): Method for meshing
        aneurysm_present (bool): Determines if aneurysm is present
        create_flow_extensions (bool): Adds flow extensions to mesh if True
        viz (bool): Visualize resulting surface model with flow rates
        config_path (str): Path to configuration file for remote simulation
        number_of_sac_points (int): Number of sac points to evaluate (inside aneurysm)
        coarsening_factor (float): Refine or coarsen the standard mesh size with given factor
        compress_mesh (bool): Compresses finalized mesh if True
    """
    # Get paths
    abs_path = path.abspath(path.dirname(__file__))
    case_name = filename_model.rsplit(path.sep, 1)[-1].rsplit('.')[0]
    dir_path = filename_model.rsplit(path.sep, 1)[0]

    # Naming conventions
    file_name_centerlines = path.join(dir_path, case_name + "_centerlines.vtp")
    file_name_aneurysm_centerlines = path.join(dir_path, case_name + "_aneurysm_centerline.vtp")
    file_name_sac_centerlines = path.join(dir_path, case_name + "_sac_centerline_{}.vtp")
    file_name_distance_to_sphere_diam = path.join(dir_path, case_name + "_distance_to_sphere_diam.vtp")
    file_name_distance_to_sphere_curv = path.join(dir_path, case_name + "_distance_to_sphere_curv.vtp")
    file_name_probe_points = path.join(dir_path, case_name + "_probe_point")
    file_name_voronoi = path.join(dir_path, case_name + "_voronoi.vtp")
    file_name_voronoi_smooth = path.join(dir_path, case_name + "_voronoi_smooth.vtp")
    file_name_surface_smooth = path.join(dir_path, case_name + "_smooth.vtp")
    file_name_model_flow_ext = path.join(dir_path, case_name + "_flowext.vtp")
    file_name_clipped_model = path.join(dir_path, case_name + "_clippedmodel.vtp")
    file_name_flow_centerlines = path.join(dir_path, case_name + "_flow_cl.vtp")
    file_name_surface_name = path.join(dir_path, case_name + "_remeshed_surface.vtp")
    file_name_xml_mesh = path.join(dir_path, case_name + ".xml")
    file_name_vtu_mesh = path.join(dir_path, case_name + ".vtu")
    file_name_run_script = path.join(dir_path, case_name + ".sh")

    print("\n--- Working on case:", case_name, "\n")

    # Open the surface file.
    print("--- Load model file\n")
    surface = ReadPolyData(filename_model)

    if not is_surface_capped and smoothing_method != "voronoi":
        print("--- Clipping the models inlet and outlets.\n")
        if not path.isfile(file_name_clipped_model):
            # TODO: Check if this is a valid call to this method
            centerline = compute_centerlines([], [], None, surface, method="pickpoint")
            surface = uncapp_surface(surface, centerline, filename=None)
        else:
            surface = ReadPolyData(file_name_clipped_model)

    parameters = get_parameters(path.join(dir_path, case_name))

    if "check_surface" not in parameters.keys():
        surface = surface_cleaner(surface)
        surface = triangulate_surface(surface)

        # Check the mesh if there is redundant nodes or NaN triangles.
        ToolRepairSTL.surfaceOverview(surface)
        ToolRepairSTL.foundAndDeleteNaNTriangles(surface)
        surface = ToolRepairSTL.cleanTheSurface(surface)
        foundNaN = ToolRepairSTL.foundAndDeleteNaNTriangles(surface)
        if foundNaN:
            raise RuntimeError(("There is an issue with the surface. "
                                "Nan coordinates or some other shenanigans."))
        else:
            parameters["check_surface"] = True
            write_parameters(parameters, path.join(dir_path, case_name))

    # Capp surface if open
    if not compute_centers(surface, test_capped=True):
        capped_surface = capp_surface(surface)
    else:
        capped_surface = surface

    # Get centerlines
    print("--- Get centerlines\n")
    inlet, outlets = get_centers(surface, path.join(dir_path, case_name))
    centerlines = compute_centerlines(inlet, outlets, file_name_centerlines,
                                      capped_surface, resampling=0.1, end_point=0)
    tol = get_tolerance(centerlines)

    if aneurysm_present:
        aneurysms = get_aneurysm_dome(capped_surface, path.join(dir_path, case_name))
        centerlineAnu = compute_centerlines(inlet, aneurysms, file_name_aneurysm_centerlines,
                                            capped_surface, resampling=0.1)

        # Extract the aneurysm centerline
        sac_centerline = []
        info = get_parameters(path.join(dir_path, case_name))
        num_anu = info["number_of_aneurysms"]

        # Compute mean distance between points
        for i in range(num_anu):
            if not path.isfile(file_name_sac_centerlines.format(i)):
                line = ExtractSingleLine(centerlineAnu, i)
                locator = get_locator(centerlines)
                for j in range(line.GetNumberOfPoints() - 1, 0, -1):
                    point = line.GetPoints().GetPoint(j)
                    ID = locator.FindClosestPoint(point)
                    tmp_point = centerlines.GetPoints().GetPoint(ID)
                    dist = np.sqrt(np.sum((np.asarray(point) - np.asarray(tmp_point)) ** 2))
                    if dist <= tol:
                        break

                tmp = ExtractSingleLine(line, 0, startID=j)
                WritePolyData(tmp, file_name_sac_centerlines.format(i))

                # List of VtkPolyData sac(s) centerline
                sac_centerline.append(tmp)

            else:
                sac_centerline.append(ReadPolyData(file_name_sac_centerlines.format(i)))

    else:
        num_anu = 0

    # Get 'center' and 'radius' of the aneurysm(s)
    sac_center = []
    misr_max = []

    if aneurysm_present:
        # Merge the sac centerline
        sac_centerlines = merge_data(sac_centerline)

        for sac in sac_centerline:
            sac_center.append(sac.GetPoints().GetPoint(sac.GetNumberOfPoints() // 2))
            tmp_misr = get_array(radiusArrayName, sac)
            misr_max.append(tmp_misr.max())

    # Smooth surface
    if smoothing_method == "voronoi":
        print("--- Smooth surface: Voronoi smoothing\n")
        if not path.isfile(file_name_surface_smooth):
            # Get Voronoi diagram
            if not path.isfile(file_name_voronoi):
                voronoi = makeVoronoiDiagram(surface, file_name_voronoi)
                WritePolyData(voronoi, file_name_voronoi)
            else:
                voronoi = ReadPolyData(file_name_voronoi)

            # Get smooth Voronoi diagram
            if not path.isfile(file_name_voronoi_smooth):
                if aneurysm_present:
                    smooth_voronoi = SmoothVoronoiDiagram(voronoi, centerlines, smoothing_factor, sac_centerlines)
                else:
                    smooth_voronoi = SmoothVoronoiDiagram(voronoi, centerlines, smoothing_factor)

                WritePolyData(smooth_voronoi, file_name_voronoi_smooth)
            else:
                smooth_voronoi = ReadPolyData(file_name_voronoi_smooth)

            # Envelope the smooth surface
            surface = create_new_surface(smooth_voronoi)

            # Uncapp the surface
            surface_uncapped = uncapp_surface(surface, centerlines, filename=None)

            # Check if there has been added new outlets
            num_outlets = centerlines.GetNumberOfLines()
            num_outlets_after = compute_centers(surface_uncapped, test_capped=True)[1]

            if num_outlets != num_outlets_after:
                surface = vmtkSmoother(surface, "laplace", iterations=200)
                WritePolyData(surface, file_name_surface_smooth)
                print(("ERROR: Automatic clipping failed. You have to open {} and " +
                       "manually clipp the branch which still is capped. " +
                       "Overwrite the current {} and restart the script.").format(
                    file_name_surface_smooth, file_name_surface_smooth))
                sys.exit(0)

            surface = surface_uncapped

            # Smoothing to improve the quality of the elements
            # Consider to add a subdivision here as well.
            surface = vmtkSmoother(surface, "laplace", iterations=200)

            # Write surface
            WritePolyData(surface, file_name_surface_smooth)

        else:
            surface = ReadPolyData(file_name_surface_smooth)


    elif smoothing_method == "laplace":
        print("--- Smooth surface: Laplacian smoothing\n")
        if not path.isfile(file_name_surface_smooth):
            surface = vmtkSmoother(surface, smoothing_method)

            # Save the smoothed surface
            WritePolyData(surface, file_name_surface_smooth)

        else:
            surface = ReadPolyData(file_name_surface_smooth)

    elif smoothing_method == "taubin":
        print("--- Smooth surface: Taubin smoothing\n")
        if not path.isfile(file_name_surface_smooth):
            surface = vmtkSmoother(surface, smoothing_method)

            # Save the smoothed surface
            WritePolyData(surface, file_name_surface_smooth)

        else:
            surface = ReadPolyData(file_name_surface_smooth)

    elif smoothing_method == "no_smooth" or None:
        print("--- No smoothing of surface\n")

    # Add flow extensions
    if create_flow_extensions:
        if not path.isfile(file_name_model_flow_ext):
            print("--- Adding flow extensions")
            extension = 5

            extender = vmtkscripts.vmtkFlowExtensions()
            extender.Surface = surface
            extender.AdaptiveExtensionLength = 1
            extender.ExtensionRatio = extension
            extender.Centerlines = centerlines
            extender.ExtensionMode = "boundarynormal"
            extender.CenterlineNormalEstimationDistanceRatio = 1.0
            extender.Interactive = 0
            extender.Execute()

            surface = extender.Surface
            surface = vmtkSmoother(surface, "laplace", iterations=200)
            WritePolyData(surface, file_name_model_flow_ext)

        else:
            surface = ReadPolyData(file_name_model_flow_ext)

    # Capp surface with flow extensions
    capped_surface = capp_surface(surface)

    # Get new centerlines with the flow extensions
    if not path.isfile(file_name_flow_centerlines):
        print("--- Compute the model centerlines with flow extension.")
        # Compute the centerlines.
        inlet, outlets = get_centers(surface, path.join(dir_path, case_name), flowext=True)
        centerlines = compute_centerlines(inlet, outlets, file_name_flow_centerlines, capped_surface, resampling=0.5)

    else:
        centerlines = ReadPolyData(file_name_flow_centerlines)

    # Choose input for the mesh
    print("--- Computing distance to sphere")
    if meshing_method == "curvature":
        if not path.isfile(file_name_distance_to_sphere_curv):
            distance_to_sphere = dist_sphere_curv(surface, centerlines,
                                                  sac_center, misr_max,
                                                  file_name_distance_to_sphere_curv,
                                                  coarsening_factor)
        else:
            distance_to_sphere = ReadPolyData(file_name_distance_to_sphere_curv)
    else:
        if not path.isfile(file_name_distance_to_sphere_diam):
            distance_to_sphere = dist_sphere_diam(surface, centerlines,
                                                  sac_center, misr_max,
                                                  file_name_distance_to_sphere_diam,
                                                  coarsening_factor)
        else:
            distance_to_sphere = ReadPolyData(file_name_distance_to_sphere_diam)

    # Compute mesh
    if not path.isfile(file_name_vtu_mesh):
        try:
            print("--- Computing mesh\n")
            mesh, remeshed_surface = generate_mesh(distance_to_sphere)
            assert remeshed_surface.GetNumberOfPoints() > 0, \
                "No points in surface mesh, try to remesh"
            assert mesh.GetNumberOfPoints() > 0, "No points in mesh, try to remesh"

        except:
            distance_to_sphere = mesh_alternative(distance_to_sphere)
            mesh, remeshed_surface = generate_mesh(distance_to_sphere)
            assert mesh.GetNumberOfPoints() > 0, "No points in mesh, after remeshing"
            assert remeshed_surface.GetNumberOfPoints() > 0, \
                "No points in surface mesh, try to remesh"

        # Write mesh in VTU format
        WritePolyData(remeshed_surface, file_name_surface_name)
        WritePolyData(mesh, file_name_vtu_mesh)

        # Write mesh to FEniCS to format
        meshWriter = vmtkscripts.vmtkMeshWriter()
        meshWriter.CellEntityIdsArrayName = "CellEntityIds"
        meshWriter.Mesh = mesh
        meshWriter.Mode = "ascii"
        meshWriter.Compressed = compress_mesh
        meshWriter.OutputFileName = file_name_xml_mesh
        meshWriter.Execute()
        polyDataVolMesh = mesh

    else:
        polyDataVolMesh = ReadPolyData(file_name_vtu_mesh)

    # Set the network object used in the scripts for 
    # boundary conditions and probes.
    network = ImportData.Network()
    centerlinesBranches = ImportData.SetNetworkStructure(centerlines, network, verbose_print)

    if not path.isfile(file_name_probe_points):
        # Get the list of coordinates for the probe points along the network centerline.
        listProbePoints = ImportData.GetListProbePoints(centerlinesBranches, network, verbose_print)
        listProbePoints += sac_center

        # Add points randomly in the sac.
        # FIXME: This is not robust enough. Suggestion to fix: Extract the
        # second half of the sac centerline, then get all points from the
        # voronoi diagram which is closest to that part compared to any ther
        # centerlines. Then randomly chose among those points. For now, simply
        # add just one point (sac_center).
        # numberOfPoints = numberOfSacPoints
        # for k in range(num_anu):
        #    u = np.random.uniform(0.0, 1.0, (numberOfPoints, 1))
        #    theta = np.random.uniform(0., 1., (numberOfPoints, 1)) * np.pi
        #    phi = np.arccos(1 - 2 * np.random.uniform(0.0, 1., (numberOfPoints, 1)))
        #    radius = misr_max[k] * u**(0.3333)
        #    x = radius * np.sin(theta) * np.cos(phi)
        #    y = radius * np.sin(theta) * np.sin(phi)
        #    z = radius * np.cos(theta)
        #    for i in range(len(x)):
        #        listProbePoints.append([np.array(misr_max_center[k][0] + x[i]).tolist()[0],
        #                                np.array(misr_max_center[k][1] + y[i]).tolist()[0],
        #                                np.array(misr_max_center[k][2] + z[i]).tolist()[0]])

        print("--- Saving probes points in: ", file_name_probe_points)
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

    # BSL method for mean inlet flow rate.
    parameters = get_parameters(path.join(dir_path, case_name))
    mean_inflow_rate = 0.27 * parameters["inlet_area"]

    # Extract the surface mesh of the wall
    wallMesh = threshold(polyDataVolMesh, "CellEntityIds", lower=0.5, upper=1.5)

    boundaryReferenceSystems = vmtkscripts.vmtkBoundaryReferenceSystems()
    boundaryReferenceSystems.Surface = wallMesh
    boundaryReferenceSystems.Execute()
    refSystem = boundaryReferenceSystems.ReferenceSystems
    cellEntityIdsArray = get_vtk_array('CellEntityIds', 0, refSystem.GetNumberOfPoints())
    refSystem.GetPointData().AddArray(cellEntityIdsArray)

    # Extract the surface mesh of the end caps
    boundarySurface = threshold(polyDataVolMesh, "CellEntityIds", upper=1.5, type="upper")

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
    idFileLine = case_name + ' ' + repr(ids[0][0] - 1) + ' '
    areaRatioLine = case_name + ' '
    for k in range(1, refSystem.GetNumberOfPoints() - 1):
        idFileLine += repr(ids[k][0] - 1) + ','
        areaRatioLine += repr(ids[k][1]) + ','
    idFileLine += repr(ids[-1][0] - 1) + ' ' + repr(ids[0][1])
    areaRatioLine += repr(ids[-1][1])

    with open(path.join(dir_path, case_name + '.txt'), 'w') as outfile:
        outfile.write("\ninlet_area: " + str(parameters["inlet_area"]))
        outfile.write('\nidFileLine: ' + str(idFileLine))
        outfile.write('\nareaRatioLine: ' + str(areaRatioLine))

    # Display the flow split at the outlets, inlet flow rate, and probes.
    if viz:
        visualize(network.elements, probe_points, surface, mean_inflow_rate)

    # Start simulation though ssh, without password
    if config_path is not None:

        # Set up simulation script
        if not path.exists(file_name_run_script):
            run_script_sample = open(path.join(abs_path, "run_script.sh"), "r").read()
            config = json.load(open(config_path))
            run_dict = dict(mesh_name=case_name,
                            num_nodes=1,
                            hours=120,
                            account="nn9249k",
                            remoteFolder=config["remoteFolder"],
                            results_folder="results")
            run_script = run_script_sample.format(**run_dict)

            # Write script
            script_file = open(file_name_run_script, "w")
            script_file.write(run_script)
            script_file.close()

        run_simulation(config_path, dir_path, case_name)


def str2bool(arg):
    """
    Convert a string to boolean.

    Args:
        arg (str): Input string.

    Returns:
        return (bool): Converted string.
    """
    if arg.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif arg.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def read_command_line():
    """
    Read arguments from commandline and return all values in a dictionary.
    """
    '''Command-line arguments.'''
    parser = argparse.ArgumentParser(
        description="Automatic pre-processing for FEniCS.")

    parser.add_argument('-v', '--verbosity',
                        dest='verbosity',
                        type=str2bool,
                        default=False,
                        help="Activates the verbose mode.")

    parser.add_argument('-i', '--inputModel',
                        type=str,
                        required=False,
                        dest='fileNameModel',
                        default='example/surface.vtp',
                        help="Input file containing the 3D model.")

    parser.add_argument('-sM', '--smoothingMethod',
                        type=str,
                        required=False,
                        dest='smoothingMethod',
                        default="no_smooth",
                        choices=["voronoi", "no_smooth", "laplace", "taubin"],
                        help="Smoothing method, for now only Voronoi smoothing is available." +
                             " For Voronoi smoothing you can also control smoothingFactor" +
                             " (default = 0.25)  and smoothingAneurysm (default = False).")

    parser.add_argument('-c', '--coarseningFactor',
                        type=float,
                        required=False,
                        dest='coarseningFactor',
                        default=1.0,
                        help="Refine or coarsen the standard mesh size")

    parser.add_argument('-sF', '--smoothingFactor',
                        type=float,
                        required=False,
                        dest='smoothingFactor',
                        default=0.25,
                        help="smoothingFactor for VoronoiSmoothing, removes all spheres which" +
                             " has a radius < MISR*(1-0.25), where MISR varying along the centerline.")

    parser.add_argument('-sA', '--smoothAneurysm',
                        type=str2bool,
                        required=False,
                        dest="smoothingAneurysm",
                        help="When using Voronoi smoothing one can choose not to smooth the" +
                             " aneurysm as this method often removes to much of the aneurysm")

    parser.add_argument('-m', '--meshingMethod',
                        dest="meshingMethod",
                        type=str,
                        choices=["diameter", "curvature"],
                        default="diameter")

    parser.add_argument('-a', '--aneurysm',
                        dest="aneu",
                        type=str2bool,
                        default=True,
                        help="Determine weather or not the model has a aneurysm. Default is False.")

    parser.add_argument('-f', '--flowext',
                        dest="fext",
                        default=True,
                        type=str2bool,
                        help="Add flow extensions to to the model.")

    parser.add_argument('-vz', '--visualize',
                        dest="viz",
                        default=True,
                        type=str2bool,
                        help="Visualize surface, inlet, outlet and probes after meshing.")

    parser.add_argument('-sp', '--sacpoints',
                        type=int,
                        help='Number of sac points to add',
                        default=20,
                        dest="sacpts")

    parser.add_argument('--simulationConfig',
                        type=str,
                        dest="config",
                        default=None,
                        help='Path to configuration file for remote simulation. ' +
                             'See example/ssh_config.json for details')

    args, _ = parser.parse_known_args()

    if args.verbosity:
        print()
        print("--- VERBOSE MODE ACTIVATED ---")

        def verbose_print(*args):
            for arg in args:
                print(arg, end=' ')
                print()
    else:
        verbose_print = lambda *a: None

    verbose_print(args)

    return dict(filename_model=args.fileNameModel, verbose_print=verbose_print, smoothing_method=args.smoothingMethod,
                smoothing_factor=args.smoothingFactor, smooth_aneurysm=args.smoothingAneurysm,
                meshing_method=args.meshingMethod, aneurysm_present=args.aneu, create_flow_extensions=args.fext,
                viz=args.viz, config_path=args.config, number_of_sac_points=args.sacpts,
                coarsening_factor=args.coarseningFactor)


if __name__ == "__main__":
    run_pre_processing(**read_command_line())
