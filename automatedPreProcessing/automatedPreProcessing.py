###!/usr/bin/python

from __future__ import print_function # Python 2 / 3 support

# Local imports
from common import *
from set_aneurysm import *
from visualize import visualize
from DisplayData import DisplayModel, VtkPointCloud
from NetworkBoundaryConditions import FlowSplitting
from simulate import run_simulation
import ImportData
import ToolRepairSTL

from os import path
from IPython import embed
from vmtk import vmtkscripts

import vtk
import json
import numpy as np
import argparse


def str2bool(arg):
    if arg.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif arg.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def Program(fileNameModel, verboseprint, smoothing, smoothing_factor,
            smooth_aneurysm, meshingMethod, aneurysm,
            createFlowext, viz, configPath, numberOfSacPoints, coarsening_factor):
    # Get paths
    abs_path = path.abspath(path.dirname(__file__))
    caseName = fileNameModel.rsplit(path.sep, 1)[-1].rsplit('.')[0]
    dir_path = fileNameModel.rsplit(path.sep, 1)[0]

    # Naming conventions
    fileNameInitialCenterlines = path.join(dir_path, caseName + "_initial_centerlines.vtp")
    fileNameCenterlines = path.join(dir_path, caseName + "_centerlines.vtp")
    filenameAneurysmCenterlines = path.join(dir_path, caseName + "_aneurysm_centerline.vtp")
    filenameSacCenterlines = path.join(dir_path, caseName + "_sac_centerline_{}.vtp")
    fileNameDistanceToSphereDiam = path.join(dir_path, caseName + "_distance_to_sphere_diam.vtp")
    fileNameDistanceToSphereCurv = path.join(dir_path, caseName + "_distance_to_sphere_curv.vtp")
    fileNameProbePoints = path.join(dir_path, caseName + "_probe_point")
    fileNameVoronoi = path.join(dir_path, caseName + "_voronoi.vtp")
    fileNameVoronoiSmooth = path.join(dir_path, caseName + "_voronoi_smooth.vtp")
    fileNameSurfaceSmooth = path.join(dir_path, caseName + "_smooth.vtp")
    fileNameModelFlowExt = path.join(dir_path, caseName + "_flowext.vtp")
    fileNameClippedModel = path.join(dir_path, caseName + "_clippedmodel.vtp")
    fileNameFlowCenterlines = path.join(dir_path, caseName + "_flow_cl.vtp")
    fileNameCappedFlowExt = path.join(dir_path, caseName + "_capped_flowext.vtp")
    fileNameSurfaceName = path.join(dir_path, caseName + "_remeshed_surface.vtp")
    fileNameOutput = path.join(dir_path, caseName + "_output.txt")
    fileNameXMLMesh = path.join(dir_path, caseName + ".xml")
    fileNameVTUMesh = path.join(dir_path, caseName + ".vtu")
    fileNameRunScript = path.join(dir_path, caseName + ".sh")

    print("\n--- Working on case:", caseName, "\n")

    # Open the surface file.
    print("--- Load model file\n")
    surface = ReadPolyData(fileNameModel)

    if not is_surface_capped and smoothing != "voronoi":
        print("--- Clipping the models inlet and outlets.\n")
        if not path.isfile(fileNameClippedModel):
                # TODO: Check if this is a valid call to this method
                centerline = compute_centerlines([], [], None, surface,
                                                 method="pickpoint")
                surface = uncapp_surface(surface, centerline, filename=None, clipspheres=0)
        else:
            surface = ReadPolyData(fileNameClippedModel)

    parameters = get_parameters(dir_path)

    if not "check_surface" in parameters.keys():
        surface = surface_cleaner(surface)
        surface = triangulate_surface(surface)

        # Check the mesh if there is redundant nodes or NaN triangles.
        ToolRepairSTL.surfaceOverview(surface)
        ToolRepairSTL.foundAndDeleteNaNTriangles(surface)
        surface = ToolRepairSTL.cleanTheSurface(surface)
        foundNaN = ToolRepairSTL.foundAndDeleteNaNTriangles(surface)
        if foundNaN == True:
            raise RuntimeError(("There is an issue with the surface. "
                                "Nan coordinates or some other shenanigans."))
        else:
            parameters["check_surface"] = True
            write_parameters(parameters, dir_path)

    # Capp surface if open
    if not compute_centers(surface, test_capped=True):
        capped_surface = capp_surface(surface)
    else:
        capped_surface = surface

    # Get centerlines
    print("--- Get centerlines\n")
    inlet, outlets = get_centers(surface, dir_path)
    centerlines = compute_centerlines(inlet, outlets, fileNameCenterlines,
                                          capped_surface, resampling=0.1, endPoint=0)
    tol = get_tolerance(centerlines)

    if aneurysm:
        aneurysms = get_aneurysm_dome(capped_surface, dir_path)
        centerlineAnu = compute_centerlines(inlet, aneurysms, filenameAneurysmCenterlines,
                                  capped_surface, resampling=0.1)

        # Extract the aneurysm centerline
        sac_centerline = []
        info = get_parameters(dir_path)
        num_anu = info["number_of_aneurysms"]

        # Compute mean distance between points
        for i in range(num_anu):
            if not path.isfile(filenameSacCenterlines.format(i)):
                line = ExtractSingleLine(centerlineAnu, i)
                locator = get_locator(centerlines)
                for j in range(line.GetNumberOfPoints() -1, 0, -1):
                    point = line.GetPoints().GetPoint(j)
                    ID = locator.FindClosestPoint(point)
                    tmp_point = centerlines.GetPoints().GetPoint(ID)
                    dist = np.sqrt(np.sum((np.asarray(point) - np.asarray(tmp_point))**2))
                    if dist <= tol:
                        break

                tmp = ExtractSingleLine(line, 0, startID=j)
                WritePolyData(tmp, filenameSacCenterlines.format(i))

                # List of VtkPolyData sac(s) centerline
                sac_centerline.append(tmp)

            else:
                sac_centerline.append(ReadPolyData(filenameSacCenterlines.format(i)))

    else:
        num_anu = 0

    # Get 'center' and 'radius' of the aneurysm(s)
    sac_center = []
    misr_max = []
    #misr_max_center = []

    if aneurysm:
        # Merge the sac centerliens
        sac_centerlines = merge_data(sac_centerline)

        for sac in sac_centerline:
            sac_center.append(sac.GetPoints().GetPoint(sac.GetNumberOfPoints() // 2))
            tmp_misr = get_array(radiusArrayName, sac)
            misr_max.append(tmp_misr.max())
            #misr_max_center.append(sac.GetPoint(np.argsort(tmp_misr)[0][-1]))

    # Smooth surface
    if smoothing == "voronoi":
        print("--- Smooth surface: Voronoi smoothing\n")
        if not path.isfile(fileNameSurfaceSmooth):
            # Get Voronoi diagram
            if not path.isfile(fileNameVoronoi):
                voronoi = makeVoronoiDiagram(surface, fileNameVoronoi)
                WritePolyData(voronoi, fileNameVoronoi)
            else:
                voronoi = ReadPolyData(fileNameVoronoi)

            # Get smooth Voronoi diagram
            if not path.isfile(fileNameVoronoiSmooth):
                if aneurysm:
                    smooth_voronoi = SmoothVoronoiDiagram(voronoi, centerlines,
                                                          smoothing_factor,
                                                          sac_centerlines)
                else:
                    smooth_voronoi = SmoothVoronoiDiagram(voronoi, centerlines,
                                                          smoothing_factor)

                WritePolyData(smooth_voronoi, fileNameVoronoiSmooth)
            else:
                smooth_voronoi = ReadPolyData(fileNameVoronoiSmooth)

            # Envelope the smooth surface
            surface = create_new_surface(smooth_voronoi)

            # Uncapp the surface
            surface_uncapped = uncapp_surface(surface, centerlines, filename=None)

            # Check if there has been added new outlets
            num_outlets = centerlines.GetNumberOfLines()
            num_outlets_after = compute_centers(surface_uncapped, test_capped=True)[1]

            if num_outlets != num_outlets_after:
                surface = vmtkSmoother(surface, "laplace", iterations=200)
                WritePolyData(surface, fileNameSurfaceSmooth)
                print(("ERROR: Automatic clipping failed. You have to open {} and " + \
                        "manually clipp the branch which still is capped. " + \
                        "Overwrite the current {} and restart the script.").format(
                        fileNameSurfaceSmooth, fileNameSurfaceSmooth))
                sys.exit(0)

                #for i in range(num_outlets):
                #    lines = [ExtractSingleLine(centerlines, j) for j in range(num_outlets) if j != i]
                #    tmp_centerlines = merge_data(lines)
                #    tmp_uncapped = uncapp_surface(surface, tmp_centerlines, filename=None)

                #    WritePolyData(tmp_uncapped, "tmp_uncapped{}.vtp".format(i))
                #    WritePolyData(tmp_centerlines, "tmp_centerline{}.vtp".format(i))
                #    if num_outlets == compute_centers(tmp_uncapped, test_capped=True)[1]:
                #        surface = vmtkSmoother(tmp_uncapped, "laplace", iterations=200)
                #        WritePolyData(surface, fileNameSurfaceSmooth)
                #        print(("ERROR: Automatic clipping failed. You have to open {} and " + \
                #              "manually clipp the branch which still is capped. " + \
                #              "Overwrite the current {} and restart the script.").format(
                #                  fileNameSurfaceSmooth, fileNameSurfaceSmooth))
                #        break
                #else:
                #    print("ERROR: Something went wrong with the capping. Please overwrite" + \
                #          " the file {} with a manually capped version.".format(fileNameSurfaceSmooth))

            surface = surface_uncapped

            # Smoothing to improve the quality of the elements
            # Consider to add a subdivition here as well.
            surface = vmtkSmoother(surface, "laplace", iterations=200)

            # Write surface
            WritePolyData(surface, fileNameSurfaceSmooth)

        else:
            surface = ReadPolyData(fileNameSurfaceSmooth)


    elif smoothing == "laplace":
        print("--- Smooth surface: Laplacian smoothing\n")
        if not path.isfile(fileNameSurfaceSmooth):
            surface = vmtkSmoother(surface, smoothing)

            # Save the smoothed surface
            WritePolyData(surface, fileNameSurfaceSmooth)

        else:
            surface = ReadPolyData(fileNameSurfaceSmooth)

    elif smoothing == "taubin":
        print("--- Smooth surface: Taubin smoothing\n")
        if not path.isfile(fileNameSurfaceSmooth):
            surface = vmtkSmoother(surface, smoothing)

            # Save the smoothed surface
            WritePolyData(surface, fileNameSurfaceSmooth)

        else:
            surface = ReadPolyData(fileNameSurfaceSmooth)

    elif smoothing == "no_smooth" or None:
        print("--- No smoothing of surface\n")

    # Add flow extensions
    if createFlowext:
        if not path.isfile(fileNameModelFlowExt):
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
            surface = vmtkSmoother(surface, "laplace", iterations=100)
            WritePolyData(surface, fileNameModelFlowExt)

        else:
            surface = ReadPolyData(fileNameModelFlowExt)

    # Capp surface with flow extensions
    capped_surface = capp_surface(surface)

    # Get new centerlines with the flow extensions
    if not path.isfile(fileNameFlowCenterlines):
        print("--- Compute the model centerlines with flow extension.")
        # Compute the centerlines.
        inlet, outlets = get_centers(surface, dir_path, flowext=True)
        centerlines = compute_centerlines(inlet, outlets,
                                             fileNameFlowCenterlines,
                                             capped_surface, resampling = 0.5)

    else:
        centerlines = ReadPolyData(fileNameFlowCenterlines)

    # Choose input for the mesh
    print("--- Computing distance to sphere")
    if meshingMethod == "curvature":
        if not path.isfile(fileNameDistanceToSphereCurv):
            distance_to_sphere = dist_sphere_curv(surface, centerlines,
                                                  sac_center, misr_max,
                                                  fileNameDistanceToSphereCurv,
                                                  coarsening_factor)
        else:
            distance_to_sphere = ReadPolyData(fileNameDistanceToSphereCurv)
    else:
        if not path.isfile(fileNameDistanceToSphereDiam):
            distance_to_sphere = dist_sphere_diam(surface, centerlines,
                                                  sac_center, misr_max,
                                                  fileNameDistanceToSphereDiam,
                                                  coarsening_factor)
        else:
            distance_to_sphere = ReadPolyData(fileNameDistanceToSphereDiam)

    # Compute mesh
    if not path.isfile(fileNameVTUMesh):
        try:
            print("--- Computing mesh\n")
            mesh, remeshSurface = generate_mesh(distance_to_sphere)
            assert remeshSurface.GetNumberOfPoints() > 0, \
                    "No points in surface mesh, try to remesh"
            assert mesh.GetNumberOfPoints() > 0, "No points in mesh, try to remesh"

        except:
            distance_to_sphere = mesh_alternative(distance_to_sphere)
            mesh, remeshSurface = generate_mesh(distance_to_sphere)
            assert mesh.GetNumberOfPoints() > 0, "No points in mesh, after remeshing"
            assert remeshSurface.GetNumberOfPoints() > 0, \
                    "No points in surface mesh, try to remesh"

        # Write mesh in VTU format
        WritePolyData(remeshSurface, fileNameSurfaceName)
        WritePolyData(mesh, fileNameVTUMesh)

        # Write mesh to FEniCS to format
        meshWriter = vmtkscripts.vmtkMeshWriter()
        meshWriter.CellEntityIdsArrayName = "CellEntityIds"
        meshWriter.Mesh = mesh
        meshWriter.OutputFileName = fileNameXMLMesh
        meshWriter.Execute()
        polyDataVolMesh = mesh

    else:
        polyDataVolMesh = ReadPolyData(fileNameVTUMesh)

    # Set the network object used in the scripts for 
    # boundary conditions and probes.
    network = ImportData.Network()
    centerlinesBranches = ImportData.SetNetworkStructure(centerlines, network, verboseprint,
                                                         isConnectivityNeeded=True)

    if not path.isfile(fileNameProbePoints):
        # Get the list of coordinates for the probe points along the network centerline.
        listProbePoints = ImportData.GetListProbePoints(centerlinesBranches, network, verboseprint)
        listProbePoints += sac_center

        # Add points randomly in the sac.
        # FIXME: This is not robust enough. Suggestion to fix: Extract the
        # second half of the sac centerline, then get all points from the
        # voronoi diagram which is closest to that part compared to any ther
        # centerlines. Then randomly chose among those points. For now, simply
        # add just one point (sac_center).
        #embed()
        #numberOfPoints = numberOfSacPoints
        #for k in range(num_anu):
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

        print("--- Saving probes points in: ", fileNameProbePoints)
        dataNumpy = np.array(listProbePoints)
        dataNumpy.dump(fileNameProbePoints)
    else:
        dataNumpy = np.load(fileNameProbePoints)

    # Set the flow split and inlet boundary condition
    # Compute the outlet boundary condition percentages.
    flowSplitting = FlowSplitting()
    flowSplitting.ComputeAlphas(network, verboseprint)
    flowSplitting.ComputeBetas(network, verboseprint)
    flowSplitting.CheckTotalFlowRate(network, verboseprint)

    # BSL method for mean inlet flow rate.
    parameters = get_parameters(dir_path)
    meanInflow = 0.27 * parameters["inlet_area"]

    # Extract the surface mesh of the wall
    wallMesh = threshold(polyDataVolMesh, "CellEntityIds", lower=0.5, upper=1.5,
                         type="between", source=1)

    boundaryReferenceSystems = vmtkscripts.vmtkBoundaryReferenceSystems()
    boundaryReferenceSystems.Surface = wallMesh
    boundaryReferenceSystems.Execute()
    refSystem = boundaryReferenceSystems.ReferenceSystems
    cellEntityIdsArray = get_vtk_array('CellEntityIds', 0, refSystem.GetNumberOfPoints())
    refSystem.GetPointData().AddArray(cellEntityIdsArray)

    # Extract the surface mesh of the end caps
    boundarySurface = threshold(polyDataVolMesh, "CellEntityIds", upper=1.5,
                                type="upper", source=1)

    pointCells = vtk.vtkIdList()
    surfaceCellEntityIdsArray = vtk.vtkIntArray()
    surfaceCellEntityIdsArray.DeepCopy(boundarySurface.GetCellData().GetArray('CellEntityIds'))

    # Find the corresponding couple (mesh outlet ID, network ID).
    ids = []
    for i in range(refSystem.GetNumberOfPoints()):
        distancePoints = 10000000
        pointId = boundarySurface.FindPoint(refSystem.GetPoint(i))
        boundarySurface.GetPointCells(pointId,pointCells)
        cellId = pointCells.GetId(0)
        cellEntityId = surfaceCellEntityIdsArray.GetValue(cellId)
        cellEntityIdsArray.SetValue(i,cellEntityId)

        meshPoint = refSystem.GetPoint(i)
        for element in network.elements:
            if element.IsAnOutlet():
                networkPoint = element.GetOutPointsx1()[0]
            if element.IsAnInlet():
                networkPoint = element.GetInPointsx0()[0]
            if vtk.vtkMath.Distance2BetweenPoints(meshPoint,networkPoint) < distancePoints:
                distancePoints = vtk.vtkMath.Distance2BetweenPoints(meshPoint,networkPoint)
                closest = element.GetId()
        if network.elements[closest].IsAnInlet():
            verboseprint('I am the inlet, Sup?')
            verboseprint(network.elements[closest].GetInPointsx0()[0])
            ids.insert(0, [cellEntityId, meanInflow])
        else:
            beta = network.elements[closest].GetBeta()
            ids.append([cellEntityId, beta])
            verboseprint(beta)
            verboseprint(network.elements[closest].GetOutPointsx1()[0])
        verboseprint('CellEntityId: %d\n' % cellEntityId)
        verboseprint('meshPoint: %f, %f, %f\n' % (meshPoint[0],meshPoint[1],meshPoint[2]))
        verboseprint(ids)

    # Store information for the solver.
    idFileLine = caseName + ' ' + repr(ids[0][0] - 1) + ' '
    areaRatioLine = caseName + ' '
    for k in range(1, refSystem.GetNumberOfPoints() - 1):
        idFileLine += repr(ids[k][0] - 1) + ','
        areaRatioLine += repr(ids[k][1]) + ','
    idFileLine += repr(ids[-1][0] - 1) + ' ' + repr(ids[0][1])
    areaRatioLine += repr(ids[-1][1])

    with open(path.join(dir_path, caseName + '.txt'), 'w') as outfile:
        outfile.write("\ninlet_area: " + str(parameters["inlet_area"]))
        outfile.write('\nidFileLine: ' + str(idFileLine))
        outfile.write('\nareaRatioLine: ' + str(areaRatioLine))

    # Display the flow split at the outlets, inlet flowrate, and probes.
    if viz:
        visualize(network.elements, dataNumpy, surface, meanInflow)

    # Start simulation though ssh, without password
    if configPath is not None:

        # Set up simulation script
        if not path.exists(fileNameRunScript):
            run_script_sample = open(path.join(abs_path, "run_script.sh"), "r").read()
            config = json.load(open(configPath))
            run_dict = dict(mesh_name=caseName,
                            num_nodes=1,
                            hours=120,
                            account="nn9249k",
                            remoteFolder=config["remoteFolder"],
                            results_folder="results")
            run_script = run_script_sample.format(**run_dict)

            # Write script
            script_file = open(fileNameRunScript, "w")
            script_file.write(run_script)
            script_file.close()

        run_simulation(configPath, dir_path, caseName)



if __name__ == "__main__":
    '''Command-line arguments.'''
    parser = argparse.ArgumentParser(
                description="Autom. Pre-processing for FEniCS.")

    parser.add_argument('-v', '--verbosity',
        dest = 'verbosity',
        type = str2bool,
        default = False,
        help = "Activates the verbose mode.")

    parser.add_argument('-i', '--inputModel',
        type = str,
        required = False,
        dest = 'fileNameModel',
        default = 'example/surface.vtp',
        help = "Input file containing the 3D model.")

    parser.add_argument('-sM', '--smoothingMethod',
        type = str,
        required = False,
        dest = 'smoothingMethod',
        default = "no_smooth",
        choices = ["voronoi", "no_smooth", "laplace", "taubin"],
        help = "Smoothing method, for now only Voronoi smoothing is avaleble." +\
               " For Voronoi smoothing you can also controll smoothingFactor" + \
               " (default = 0.25)  and smoothingAneurysm (default = False).")

    parser.add_argument('-c', '--coarseningFactor',
        type = float,
        required = False,
        dest = 'coarseningFactor',
        default = 1.0,
        help = "Refine or coarsen the standard mesh size")

    parser.add_argument('-sF', '--smoothingFactor',
        type = float,
        required = False,
        dest = 'smoothingFactor',
        default = 0.25,
        help = "smoothingFactor for VoronoiSmoothing, removes all spheres which" + \
               " has a radious < MISR*(1-0.25), where MISR varying along the centerline.")

    parser.add_argument('-sA', '--smoothAneurysm',
        type = str2bool,
        required = False,
        dest = "smoothingAneurysm",
        help = "When using Voronoi smoothing one can choose not to smooth the" + \
               " aneurysm as this method often removes to much of the aneurysm")

    parser.add_argument('-m', '--meshingMethod',
        dest = "meshingMethod",
        type = str,
        required = True,
        choices = ["diameter", "curvature"],
        default = False)

    parser.add_argument('-a', '--aneurysm',
        dest = "aneu",
        type = str2bool,
        default = True,
        help = "Determine wether or not the model has a aneurysm. Default is False.")

    parser.add_argument('-f', '--flowext',
        dest = "fext",
        default = True,
        type = str2bool,
        help = "Add flow extensions to to the model.")

    parser.add_argument('-vz', '--visualize',
        dest = "viz",
        default = True,
        type = str2bool,
        help = "Visualize surface, inlet, outlet and propes after meshing.")

    parser.add_argument('-sp', '--sacpoints',
        type = int,
        help = 'Number of sac points to add',
        default = 20,
        dest = "sacpts")

    parser.add_argument('--simulationConfig',
        type = str,
        dest = "config",
        default = None,
        help = 'Path to configuration file for remote simulation. ' + \
               'See example/ssh_config.json for details')

    args = parser.parse_args()

    if args.verbosity:
        print()
        print("--- VERBOSE MODE ACTIVATED ---")
        def verboseprint(*args):
            for arg in args:
                print(arg, end=' ')
                print()
    else:
        verboseprint = lambda *a: None

    verboseprint(args)

    # Start the script.
    print(args.config)
    Program(args.fileNameModel, verboseprint, args.smoothingMethod,
            args.smoothingFactor, args.smoothingAneurysm, args.meshingMethod,
            args.aneu, args.fext, args.viz, args.config, args.sacpts,
            args.coarseningFactor)
