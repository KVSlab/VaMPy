import vtk
from vmtk import vtkvmtk, vmtkscripts

try:
    from vmtkpointselector import *
except ImportError:
    pass
import numpy as np
import sys
import re
from os import path, listdir
import math

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


def ReadPolyData(filename):
    """
    Load the given file, and return a vtkPolyData object for it.

    Args:
        filename (str): Path to input file.

    Returns:
        polyData (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Output data.
    """

    # Check if file exists
    if not path.exists(filename):
        raise RuntimeError("Could not find file: %s" % filename)

    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get reader
    if fileType == 'stl':
        reader = vtk.vtkSTLReader()
        reader.MergingOn()
    elif fileType == 'vtk':
        reader = vtk.vtkPolyDataReader()
    elif fileType == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif fileType == 'vts':
        reader = vtk.vtkXMLStructuredGridReader()
    elif fileType == 'vtr':
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fileType == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fileType == "vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fileType == "tec":
        polyData = ReadTecplotSurfaceFile(filename)
        return polyData
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Read
    reader.SetFileName(filename)
    reader.Update()
    polyData = reader.GetOutput()

    return polyData


def WritePolyData(input_data, filename):
    """
    Write the given input data based on the file name extension.

    Args:
        input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Input data.
        filename (str): Save path location.
    """

    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get writer
    if fileType == 'stl':
        writer = vtk.vtkSTLWriter()
    elif fileType == 'vtk':
        writer = vtk.vtkPolyDataWriter()
    elif fileType == 'vts':
        writer = vtk.vtkXMLStructuredGridWriter()
    elif fileType == 'vtr':
        writer = vtk.vtkXMLRectilinearGridWriter()
    elif fileType == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif fileType == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif fileType == "vti":
        writer = vtk.vtkXMLImageDataWriter()
    elif fileType == "tec":
        WriteTecplotSurfaceFile(input_data, filename)
        return
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Set filename and input
    writer.SetFileName(filename)
    if version < 6:
        writer.SetInput(input_data)
    else:
        writer.SetInputData(input_data)
    writer.Update()

    # Write
    writer.Write()


def ReadTecplotSurfaceFile(filename):
    f = open(filename, 'r')
    line = f.readline()

    # Read title
    if line.split()[0] == 'TITLE':
        line = f.readline()

    # Read array names
    if line.split()[0] == 'VARIABLES' or line.split('=')[0] == 'VARIABLES':
        arrayNames = line.split('=')[1].strip().split(',')
        arrayNames[0:3] = []
        line = f.readline()

    # Read ZONE variables
    if line.split()[0] == 'ZONE':
        lineNid = line.find('N=')
        lineN = line[lineNid:lineNid + line[lineNid:].find(',')].split('=')[1]
        numberOfNodes = int(lineN)
        lineEid = line.find('E=')
        lineE = line[lineEid:lineEid + line[lineEid:].find(',')].split('=')[1]
        numberOfElements = int(lineE)
        elementType = 'TRIANGLE'
        if line.find('ET=') != -1:
            if 'TRIANGLE' in line:
                elementType = 'TRIANGLE'
            elif 'QUADRILATERAL' in line:
                elementType = 'QUADRILATERAL'

    # Initialize vtk objects
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    points.SetNumberOfPoints(numberOfNodes)
    surface = vtk.vtkPolyData()
    surface.SetPoints(points)
    surface.SetPolys(cells)

    # Fill array object
    for arrayName in arrayNames:
        array = vtk.vtkDoubleArray()
        array.SetName(arrayName)
        array.SetNumberOfTuples(numberOfNodes)
        surface.GetPointData().AddArray(array)

    # Get rest of data
    data = f.read().split()

    # Insert points
    dataCounter = 0
    for i in range(numberOfNodes):
        point = [float(data[dataCounter]), float(data[dataCounter + 1]), float(data[dataCounter + 2])]
        dataCounter += 3
        points.SetPoint(i, point)
        for j in range(len(arrayNames)):
            surface.GetPointData().GetArray(arrayNames[j]).SetComponent(i, 0, float(data[dataCounter]))
            dataCounter += 1

    # TODO: Is this necessary? It is not inserted into the surface
    # Insert ids
    cellIds = vtk.vtkIdList()
    for i in range(numberOfElements):
        cellIds.Initialize()
        cellIds.InsertNextId(int(data[dataCounter]) - 1)
        dataCounter += 1
        cellIds.InsertNextId(int(data[dataCounter]) - 1)
        dataCounter += 1
        cellIds.InsertNextId(int(data[dataCounter]) - 1)
        dataCounter += 1

        if elementType == "QUADRILATERAL":
            cellIds.InsertNextId(int(data[dataCounter]) - 1)
            dataCounter += 1

        cells.InsertNextCell(cellIds)

    return surface


def WriteTecplotSurfaceFile(surface, filename):
    # Triangulate surface
    surface = triangulate_surface(surface)

    # Open file
    f = open(filename, 'w')
    line = "VARIABLES = X,Y,Z"

    # Get array
    arrayNames = []
    for i in range(surface.GetPointData().GetNumberOfArrays()):
        array = surface.GetPointData().GetArray(i)
        arrayName = array.GetName()

        if arrayName is None:
            continue
        if arrayName[-1] == '_':
            continue

        arrayNames.append(arrayName)

        if array.GetNumberOfComponents() == 1:
            line = line + ',' + arrayName
        else:
            for j in range(array.GetNumberOfComponents()):
                line = line + ',' + arrayName + str(j)

    line = line + '\n'
    f.write(line)

    # Write header
    line = "ZONE N=%s,E=%s,F=FEPOINT,ET=TRIANGLE\n" % \
           (str(surface.GetNumberOfPoints()),
            str(surface.GetNumberOfCells()))
    f.write(line)

    # Write coordinates and array
    for i in range(surface.GetNumberOfPoints()):
        point = surface.GetPoint(i)
        line = " ".join(point)
        for arrayName in arrayNames:
            array = surface.GetPointData().GetArray(arrayName)
            for j in range(array.GetNumberOfComponents()):
                line = line + ' ' + str(array.GetComponent(i, j))
        line = line + '\n'
        f.write(line)

    # Write connectivity and ids
    for i in range(surface.GetNumberOfCells()):
        cellPointIds = surface.GetCell(i).GetPointIds()
        line = ''
        for j in range(cellPointIds.GetNumberOfIds()):
            if j > 0:
                line = line + ' '
            line = line + str(cellPointIds.GetId(j) + 1)
        line = line + '\n'
        f.write(line)


def uncapp_surface(surface, centerlines, filename, clipspheres=0):
    """
    A method for removing endcapps on a surface. The method considers the centerline
    of the surface model and clips the endcaps a distance of clipspheres * MISR away from the edge.

    Args:
        surface (vtkPolyData): Surface to be uncapped.
        centerlines (vtkPolyData): Centerlines in surface model
        filename (str): Filename to write uncapped model to
        clipspheres (int): Number of MISR to clip model

    Returns:
        surface (vtkPolyData): The uncapped surface.
    """
    extractor = vmtkscripts.vmtkEndpointExtractor()
    extractor.Centerlines = centerlines
    extractor.RadiusArrayName = radiusArrayName
    extractor.GroupIdsArrayName = groupIDsArrayName
    extractor.BlankingArrayName = branchClippingArrayName
    extractor.NumberOfEndPointSpheres = clipspheres
    extractor.Execute()
    clipped_centerlines = extractor.Centerlines

    clipper = vmtkscripts.vmtkBranchClipper()
    clipper.Surface = surface
    clipper.Centerlines = clipped_centerlines
    clipper.RadiusArrayName = radiusArrayName
    clipper.GroupIdsArrayName = groupIDsArrayName
    clipper.BlankingArrayName = branchClippingArrayName
    clipper.Execute()
    surface = clipper.Surface

    connector = vmtkscripts.vmtkSurfaceConnectivity()
    connector.Surface = surface
    connector.CleanOutput = 1
    connector.Execute()
    surface = connector.Surface

    if filename is not None:
        WritePolyData(surface, filename)

    return surface


def get_teams(dirpath):
    files = [f for f in listdir(dirpath) if path.isdir(path.join(dirpath, f))]
    teams = []
    for folder in files:
        try:
            teams.append(path.join(dirpath, folder))
        except:
            pass

    teams.sort(key=lambda x: int(x.split(path.sep)[-1]))
    return teams


def makeVoronoiDiagram(surface, file_path):
    """
    Compute the voronoi diagram of surface model.

    Args:
        surface (polydata): Capped surface model to create a Voronoi diagram of.
        file_path (str): Absolute path to surface model path.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    if path.isfile(file_path):
        return ReadPolyData(file_path)

    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = 0
    voronoi.Execute()

    WritePolyData(voronoi.VoronoiDiagram, file_path)

    return voronoi.VoronoiDiagram


def write_spheres(points, dirpath, radius=None, name="sphere%s.vtp", base=0.2):
    radius = [base] * len(points) if radius is None else radius
    for counter, point in enumerate(points):
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(point)
        sphere.SetPhiResolution(100)
        sphere.SetThetaResolution(100)
        sphere.SetRadius(radius[counter])
        sphere_ = sphere.GetOutput()

        WritePolyData(sphere_, path.join(dirpath, name % counter))


def get_tolerance(centerline):
    """
    Finds tolerance based on average length between first N points along the input centerline.

    Args:
        centerline (vtkPolyData): Centerline data.

    Returns:
        tolerance (float): Tolerance value.
    """
    line = ExtractSingleLine(centerline, 0)
    length = get_curvilinear_coordinate(line)

    return np.mean(length[1:] - length[:-1]) / divergingRatioToSpacingTolerance


def getRelevantOutlets(aneurysm_type, centerline, centerline_aneurysm, dir_path):
    """
    Extract relevant outlets of the
    input surface model.

    Args:
        aneurysm_type (bool): Type of aneurysm, saccular of terminal
        centerline_aneurysm (vtkPolyData): Centerline into aneurysm
        centerline (vtkPolyData): Surface model centerline.
        dir_path (str): Location of info-file.

    Returns:
        relevant_outlets (list): List of relevant outlet IDs.
    """

    parameters = get_parameters(dir_path)
    if "relevant_outlet_0" in parameters.keys():
        return parameters["relevant_outlet_0"], parameters["relevant_outlet_1"]

    n_aneurysm = centerline_aneurysm.GetNumberOfPoints()
    tol = get_tolerance(centerline)
    div_anu = []
    aneurysm = centerline_aneurysm.GetCell(0)
    get_point = aneurysm.GetPoints().GetPoint

    for j in range(centerline.GetNumberOfCells()):
        line = centerline.GetCell(j)
        get_point_line = line.GetPoints().GetPoint
        for i in range(min(n_aneurysm, line.GetNumberOfPoints())):
            point0 = get_point(i)
            point1 = get_point_line(i)

            dist = distance(point0, point1)
            if dist > tol:
                div_anu.append(i)
                break

    index = np.argsort(np.array(div_anu))
    line = ExtractSingleLine(centerline, index[-1])
    endpoint1 = list(line.GetPoints().GetPoint(line.GetNumberOfPoints() - 1))
    inlet = list(line.GetPoints().GetPoint(0))

    if aneurysm_type == "terminal":
        diff = True
        line_id = -2
        line_div_point = line.GetPoints().GetPoint(div_anu[index[1]] + 10)
        while diff:
            line2 = ExtractSingleLine(centerline, index[line_id])
            line2_div_point = line2.GetPoints().GetPoint(div_anu[index[line_id]])
            if distance(line2_div_point, line_div_point) > 3 * tol:
                diff = False
            endpoint2 = list(line2.GetPoints().GetPoint(line2.GetNumberOfPoints() - 1))
            line_id -= 1
    else:
        endpoint2 = endpoint1
        endpoint1 = inlet

    data = {}
    for i, e in enumerate([endpoint1, endpoint2]):
        data["relevant_outlet_%d" % i] = e

    write_parameters(data, dir_path)

    return endpoint1, endpoint2


def SmoothVoronoiDiagram(voronoi, centerlines, smoothing_factor, no_smooth_cl=None):
    """
    Smooth voronoi diagram based on a given smoothing factor. Each voronoi point
    that has a radius less then MISR*(1-smoothingFactor) at the closest centerline point is removed.

    Args:
        voronoi (vtkPolyData): Voronoi diagram to be smoothed.
        centerlines (vtkPolyData): Centerline data.
        smoothing_factor (float): Smoothing factor: remove points with radius below (1-smoothing_factor)*MISR
        no_smooth_cl (vktPolyData): Section of centerline not to smooth along.

    Returns:
        smoothedDiagram (vtkPolyData): Smoothed voronoi diagram.
    """

    numberOfPoints = voronoi.GetNumberOfPoints()

    threshold = get_array(radiusArrayName, centerlines) * (1 - smoothing_factor)
    locator = get_locator(centerlines)
    if no_smooth_cl is not None:
        no_locator = get_locator(no_smooth_cl)

    smoothedDiagram = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArrayNumpy = np.zeros(numberOfPoints)

    count = 0
    for i in range(numberOfPoints):
        point = voronoi.GetPoint(i)
        radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        id_ = locator.FindClosestPoint(point)
        cl_point = centerlines.GetPoint(id_)

        if distance(point, cl_point) > 3 * threshold[id_]:
            points.InsertNextPoint(point)
            cellArray.InsertNextCell(1)
            cellArray.InsertCellPoint(count)
            radiusArrayNumpy[count] = radius
            count += 1

        elif no_smooth_cl is not None:
            dist1 = distance(point, centerlines.GetPoint(id_))
            id_1 = no_locator.FindClosestPoint(point)
            dist2 = distance(point, no_smooth_cl.GetPoint(id_1))

            if dist2 < dist1:
                points.InsertNextPoint(point)
                cellArray.InsertNextCell(1)
                cellArray.InsertCellPoint(count)
                radiusArrayNumpy[count] = radius
                count += 1
            else:
                if radius >= threshold[id_]:
                    points.InsertNextPoint(point)
                    cellArray.InsertNextCell(1)
                    cellArray.InsertCellPoint(count)
                    radiusArrayNumpy[count] = radius
                    count += 1
        else:
            if radius >= threshold[id_]:
                points.InsertNextPoint(point)
                cellArray.InsertNextCell(1)
                cellArray.InsertCellPoint(count)
                radiusArrayNumpy[count] = radius
                count += 1

        radiusArray = get_vtk_array(radiusArrayName, 1, count)

    for i in range(count):
        radiusArray.SetTuple1(i, radiusArrayNumpy[i])

    smoothedDiagram.SetPoints(points)
    smoothedDiagram.SetVerts(cellArray)
    smoothedDiagram.GetPointData().AddArray(radiusArray)

    return smoothedDiagram


def get_curvilinear_coordinate(line):
    """
    Get curvilinear coordinates along
    an input centerline.

    Args:
        line (vtkPolyData): Input centerline

    Returns:
        curv_coor (ndarray): Array of abscissa points.
    """

    curv_coor = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints() - 1):
        pnt1 = np.asarray(line.GetPoints().GetPoint(i))
        pnt2 = np.asarray(line.GetPoints().GetPoint(i + 1))
        curv_coor[i + 1] = np.sqrt(np.sum((pnt1 - pnt2) ** 2)) + curv_coor[i]

    return curv_coor


def merge_data(inputs):
    """
    Appends one or more polygonal dataset together into a single
    polygonal dataset.

    Args:
        inputs (list): List of vtkPolyData objects.

    Returns:
        merged_data (vtkPolyData): Single polygonal dataset.
    """
    appendFilter = vtk.vtkAppendPolyData()
    for input_ in inputs:
        appendFilter.AddInputData(input_)
    appendFilter.Update()

    return appendFilter.GetOutput()


def get_array(array_name, line, k=1):
    """
    Get data array from polydata object (PointData).

    Args:
        array_name (str): Name of array.
        line (vtkPolyData): Centerline object.
        k (int): Dimension.

    Returns:
        array (ndarray): Array containing data points.
    """

    array = np.zeros((line.GetNumberOfPoints(), k))
    if k == 1:
        getData = line.GetPointData().GetArray(array_name).GetTuple1
    elif k == 2:
        getData = line.GetPointData().GetArray(array_name).GetTuple2
    elif k == 3:
        getData = line.GetPointData().GetArray(array_name).GetTuple3

    for i in range(line.GetNumberOfPoints()):
        array[i, :] = getData(i)

    return array


def get_array_cell(array_name, line, k=1):
    """
    Get cell data array from polydata object (CellData).

    Args:
        array_name (str): Name of array.
        line (vtkPolyData): Centerline object.
        k (int): Dimension.

    Returns:
        array (ndarray): Array containing data points.
    """
    array = np.zeros((line.GetNumberOfCells(), k))
    if k == 1:
        getData = line.GetCellData().GetArray(array_name).GetTuple1
    elif k == 2:
        getData = line.GetCellData().GetArray(array_name).GetTuple2
    elif k == 3:
        getData = line.GetCellData().GetArray(array_name).GetTuple3
    elif k == 9:
        getData = line.GetCellData().GetArray(array_name).GetTuple9

    for i in range(line.GetNumberOfCells()):
        array[i, :] = getData(i)

    return array


# TODO: Get bounds and compute polyballs based on that
#       bounds = surface.GetBounds()
def create_new_surface(complete_voronoi_diagram, poly_ball_image_size=[280, 280, 280]):
    """
    Envelops an input voronoi diagram into a new surface model at a
    given resolution determined by the poly_ball_size.

    Args:
        complete_voronoi_diagram (vtkPolyData): Voronoi diagram
        poly_ball_image_size (list): List of dimensional resolution of output model

    Returns:
        envelope (vtkPolyData): Enveloped surface model.
    """

    # Get the x,y, and z range of the completeVoronoiDiagram
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    if version < 6:
        modeller.SetInput(complete_voronoi_diagram)
    else:
        modeller.SetInputData(complete_voronoi_diagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(poly_ball_image_size)
    modeller.Update()

    # Write the new surface
    marchingCube = vtk.vtkMarchingCubes()
    if version < 6:
        marchingCube.SetInput(modeller.GetOutput())
    else:
        marchingCube.SetInputData(modeller.GetOutput())
    marchingCube.SetValue(0, 0.0)
    marchingCube.Update()
    envelope = marchingCube.GetOutput()

    return envelope


def get_aneurysm_dome(surface, dir_path):
    # Check if info exists
    if not path.isfile(path.join(dir_path, dir_path + ".txt")):
        provide_aneurysm_points(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    dome = []
    for key, value in parameters.items():
        if key.startswith("aneurysm_"):
            dome.append(value)

    if dome == []:
        dome = provide_aneurysm_points(surface, dir_path)

    # Flatten list
    return [item for sublist in dome for item in sublist]


def centerline_div(centerline1, centerline2, tol):
    # Find clipping points
    N = min(centerline1.GetNumberOfPoints(), centerline2.GetNumberOfPoints())
    get_point1 = centerline1.GetPoints().GetPoint
    get_point2 = centerline2.GetPoints().GetPoint

    for i in range(0, N):
        distance_between_points = distance(get_point1(i), get_point2(i))
        if distance_between_points > tol:
            break

    return i


def provide_aneurysm_points(surface, dir_path=None):
    """
    Get relevant aneurysm points from user selected points on a input surface.

    Args:
        surface (vtkPolyData): Surface model.
        dir_path (str): Location of info.json file

    Returns:
        points (list): List of relevant outlet IDs
    """
    # Fix surface
    cleaned_surface = surface_cleaner(surface)
    triangulated_surface = triangulate_surface(cleaned_surface)

    # Select seeds
    SeedSelector = vmtkPickPointSeedSelector()
    SeedSelector.SetSurface(triangulated_surface)
    SeedSelector.Execute()

    aneurysmSeedIds = SeedSelector.GetTargetSeedIds()
    get_point = surface.GetPoints().GetPoint
    points = [list(get_point(aneurysmSeedIds.GetId(i))) for i in range(aneurysmSeedIds.GetNumberOfIds())]

    if dir_path is not None:
        info = {"number_of_aneurysms": len(points)}

        for i in range(len(points)):
            info["aneurysm_%d" % i] = points[i]
            write_parameters(info, dir_path)

    return points


def getData(centerline_par, centerline_dau1, centerline_dau2, tol, aneurysm_type):
    """
    Get info about bifurcating points and diverging points within a bifurcation.
    End points are set based on the MISR at the selected points.

    Args:
        aneurysm_type (str): Type of aneurysm, saccular or terminal
        centerline_dau1 (vtkPolyData): Centerline from inlet to relevant outlet 1.
        centerline_dau2 (vtkPolyData): Centerline from inlet to relevant outlet 2.
        centerline_par(vtkPolyData): Centerline through bifurcation.
        tol (float): Tolerance parameter.

    Returns:
        data (dict): Contains info about diverging point locations.
    """
    # Declare variables before loop if values are not found
    data = {"dau1": {}, "dau2": {}}

    # List of points conected to ID
    points_ids_0 = vtk.vtkIdList()
    points_ids_1 = vtk.vtkIdList()

    data_list = [("dau1", centerline_dau1), ("dau2", centerline_dau2)]

    # Additional point for terminal
    if aneurysm_type == "terminal":
        data_list += [("par", centerline_par)]
        data["par"] = {}

    for key, centerline in data_list:
        if centerline is None: continue
        # One is the branch to the left and the other is the one to the right
        centerline.GetCellPoints(0, points_ids_0)
        centerline.GetCellPoints(1, points_ids_1)

        # Find clipping points
        N = min(points_ids_0.GetNumberOfIds(), points_ids_1.GetNumberOfIds())
        for i in range(0, N):
            cell_point_0 = centerline.GetPoint(points_ids_0.GetId(i))
            cell_point_1 = centerline.GetPoint(points_ids_1.GetId(i))

            distance_between_points = distance(cell_point_0, cell_point_1)
            if distance_between_points > tol:
                tmpI = i
                point_ID_0 = points_ids_0.GetId(i)
                point_ID_1 = points_ids_1.GetId(i)
                center = centerline.GetPoint(point_ID_0)
                r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID_0)
                break

        end, r_end = move_past_sphere(centerline, center, r, point_ID_0,
                                      stop=point_ID_0 * 100, step=1, X=1)

        data[key]["end_point"] = end
        data[key]["r_end"] = r_end
        data[key]["div_point"] = center
        data[key]["ID_div"] = point_ID_0
        data[key]["i_div"] = tmpI
        data[key]["r_div"] = r

    return data


def write_points(points, filename):
    """
    Writes input points to file.

    Args:
        points (vtkPolyData): Point data.
        filename (str): Save location.
    """
    pointSet = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()

    for i in range(points.GetNumberOfPoints()):
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    pointSet.SetPoints(points)
    pointSet.SetVerts(cellArray)

    WritePolyData(pointSet, filename)


def surface_cleaner(surface):
    """
    Clean surface by merging duplicate points, and/or
    removing unused points and/or removing degenerate cells.

    Args:
        surface (vtkPolyData): Surface model.

    Returns:
        cleanedSurface (vtkPolyData): Cleaned surface model.
    """

    surfaceCleaner = vtk.vtkCleanPolyData()
    if version < 6:
        surfaceCleaner.SetInput(surface)
    else:
        surfaceCleaner.SetInputData(surface)
    surfaceCleaner.Update()

    return surfaceCleaner.GetOutput()


def get_centers(surface, dir_path, flowext=False):
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
    if flowext or not path.isfile(path.join(dir_path, dir_path + ".txt")):
        compute_centers(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    outlets = []
    inlet = []
    for key, value in parameters.items():
        if key == "inlet":
            inlet = value
        elif "outlet" in key and "area" not in key and "relevant" not in key:
            outlets += value

    num_outlets = len(outlets) // 3
    if num_outlets != 0:
        outlets = []
        for i in range(num_outlets):
            outlets += parameters["outlet%d" % i]

    if inlet == [] and outlets == []:
        inlet, outlets = compute_centers(surface, dir_path)

    return inlet, outlets


def triangulate_surface(surface):
    """
    Wrapper for vtkTriangleFilter.

    Args:
        surface (vtkPolyData): Surface to triangulate.
    Returns:
        surface (vtkPolyData): Triangulated surface.
    """

    surfaceTriangulator = vtk.vtkTriangleFilter()
    if version < 6:
        surfaceTriangulator.SetInput(surface)
    else:
        surfaceTriangulator.SetInputData(surface)
    surfaceTriangulator.PassLinesOff()
    surfaceTriangulator.PassVertsOff()
    surfaceTriangulator.Update()

    return surfaceTriangulator.GetOutput()


def geometryFilter(unstructured_grid):
    # Convert unstructured grid to polydata
    filter = vtk.vtkGeometryFilter()
    if version < 6:
        filter.SetInput(unstructured_grid)
    else:
        filter.SetInputData(unstructured_grid)
    filter.Update()
    polydata = filter.GetOutput()

    return polydata


def threshold(surface, name, lower=0, upper=1, type="between", source=1):
    """
    Wrapper for vtkThreshold. Extract a section of a surface given a criteria.

    Args:
        surface (vtkPolyData): The input data to be extracted.
        name (str): Name of scalar array.
        lower (float): Lower bound.
        upper (float): Upper bound.
        type (str): Type of threshold (lower, upper, between)
        source (int): PointData or CellData.

    Returns:
        surface (vtkPolyData): The extracted surface based on the lower and upper limit.
    """

    # source = 1 uses cell data as input
    # source = 0 uses point data as input

    # Apply threshold
    threshold = vtk.vtkThreshold()
    if version < 6:
        threshold.SetInput(surface)
    else:
        threshold.SetInputData(surface)
    if type == "between":
        threshold.ThresholdBetween(lower, upper)
    elif type == "lower":
        threshold.ThresholdByLower(lower)
    elif type == "upper":
        threshold.ThresholdByUpper(upper)
    else:
        print(("%s is not a threshold type. Pleace chose from: upper, lower" + \
               ", or between") % type)
        sys.exit(0)

    threshold.SetInputArrayToProcess(0, 0, 0, source, name)
    threshold.Update()
    surface = threshold.GetOutput()

    # Convert to polydata
    surface = geometryFilter(surface)

    return surface


def compute_area(surface):
    """
    Calculate the area from the given polydata

    Args:
        surface (vtkPolyData): Surface to compute are off

    Returns:
        area (float): Area of the input surface
    """
    mass = vtk.vtkMassProperties()
    if version < 6:
        mass.SetInput(surface)
    else:
        mass.SetInputData(surface)

    return mass.GetSurfaceArea()


def capp_surface(surface):
    """
    Wrapper for vmtkCapPolyData.
    Close holes in a surface model.

    Args:
        surface (vtkPolyData): Surface to be capped.

    Returns:
        surface (vtkPolyData): Capped surface.
    """

    surfaceCapper = vtkvmtk.vtkvmtkCapPolyData()
    if version < 6:
        surfaceCapper.SetInput(surface)
    else:
        surfaceCapper.SetInputData(surface)
    surfaceCapper.SetDisplacement(0.0)
    surfaceCapper.SetInPlaneDisplacement(0.0)
    surfaceCapper.Update()

    return surfaceCapper.GetOutput()


def compute_distance_to_sphere(surface, centerSphere, radiusSphere=0.0,
                               distanceOffset=0.0, distanceScale=0.01,
                               minDistance=0.2, maxDistance=0.3,
                               distanceToSpheresArrayName="DistanceToSpheres"):
    # Check if there allready exists a distance to spheres
    N = surface.GetNumberOfPoints()
    number, names = get_number_of_arrays(surface)
    add = False
    if distanceToSpheresArrayName not in names: add = True

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
        newDist = distance(centerSphere, surface.GetPoints().GetPoint(i)) - radiusSphere

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


def is_surface_capped(surface):
    """
    Checks if the surface is closed, and how many openings there are.

    Args:
        surface (vtkPolyData): Surface to be checked

    Returns:
        open (boolean): Open or closed surface
        number (int): Number of integer
    """

    return compute_centers(surface, test_capped=True)[0]


def getConnectivity(surface):
    """
    Wrapper of vtkPolyDataConnectivityFilter. Compute connectivity.

    Args:
        surface (vtkPolyData): Input surface data.
    """

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    # Backwards compatibility
    if version < 6:
        connectivity.SetInput(surface)
    else:
        connectivity.SetInputData(surface)

    # Mark each region with "RegionId"
    connectivity.SetExtractionModeToAllRegions()
    connectivity.ColorRegionsOn()
    connectivity.Update()
    output = connectivity.GetOutput()

    return output


def computeCircleness(surface):
    """
    Compute the area ratio between minimum circle and the maximum circle.

    Args:
        surface (vtkPolyData): Boundary edges of an opening

    Returns:
        circleness (float): Area ratio
        center (list): Center of the opening.
    """

    edges = getFeatureEdges(surface)

    # Get points
    points = []
    for i in range(edges.GetNumberOfPoints()):
        points.append(edges.GetPoint(i))

    # Compute center
    points = np.array(points)
    center = np.mean(np.array(points), axis=0)

    # Compute ratio between max inscribed sphere, and min inscribed "area"
    point_radius = np.sqrt(np.sum((points - center) ** 2, axis=1))
    argsort = np.argsort(point_radius)
    if point_radius[argsort[1]] / point_radius[argsort[0]] > 15:
        radius_min = point_radius[argsort[1]]
    else:
        radius_min = point_radius.min()

    min_area = math.pi * radius_min ** 2
    max_area = math.pi * point_radius.max() ** 2

    return max_area / min_area, center


def getFeatureEdges(polyData):
    """
    Wrapper for vtkFeatureedges. Extracts the edges of the cells that are open.

    Args:
        polydata (vtkPolyData): surface to extract the openings from.

    Returns:
        feature_edges (vtkPolyData): The boundary edges of the surface.
    """

    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    if version < 6:
        featureEdges.SetInput(polyData)
    else:
        featureEdges.SetInputData(polyData)
    featureEdges.Update()

    return featureEdges.GetOutput()


def compute_centers(polyData, case_path=None, test_capped=False):
    """
    Compute the center of all the openings in the surface. The inlet is chosen based on
    the largest area.

    Args:
        test_capped (bool): Check if surface is capped
        polyData (vtkPolyData): centers of the openings
        case_path (str): path to case directory.

    Returns:
        inlet_center (list): Inlet center.
        outlet_centers (list): A flattened list with all the outlet centers.
    """

    # Get cells which are open
    cells = getFeatureEdges(polyData)

    if cells.GetNumberOfCells() == 0 and not test_capped:
        print("WARNING: The model is capped, so it is uncapped, but the method is experimental.")
        uncapped_surface = uncapp_surface(polyData)
        compute_centers(uncapped_surface, case_path, test_capped)
    elif cells.GetNumberOfCells() == 0 and test_capped:
        return False, 0

    # Compute connectivity of the cells
    outputs = getConnectivity(cells)

    # Get connectivity array
    region_array = get_array("RegionId", outputs)

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
        boundary = threshold(outputs, "RegionId", lower=i - 0.1, upper=i + 0.1,
                             type="between", source=0)

        delaunay_filter = vtk.vtkDelaunay2D()
        delaunay_filter.SetInputData(boundary)
        delaunay_filter.Update()
        area.append(compute_area(delaunay_filter.GetOutput()))

        # Get center
        center.append(np.mean(points[(region_array == i).nonzero()[0]], axis=0))

    # Store the center and area
    inlet_ind = area.index(max(area))
    if case_path is not None:
        info = {"inlet": center[inlet_ind].tolist(), "inlet_area": area[inlet_ind]}
        p = 0
        for i in range(len(area)):
            if i == inlet_ind: p = -1; continue
            info["outlet%d" % (i + p)] = center[i].tolist()
            info["outlet%s_area" % (i + p)] = area[i]

        write_parameters(info, case_path)

    inlet_center = center[inlet_ind].tolist()
    center.pop(inlet_ind)

    center_ = [item for sublist in center for item in sublist]

    return inlet_center, center_


def compute_bary_center(points):
    # Get i+1
    shifted = np.zeros(points.shape)
    shifted[1:, :] = points[:-1, :]
    shifted[0, :] = points[-1, :]

    # Compute weights
    weight = np.sqrt(np.sum((points - shifted) ** 2, axis=1))
    weight_sum = np.sum(weight)

    # Compute center
    center_x = np.sum((points[:, 0] + shifted[:, 0]) / 2 * weight) / weight_sum
    center_y = np.sum((points[:, 1] + shifted[:, 1]) / 2 * weight) / weight_sum
    center_z = np.sum((points[:, 2] + shifted[:, 2]) / 2 * weight) / weight_sum

    return [center_x, center_y, center_z]


def get_vtk_array(name, comp, num):
    """An empty vtkDoubleArray.

    Args:
        name (str): Name of array.
        comp (int): Number of components
        num (int): Number of tuples.

    Returns:
        array (vtkDoubleArray): An empty vtk array.
    """
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(comp)
    array.SetNumberOfTuples(num)
    for i in range(comp):
        array.FillComponent(i, 0.0)
    array.SetName(name)

    return array


def get_locator_cell(surface):
    """Wrapper for vtkCellLocator

    Args:
        surface (vtkPolyData): input surface

    Returns:
        return (vtkCellLocator): Cell locator of the input surface.
    """
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()

    return locator


def get_locator(centerline):
    """Wrapper for vtkStaticPointLocator.

    Args:
        centerline (vtkPolyData): Input vtkPolyData.

    Returns:
        locator (vtkStaticPointLocator): Point locator of the input surface.
    """
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(centerline)
    locator.BuildLocator()

    return locator


def distance(point1, point2):
    """Distance between two points.

    Args:
        point1 (ndarray): A point
        point2 (ndarray): A point

    Returns:
        distance (float): Distance between point1 and point2
    """

    return np.sqrt(np.sum((np.asarray(point1) - np.asarray(point2)) ** 2))


def remove_distant_points(voronoi, centerline):
    """Take a voronoi diagram and a centerline remove points that are far away.

    Args:
        voronoi (vtkPolyData): Voronoi data.
        centerline (vtkPolyData): centerline.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram without the extreme points
    """

    N = voronoi.GetNumberOfPoints()
    newVoronoi = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radius = np.zeros(N)

    locator = get_locator(centerline)
    get_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1
    limit = get_data(0)
    limit = limit * 10

    count = 0
    for i in range(N):
        point = voronoi.GetPoint(i)
        ID = locator.FindClosestPoint(point)
        cl_point = centerline.GetPoint(ID)
        dist = distance(point, cl_point)
        if dist / 3 > get_data(i) or get_data(i) > limit:
            count += 1
            continue

        points.InsertNextPoint(point)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i - count)
        value = get_data(i)
        radius[i - count] = value

    print("Removed %s points from the voronoi diagram" % count)

    radiusArray = get_vtk_array(radiusArrayName, 1, N - count)
    for i in range(N - count):
        radiusArray.SetTuple(i, [float(radius[i])])

    newVoronoi.SetPoints(points)
    newVoronoi.SetVerts(cellArray)
    newVoronoi.GetPointData().AddArray(radiusArray)

    return newVoronoi


def success(text):
    if not "error: " in text.lower():
        return True, ""
    else:
        error_message = re.search(r'error: (.*)', text.lower()).groups()[0]
        return False, error_message


def compute_centerline_sections(surface, centerline):
    """
    Wrapper for vmtk centerline sections.

    Args:
        surface (vtkPolyData): Surface to measure area.
        centerline (vtkPolyData): centerline to measure along.

    Returns:
        line (vtkPolyData): centerline with the attributes
        centerline_sections_area (vtkPolyData): sections along the centerline
    """
    centerlineSections = vtkvmtk.vtkvmtkPolyDataCenterlineSections()
    centerlineSections.SetInputData(surface)
    centerlineSections.SetCenterlines(centerline)
    centerlineSections.SetCenterlineSectionAreaArrayName('CenterlineSectionArea')
    centerlineSections.SetCenterlineSectionMinSizeArrayName('CenterlineSectionMinSize')
    centerlineSections.SetCenterlineSectionMaxSizeArrayName('CenterlineSectionMaxSize')
    centerlineSections.SetCenterlineSectionShapeArrayName('CenterlineSectionShape')
    centerlineSections.SetCenterlineSectionClosedArrayName('CenterlineSectionClosed')
    centerlineSections.Update()

    CenterlineSections = centerlineSections.GetOutput()
    line = centerlineSections.GetCenterlines()

    return line, CenterlineSections


def compute_centerlines(inlet, outlet, filepath, surface, resampling=1, smooth=False, num_iter=100, smooth_factor=0.1,
                        end_point=1, method="pointlist"):
    """Wrapper for vmtkcenterlines and vmtkcenterlinesmoothing.

    Args:
        inlet (list): point of the inlet
        outlet (list): flatt list of the outlet points
        filepath (str): path to where to store the centerline
        surface (vtkPolyData): surface to get the centerline from.
        resampling (float): resampling step length.
        smooth (bool): smooth centerline or not.
        num_iter (int): number of iterations in smooth.
        smooth_factor (float): smoothing factor.
        end_point (int): 0 or 1, include end point in centerline.
        method (str): method for setting the inlet and outlet location

    Returns:
        centerline (vtkPolyData): centerline of the surface.
    """

    if path.isfile(filepath):
        return ReadPolyData(filepath)

    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = 'pointlist'
    centerlines.AppendEndPoints = end_point
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet
    centerlines.Execute()
    centerlines = centerlines.Centerlines

    if smooth:
        centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
        centerlineSmoothing.SetInputData(centerlines)
        centerlineSmoothing.SetNumberOfSmoothingIterations(num_iter)
        centerlineSmoothing.SetSmoothingFactor(smooth_factor)
        centerlineSmoothing.Update()

        centerlines = centerlineSmoothing.GetOutput()

    # Save the computed centerline.
    if filepath is not None:
        WritePolyData(centerlines, filepath)

    return centerlines


def CenterlineAttribiutes(line, smooth=False, iterations=300, factor=0.1):
    """
    Wrapper for centerline attributes.

    Args:
        line (vtkPolyData): Line to investigate.
        smooth (bool): Turn on and off smoothing before computing the attributes features.
        iterations (int): Number of iterations.
        factor (float): Smoothing factor.

    Returns:
        line (vtkPolyData): Line with centerline attributes.
    """
    attributes = vmtkscripts.vmtkCenterlineAttributes()
    attributes.Centerlines = line
    attributes.NormalsArrayName = parallelTransportNormalsArrayName
    attributes.AbscissaArrayName = abscissasArrayName
    attributes.Execute()
    centerline = attributes.Centerlines

    geometry = vmtkscripts.vmtkCenterlineGeometry()
    geometry.Centerlines = centerline

    if smooth:
        geometry.LineSmoothing = 1
        geometry.OutputSmoothedLines = smooth
        geometry.SmoothingFactor = factor
        geometry.NumberOfSmoothingIterations = iterations
    geometry.FernetTangentArrayName = "FernetTangent"
    geometry.FernetNormalArrayName = "FernetNormal"
    geometry.FrenetBinormalArrayName = "FernetBiNormal"
    geometry.CurvatureArrayName = "Curvature"
    geometry.TorsionArrayName = "Torsion"
    geometry.TortuosityArrayName = "Tortuosity"
    geometry.Execute()

    return geometry.Centerlines


def generate_mesh(surface):
    # Compute the mesh.
    meshGenerator = vmtkscripts.vmtkMeshGenerator()
    meshGenerator.Surface = surface
    meshGenerator.ElementSizeMode = "edgelengtharray"
    meshGenerator.TargetEdgeLengthArrayName = "Size"
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


def create_vtk_array(values, name, k=1):
    vtkArray = get_vtk_array(name, k, values.shape[0])

    if k == 1:
        for i in range(values.shape[0]):
            vtkArray.SetTuple1(i, values[i])
    elif k == 2:
        for i in range(values.shape[0]):
            vtkArray.SetTuple2(i, values[i, 0], values[i, 1])
    elif k == 3:
        for i in range(values.shape[0]):
            vtkArray.SetTuple3(i, values[i, 0], values[i, 1], values[i, 2])
    elif k == 9:
        for i in range(values.shape[0]):
            vtkArray.SetTuple9(i, values[i, 0], values[i, 1], values[i, 2],
                               values[i, 3], values[i, 4], values[i, 5],
                               values[i, 6], values[i, 7], values[i, 8])

    return vtkArray


def GramSchmidt(V):
    V = 1.0 * V
    U = np.copy(V)

    def proj(u, v):
        return u * np.dot(v, u) / np.dot(u, u)

    for i in range(1, V.shape[1]):
        for j in range(i):
            U[:, i] -= proj(U[:, j], V[:, i])

    # normalize column
    den = (U ** 2).sum(axis=0) ** 0.5
    E = U / den
    return E


def get_parameters(dir_path):
    # If info.txt file, return an empty dict
    if not path.isfile(dir_path + ".txt"): return {}

    # Get text
    f = open(dir_path + ".txt", "r")
    text = f.read()
    f.close()
    text = text.split("\n")

    # Get values
    data = {}
    for par in text:
        if par != "":
            split = par.split(": ")
            if len(split) == 2:
                key, value = split
            else:
                key = split[0]
                value = ": ".join(split[1:])
            try:
                data[key] = eval(value)
            except:
                data[key] = value

    return data


def write_parameters(data, dir_path):
    # Get old parameters
    parameters = get_parameters(dir_path)

    # Add new parameters (can overwrite old as well)
    for key, value in data.items():
        parameters[key] = value

    # Get new text
    text = ["%s: %s" % (k, v) for k, v in parameters.items()]
    text = "\n".join(text)

    # Write text
    f = open(dir_path + ".txt", "w")
    f.write(text)
    f.close()


def data_to_vtkPolyData(data, header, TNB=None, PT=None):
    line = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    cellArray.InsertNextCell(data.shape[0])
    linePoints = vtk.vtkPoints()

    info_array = []
    for i in range(3, data.shape[1]):
        radiusArray = get_vtk_array(header[i], 1, data.shape[0])
        info_array.append(radiusArray)

    if TNB is not None:
        for i in range(3):
            radiusArray = get_vtk_array(header[i + data.shape[1]], 3, data.shape[0])
            info_array.append(radiusArray)

    if PT is not None:
        start = data.shape[1] if TNB is None else data.shape[1] + 3
        for i in range(2):
            radiusArray = get_vtk_array(header[i + start], 3, PT[0].shape[0])
            info_array.append(radiusArray)

    for i in range(data.shape[0]):
        cellArray.InsertCellPoint(i)
        linePoints.InsertNextPoint(data[i, :3])
        for j in range(3, data.shape[1]):
            info_array[j - 3].SetTuple1(i, data[i, j])

    if TNB is not None:
        for i in range(data.shape[0]):
            for j in range(data.shape[1] - 3, data.shape[1], 1):
                tnb_ = TNB[j - data.shape[1]][i, :]
                info_array[j].SetTuple3(i, tnb_[0], tnb_[1], tnb_[2])

    if PT is not None:
        start = data.shape[1] - 3 if TNB is None else data.shape[1]
        for i in range(PT[-1].shape[0]):
            for j in range(start, start + 2, 1):
                pt_ = PT[j - start][i, :]
                info_array[j].SetTuple3(i, pt_[0], pt_[1], pt_[2])

    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for i in range(len(header) - 3):
        line.GetPointData().AddArray(info_array[i])

    return line


def get_number_of_arrays(line):
    count = 0
    names = []
    name = 0
    while name is not None:
        name = line.GetPointData().GetArrayName(count)
        if name is not None:
            names.append(name)
            count += 1

    return count, names


def ExtractSingleLine(centerlines, id, startID=0, endID=None):
    cell = vtk.vtkGenericCell()
    centerlines.GetCell(id, cell)
    N = cell.GetNumberOfPoints() if endID is None else endID + 1

    line = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    cellArray.InsertNextCell(N - startID)
    linePoints = vtk.vtkPoints()

    arrays = []
    N_, names = get_number_of_arrays(centerlines)
    for i in range(N_):
        tmp = centerlines.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        radiusArray = get_vtk_array(names[i], tmp_comp, N - startID)
        arrays.append(radiusArray)

    getArray = []
    for i in range(N_):
        getArray.append(centerlines.GetPointData().GetArray(names[i]))

    count = 0
    for i in range(startID, N):
        cellArray.InsertCellPoint(count)
        linePoints.InsertNextPoint(cell.GetPoints().GetPoint(i))

        for j in range(N_):
            num = getArray[j].GetNumberOfComponents()
            if num == 1:
                tmp = getArray[j].GetTuple1(i)
                arrays[j].SetTuple1(count, tmp)
            elif num == 2:
                tmp = getArray[j].GetTuple2(i)
                arrays[j].SetTuple2(count, tmp[0], tmp[1])
            elif num == 3:
                tmp = getArray[j].GetTuple3(i)
                arrays[j].SetTuple3(count, tmp[0], tmp[1], tmp[2])
            elif num == 9:
                tmp = getArray[j].GetTuple9(i)
                arrays[j].SetTuple9(count, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                    tmp[5], tmp[6], tmp[7], tmp[8])
        count += 1

    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for j in range(N_):
        line.GetPointData().AddArray(arrays[j])
    return line


def move_past_sphere(centerline, center, r, start, step=-1, stop=0, X=0.8):
    """Moves a point along the centerline until it as outside MIS"""
    # Create the minimal inscribed sphere
    MISphere = vtk.vtkSphere()
    MISphere.SetCenter(center)
    MISphere.SetRadius(r * X)
    tempPoint = [0.0, 0.0, 0.0]

    # Go the length of one MISR backwards
    for i in range(start, stop, step):
        value = MISphere.EvaluateFunction(centerline.GetPoint(i))
        if (value >= 0.0):
            tempPoint = centerline.GetPoint(i)
            break

    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

    return tempPoint, r


def dist_sphere_curv(surface, centerlines, sac_center, misr_max, fileName, factor):
    # Get longest centerline
    length = []
    for i in range(centerlines.GetNumberOfLines()):
        line = ExtractSingleLine(centerlines, i)
    length.append(get_curvilinear_coordinate(line)[-1])
    ind_longest = length.index(max(length))

    # Get all bifurcations along the longest centerline
    bif = []
    bif_id = []
    longest_line = ExtractSingleLine(centerlines, ind_longest)
    tol = get_tolerance(centerlines)

    for i in range(centerlines.GetNumberOfLines()):
        if i == ind_longest: continue
        comp_line = ExtractSingleLine(centerlines, i)
        for j in range(comp_line.GetNumberOfPoints()):
            pnt1 = longest_line.GetPoints().GetPoint(j)
            pnt2 = comp_line.GetPoints().GetPoint(j)
            if distance(pnt1, pnt2) > tol:
                bif.append(pnt1)
                bif_id.append(j)
                break

    # Remove bifurcations detected twice
    pop = []
    for i in range(len(bif)):
        for j in range(i + 1, len(bif)):
            dist = distance(bif[i], bif[j])
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
    curvatureArray = get_array("Curvature", curvatureSurface)
    distance_to_sphere_array = get_array("DistanceToSpheres", distance_to_sphere)
    size_array = curvatureArray * distance_to_sphere_array * factor

    size_vtk_array = create_vtk_array(size_array, "Size")
    curvatureSurface.GetPointData().AddArray(size_vtk_array)

    WritePolyData(curvatureSurface, fileName)

    return distance_to_sphere


def dist_sphere_diam(surface, centerlines, sac_center, misr_max, fileName, factor):
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
    diameter_array = 2 * get_array("DistanceToCenterlines", distance_to_sphere)
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
        distance_to_spheres_array = get_array("DistanceToSpheres", distance_to_sphere)
        element_size = np.minimum(element_size, distance_to_spheres_array) * factor

    vtk_array = create_vtk_array(element_size, "Size")
    distance_to_sphere.GetPointData().AddArray(vtk_array)
    WritePolyData(distance_to_sphere, fileName)

    return distance_to_sphere


def mesh_alternative(surface):
    print("--- Meshing failed.")
    print("--- Proceeding with surface smooting and meshing.")
    surface = vmtkSmoother(surface, "laplace", iterations=500)

    subdiv = vmtkscripts.vmtkSurfaceSubdivision()
    subdiv.Surface = surface
    subdiv.Method = "butterfly"
    subdiv.Execute()
    surface = subdiv.Surface

    return vmtkSmoother(surface, "laplace", iterations=500)


def vmtkSmoother(surface, method, iterations=600):
    smoother = vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.NumberOfIterations = iterations
    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface
