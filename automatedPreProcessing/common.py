import vtk
from vmtk import vtkvmtk, vmtkscripts
try:
    from vmtkpointselector import *
except:
    pass
import numpy as np
import sys
import re
from os import path, listdir
from IPython import embed
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

# Options not availeble from commandline
divergingRatioToSpacingTolerance = 2.0
interpolationHalfSize = 3
voronoiCoreCutOffThreshold = 0.75
numberOfSplineAnalyzedPoints = 40
phiValues = [float(i) for i in range(2, 43, 2)]
thetaStep = 2.0

# Shortcuts
#distance = vtk.vtkMath.Distance2BetweenPoints
version = vtk.vtkVersion().GetVTKMajorVersion()

def ReadPolyData(filename):
    '''Load the given file, and return a vtkPolyData object for it. '''

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
    """Write the given input data based on the file name extension"""
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
        lineN = line[lineNid:lineNid + line[lineNid:].find(',') ].split('=')[1]
        numberOfNodes = int(lineN)
        lineEid = line.find('E=')
        lineE = line[lineEid:lineEid + line[lineEid:].find(',') ].split('=')[1]
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
        point = [float(data[dataCounter]), float(data[dataCounter+1]), float(data[dataCounter+2])]
        dataCounter += 3
        points.SetPoint(i,point)
        for j in range(len(arrayNames)):
            surface.GetPointData().GetArray(arrayNames[j]).SetComponent(i, 0, float(data[dataCounter]))
            dataCounter += 1

    # TODO: Is this necessary? It is not inserted into the surface
    # Insert ids
    cellIds = vtk.vtkIdList()
    for i in range(numberOfElements):
        cellIds.Initialize()
        cellIds.InsertNextId(int(data[dataCounter])-1)
        dataCounter += 1
        cellIds.InsertNextId(int(data[dataCounter])-1)
        dataCounter += 1
        cellIds.InsertNextId(int(data[dataCounter])-1)
        dataCounter += 1

        if elementType == "QUADRILATERAL":
            cellIds.InsertNextId(int(data[dataCounter])-1)
            dataCounter += 1

        cells.InsertNextCell(cellIds)

    return surface


def WriteTecplotSurfaceFile(surface, filename):
    # Triangulate surface
    surface = triangulate_surface(surface)

    # Open file
    f = open(fileName, 'w')
    line = "VARIABLES = X,Y,Z"

    # Get array
    arrayNames = []
    for i in range(self.Surface.GetPointData().GetNumberOfArrays()):
        array = self.Surface.GetPointData().GetArray(i)
        arrayName = array.GetName()

        if arrayName == None:
            continue
        if (arrayName[-1] == '_'):
            continue

        arrayNames.append(arrayName)

        if (array.GetNumberOfComponents() == 1):
            line = line + ',' + arrayName
        else:
            for j in range(array.GetNumberOfComponents()):
                line = line + ',' + arrayName + str(j)

    line = line + '\n'
    f.write(line)

    # Write header
    line = "ZONE N=%s,E=%s,F=FEPOINT,ET=TRIANGLE\n" % \
            (str(self.Surface.GetNumberOfPoints()),
             str(self.Surface.GetNumberOfCells()))
    f.write(line)

    # Write coordintes and array
    for i in range(self.Surface.GetNumberOfPoints()):
        point = self.Surface.GetPoint(i)
        line = " ".join(point)
        for arrayName in arrayNames:
            array = self.Surface.GetPointData().GetArray(arrayName)
            for j in range(array.GetNumberOfComponents()):
                line = line + ' ' + str(array.GetComponent(i,j))
        line = line + '\n'
        f.write(line)

    # Write conectivity(?) and ids
    for i in range(self.Surface.GetNumberOfCells()):
        cellPointIds = self.Surface.GetCell(i).GetPointIds()
        line = ''
        for j in range(cellPointIds.GetNumberOfIds()):
            if (j > 0):
                line = line + ' '
            line = line + str(cellPointIds.GetId(j) + 1)
        line = line + '\n'
        f.write(line)


def uncapp_surface(surface, centerlines, filename, clipspheres=0):
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
            a = int(folder)
            teams.append(path.join(dirpath, folder))
        except:
            pass

    teams.sort(key=lambda x: int(x.split(path.sep)[-1]))
    return teams


def makeVoronoiDiagram(surface, file_path):
    if path.isfile(file_path):
        return ReadPolyData(file_path)

    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = 0
    voronoi.Execute()

    WritePolyData(voronoi.VoronoiDiagram, file_path)

    return voronoi.VoronoiDiagram


def write_spheres(points, dirpath, radius=None, name="sphere%s.vtp", base=0.2):
    radius = [base]*len(points) if radius is None else radius
    for counter, point in enumerate(points):
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(point)
        sphere.SetPhiResolution(100)
        sphere.SetThetaResolution(100)
        sphere.SetRadius(radius[counter])
        sphere_ = sphere.GetOutput()

        WritePolyData(sphere_, path.join(dirpath, name % counter))


def get_tolerance(centerline):
    line = ExtractSingleLine(centerline, 0)
    length = get_curvilinear_coordinate(line)

    return np.mean(length[1:] - length[:-1]) / divergingRatioToSpacingTolerance


def getRelevantOutlets(aneurysm_type, centerline, centerline_aneurysm, dirpath):
    parameters = get_parameters(dirpath)
    if "relevant_outlet_0" in parameters.keys():
        return parameters["relevant_outlet_0"], parameters["relevant_outlet_1"]

    N_aneurysm = centerline_aneurysm.GetNumberOfPoints()
    tol = get_tolerance(centerline)
    div_anu = []
    aneurysm = centerline_aneurysm.GetCell(0)
    get_point = aneurysm.GetPoints().GetPoint

    for j in range(centerline.GetNumberOfCells()):
        line = centerline.GetCell(j)
        get_point_line = line.GetPoints().GetPoint
        for i in range(min(N_aneurysm, line.GetNumberOfPoints())):
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
            if distance(line2_div_point, line_div_point) > 3*tol:
                diff = False
            endpoint2 = list(line2.GetPoints().GetPoint(line2.GetNumberOfPoints() - 1))
            line_id -= 1
    else:
        endpoint2 = endpoint1
        endpoint1 = inlet

    data = {}
    for i, e in enumerate([endpoint1, endpoint2]):
        data["relevant_outlet_%d" % i] = e

    write_parameters(data, dirpath)

    return endpoint1, endpoint2


def SmoothVoronoiDiagram(voronoi, centerlines, smoothingFactor,
                         no_smooth_cl=None):
    numberOfPoints = voronoi.GetNumberOfPoints()

    threshold = get_array(radiusArrayName, centerlines) * (1 - smoothingFactor)
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

        if distance(point, cl_point) > 3*threshold[id_]:
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
    curv_coor = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints() - 1):
        pnt1 = np.asarray(line.GetPoints().GetPoint(i))
        pnt2 = np.asarray(line.GetPoints().GetPoint(i+1))
        curv_coor[i+1] = np.sqrt(np.sum((pnt1 - pnt2)**2)) + curv_coor[i]

    return curv_coor


def merge_data(inputs):
    appendFilter = vtk.vtkAppendPolyData()
    for input_ in inputs:
        appendFilter.AddInputData(input_)
    appendFilter.Update()

    return appendFilter.GetOutput()


def get_array(arrayName, line, k=1):
    array = np.zeros((line.GetNumberOfPoints(), k))
    if k == 1:
        getData = line.GetPointData().GetArray(arrayName).GetTuple1
    elif k == 2:
        getData = line.GetPointData().GetArray(arrayName).GetTuple2
    elif k ==3:
        getData = line.GetPointData().GetArray(arrayName).GetTuple3

    for i in range(line.GetNumberOfPoints()):
        array[i,:] = getData(i)

    return array


def get_array_cell(arrayName, line, k=1):
    array = np.zeros((line.GetNumberOfCells(), k))
    if k == 1:
        getData = line.GetCellData().GetArray(arrayName).GetTuple1
    elif k == 2:
        getData = line.GetCellData().GetArray(arrayName).GetTuple2
    elif k ==3:
        getData = line.GetCellData().GetArray(arrayName).GetTuple3
    elif k == 9:
        getData = line.GetCellData().GetArray(arrayName).GetTuple9

    for i in range(line.GetNumberOfCells()):
        array[i,:] = getData(i)

    return array


# TODO: Get bounds and compute polyballs based on that
#       bounds = surface.GetBounds()
def create_new_surface(completeVoronoiDiagram, polyBallImageSize=[280, 280, 280]):
    # Get the x,y, and z range of the completeVoronoiDiagram
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    if version < 6:
        modeller.SetInput(completeVoronoiDiagram)
    else:
        modeller.SetInputData(completeVoronoiDiagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(polyBallImageSize)
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


def get_relevant_outlets_tavel(surface, dir_path):
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
        provide_relevant_outlets_tawss(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    relevant_outlets = []
    inlet = []
    for key, value in parameters.items():
        if key.startswith("tavel_relevant_outlet_"):
            relevant_outlets.append(value)
        elif key.startswith("tavel_inlet"):
            inlet = value

    if relevant_outlets == []:
        inlet, relevant_outlets = provide_relevant_outlets_inlet_tavel(surface, dir_path)

    return inlet, relevant_outlets


def get_relevant_outlets_tawss(surface, dir_path):
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
        provide_relevant_outlets_tawss(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    relevant_outlets = []
    inlet = []
    for key, value in parameters.items():
        if key.startswith("tawss_relevant_outlet_"):
            relevant_outlets.append(value)
        elif key.startswith("tawss_inlet"):
            inlet = value

    if relevant_outlets == []:
        inlet, relevant_outlets = provide_relevant_outlets_inlet_tawss(surface, dir_path)

    return inlet, relevant_outlets


def get_relevant_outlets(surface, dir_path):
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
        provide_relevant_outlets(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    relevant_outlets = []
    for key, value in parameters.items():
        if key.startswith("relevant_outlet_"):
            relevant_outlets.append(value)

    if relevant_outlets == []:
        relevant_outlets = provide_relevant_outlets(surface, dir_path)

    return relevant_outlets


def get_aneurysm_dome(surface, dir_path):
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
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


def provide_relevant_outlets_tavel(surface, dir_path=None):
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
    info = {}

    if dir_path is not None:
        info["tavel_inlet"] = points[0]
        for i in range(1, len(points)):
            info["tavel_relevant_outlet_%d" % (i-1)] = points[i]
            write_parameters(info, dir_path)

    return points[0], points[1:]


def provide_relevant_outlets_tawss(surface, dir_path=None):
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
    info = {}

    if dir_path is not None:
        info["tawss_inlet"] = points[0]
        for i in range(1, len(points)):
            info["tawss_relevant_outlet_%d" % (i-1)] = points[i]
            write_parameters(info, dir_path)

    return points[0], points[1:]


def provide_relevant_outlets(surface, dir_path=None):
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
    info = {}

    if dir_path is not None:
        for i in range(len(points)):
            info["relevant_outlet_%d" % i] = points[i]
            write_parameters(info, dir_path)

    return points


def provide_aneurysm_points(surface, dir_path=None):
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
    # Declear variables before loop if values are not found
    data = {"dau1":{}, "dau2":{}}

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
                                        stop=point_ID_0*100, step=1, X=1)

        data[key]["end_point"] = end
        data[key]["r_end"] = r_end
        data[key]["div_point"] = center
        data[key]["ID_div"] = point_ID_0
        data[key]["i_div"] = tmpI
        data[key]["r_div"] = r

    return data


def write_points(points, filename):
    pointSet = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()

    for i in range(points.GetNumberOfPoints()):
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    pointSet.SetPoints(points)
    pointSet.SetVerts(cellArray)

    WritePolyData(pointSet, filename)


def surface_cleaner(surface):
    "Clean surface"
    surfaceCleaner = vtk.vtkCleanPolyData()
    if version < 6:
        surfaceCleaner.SetInput(surface)
    else:
        surfaceCleaner.SetInputData(surface)
    surfaceCleaner.Update()

    return surfaceCleaner.GetOutput()


def get_area(dir_path):
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
        compute_centers(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    outlets_area = []
    for key, value in parameters.items():
        if key == "inlet_area":
            inlet_area = value
        elif "outlet" in key and "area" in key and "relevant" not in key:
            outlets_area.append(value)

    if len(outlets_area) != 0:
        outlet_area = []
        for i in range(len(outlets_area)):
            outlet_area.append(parameters["outlet%d_area" % i])

    return inlet_area, outlet_area


def get_centers(surface, dir_path, flowext=False):
    # Check if info exists
    if flowext or not path.isfile(path.join(dir_path, "info.txt")):
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
    """Triangulate surface or polygon(?)"""
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
    # source = 1 uses cell data as input
    # source = 0 uses point data as input

    # Apply threshold
    threshold = vtk.vtkThreshold()
    if version < 6:
        threshold.SetInput(surface)
    else:
        threshold.SetInputData(surface)
    if type=="between":
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
    "Compute area of polydata"
    mass = vtk.vtkMassProperties()
    if version < 6:
        mass.SetInput(surface)
    else:
        mass.SetInputData(surface)

    return mass.GetSurfaceArea()


def uncapp_surface_old(surface):
    # Add-hoc method for removing capps on surfaces
    # This could proboly be highly improved, but is sufficient for now.

    # Get cell normals
    normal_generator = vtk.vtkPolyDataNormals()
    if version < 6:
        normal_generator.SetInput(surface)
    else:
        normal_generator.SetInputData(surface)
    normal_generator.ComputePointNormalsOff()
    normal_generator.ComputeCellNormalsOn()
    normal_generator.Update()
    cell_normals = normal_generator.GetOutput()

    # Compute gradients of the normals
    gradients_generator = vtk.vtkGradientFilter()
    if version < 6:
        gradients_generator.SetInput(cell_normals)
    else:
        gradients_generator.SetInputData(cell_normals)
    gradients_generator.SetInputArrayToProcess(0, 0, 0, 1, "Normals")
    gradients_generator.Update()
    gradients = gradients_generator.GetOutput()

    # Compute the magnitude of the gradient
    gradients_array = get_array_cell("Gradients", gradients, k=9)
    gradients_magnitude = np.sqrt(np.sum(gradients_array**2, axis=1))

    # Mark all cells with a gradient magnitude less then 0.1
    end_capp_array = gradients_magnitude < 0.08
    end_capp_vtk = get_vtk_array("Gradients_mag", 1, end_capp_array.shape[0])
    for i, p in enumerate(end_capp_array):
        end_capp_vtk.SetTuple(i, [p])
    gradients.GetCellData().AddArray(end_capp_vtk)

    # Extract capps
    end_capps = threshold(gradients, "Gradients_mag", lower=0.5, upper=1.5,
                            type="between", source=1)

    # Get connectivity
    end_capps_connectivity = getConnectivity(end_capps)
    region_array = get_array("RegionId", end_capps_connectivity)

    # Compute area for each region
    area = []
    circleness = []
    regions = []
    centers_edge = []
    limit = 0.1
    for i in range(int(region_array.max())+1):
        regions.append(threshold(end_capps_connectivity, "RegionId",  lower=(i-limit),
                            upper=(i+limit), type="between", source=0))
        circ, center = computeCircleness(regions[-1])
        circleness.append(circ)
        centers_edge.append(center)
        area.append(compute_area(regions[-1]))

    # Only keep outlets with circleness < 3 and area > 0.3 mm^2
    circleness_IDs = np.where(np.array(circleness) < 3)
    region_IDs = np.where(np.array(area) > 0.3)
    regions = [regions[i] for i in region_IDs[0] if i in circleness_IDs[0]]
    centers_edge = [centers_edge[i] for i in region_IDs[0] if i in circleness_IDs[0]]

    # Mark the outlets on the original surface
    mark_outlets = create_vtk_array(np.zeros(surface.GetNumberOfCells()), "outlets", k=1)
    locator = get_locator_cell(surface)
    tmp_center = [0, 0, 0]
    for region in regions:
        centers_filter = vtk.vtkCellCenters()
        if version < 6:
            centers_filter.SetInput(region)
        else:
            centers_filter.SetInputData(region)
        centers_filter.VertexCellsOn()
        centers_filter.Update()
        centers = centers_filter.GetOutput()

        for i in range(centers.GetNumberOfPoints()):
            centers.GetPoint(i, tmp_center)
            p = [0, 0, 0]
            cellId = vtk.mutable(0)
            subId = vtk.mutable(0)
            dist = vtk.mutable(0)
            locator.FindClosestPoint(tmp_center, p, cellId, subId, dist)
            mark_outlets.SetTuple(cellId, [1])

    surface.GetCellData().AddArray(mark_outlets)

    # Remove the outlets from the original surface
    uncapped_surface = threshold(surface, "outlets", lower=0, upper=0.5,
                                type="between", source=1)

    # Check if some cells where not marked
    remove = True
    while remove:
        locator = get_locator_cell(uncapped_surface)
        mark_outlets = create_vtk_array(np.zeros(uncapped_surface.GetNumberOfCells()), "outlets", k=1)
        remove = False
        for center in centers_edge:
            locator.FindClosestPoint(center, p, cellId, subId, dist)
            if dist < 0.01:
                remove = True
                mark_outlets.SetTuple(cellId, [1])

        uncapped_surface.GetCellData().AddArray(mark_outlets)

        if remove:
            uncapped_surface = threshold(uncapped_surface, "outlets", lower=0,
                                            upper=0.5, type="between", source=1)

    return uncapped_surface


def capp_surface(surface):
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
    return compute_centers(surface, test_capped=True)[0]


def getConnectivity(surface):
    """Compute connectivity of the cells"""
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
    edges = getFeatureEdges(surface)

    # Get points
    points = []
    for i in range(edges.GetNumberOfPoints()):
        points.append(edges.GetPoint(i))

    # Compute center
    points = np.array(points)
    center = np.mean(np.array(points), axis=0)

    # Compute ratio between max inscribed sphere, and min inscribed "area"
    point_radius = np.sqrt(np.sum((points-center)**2, axis=1))
    argsort = np.argsort(point_radius)
    if point_radius[argsort[1]] / point_radius[argsort[0]] > 15:
        radius_min = point_radius[argsort[1]]
    else:
        radius_min = point_radius.min()

    min_area = math.pi * radius_min**2
    max_area = math.pi * point_radius.max()**2

    return max_area / min_area, center


def getFeatureEdges(polyData):
    """Extracts the edges of the cells that are open"""
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
        boundary = threshold(outputs, "RegionId", lower=i-0.1, upper=i+0.1,
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
    shifted[1:,:] = points[:-1,:]
    shifted[0,:] = points[-1,:]

    # Compute weights
    weight = np.sqrt(np.sum((points - shifted)**2, axis=1))
    weight_sum = np.sum(weight)

    # Compute center
    center_x = np.sum((points[:,0] + shifted[:,0])/2 * weight) / weight_sum
    center_y = np.sum((points[:,1] + shifted[:,1])/2 * weight) / weight_sum
    center_z = np.sum((points[:,2] + shifted[:,2])/2 * weight) / weight_sum

    return [center_x, center_y, center_z]


def get_vtk_array(name, comp, num):
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(comp)
    array.SetNumberOfTuples(num)
    for i in range(comp):
        array.FillComponent(i, 0.0)
    array.SetName(name)

    return array


def get_locator_cell(surface):
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()

    return locator


def get_locator(centerline):
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(centerline)
    locator.BuildLocator()

    return locator


def distance(point1, point2):
    return np.sqrt(np.sum((np.asarray(point1) - np.asarray(point2))**2))


def remove_distant_points(voronoi, centerline):
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
        #comp = (47.424041748046875, 43.039527893066406, 41.241416931152344)
        if dist/3 > get_data(i) or get_data(i) > limit:
            #print(point)
            count += 1
            continue

        points.InsertNextPoint(point)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i-count)
        value = get_data(i)
        radius[i-count] = value

    print("Removed %s points from the voronoi diagram" % count)

    radiusArray = get_vtk_array(radiusArrayName, 1, N-count)
    for i in range(N-count):
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


def compute_centerlines(inlet, outlet, filepath, surface, resampling=1,
                        smooth=False, num_iter=100, smooth_factor=0.1,
                        endPoint=1, method="pointlist"):
    if path.isfile(filepath):
        return ReadPolyData(filepath)

    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = 'pointlist'
    centerlines.AppendEndPoints = endPoint
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet
    centerlines.Execute()
    centerlines = centerlines.Centerlines

    if smooth:
        centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
        centerlineSmoothing.SetInputData(self.Centerlines)
        centerlineSmoothing.SetNumberOfSmoothingIterations(num_iter)
        centerlineSmoothing.SetSmoothingFactor(smooth_factor)
        centerlineSmoothing.Update()

        centerlines = centerlinesSmooth.GetOutput()

    # Save the computed centerline.
    if filepath is not None:
        WritePolyData(centerlines, filepath)

    return centerlines


# TODO: Replace with vmtkscript
def CenterlineAttribiutes(line, remove=True, filename=None, smooth=False,
                         it=300, factor=0.1):
    if filename is None:
        filename = "tmp_cl.vtp"
        WritePolyData(line, filename)

    command = ('vmtkcenterlineattributes -ifile %s --pipe vmtkcenterlinegeometry ' + \
               '-ofile %s') % (filename, filename)
    if smooth:
        command += ' -smoothing 1 iterations %s -factor %s -outputsmoothd 1' % \
                   (it, factor)
    else:
        command += ' -smoothing 0'
    a = check_output(command, stderr=STDOUT, shell=True)
    status, text = success(a)

    if not status:
        print("smoething went wront when finding the attributes for the centerline")
        print(text)
        sys.exit(0)

    line = ReadPolyData(filename)
    if remove:
        check_output('rm ' + filename, stderr=STDOUT, shell=True)

    return line


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
            vtkArray.SetTuple2(i, values[i,0], values[i,1])
    elif k == 3:
        for i in range(values.shape[0]):
            vtkArray.SetTuple3(i, values[i,0], values[i,1], values[i,2])
    elif k == 9:
        for i in range(values.shape[0]):
            vtkArray.SetTuple9(i, values[i,0], values[i,1], values[i,2],
                                  values[i,3], values[i,4], values[i,5],
                                  values[i,6], values[i,7], values[i,8])

    return vtkArray


def GramSchmidt(V):
    V = 1.0 * V
    U = np.copy(V)

    def proj(u, v):
        return u * np.dot(v,u) / np.dot(u,u)

    for i in xrange(1, V.shape[1]):
        for j in xrange(i):
            U[:,i] -= proj(U[:,j], V[:,i])

    # normalize column
    den=(U**2).sum(axis=0)**0.5
    E = U/den
    return E


def get_parameters(folder):
    # If info.txt file, return an empty dict
    if not path.isfile(path.join(folder, "info.txt")): return {}

    # Get text
    f = open(path.join(folder, "info.txt"), "r")
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


def write_parameters(data, folder):
    # Get old parameters
    parameters = get_parameters(folder)

    # Add new parameters (can overwrite old as well)
    for key, value in data.items():
        parameters[key] = value

    # Get new text
    text = ["%s: %s" % (k, v) for k, v in parameters.items()]
    text = "\n".join(text)

    # Write text
    f = open(path.join(folder, "info.txt"), "w")
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
            radiusArray = get_vtk_array(header[i+data.shape[1]], 3, data.shape[0])
            info_array.append(radiusArray)

    if PT is not None:
        start = data.shape[1] if TNB is None else data.shape[1] + 3
        for i in range(2):
            radiusArray = get_vtk_array(header[i+start], 3, PT[0].shape[0])
            info_array.append(radiusArray)

    for i in range(data.shape[0]):
        cellArray.InsertCellPoint(i)
        linePoints.InsertNextPoint(data[i,:3])
        for j in range(3, data.shape[1]):
            info_array[j-3].SetTuple1(i, data[i, j])

    if TNB is not None:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]-3, data.shape[1], 1):
                tnb_ = TNB[j - data.shape[1]][i,:]
                info_array[j].SetTuple3(i, tnb_[0], tnb_[1], tnb_[2])

    if PT is not None:
        start = data.shape[1]-3 if TNB is None else data.shape[1]
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
        if (value>=0.0):
            tempPoint = centerline.GetPoint(i)
            break

    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

    return tempPoint, r


def dist_sphere_curv(surface, centerlines, sac_center, misr_max, fileName):
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
    size_array = curvatureArray * distance_to_sphere_array*factor

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
    element_size = 13. / 35 * diameter_array**2 + lower
    element_size[element_size > upper] = upper
    element_size[element_size < lower] = lower
    elements_vtk = create_vtk_array(element_size, "Num elements")
    distance_to_sphere.GetPointData().AddArray(elements_vtk)
    element_size = diameter_array / element_size
    #element_size[element_size < 0.12] = 0.12

    # Reduce element size in aneurysm
    for i in range(len(sac_center)):
        distance_to_sphere = compute_distance_to_sphere(distance_to_sphere,
                                                        sac_center[i],
                                                        maxDistance=100,
                                                        distanceScale=0.2 / (misr_max[i] * 2.))

    distance_to_spheres_array = get_array("DistanceToSpheres", distance_to_sphere)
    element_size = np.minimum(element_size, distance_to_spheres_array)*factor
    vtk_array = create_vtk_array(element_size, "Size")
    distance_to_sphere.GetPointData().AddArray(vtk_array)
    #distance_to_sphere.GetPointData().RemoveArray("DistanceToCenterlines")
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
    smoother= vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.NumberOfIterations = iterations
    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface



