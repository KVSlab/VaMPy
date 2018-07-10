from common import *
import sys
import pandas as pn


def save_data(relative_path, case, data):
    dataFrame = pn.DataFrame(data, columns=data.keys())
    dataFrame.to_csv(path.join(relative_path, "morphology", "results", "Case%s.csv" % case),
                     index=False)


def get_data(relative_path, case):
    if not path.exists(path.join(relative_path, "morphology", "results", "Case%s.csv" % case)):
        data = create_empty_data(relative_path, case)
        teams = get_teams(relative_path)
        for i, team in enumerate(teams):
            data["teamId"][i] = team.split(path.sep)[-1]
        collect_precomputed(teams, case, data)
    else:
        dataFrame = pn.read_csv(path.join(relative_path, "morphology", "results", "Case%s.csv" % case))
        data = dataFrame.to_dict("list")

    indicies = get_indices(relative_path)
    keys = data.keys()
    for index in indicies:
        if index not in keys:
            data[index] = ["na"]*28

    return data


def collect_precomputed(teams, case, data):
    # Loop over all teams
    for i, team in enumerate(teams):
        files = listdir(path.join(team, case, "csvfile"))

        # Lopp over all files
        for f in files:
            # Read File
            text = open(path.join(team, case, "csvfile", f), "r").readlines()
            indicies = text[0].split()
            values = text[1].split()

            # Insert values into dict
            for index, value in zip(indicies, values):
                # Skipp indicies not in data
                if not index in data.keys(): continue
                try:
                    data[index][i] = float(value)
                except:
                    data[index][i] = value


def get_indices(relative_path):
    indicies = open(path.join(relative_path, "morphology", "indices.txt"), "r").readlines()
    indicies = [ind.strip() for ind in indicies if not ind.startswith("#")] # or ind != ""]
    indicies = [ind.strip() for ind in indicies if ind != ""]
    return indicies


def create_empty_data(relative_path, case):
    indicies = get_indices(relative_path)
    data = {}
    for ind in indicies:
        data[ind] = ["na"]*28

    return data


def rotate_surface(centeroid, normal, surface):
    translate = vtk.vtkTransform()
    translate.Translate(-centeroid[0], -centeroid[1], -centeroid[2])
    translateFilter = vtk.vtkTransformPolyDataFilter()
    translateFilter.SetTransform(translate)
    translateFilter.SetInputData(surface)
    translateFilter.Update()

    tmp_surf = translateFilter.GetOutput()
    WritePolyData(tmp_surf, "test_tmp.vtp")

    # Rotate around x
    rotate = vtk.vtkTransform()
    rotate_vector = np.cross([0, 0, 1], normal)
    rotate_vector = rotate_vector / np.linalg.norm(rotate_vector)
    rotate.RotateWXYZ(-180*math.acos(normal[2]) / math.pi, rotate_vector)
    rotateFilter = vtk.vtkTransformPolyDataFilter()
    rotateFilter.SetTransform(rotate)
    rotateFilter.SetInputData(translateFilter.GetOutput())
    rotateFilter.Update()

    return rotateFilter.GetOutput()


def compute_morphology(team, case, data, i):
    if not path.exists(path.join(team, case, "Case%s_sac.vtp" % case)):
        print("Skipp team", team.split(path.sep)[-1], "case", case)
        return
    centerline_name = "Case%s_centerline_relevant_outlets.vtp" % case
    surface = ReadPolyData(path.join(team, case, "Case%s_sac.vtp" % case))
    plane = ReadPolyData(path.join(team, case, "Case%s_firstclosedsection.vtp" % case))
    parameters = get_parameters(path.join(team, case))
    model = ReadPolyData(path.join(team, case, "Case%s.vtp" % case))
    closed_model = ReadPolyData(path.join(team, case, "Case%s_closed.vtp" % case))
    relevant_centerline = ReadPolyData(path.join(team, case, centerline_name))
    normal = [data["firstClosedN1"][i], data["firstClosedN2"][i], data["firstClosedN3"][i]]
    plane_points = np.array([plane.GetPoint(i) for i in range(plane.GetNumberOfPoints())])
    centeroid = compute_bary_center(plane_points)

    # Compute average neck diameter
    if data["averageNeckDiameter"][i] == "na":
        print("Compute neck diam")
        edgeFilter = vtk.vtkFeatureEdges()
        edgeFilter.SetInputData(plane)
        edgeFilter.BoundaryEdgesOn()
        edgeFilter.ManifoldEdgesOff()
        edgeFilter.FeatureEdgesOff()
        edgeFilter.NonManifoldEdgesOff()
        edgeFilter.Update()
        edge = edgeFilter.GetOutput()
        averageNeckDiameter = 2*np.sum(np.array([distance(edge.GetPoint(i), centeroid) for i in
                                        range(edge.GetNumberOfPoints())])) / edge.GetNumberOfPoints()
        data["averageNeckDiameter"][i] = averageNeckDiameter
    else:
        averageNeckDiameter = data["averageNeckDiameter"][i]

    # Compute prependicular hight, rotate the clipp plane to be in the x-y
    # axsis
    if data["prependicularHeight"][i] == "na":
        print("Compute prepen height")
        rotated_surface = rotate_surface(centeroid, normal, surface)
        WritePolyData(rotated_surface, path.join(team, case, "Case%s_rotated_surface.vtp" % case))
        rotated_plane = rotate_surface(centeroid, normal, plane)

        # Extrude plane
        extrude = vtk.vtkLinearExtrusionFilter()
        extrude.SetInputData(rotated_plane)
        extrude.SetExtrusionTypeToNormalExtrusion()
        extrude.SetVector(0, 0, 1) #Check if some are: [0, 0, -1] instead
        extrude.SetScaleFactor(50)
        extrude.Update()
        inclusion_box = extrude.GetOutput()

        # Check if point is inside plane
        inside = vtk.vtkSelectEnclosedPoints()
        inside.SetSurfaceData(inclusion_box)
        inside.SetInputData(rotated_surface)
        inside.Update()

        max_point = [0, 0, 0]
        for j in range(rotated_surface.GetNumberOfPoints()):
            point = rotated_surface.GetPoint(j)
            if point[2] > max_point[2] and inside.IsInside(j):
                max_point = point

        prependicularHight = max_point[2]
        data["prependicularHeight"][i] = prependicularHight
    else:
        prependicularHight = data["prependicularHeight"][i]

    # Compute maximum height
    if data["maximumHeight"][i] == "na":
        print("Maximum height")
        locator = get_locator(surface)
        ids = vtk.vtkIdList()
        locator.FindClosestNPoints(surface.GetNumberOfPoints(), centeroid, ids)
        point = surface.GetPoint(ids.GetId(surface.GetNumberOfPoints()-1))
        maximumHeight = distance(point, centeroid)
        data["maximumHeight"][i] = maximumHeight
    else:
        maximumHeight = data["maximumHeight"][i]

    # Compute Aspect ratio
    if data["aspectRatio"][i] == "na":
        print("aspect ratio")
        aspectRatio = prependicularHight / averageNeckDiameter
        data["aspectRatio"][i] = aspectRatio
    else:
        aspectRatio = data["aspectRatio"][i]

    # Size ratio
    if ("14" in team and "3" in case) or ("13" in team and "4" in case):
        pass
    elif data["sizeRatio"][i] == "na":
        print("size ratio")
        attri_path = path.join(team, case, "Case%s_centerline_area.vtp" % case)
        slices_path = path.join(team, case, "Case%s_centerline_slices.vtp" % case)

        # Get the location to meassure the area
        points = [parameters["par"]["div_point"], parameters["dau1"]["end_point"],
                  parameters["dau2"]["end_point"]]

        # Merge the two centerlines
        section1 = ExtractSingleLine(relevant_centerline, 0, startID=5)
        section2 = ExtractSingleLine(relevant_centerline, 1, startID=5)
        section1 = ExtractSingleLine(section1, 0, endID=section1.GetNumberOfPoints() - 10)
        section2 = ExtractSingleLine(section2, 0, endID=section2.GetNumberOfPoints() - 10)
        centerline_no_end = merge_data([section1, section2])

        # Compute area along the centerline
        print("Compute centerline sections")
        #if "19b" in team:# and "3" in case:
        #    embed()
        attributes, slices = compute_centerline_sections(model, centerline_no_end)
        inlet_area = get_array("CenterlineSectionArea", attributes)[0]*1.5
        WritePolyData(attributes, attri_path)
        WritePolyData(slices, slices_path)

        # Remove non-interesting slices
        slices = threshold(slices, "CenterlineSectionArea", upper=inlet_area, source=1)

        print("The rest")
        # Get all cells, NB! GetCell(int) returns a shallow copy
        sections = []
        for p in range(slices.GetNumberOfCells()):
            tmp_cell = vtk.vtkGenericCell()
            slices.GetCell(p, tmp_cell)
            sections.append(tmp_cell)

        # Compute all centers
        centers = []
        for section in sections:
            x_sum = 0; y_sum = 0; z_sum = 0
            for k in range(section.GetNumberOfPoints()):
                tmp_point = section.GetPoints().GetPoint(k)
                x_sum += tmp_point[0]
                y_sum += tmp_point[1]
                z_sum += tmp_point[2]
            centers.append([x_sum / section.GetNumberOfPoints(),
                            y_sum / section.GetNumberOfPoints(),
                            z_sum / section.GetNumberOfPoints()])

        # Get slices closest to end points
        dist1 = np.sqrt(np.sum((np.array(centers) - np.array(points[0]))**2, axis=1))
        dist2 = np.sqrt(np.sum((np.array(centers) - np.array(points[1]))**2, axis=1))
        dist3 = np.sqrt(np.sum((np.array(centers) - np.array(points[2]))**2, axis=1))

        ids = [list(dist1).index(np.min(dist1)),
               list(dist2).index(np.min(dist2)),
               list(dist3).index(np.min(dist3))]
        sections = [sections[ids[0]], sections[ids[1]], sections[ids[2]]]
        centers = [centers[ids[0]], centers[ids[1]], centers[ids[2]]]

        # FIXME: Store the selected slices
        #poly_sections = []
        text = ""
        for section in sections:
            for p in range(section.GetNumberOfPoints()):
                text += str(section.GetPoints().GetPoint(p))[1:-1] + "\n"
        f = open(path.join(team, case, "Case%s_sections.particles" % case), "w")
        f.write(text)
        f.close()
        #    tmp = vtk.vtkUnstructuredGrid()
        #    tmp.InsertNextCell(sections[k].GetCellType(), sections[k].GetPointIds())
        #    tmp.SetPoints(sections[k].GetPoints())
        #    geometryFilter = vtk.vtkGeometryFilter()
        #    geometryFilter.SetInputData(tmp)
        #    geometryFilter.Update()
        #    poly_sections.append(geometryFilter.GetOutput())
        #merge = merge_data(poly_sections)
        #WritePolyData(merge, path.join(team, case, "Case%s_sections.vtp" % case))

        # Compute mean parent artert diameter
        mean_diameter = []
        for j, section in enumerate(sections):
            sum_dist = 0
            for k in range(section.GetNumberOfPoints()):
                sum_dist += distance(centers[j], section.GetPoints().GetPoint(k))
            mean_diameter.append(2*sum_dist/section.GetNumberOfPoints())
        parentArteryDiameter = np.mean(mean_diameter)

        sizeRatio = maximumHeight / parentArteryDiameter
        data["parentArteryDiameter"][i] = parentArteryDiameter
        data["sizeRatio"][i] = sizeRatio
    else:
        sizeRatio = data["sizeRatio"][i]
        parentArteryArea = data["parentArteryDiameter"][i]

    # Compute the convex hull
    if data["convexHullVolume"][i] == "na":
        print("Convex hull")
        delaunay = vtk.vtkDelaunay3D()
        delaunay.SetInputData(surface)
        delaunay.Update()
        surfacefilter = vtk.vtkDataSetSurfaceFilter()
        surfacefilter.SetInputData(delaunay.GetOutput())
        surfacefilter.Update()
        convexHull = surfacefilter.GetOutput()
        WritePolyData(convexHull, path.join(team, case, "Case%s_convex_hull.vtp" % case))

        # Compute measures of the convex hull
        mass = vtk.vtkMassProperties()
        mass.SetInputData(convexHull)
        mass.Update()
        convexHullSurface = mass.GetSurfaceArea()
        convexHullVolume = mass.GetVolume()
        data["convexHullSurface"][i] = convexHullSurface
        data["convexHullVolume"][i] = convexHullVolume
    else:
        convexHullSurface = data["convexHullSurface"][i]
        convexHullVolume = data["convexHullVolume"][i]

    # Convex hull indicies
    if data["nonsphericityIndex"][i] == "na":
        print("convex indesies")
        undulationIndex = 1 - (data["sacVolume"][i] / convexHullVolume)
        ellipticityIndex = 1 - ((18*math.pi)**(1./3)*convexHullVolume**(2./3) / convexHullSurface)
        nonsphericityIndex = 1 - (18*math.pi)**(1./3)*data["sacVolume"][i]**(2./3) / data["sacSurface"][i]
        data["undulationIndex"][i] = undulationIndex
        data["ellipticityIndex"][i] = ellipticityIndex
        data["nonsphericityIndex"][i] = nonsphericityIndex
    else:
        undulationIndex = data["undulationIndex"][i]
        ellipticityIndex = data["ellipticityIndex"][i]
        nonsphericityIndex = data["nonsphericityIndex"][i]


if __name__ == "__main__":
    relative_path = path.join(path.dirname(path.abspath(__file__)), "..")
    teams = get_teams(relative_path)

    for case in ["1", "2", "3", "4", "5"]:
        data = get_data(relative_path, case)
        for i, team in enumerate(teams):
            print("Work on team", team.split(path.sep)[-1])
            compute_morphology(team, case, data, i)
        save_data(relative_path, case, data)
