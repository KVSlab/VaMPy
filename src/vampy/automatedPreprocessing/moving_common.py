from os import listdir, mkdir

from morphman.common import *
from vtk.numpy_interface import dataset_adapter as dsa

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


##############################################################
# A Collection of utility scripts for moving mesh generation #
##############################################################

#######################
# MOVING MESH SCRIPTS #
#######################

def get_point_map(remeshed, remeshed_extended):
    """
    Create a map between original surface model and model with
    flow extensions; remeshed_extended, for clamping the inlet(s)
    and outlet. A linear profile for movement between endpoints
    and moving domain is prescribed.
    Args:
        remeshed (vtkPolyData): Input surface
        remeshed_extended (vtkPolyData): Input surface with flow extensions
    Returns:
        distances (ndarray): Array of linear profile between boundary and moving mesh
        point_map (ndarray): Mapping between surface IDs and points along extension
    """
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
        tmp_id = -1
        tmp_cross_norm = 1e16
        # Some hacks to find the correct corresponding points
        for region_id in range(len(outer_features)):
            region_id_out = boundary_map[region_id]
            id_i = inner_locators[region_id].FindClosestPoint(point)
            id_o = outer_locators[region_id_out].FindClosestPoint(point)

            p_i = inner_features[region_id].GetPoint(id_i)
            p_o = outer_features[region_id_out].GetPoint(id_o)

            # Compute norm of cross product between vectors pointing to surface point (on cylinder)
            p_a = np.array(p_i) - np.array(point)
            p_b = np.array(p_o) - np.array(point)
            p_cross = np.cross(p_a, p_b)
            p_cross_norm = np.linalg.norm(p_cross)
            if p_cross_norm < tmp_cross_norm:
                tmp_cross_norm = p_cross_norm
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


def move_surface_model(surface, original, remeshed, remeshed_extended, distance, point_map, file_path, i, points,
                       clamp_boundaries):
    """
    Args:
        surface:
        original:
        remeshed:
        remeshed_extended:
        distance:
        point_map:
        file_path:
        i:
        points:
        clamp_boundaries:
    Returns:
    """
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
    if clamp_boundaries:
        displacement[remeshed.Points.shape[0]:] = distance * displacement[point_map]
    else:
        displacement[remeshed.Points.shape[0]:] = displacement[point_map]

    # Move the mesh points
    new_surface.Points += displacement
    write_polydata(new_surface.VTKObject, file_path)
    points[:, :, i] = new_surface.Points.copy()
    new_surface.Points -= displacement


def save_displacement(file_name_displacement_points, points):
    """
    Resample displacement points and write them
    to numpy data array
    Args:
        file_name_displacement_points (str): Path to numpy point path
        points (ndarray): Numpy array consisting of displacement surface points
    """
    N = 200

    points[:, :, -1] = points[:, :, 0]
    time = np.linspace(0, 1, points.shape[2])
    N2 = N + N // (time.shape[0] - 1)

    move = np.zeros((points.shape[0], points.shape[1], N + 1))
    move[:, 0, :] = resample(points[:, 0, :], N2, time, axis=1)[0][:, :N - N2 + 1]
    move[:, 1, :] = resample(points[:, 1, :], N2, time, axis=1)[0][:, :N - N2 + 1]
    move[:, 2, :] = resample(points[:, 2, :], N2, time, axis=1)[0][:, :N - N2 + 1]

    points = move
    points.dump(file_name_displacement_points)


def project_displacement(clamp_boundaries, distance, folder_extended_surfaces, folder_moved_surfaces, point_map,
                         surface, surface_extended, remeshed):
    """
    Args:
        clamp_boundaries:
        distance:
        folder_extended_surfaces:
        folder_moved_surfaces:
        point_map:
        surface:
        surface_extended:
    Returns:
    """
    # Add extents to all surfaces
    extended_surfaces = sorted(
        [f for f in listdir(folder_moved_surfaces) if f[:2] in ["LA", "Co", "Fu", "Ge", "Sm", "AF", "an"]])
    if not path.exists(folder_extended_surfaces):
        mkdir(folder_extended_surfaces)

    n_surfaces = len(extended_surfaces)

    print("--- Projecting surfaces ---")
    points = np.zeros((surface_extended.GetNumberOfPoints(), 3, n_surfaces))
    for i in range(n_surfaces):
        model_path = path.join(folder_moved_surfaces, extended_surfaces[i])
        if i == 0:
            points[:, :, i] = dsa.WrapDataObject(surface_extended).Points
            continue

        tmp_surface = read_polydata(model_path)
        new_path = path.join(folder_extended_surfaces, model_path.split("/")[-1])
        if not path.exists(new_path):
            move_surface_model(tmp_surface, surface, remeshed, surface_extended, distance, point_map, new_path, i,
                               points, clamp_boundaries)
    return points


def create_funnel(surface, cl, cl_ext, radius):
    surface = dsa.WrapDataObject(surface)
    points = surface.Points
    p0 = np.array(cl.GetPoint(0))  # At MV
    p1 = np.array(cl_ext.GetPoint(0))  # At flow extension boundary
    n = n_z = (p1 - p0) / la.norm(p1 - p0)
    n_x = np.array([0, -n[2], n[1]])
    n_y = np.cross(n_z, n_x)
    tol = 1E-1
    max_dist = get_distance_between_points(p0, p1) + tol

    for j in range(len(points)):
        p = points[j]

        current_dist_p0 = get_distance_to_plane(n, p0, p)
        current_dist_p1 = get_distance_to_plane(n, p1, p)
        if current_dist_p0 <= max_dist and current_dist_p1 <= max_dist:
            x_mark = (p - p0).dot(n_x)
            y_mark = (p - p0).dot(n_y)
            z_mark = (p - p0).dot(n_z)

            def get_scale(d):
                # TODO: Find choice for "1.2" based on MISR
                return -d / max_dist * 0.95

            scale = get_scale(z_mark)
            x_new = x_mark * scale
            y_new = y_mark * scale
            z_new = 0

            X_mark = np.array([x_new, y_new, z_new])
            R_inv = np.array([n_x, n_y, n_z]).T

            p_new = R_inv.dot(X_mark)

        else:
            p_new = np.array([0, 0, 0])

        points[j] += p_new

    surface.SetPoints(points)

    return surface.VTKObject


def get_distance_between_points(p0, p1):
    p0 = np.array(p0)
    p1 = np.array(p1)

    return la.norm(p1 - p0)


def get_distance_to_plane(n, P, Q):
    return np.abs(n.dot(P - Q)) / la.norm(n)
