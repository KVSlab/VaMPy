# Copyright (c) 2025 Simula Research Laboratory
# SPDX-License-Identifier: GPL-3.0-or-later
from os import listdir, mkdir

from morphman.common import *
from scipy.interpolate import interp1d
from vtk.numpy_interface import dataset_adapter as dsa

from vampy.automatedPreprocessing.preprocessing_common import scale_surface


##############################################################
# A Collection of utility scripts for moving mesh generation #
##############################################################

def get_point_map(remeshed, remeshed_extended, profile="linear"):
    """
    Create a map between original surface model and model with
    flow extensions; remeshed_extended, for clamping the inlet(s)
    and outlet. A linear or sinusodial profile for movement between endpoints
    and moving domain is prescribed.
    Args:
        remeshed (vtkPolyData): Input surface
        remeshed_extended (vtkPolyData): Input surface with flow extensions
        profile (str): 'linear' or 'sine' for respective profiles
    Returns:
        distances (ndarray): Array of linear or sinusodial profile between boundary and moving mesh
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
        tmp_p_norm = 1e16
        # Some hacks to find the correct corresponding points
        for region_id in range(len(outer_features)):
            region_id_out = boundary_map[region_id]
            id_i = inner_locators[region_id].FindClosestPoint(point)
            id_o = outer_locators[region_id_out].FindClosestPoint(point)

            p_i = inner_features[region_id].GetPoint(id_i)
            p_o = outer_features[region_id_out].GetPoint(id_o)

            # Compute distances from current point to boundaries
            p_a = np.array(p_i) - np.array(point)
            p_b = np.array(p_o) - np.array(point)
            p_norm = np.linalg.norm(p_b) + np.linalg.norm(p_a)
            if p_norm < tmp_p_norm:
                tmp_id = region_id
                tmp_p_norm = p_norm

        regionId = tmp_id
        regionId_out = boundary_map[regionId]
        inner_id = inner_locators[regionId].FindClosestPoint(point)
        outer_id = outer_locators[regionId_out].FindClosestPoint(point)

        dist_to_boundary = get_distance(point, outer_regions[regionId_out].Points[outer_id])
        dist_between_boundaries = get_distance(inner_regions[regionId].Points[inner_id],
                                               outer_regions[regionId_out].Points[outer_id])
        if profile == "linear":
            distances[i] = (dist_to_boundary / dist_between_boundaries)
        elif profile == "sine":
            distances[i] = easing_cos(dist_to_boundary, dist_between_boundaries)
        point_map[i] = locator_remeshed.FindClosestPoint(inner_regions[regionId].Points[inner_id])

    # Let the points corresponding to the caps have distance 0
    point_map = point_map.astype(int)

    return distances, point_map


def easing_cos(dist_b, dist_bb):
    """
    A simple easing function between 0 and 1 for creating a sinusoidal profile
    Args:
        dist_b (ndarray): Distance to boundary for current point
        dist_bb (ndarray): Distance between boundaries

    Returns:
        distances (ndarray): Sinusoidal profile based on distance between point and boundary
    """
    return - 0.5 * (np.cos(np.pi * dist_b / dist_bb) - 1)


def move_surface_model(surface, original, remeshed, remeshed_extended, distance, point_map, file_path, i, points,
                       clamp_boundaries):
    """
    Computes and applies the displacement between the original and the given surface, then projects this
    displacement onto a remeshed surface. The displaced points of the remeshed surface are then saved
    to a specified file path.

    Args:
        surface (vtkPolyData): The source surface used to compute the displacement.
        original (vtkPolyData): The reference surface for computing the displacement.
        remeshed (vtkPolyData): A remeshed version of the original surface.
        remeshed_extended (vtkPolyData): An extended version of the remeshed surface.
        distance (float): Scalar factor for the displacement.
        point_map (ndarray): Mapping of the points from remeshed to remeshed_extended.
        file_path (str): Path to save the displaced points of the remeshed surface.
        i (int): Index to specify where in the 'points' array the data should be saved.
        points (ndarray): An array where the displaced points will be saved.
        clamp_boundaries (bool): If True, the displacement will be clamped by the distance in the extensions.

    Returns:
        None: The function writes results to a file and modifies the points array in place.
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

    # Get surface projection
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


def save_displacement(file_name_displacement_points, points, number_of_points=100):
    """
    Resample displacement points and write them
    to numpy data array
    Args:
        file_name_displacement_points (str): Path to numpy point path
        points (ndarray): Numpy array consisting of displacement surface points
        number_of_points (int): Number of points to interpolate displacement points over
    """
    points[:, :, -1] = points[:, :, 0]

    time = np.linspace(0, 1, points.shape[2])
    move = np.zeros((points.shape[0], points.shape[1], number_of_points + 1))
    time_r = np.linspace(0, 1, number_of_points + 1)
    # Use interp1d if smooth displacement
    x = interp1d(time, points[:, 0, :], axis=1)
    y = interp1d(time, points[:, 1, :], axis=1)
    z = interp1d(time, points[:, 2, :], axis=1)

    move[:, 0, :] = x(time_r)
    move[:, 1, :] = y(time_r)
    move[:, 2, :] = z(time_r)

    points = move
    points.dump(file_name_displacement_points)


def project_displacement(clamp_boundaries, distance, folder_extended_surfaces, folder_moved_surfaces, point_map,
                         surface, surface_extended, remeshed, scale_factor):
    """
    Projects the displacement of moved surfaces located in a specified folder onto an extended surface. The projected
    surfaces are then saved to a designated output directory.
    Only processes files in "folder_moved_surfaces" that end with ".vtp" or ".stl".

    Args:
        clamp_boundaries (bool): If True, the displacement will be clamped by the distance in the extensions.
        distance (float): Scalar factor for the displacement.
        folder_extended_surfaces (str): Path to the directory where extended surfaces will be saved.
        folder_moved_surfaces (str): Path to the directory containing the moved surfaces to be projected.
        point_map (ndarray): Mapping of the points from remeshed to surface_extended.
        surface (vtkPolyData): The reference surface for computing the displacement.
        surface_extended (vtkPolyData): An extended version of the reference surface.
        remeshed (vtkPolyData): A remeshed version of the original surface.
        scale_factor (float): Factor to scale the surfaces.

    Returns:
        array: A 3D numpy array where the displaced points will be saved. The shape is
        (num_points, 3, num_surfaces).
    """
    # Add extents to all surfaces
    extended_surfaces = sorted([f for f in listdir(folder_moved_surfaces) if f.endswith(".vtp") or f.endswith(".stl")])
    if not path.exists(folder_extended_surfaces):
        mkdir(folder_extended_surfaces)

    n_surfaces = len(extended_surfaces)

    print("--- Projecting surfaces\n")
    points = np.zeros((surface_extended.GetNumberOfPoints(), 3, n_surfaces))
    for i in range(n_surfaces):
        model_path = path.join(folder_moved_surfaces, extended_surfaces[i])
        if i == 0:
            points[:, :, i] = dsa.WrapDataObject(surface_extended).Points
            continue

        tmp_surface = read_polydata(model_path)
        if scale_factor is not None:
            tmp_surface = scale_surface(tmp_surface, scale_factor)

        new_path = path.join(folder_extended_surfaces, model_path.split("/")[-1])
        if not path.exists(new_path):
            move_surface_model(tmp_surface, surface, remeshed, surface_extended, distance, point_map, new_path, i,
                               points, clamp_boundaries)
    return points
