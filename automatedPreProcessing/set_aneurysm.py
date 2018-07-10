from os import path
from common import *
from vmtk import vmtkscripts
from vmtkpointselector import *

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
    points = [get_point(aneurysmSeedIds.GetId(i)) for i in range(aneurysmSeedIds.GetNumberOfIds())]

    if dir_path is not None:
        info = {"number_of_aneurysms": len(points)}

    for i in range(len(points)):
        info["aneurysm_%d" % i] = points[i]

    write_parameters(info, dir_path)

    return points


if __name__ == "__main__":
    # TODO add argpars for folder where case folders are located
    # then do listfiles(), and get all folders or similar.

    for c in cases:
        provide_outlets(surface, path.join(dir_path, c))
