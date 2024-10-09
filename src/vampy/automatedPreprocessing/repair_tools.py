#!/usr/bin/env python

# Program:   AneuTool
# Module:    ToolRepairSTL.py
# Language:  Python
# Date:      $Date: 2016/17/04 00:00:00 $
# Version:   $Revision: 0.0.1 $
# Author:    Christophe Chnafa

#   Copyright (c) Christophe Chnafa. All rights reserved.

import vtk


def is_nan(num):
    return num != num


def print_surface_info(surface):
    nTriangles = surface.GetNumberOfCells()
    nPoints = surface.GetNumberOfPoints()
    print("> --- Surface overview:")
    print(("> Total number of triangles: %s." % nTriangles))
    print(("> Total number of points: %s." % nPoints))
    print(">")


def find_and_delete_nan_triangles(surface):
    ctrNaN = 0
    foundNaN = False
    nTriangles = surface.GetNumberOfCells()
    print("> --- Check the surface.")
    # The links from points to cells need to be build.
    surface.BuildLinks()
    for i in range(0, nTriangles):
        killThisTriangle = False
        nPointsForThisCell = surface.GetCell(i).GetPoints().GetNumberOfPoints()
        if nPointsForThisCell > 3:
            print(
                "> WARNING: found Cell with more than 3 points: there is more than triangles."
            )
        for j in range(0, nPointsForThisCell):
            x = [0.0, 0.0, 0.0]
            surface.GetCell(i).GetPoints().GetPoint(j, x)
            if is_nan(x[0]) | is_nan(x[1]) | is_nan(x[2]):
                ctrNaN += 1
                killThisTriangle = True
        if killThisTriangle:
            surface.DeleteCell(i)
    surface.RemoveDeletedCells()
    print(("> Found %s NaN cells." % ctrNaN))
    print(">")
    if ctrNaN > 0:
        foundNaN = True

    return foundNaN


def clean_surface(surface):
    print("> --- Cleaning the surface.")
    cleanPolyData = vtk.vtkCleanPolyData()
    if vtk.VTK_MAJOR_VERSION <= 5:
        cleanPolyData.SetInput(surface)
    else:
        cleanPolyData.SetInputData(surface)
    cleanPolyData.PointMergingOff()  # point locator will not be used, and points
    # that are not used by any cells will be eliminated,
    # but never merged.

    # In VTK the tolerance is defined as a fraction
    # of the bounding box length.
    tol = 0.0  # 0.0005
    cleanPolyData.SetTolerance(tol)
    cleanPolyData.Update()

    cleanPolyData.Update()
    outputPolyData = cleanPolyData.GetOutput()

    print("> Done.")
    print("> ")

    return outputPolyData
