#!/usr/bin/env python

# Program:   AneuTools
# Module:    DisplayData.py
# Language:  Python
# Date:      $Date: 2016/17/04 00:00:00 $
# Version:   $Revision: 0.0.1 $
# Author:    Christophe Chnafa

#   Copyright (c) Christophe Chnafa. All rights reserved.

import vtk


class VtkPointCloud:  # pragma: no cover
    def __init__(self, maxNumPoints=1e6):
        self.maxNumPoints = maxNumPoints
        self.vtkPolyData = vtk.vtkPolyData()
        self.clearPoints()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.vtkPolyData)
        self.vtkActor = vtk.vtkActor()
        self.vtkActor.SetMapper(mapper)
        self.vtkActor.GetProperty().SetPointSize(9.)
        self.vtkActor.GetProperty().SetColor(1., .0, .0)

    def addPoint(self, point):
        if self.vtkPoints.GetNumberOfPoints() < self.maxNumPoints:
            pointId = self.vtkPoints.InsertNextPoint(point[:])
            self.vtkCells.InsertNextCell(1)
            self.vtkCells.InsertCellPoint(pointId)
        else:
            print(">>> Warning: number of points > 1e6.")
        self.vtkCells.Modified()
        self.vtkPoints.Modified()

    def clearPoints(self):
        self.vtkPoints = vtk.vtkPoints()
        self.vtkCells = vtk.vtkCellArray()
        self.vtkPolyData.SetPoints(self.vtkPoints)
        self.vtkPolyData.SetVerts(self.vtkCells)


class DisplayModel(object):  # pragma: no cover
    def polyDataToActor(self, polyData, opacity=1.0):
        """Wrap the provided vtkPolyData object in a mapper and an actor,
        returning the actor. """
        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION > 5:
            mapper.SetInputData(polyData)
        else:
            mapper.SetInputConnection(polyData.GetProducerPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(opacity)

        return (actor)

    def setLight(self, renderer):
        lightKit = vtk.vtkLightKit()
        lightKit.MaintainLuminanceOn()
        lightKit.SetKeyLightIntensity(0.8)
        # 0 cold blue, 0.5 neutral white, 1 is reddish
        lightKit.SetKeyLightWarmth(0.5)
        lightKit.SetFillLightWarmth(0.5)

        # The function is called SetHeadLightWarmth starting from VTK 5.0
        try:
            lightKit.SetHeadLightWarmth(0.5)
        except Exception:
            lightKit.SetHeadlightWarmth(0.5)

        # intensity ratios
        lightKit.SetKeyToFillRatio(2.)
        lightKit.SetKeyToHeadRatio(7.)
        lightKit.SetKeyToBackRatio(1000.)
        lightKit.AddLightsToRenderer(renderer)

    def renderWindow(self, renderer, titleWindow):
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.SetSize(700, 700)
        renderWindow.AddRenderer(renderer)
        renderWindow.SetWindowName(titleWindow)

        # Assign a control style.
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderWindow)
        interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
        interactor.SetInteractorStyle(interactorStyle)

        # Initialize and start the interactor.
        interactor.Initialize()
        interactor.Start()

    def DisplayProbesAndModel(self, centerline, fileNameCenterline,
                              listProbePoints, model=None):
        """Displays a model and the corresponding probe points along
        the centerline. """

        if model is None:
            isDisplayingModel = False
        else:
            isDisplayingModel = True

        # Create a cloud of points from the list of probe points.
        pointCloud = VtkPointCloud()
        nPoints = len(listProbePoints)
        for i in range(0, nPoints - 1):
            pointCloud.addPoint(listProbePoints[i])

        # Create a rendering window and renderer.
        ren = vtk.vtkRenderer()
        renWindows = vtk.vtkRenderWindow()
        renWindows.AddRenderer(ren)

        # Create a renderwindowinteractor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWindows)

        # Assign a control style.
        style = vtk.vtkInteractorStyleTrackballCamera()
        iren.SetInteractorStyle(style)

        # Assign actor to the renderer.
        ren.AddActor(pointCloud.vtkActor)
        if isDisplayingModel:
            opacity = 0.25
            ren.AddActor(self.polyDataToActor(model, opacity))
        ren.AddActor(self.polyDataToActor(centerline))
        ren.SetBackground(.2, .3, .4)
        renWindows.SetSize(700, 700)

        # Create a text actor.
        txt = vtk.vtkTextActor()
        guiText = ("Centerline file name: "
                   + repr(fileNameCenterline.rsplit('/', 1)[-1]) + "\n"
                   + "Number of probes: " + repr(len(listProbePoints)) + "\n"
                   + "Q to exit.")
        txt.SetInputData(guiText)
        txtprop = txt.GetTextProperty()
        txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(15)
        txtprop.SetColor(1, 1, 1)
        txt.SetDisplayPosition(20, 30)

        # Assign actor to the renderer.
        ren.AddActor(txt)

        # Enable user interface interactor.
        iren.Initialize()
        renWindows.Render()
        renWindows.SetWindowName("Probe Points.")
        iren.Start()
