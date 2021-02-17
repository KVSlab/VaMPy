#!/usr/bin/env python

import sys
import math
import os
import vtk
from common import *

import pythonprofileanalysis
from vmtk import vtkvmtk


def ExtractIds(points, centerlines, aneurysmType):
   ids = []
   numberOfCells = centerlines.GetNumberOfCells()

   if (aneurysmType == 'lateral'):
      if (points.GetNumberOfPoints()<=1):
         print('Error: number of clipping points not 2')
         sys.exit(0)

      point0 = points.GetPoint(0)
      point1 = points.GetPoint(1)

      line = ExtractCenterlineCell(0,centerlines)
      locator = vtk.vtkPointLocator()
      locator.SetDataSet(line)
      locator.BuildLocator()

      id0 = locator.FindClosestPoint(point0)
      id1 = locator.FindClosestPoint(point1)
      ids.append(id0)
      ids.append(id1)

   if (aneurysmType == 'terminal'):
       if (points.GetNumberOfPoints() != 3):
          print('Error: number of clipping points not 3')
          sys.exit(0)

       for i in range(1, numberOfCells+1):
          point0 = points.GetPoint(0)
          point1 = points.GetPoint(i)

          line = ExtractSinleLine(i-1, centerlines)

          locator = vtk.vtkPointLocator()
          locator.SetDataSet(line)
          locator.BuildLocator()

          tempid0 = locator.FindClosestPoint(point0)
          tempid1 = locator.FindClosestPoint(point1)
          ids.append(tempid0)
          ids.append(tempid1)

   return ids


def CreateCenterlines(points, centerlines, aneurysmType):
    # This function assumes that the centerlines are "sorted", and thus created
    # by remove_aneurysm.py (cell 0 and 1 are the centerlines corresponding to
    # the relevant outlets). The output is one/two centerlines between the clipping
    # points

    if points.GetNumberOfPoints() == 2:
        # Extract relevant centerline
        line = ExtractSingleLine(centerlines, 0)

        # Get locator
        locator = get_locator(line)

        # Get start and end ID
        id_start = locator.FindClosestPoint(points.GetPoint(0))
        id_end = locator.FindClosestPoint(points.GetPoint(1))

        # Extract relevant section
        section = ExtractSingleLine(line1, 0, start_id=id1_start, end_id=id1_end)

    # For Terminal Aneurysms
    else:
        # Extract relevant centerline
        line1 = ExtractSingleLine(centerlines, 0)
        line2 = ExtractSingleLine(centerlines, 1)

        # Get locator
        locator1 = get_locator(line1)
        locator2 = get_locator(line2)

        # Get start and end IDs
        id1_start = locator1.FindClosestPoint(points.GetPoint(0))
        id2_start = locator2.FindClosestPoint(points.GetPoint(0))
        id1_end1 = locator1.FindClosestPoint(points.GetPoint(1))
        id1_end2 = locator1.FindClosestPoint(points.GetPoint(2))
        id2_end1 = locator2.FindClosestPoint(points.GetPoint(1))
        id2_end2 = locator2.FindClosestPoint(points.GetPoint(2))

        # Chose end point for each centerline
        dist1 = distance(line1.GetPoint(id1_end1), points.GetPoint(1))
        dist2 = distance(line1.GetPoint(id1_end2), points.GetPoint(2))

        if dist1 < dist2:
            id1_end = id1_end1
            id2_end = id2_end2
        else:
            id1_end = id1_end2
            id2_end = id2_end1

        # Extract relevant section
        section1 = ExtractSingleLine(line1, 0, start_id=id1_start, end_id=id1_end)
        section2 = ExtractSingleLine(line2, 0, start_id=id2_start, end_id=id2_end)

        # Merge the two centerlines
        section = merge_data([section1, section2])

        assert section.GetNumberOfCells() == 2, "There are not 2 cells"
        assert section1.GetNumberOfPoints() + section2.GetNumberOfPoints() == section.GetNumberOfPoints()

    return section


def ExtractAneurysmSacVoronoiDiagram(voronoi, centerlines):
    numberOfPoints = voronoi.GetNumberOfPoints()

    clippingArray = vtk.vtkDoubleArray()
    clippingArray.SetNumberOfComponents(1)
    clippingArray.SetName(clippingArrayName)
    clippingArray.SetNumberOfTuples(numberOfPoints)

    voronoi.GetPointData().AddArray(clippingArray)

    tube = vtkvmtk.vtkvmtkPolyBallLine()
    tube.SetInput(centerlines)
    tube.SetPolyBallRadiusArrayName(radiusArrayName)
    tube.UseRadiusInformationOn()

    for i in range(numberOfPoints):
        point = voronoi.GetPoint(i)
        tubeValue = tube.EvaluateFunction(point)
        clippingArray.SetValue(i,tubeValue)

    voronoi.GetPointData().SetActiveScalars(clippingArrayName)

    extractor = vtk.vtkClipPolyData()
    extractor.SetInputData(voronoi)
    extractor.SetValue(0.0)
    extractor.InsideOutOff()
    extractor.Update()

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(extractor.GetOutput())
    connectivity.SetExtractionModeToLargestRegion()
    connectivity.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(connectivity.GetOutput())
    cleaner.Update()
    sacvoronoi = cleaner.GetOutput()

    return sacvoronoi


def ExtractAneurysmSacVoronoiCore(voronoi):
   voronoi.GetPointData().SetActiveScalars(radiusArrayName)
   radiusArray = voronoi.GetPointData().GetScalars()
   voronoiCoreCutOffValue = voronoiCoreCutOffThreshold*radiusArray.GetRange()[1]

   extractor = vtk.vtkClipPolyData()
   extractor.SetInputData(voronoi)
   extractor.SetValue(voronoiCoreCutOffValue)
   extractor.Update()

   connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
   connectivityFilter.SetInputData(extractor.GetOutput())
   connectivityFilter.ColorRegionsOff()
   connectivityFilter.SetExtractionModeToLargestRegion()
   connectivityFilter.Update()

   cleaner = vtk.vtkCleanPolyData()
   cleaner.SetInputData(connectivityFilter.GetOutput())
   cleaner.Update()

   voronoicore = cleaner.GetOutput()
   return voronoicore

def ClipVoronoiDiagramWithCore(voronoi,core,flag):
   polyBall = vtkvmtk.vtkvmtkPolyBall()
   polyBall.SetInputData(core)
   polyBall.SetPolyBallRadiusArrayName(radiusArrayName)

   clippingArray = vtk.vtkDoubleArray()
   clippingArray.SetName(clippingArrayName)
   clippingArray.SetNumberOfComponents(1)
   clippingArray.SetNumberOfTuples(voronoi.GetNumberOfPoints())

   voronoi.GetPointData().AddArray(clippingArray)

   for i in range(voronoi.GetNumberOfPoints()):
      point = voronoi.GetPoint(i)
      polyBallValue = polyBall.EvaluateFunction(point)
      clippingArray.SetValue(i,polyBallValue)

   voronoi.GetPointData().SetActiveScalars(clippingArrayName)

   voronoiExtractor = vtk.vtkClipPolyData()
   voronoiExtractor.SetInputData(voronoi)
   voronoiExtractor.SetValue(0.0)
   voronoiExtractor.SetInsideOut(flag)
   voronoiExtractor.Update()

   connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
   connectivityFilter.SetInputData(voronoiExtractor.GetOutput())
   connectivityFilter.ColorRegionsOff()
   connectivityFilter.SetExtractionModeToLargestRegion()
   connectivityFilter.Update()

   voronoiCleaner = vtk.vtkCleanPolyData()
   voronoiCleaner.SetInputData(connectivityFilter.GetOutput())
   voronoiCleaner.Update()

   clippedvoronoi = voronoiCleaner.GetOutput()
   return clippedvoronoi

def ExtractAneurysmSacSurface(surface, voronoisac, bjrvoronoi):
   numberOfPoints = surface.GetNumberOfPoints()

   clippingArray = vtk.vtkDoubleArray()
   clippingArray.SetNumberOfComponents(1)
   clippingArray.SetName(clippingArrayName)
   clippingArray.SetNumberOfTuples(numberOfPoints)

   surface.GetPointData().AddArray(clippingArray)

   polyBall = vtkvmtk.vtkvmtkPolyBall()
   polyBall.SetInputData(voronoisac)
   polyBall.SetPolyBallRadiusArrayName(radiusArrayName)

   tube = vtkvmtk.vtkvmtkPolyBall()
   tube.SetInputData(bjrvoronoi)
   tube.SetPolyBallRadiusArrayName(radiusArrayName)

   for i in range(numberOfPoints):
      point = surface.GetPoint(i)
      polyBallValue = polyBall.EvaluateFunction(point)
      tubeValue = tube.EvaluateFunction(point)
      clippingValue = polyBallValue-tubeValue
      clippingArray.SetValue(i,clippingValue)

   surface.GetPointData().SetActiveScalars(clippingArrayName)

   surfaceClipper = vtk.vtkClipPolyData()
   surfaceClipper.SetInputData(surface)
   surfaceClipper.SetValue(1E-4)
   surfaceClipper.InsideOutOn()
   surfaceClipper.Update()

   connectivity = vtk.vtkPolyDataConnectivityFilter()
   connectivity.SetInputData(surfaceClipper.GetOutput())
   connectivity.SetExtractionModeToLargestRegion()
   connectivity.Update()

   cleaner = vtk.vtkCleanPolyData()
   cleaner.SetInputData(connectivity.GetOutput())
   cleaner.Update()

   sacsurface = cleaner.GetOutput()
   return sacsurface

def RemoveBranchVesselFromSac(surface,notvoronoi,voronoi):
   clippingArray = vtk.vtkDoubleArray()
   clippingArray.SetNumberOfComponents(1)
   clippingArray.SetName(branchClippingArrayName)
   clippingArray.SetNumberOfTuples(surface.GetNumberOfPoints()) 
   surface.GetPointData().AddArray(clippingArray)

   polyBall1 = vtkvmtk.vtkvmtkPolyBall()
   polyBall1.SetInputData(notvoronoi)
   polyBall1.SetPolyBallRadiusArrayName(radiusArrayName)

   polyBall2 = vtkvmtk.vtkvmtkPolyBall()
   polyBall2.SetInputData(voronoi)
   polyBall2.SetPolyBallRadiusArrayName(radiusArrayName)

   for i in range(surface.GetNumberOfPoints()):
     point = surface.GetPoint(i)
     value1 = polyBall1.EvaluateFunction(point)
     value2 = polyBall2.EvaluateFunction(point)
     value =  value1 - value2
     clippingArray.SetValue(i,value)

   surface.GetPointData().SetActiveScalars(branchClippingArrayName)
   
   surfaceClipper = vtk.vtkClipPolyData()
   surfaceClipper.SetInputData(surface)
   surfaceClipper.SetValue(1E-4)
   surfaceClipper.InsideOutOff()
   surfaceClipper.Update()
   
   connectivity = vtk.vtkPolyDataConnectivityFilter()
   connectivity.SetInputData(surfaceClipper.GetOutput())
   connectivity.SetExtractionModeToLargestRegion()
   connectivity.Update()
   
   cleaner = vtk.vtkCleanPolyData()
   cleaner.SetInputData(connectivity.GetOutput())
   cleaner.Update()
   
   sacsurface = cleaner.GetOutput()
   return sacsurface

def ComputeDistanceToParentVesselCenterlines(surface,voronoi,centerline):
   tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
   tubeFunction.SetInput(centerline)
   tubeFunction.SetPolyBallRadiusArrayName(radiusArrayName) 
   tubeFunction.UseRadiusInformationOn()
   # TODO: Check if there is another way to avoid this
   tubeFunction.SetGlobalWarningDisplay(0)

   polyBallFunction = vtkvmtk.vtkvmtkPolyBall()
   polyBallFunction.SetInput(voronoi)
   polyBallFunction.SetPolyBallRadiusArrayName(radiusArrayName)

   distanceArray = vtk.vtkDoubleArray()
   distanceArray.SetNumberOfComponents(1)
   distanceArray.SetNumberOfTuples(surface.GetNumberOfPoints())
   distanceArray.SetName(distanceToTubeArrayName)
   distanceArray.FillComponent(0,0.0) 

   surface.GetPointData().AddArray(distanceArray)

   for i in range(surface.GetNumberOfPoints()):
      point = surface.GetPoint(i)
      tubeValue = tubeFunction.EvaluateFunction(point)
      polyBallValue = polyBallFunction.EvaluateFunction(point)
      distance = polyBallValue - tubeValue
      distanceArray.SetValue(i,distance)

   return surface


def ExtractSplineFromContours(sac, numberofcontours):
   #embed()
   sac.GetPointData().RemoveArray(branchClippingArrayName)
   sac.GetPointData().SetActiveScalars(distanceToTubeArrayName)

   distanceArray = sac.GetPointData().GetArray(distanceToTubeArrayName)
   distanceArrayRange = distanceArray.GetRange()

   contourStep = (distanceArrayRange[1]-distanceArrayRange[0])/float(numberofcontours-1)

   contourBarycenterPoints = vtk.vtkPoints()

   for i in range(numberofcontours):
      contourValue = distanceArrayRange[1]-i*contourStep

      contourFilter = vtk.vtkContourFilter()
      contourFilter.SetInputData(sac)
      contourFilter.SetValue(0,contourValue)
      contourFilter.Update()

      contour = contourFilter.GetOutput()

      close = IsSectionClosed(contour)
      if (close==0):
          continue
      else:
         stripper = vtk.vtkStripper()
         stripper.SetInputData(contour)
         stripper.Update()

         connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
         connectivityFilter.SetInputData(stripper.GetOutput())
         connectivityFilter.ColorRegionsOff()
         connectivityFilter.SetExtractionModeToLargestRegion()
         connectivityFilter.Update()

         numberOfExtractedRegions = connectivityFilter.GetNumberOfExtractedRegions()

         if (connectivityFilter.GetOutput().GetNumberOfPoints()==0): continue
         contourBarycenter = ComputeBarycenter(contour)
         contourBarycenterPoints.InsertNextPoint(contourBarycenter)

   zeroProfile = ExtractCutSacProfile(sac,contourBarycenterPoints)
   zeroProfileBarycenter = ComputeBarycenter(zeroProfile)

   contourBarycenterPoints.InsertNextPoint(zeroProfileBarycenter)

   numberOfSplinePoints = 400
   spline = ComputeBarycenterSpline(contourBarycenterPoints, numberOfSplinePoints)
   return spline


def IsSectionClosed(singleProfile):
   closed = -1
   numberOfProfilePoints = singleProfile.GetNumberOfPoints()
   for i in range(numberOfProfilePoints):
      cellIds = vtk.vtkIdList()
      singleProfile.GetPointCells(i,cellIds)
      if (cellIds.GetNumberOfIds()==1):
         closed = 0
         break
      closed = 1
   return closed

def ComputeBarycenter(dataset):
   bc = [0.0,0.0,0.0]
   for i in range(dataset.GetNumberOfPoints()):
      point = [0.0,0.0,0.0]
      dataset.GetPoint(i,point)
      bc[0] = bc[0] + point[0]
      bc[1] = bc[1] + point[1]
      bc[2] = bc[2] + point[2]

   bc[0] = bc[0] / dataset.GetNumberOfPoints()
   bc[1] = bc[1] / dataset.GetNumberOfPoints()
   bc[2] = bc[2] / dataset.GetNumberOfPoints()

   return bc

def ExtractCutSacProfile(surface, points):
   featureEdges = vtk.vtkFeatureEdges()
   featureEdges.SetInputData(surface)
   featureEdges.BoundaryEdgesOn()
   featureEdges.FeatureEdgesOff()
   featureEdges.NonManifoldEdgesOff()
   featureEdges.ColoringOff()
   featureEdges.Update()

   stripper = vtk.vtkStripper()
   stripper.SetInputData(featureEdges.GetOutput())
   stripper.Update()

   firstPoint = points.GetPoint(0)
   cutConnectivityFilter = vtk.vtkPolyDataConnectivityFilter()
   cutConnectivityFilter.SetInputData(stripper.GetOutput())
   cutConnectivityFilter.ColorRegionsOff()
   cutConnectivityFilter.SetExtractionModeToClosestPointRegion()
   cutConnectivityFilter.SetClosestPoint(firstPoint)
   cutConnectivityFilter.Update()

   zeroProfile = cutConnectivityFilter.GetOutput()
   return zeroProfile

def ComputeBarycenterSpline(points,numberOfSplinePts):
   numberofpts = points.GetNumberOfPoints()

   xspline = vtk.vtkCardinalSpline()
   yspline = vtk.vtkCardinalSpline()
   zspline = vtk.vtkCardinalSpline()
 
   firstPoint = points.GetPoint(numberofpts-1)
   xspline.AddPoint(float(0),firstPoint[0])
   yspline.AddPoint(float(0),firstPoint[1])
   zspline.AddPoint(float(0),firstPoint[2])
  
   for i in range(1,points.GetNumberOfPoints()-1):
      point = [0.0,0.0,0.0]
      points.GetPoint(i,point)

      xspline.AddPoint(float(i),point[0])
      yspline.AddPoint(float(i),point[1])
      zspline.AddPoint(float(i),point[2])

   xspline.Compute()
   yspline.Compute()
   zspline.Compute()

   spline = vtk.vtkPolyData()
   splinePoints = vtk.vtkPoints()
   splineCellArray = vtk.vtkCellArray()
   splineCellArray.InsertNextCell(numberOfSplinePts)
 
   for i in range(numberOfSplinePts):
      xvalue = xspline.Evaluate(float(i)*(float(numberofpts)/float(numberOfSplinePts)))
      yvalue = yspline.Evaluate(float(i)*(float(numberofpts)/float(numberOfSplinePts)))
      zvalue = zspline.Evaluate(float(i)*(float(numberofpts)/float(numberOfSplinePts)))
      point = [xvalue,yvalue,zvalue]    
     
      splinePoints.InsertNextPoint(point)
      splineCellArray.InsertCellPoint(i)

   spline.SetPoints(splinePoints)
   spline.SetLines(splineCellArray)
   return spline

def ComputeAneurysmSacSectionsAlongBarycenterSpline(spline, sac,
                                                    csvdatafilename,
                                                    numberOfPoints, ID):
   datafile = open(csvdatafilename,'w')
   dataline = 'ID profile phi theta Closed Area Length MinSize MaxSize ShapeFactor'+'\n'
   datafile.write(dataline)

   for i in range(0,numberOfPoints):
      print('profile ',i)

      tangent = spline.GetPointData().GetArray('FrenetTangent').GetTuple3(i)
      normal = spline.GetPointData().GetArray('FrenetNormal').GetTuple3(i)
      binormal = spline.GetPointData().GetArray('FrenetBinormal').GetTuple3(i)
      point = spline.GetPoint(i)

      inverseRotationMatrix = ComputeRotationMatrix(normal,binormal,tangent)

      section, closed = CutSacWithPlane(sac,tangent,point)
      barycenter = ComputeBarycenter(section)

      if section.GetNumberOfPoints() < 4:
               sectionArea = 'NA'
               sectionLength = 'NA'
               sectionSizeRange = ['NA','NA']
               sectionShape = 'NA'
      else:
         sectionArea = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionArea(section)
         sectionLength = ComputeProfileLength(section)
         sectionSizeRange = [0.0, 0.0]
         sectionShape = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionShape\
                             (section,barycenter,sectionSizeRange)

      dataline = str(ID)+' '+str(i)+' '+str(0.0)+' '+str(0.0)+' '+str(closed)+' '+str(sectionArea)+' '+str(sectionLength)+' '+str(sectionSizeRange[0])+' '+str(sectionSizeRange[1])+' '+str(sectionShape)+'\n'
      datafile.write(dataline)

      r = vtk.vtkMath.Norm(tangent)

      thetaIncrement = math.radians(thetaStep)

      for k in range(len(phiValues)):
         phi = math.radians(phiValues[k])
         numberOfThetaIncrements = int((2.0*math.pi)/thetaIncrement)

         for j in range(numberOfThetaIncrements):
            theta = 0.0 +j*thetaIncrement

            sphericalx = r*math.cos(theta)*math.sin(phi)
            sphericaly = r*math.sin(theta)*math.sin(phi)
            sphericalz = r*math.cos(phi)
            sphericalPoint = [sphericalx,sphericaly,sphericalz]
            sphericalPoint.append(1.0)

            rotatedPoint = [0.0,0.0,0.0]

            sphericalBackPoint = inverseRotationMatrix.MultiplyDoublePoint(sphericalPoint)
            rotatedPoint[0] = point[0]+sphericalBackPoint[0]
            rotatedPoint[1] = point[1]+sphericalBackPoint[1]
            rotatedPoint[2] = point[2]+sphericalBackPoint[2]

            newPlaneNormal = [0.0,0.0,0.0]
            newPlaneNormal[0] = rotatedPoint[0]-point[0]
            newPlaneNormal[1] = rotatedPoint[1]-point[1]
            newPlaneNormal[2] = rotatedPoint[2]-point[2]
            vtk.vtkMath.Normalize(newPlaneNormal)

            rotatedSection,rotatedClosed = CutSacWithPlane(sac,newPlaneNormal,point)
            rotatedBarycenter = ComputeBarycenter(rotatedSection)

            numberOfRotatedSectionPoints = rotatedSection.GetNumberOfPoints()
            if (numberOfRotatedSectionPoints<10):
               rotatedSectionArea = 'NA'
               rotatedSectionLength = 'NA'
               rotatedSectionSizeRange = ['NA','NA']
               rotatedSectionShape = 'NA'
            else:
               rotatedSectionArea = vtkvmtk.vtkvmtkPolyDataBranchSections.\
                                        ComputeBranchSectionArea(rotatedSection)
               rotatedSectionLength = ComputeProfileLength(rotatedSection)
               rotatedSectionSizeRange = [0.0,0.0]
               rotatedSectionShape = vtkvmtk.vtkvmtkPolyDataBranchSections.\
                                        ComputeBranchSectionShape(rotatedSection,
                                                                  rotatedBarycenter,
                                                                  rotatedSectionSizeRange)

            dataline = str(ID)+' '+str(i)+' '+str(math.degrees(phi))+' '+str(math.degrees(theta))+' '+str(rotatedClosed)+' '+str(rotatedSectionArea)+' '+str(rotatedSectionLength)+' '+str(rotatedSectionSizeRange[0])+' '+str(rotatedSectionSizeRange[1])+' '+str(rotatedSectionShape)+'\n'
            datafile.write(dataline)

   datafile.close()


def CutSacWithPlane(surface, normal, origin):
   plane = vtk.vtkPlane()
   plane.SetOrigin(origin)
   plane.SetNormal(normal)

   cutter = vtk.vtkCutter()
   cutter.SetInputData(surface)
   cutter.SetCutFunction(plane)
   cutter.Update()

   close = IsSectionClosed(cutter.GetOutput())

   stripper = vtk.vtkStripper()
   stripper.SetInputData(cutter.GetOutput())
   stripper.Update()

   connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
   connectivityFilter.SetInputData(stripper.GetOutput())
   connectivityFilter.SetClosestPoint(origin)
   connectivityFilter.SetExtractionModeToClosestPointRegion()
   connectivityFilter.Update()

   contour = connectivityFilter.GetOutput()
   numberOfPoints = contour.GetNumberOfPoints()

   section = vtk.vtkPolyData()
   sectionPoints = vtk.vtkPoints()
   sectionCellArray = vtk.vtkCellArray()
   sectionCellArray.InsertNextCell(numberOfPoints)

   closedArray = vtk.vtkIntArray()
   closedArray.SetNumberOfComponents(1)
   closedArray.SetNumberOfTuples(1)
   closedArray.SetName(closedArrayName)
   closedArray.SetValue(0,close)

   for i in range(numberOfPoints):
      point = contour.GetPoint(i)
      sectionPoints.InsertNextPoint(point)
      sectionCellArray.InsertCellPoint(i)

   sectionCellArray.InsertCellPoint(0)

   section.SetPoints(sectionPoints)
   section.SetPolys(sectionCellArray)
   section.GetCellData().AddArray(closedArray)

   return section, close

def ComputeProfileLength(polygon):
   length = 0.0
   numberOfPoints = polygon.GetNumberOfPoints()
   for i in range(numberOfPoints-1):
      point0 = polygon.GetPoint(i)
      point1 = polygon.GetPoint(i+1)
      distance = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
      length += distance 

   lastPoint = polygon.GetPoint(numberOfPoints-1)
   firstPoint = polygon.GetPoint(0)
   distance = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(lastPoint,firstPoint))
   length += distance
   return length

def ComputeRotationMatrix(n,b,t):
   matrix = vtk.vtkMatrix4x4()
   matrix.SetElement(0,0,n[0])
   matrix.SetElement(0,1,n[1])
   matrix.SetElement(0,2,n[2])
   matrix.SetElement(0,3,0.0)

   matrix.SetElement(1,0,b[0])
   matrix.SetElement(1,1,b[1])
   matrix.SetElement(1,2,b[2])
   matrix.SetElement(1,3,0.0)

   matrix.SetElement(2,0,t[0])
   matrix.SetElement(2,1,t[1])
   matrix.SetElement(2,2,t[2])
   matrix.SetElement(2,3,0.0)

   matrix.SetElement(3,0,0.0)
   matrix.SetElement(3,1,0.0)
   matrix.SetElement(3,2,0.0)
   matrix.SetElement(3,3,1.0)

   matrix.Invert()
   return matrix


def ComputeFirstClosedSectionParameters(pointId, phi, theta, sac, spline, filename):
   point = spline.GetPoint(pointId)

   FFtangent = spline.GetPointData().GetArray('FrenetTangent').GetTuple3(pointId)
   FFnormal = spline.GetPointData().GetArray('FrenetNormal').GetTuple3(pointId)
   FFbinormal = spline.GetPointData().GetArray('FrenetBinormal').GetTuple3(pointId)

   inverseRotationMatrix = ComputeRotationMatrix(FFnormal,FFbinormal,FFtangent)

   r = vtk.vtkMath.Norm(FFtangent)
   phi = math.radians(phi)
   theta = math.radians(theta)

   x = r*math.cos(theta)*math.sin(phi)
   y = r*math.sin(theta)*math.sin(phi)
   z = r*math.cos(phi)
   newPoint = [x,y,z]
   newPoint.append(1.0)

   rotatedPoint = [0.0,0.0,0.0]
   backPoint = inverseRotationMatrix.MultiplyDoublePoint(newPoint)
   rotatedPoint[0] = point[0]+backPoint[0]
   rotatedPoint[1] = point[1]+backPoint[1]
   rotatedPoint[2] = point[2]+backPoint[2]

   newPlaneNormal = [0.0,0.0,0.0]
   newPlaneNormal[0] = rotatedPoint[0]-point[0]
   newPlaneNormal[1] = rotatedPoint[1]-point[1]
   newPlaneNormal[2] = rotatedPoint[2]-point[2]
   vtk.vtkMath.Normalize(newPlaneNormal)

   section, closed = CutSacWithPlane(sac,newPlaneNormal,point)
   barycenter = ComputeBarycenter(section)

   sectionArea = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionArea(section)
   sectionLength = ComputeProfileLength(section)
   sectionSizeRange = [0.0,0.0]
   sectionShape = vtkvmtk.vtkvmtkPolyDataBranchSections.ComputeBranchSectionShape(section,
                                                                                  barycenter,
                                                                                  sectionSizeRange)

   areaArray = vtk.vtkDoubleArray()
   areaArray.SetNumberOfComponents(1)
   areaArray.SetNumberOfTuples(1)
   areaArray.SetName('SectionArea')
   areaArray.SetTuple1(0,sectionArea)

   lengthArray = vtk.vtkDoubleArray()
   lengthArray.SetNumberOfComponents(1)
   lengthArray.SetNumberOfTuples(1)
   lengthArray.SetName('SectionLength')
   lengthArray.SetTuple1(0,sectionLength)

   minArray = vtk.vtkDoubleArray()
   minArray.SetNumberOfComponents(1)
   minArray.SetNumberOfTuples(1)
   minArray.SetName('SectionMinSize')
   minArray.SetTuple1(0,sectionSizeRange[0])

   maxArray = vtk.vtkDoubleArray()
   maxArray.SetNumberOfComponents(1)
   maxArray.SetNumberOfTuples(1)
   maxArray.SetName('SectionMaxSize')
   maxArray.SetTuple1(0,sectionSizeRange[1])

   shapeArray = vtk.vtkDoubleArray()
   shapeArray.SetNumberOfComponents(1)
   shapeArray.SetNumberOfTuples(1)
   shapeArray.SetName('SectionShapeFactor')
   shapeArray.SetTuple1(0,sectionShape)

   section.GetCellData().AddArray(areaArray)
   section.GetCellData().AddArray(lengthArray)
   section.GetCellData().AddArray(minArray)
   section.GetCellData().AddArray(maxArray)
   section.GetCellData().AddArray(shapeArray)
   WritePolyData(section, filename)

   return point,newPlaneNormal,section



# -----------------------------------------------------------------------------------------


## Program:	extractaneurysmneckplanesection.py	
## Language:	Python
## Date:	2012/02/27
## Version:	1.0
## Application: Cerebral Aneurysms - Aneurysm Sac Morphology 
## Author:	Marina Piccinelli

## Description:	Given the surface model of the aneurysm sac and its close vasculature, the 
##		surface model of the parent artery without the aneurysm and their Voronoi
##		Diagrams, the neck planar cross section is extracted at the interface 
##		between parent artery and sac. The Voronoi Diagram and the core Voronoi 
##		Diagram of the aneurysm sac are also extracted. Data about the neck cross
##		section position and orientation are saved in a txt file.


# -----------------------------------------------------------------------------------------

# VMTK COMMON DATA ARRAY
def extractaneurysmneckplanesection_run(inputfiledirectory, ID, aneurysmType):
    radiusArrayName = 'MaximumInscribedSphereRadius'

    # ADDITIONAL DATA ARRAYS 
    clippingArrayName = 'ClippingArray'
    branchClippingArrayName = 'BranchClippingArray'
    distanceToTubeArrayName = 'DistanceToTubeFunction'
    closedArrayName = 'ClosedSection'

    # OPTIONS TO SET
    voronoiCoreCutOffThreshold = 0.75
    clipWithCore = 0
    removeBranchVessel = 0

    numberOfSplineAnalyzedPoints = 40
    phiValues = [2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,42.0]
    thetaStep = 2.0

    #print('USAGE:')
    #print('      ./extractaneurysmneckplanesection.py inputfilesDirectory caseID aneurysmType')
    #print('')

    #print('Inputfiles Directory	', inputfiledirectory)
    #print('case ID			', ID)
    #print('aneurysm Type		', aneurysmType)
    #print('')

    # input filenames
    modelfilename = inputfiledirectory + '/' + ID + '.vtp'
    if os.path.isfile(modelfilename):
        print(modelfilename)
    else:
        print('couldnt get the file', modelfilename)
        exit(1)

    modelvoronoidiagramfilename = inputfiledirectory + '/' + ID + '_voronoi.vtp'
    if os.path.isfile(modelvoronoidiagramfilename):
        print(modelvoronoidiagramfilename)
    else:
        print('couldnt get the file', modelvoronoidiagramfilename)
        exit(1)

    clippingpointsfilename = inputfiledirectory + '/' + ID + '_clippingpoints.vtp'
    if os.path.isfile(clippingpointsfilename):
        print(clippingpointsfilename)
    else:
        print('couldnt get the file', clippingpointsfilename)
        exit(1)

    parentvesselcenterlinesfilename	= inputfiledirectory + '/' + ID + '_reconstructedmodel_cl.vtp'
    if os.path.isfile(parentvesselcenterlinesfilename):
        print(parentvesselcenterlinesfilename)
    else:
        print('couldnt get the file', parentvesselcenterlinesfilename)
        exit(1)

    parentvesselvoronoidiagramfilename = inputfiledirectory + '/' + ID + '_reconstructedmodel_vor.vtp'
    if os.path.isfile(parentvesselvoronoidiagramfilename):
        print(parentvesselvoronoidiagramfilename)
    else:
        print('couldnt get the file', parentvesselvoronoidiagramfilename)
        exit(1)

    # output filenames
    tubecenterlinesfilename = inputfiledirectory + '/' + ID + '_tubecl.vtp'
    smoothedsplinefilename  = inputfiledirectory + '/' + ID + '_spline.vtp'

    aneurysmsacsurfacefilename = inputfiledirectory + '/' + ID + '_aneurysmsurface.vtp'

    aneurysmsacvoronoifilename = inputfiledirectory + '/' + ID + '_sacvoronoi.vtp'
    aneurysmnotsacvoronoifilename = inputfiledirectory + '/' + ID + '_notsacvoronoi.vtp'
    aneurysmsacvoronoicorefilename = inputfiledirectory + '/' + ID + '_saccorevoronoi.vtp'

    firstclosedsectionfilename = inputfiledirectory + '/' + ID + '_firstclosedsection.vtp'
    minsectionfilename = inputfiledirectory + '/' + ID + '_minsection.vtp'

    # creation of directory for csv files  
    csvdirectory = inputfiledirectory + '/' + 'csvfile'
    if not (os.path.isdir(csvdirectory)):
        os.mkdir(csvdirectory)

    print('Reading input files')
    model          = ReadPolyData(modelfilename)
    voronoiDiagram = ReadPolyData(modelvoronoidiagramfilename)
    clippingPoints = ReadPolyData(clippingpointsfilename)

    BJRCenterlines = ReadPolyData(parentvesselcenterlinesfilename)
    BJRVoronoi     = ReadPolyData(parentvesselvoronoidiagramfilename)

    print('Extracting centerline between clipping points')

    # Check if the new clipping ids match with previous
    #clippingIds = ExtractIds(clippingPoints,BJRCenterlines, aneurysmType)
    #test_points = [str(clippingPoints.GetPoint(i)) for i in range(clippingPoints.GetNumberOfPoints())]
    #text_points = "\n".join(test_points)
    #print("Original clipping points\n", text_points)
    #test_points = [str(BJRCenterlines.GetPoint(i)) for i in clippingIds]
    #text_points = "\n".join(test_points)
    #print("New clipping points\n", text_points)
    #sys.exit(0)
    tubeCenterlines = CreateCenterlines(clippingPoints, BJRCenterlines,
                                        aneurysmType)
    WritePolyData(tubeCenterlines, tubecenterlinesfilename)

    print('Extract aneurysm sac voronoi diagram')
    sacVoronoiDiagram = ExtractAneurysmSacVoronoiDiagram(voronoiDiagram, BJRCenterlines)
    WritePolyData(sacVoronoiDiagram, aneurysmsacvoronoifilename)

    print('Extract aneurysm sac core voronoi diagram')
    sacCoreVoronoiDiagram = ExtractAneurysmSacVoronoiCore(sacVoronoiDiagram)
    WritePolyData(sacCoreVoronoiDiagram,aneurysmsacvoronoicorefilename)

    #sys.exit()

    if (clipWithCore == 1):
        print('Clip aneurysm sac voronoi diagram with core')
        clippedSacVoronoiDiagram = ClipVoronoiDiagramWithCore(sacVoronoiDiagram,sacCoreVoronoiDiagram,1)

        notSacVoronoiDiagram = ClipVoronoiDiagramWithCore(sacVoronoiDiagram,sacCoreVoronoiDiagram,0)
        WritePolyData(notSacVoronoiDiagram,aneurysmnotsacvoronoifilename)

        sacVoronoiDiagram = ClipVoronoiDiagramWithCore(sacVoronoiDiagram,sacCoreVoronoiDiagram,1)
        WritePolyData(sacVoronoiDiagram,aneurysmsacvoronoifilename)

    print('Extract aneurysm sac surface')
    sacSurface = ExtractAneurysmSacSurface(model,sacVoronoiDiagram,BJRVoronoi)

    if (removeBranchVessel == 1):
        print('Remove Branch vessel')
        sacSurface = RemoveBranchVesselFromSac(sacSurface,notSacVoronoiDiagram,sacVoronoiDiagram)

    print('Compute aneurysm sac distance from parent vessel centerline')
    sacSurface = ComputeDistanceToParentVesselCenterlines(sacSurface,sacVoronoiDiagram,tubeCenterlines)
    WritePolyData(sacSurface,aneurysmsacsurfacefilename)

    print('Compute barycenters spline')
    sacSurface.GetPointData().SetActiveScalars(distanceToTubeArrayName)

    numberOfContours = 199
    barycenterSpline = ExtractSplineFromContours(sacSurface, numberOfContours)

    print('Smoothing spline')
    centerlineSmoothing = vtkvmtk.vtkvmtkCenterlineSmoothing()
    centerlineSmoothing.SetInputData(barycenterSpline)
    centerlineSmoothing.SetNumberOfSmoothingIterations(50)
    centerlineSmoothing.SetSmoothingFactor(1.5)
    centerlineSmoothing.Update()
    smoothedSpline = centerlineSmoothing.GetOutput()

    print('Computing spline geometry')
    centerlineGeometry = vtkvmtk.vtkvmtkCenterlineGeometry()
    centerlineGeometry.SetInputData(smoothedSpline)
    centerlineGeometry.SetLengthArrayName('Length')
    centerlineGeometry.SetCurvatureArrayName('Curvature')
    centerlineGeometry.SetTorsionArrayName('Torsion')
    centerlineGeometry.SetTortuosityArrayName('Tortuosity')
    centerlineGeometry.SetFrenetTangentArrayName('FrenetTangent')
    centerlineGeometry.SetFrenetNormalArrayName('FrenetNormal')
    centerlineGeometry.SetFrenetBinormalArrayName('FrenetBinormal')
    centerlineGeometry.SetLineSmoothing(1)
    centerlineGeometry.SetNumberOfSmoothingIterations(100)
    centerlineGeometry.SetSmoothingFactor(0.8)
    centerlineGeometry.Update()

    geometrySpline = centerlineGeometry.GetOutput()
    WritePolyData(geometrySpline,smoothedsplinefilename)

    # file for all sections along centerlines
    print('Compute sections along barycenter spline')
    datafilename = inputfiledirectory + '/csvfile/' + ID + '_profiles_datafile.csv'
    ComputeAneurysmSacSectionsAlongBarycenterSpline(geometrySpline, sacSurface, datafilename,
                                                    numberOfSplineAnalyzedPoints,
                                                    ID)

    print('Find minimum/maximum area sections')
    sectionsData = pythonprofileanalysis.ExtractMinimaAreasProfiles(datafilename,numberOfSplineAnalyzedPoints,len(phiValues),thetaStep)
    firstClosedPointID = sectionsData[0]
    firstClosedPhi = sectionsData[2]
    firstClosedTheta = sectionsData[3]
    print('First Section ', firstClosedPointID, firstClosedPhi, firstClosedTheta)

    firstPoint, firstNormal, firstClosedSection = ComputeFirstClosedSectionParameters(firstClosedPointID,
                                                                                      firstClosedPhi,
                                                                                      firstClosedTheta,
                                                                                      sacSurface,
                                                                                      geometrySpline,
                                                                                      firstclosedsectionfilename)

    if (sectionsData[4]==1):
        localMinPointID = sectionsData[5]
        localMinPhi = sectionsData[7]
        localMinTheta = sectionsData[8]
        print('Local Min Section ',localMinPointID,localMinPhi,localMinTheta)
        minPoint,minNormal,minClosedSection = ComputeFirstClosedSectionParameters(localMinPointID,
                                                                                  localMinPhi,
                                                                                  localMinTheta,
                                                                                  sacSurface,
                                                                                  geometrySpline,
                                                                                  minsectionfilename)

    # file for the neck section parameters
    parametersfilename = inputfiledirectory + '/'  + '/csvfile/' + ID +'_parameters.csv'
    parametersfile = open(parametersfilename,'w')

    print('Saving Data')
    parametersline = 'ID firstCloseId firstClosedPhi firstClosedTheta firstClosedN1 firstClosedN2 firstClosedN3 firstClosedArea firstClosedProfileLength firstClosedMinSize firstClosedMaxSize firstClosedShapeFactor minClosedId  minClosedPhi minClosedTheta minClosedN1 minClosedN2 minClosedN3 minClosedArea minClosedProfileLength minClosedMinSize minClosedMaxSize minClosedShapeFactor '+'\n'
    parametersfile.write(parametersline)

    parametersline = str(ID)+' '+str(firstClosedPointID)+' '+str(firstClosedPhi)+' '+str(firstClosedTheta)+' '+str(firstNormal[0])+' '+str(firstNormal[1])+' '+str(firstNormal[2])+' '+str(firstClosedSection.GetCellData().GetArray('SectionArea').GetTuple1(0))+' '+str(firstClosedSection.GetCellData().GetArray('SectionLength').GetTuple1(0))+' '+str(firstClosedSection.GetCellData().GetArray('SectionMinSize').GetTuple1(0))+' '+str(firstClosedSection.GetCellData().GetArray('SectionMaxSize').GetTuple1(0))+' '+str(firstClosedSection.GetCellData().GetArray('SectionShapeFactor').GetTuple1(0))+' '

    if (sectionsData[4]==1):
        parametersline = parametersline+str(localMinPointID)+' '+str(localMinPhi)+' '+str(localMinTheta)+' '+str(minNormal[0])+' '+str(minNormal[1])+' '+str(minNormal[2])+' '+str(minClosedSection.GetCellData().GetArray('SectionArea').GetTuple1(0))+' '+str(minClosedSection.GetCellData().GetArray('SectionLength').GetTuple1(0))+' '+str(minClosedSection.GetCellData().GetArray('SectionMinSize').GetTuple1(0))+' '+str(minClosedSection.GetCellData().GetArray('SectionMaxSize').GetTuple1(0))+' '+str(minClosedSection.GetCellData().GetArray('SectionShapeFactor').GetTuple1(0))+'\n'
    else:
        parametersline = parametersline+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+' '+str('NA')+'\n'

    parametersfile.write(parametersline)
    parametersfile.close()
