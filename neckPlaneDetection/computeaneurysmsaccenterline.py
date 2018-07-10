#!/usr/bin/env python

import sys,math,os
import vtk
from vmtk import vtkvmtk
from common import *


def ExtractSectionNormal(filename):
   normal = [0.0,0.0,0.0]
   file_ = open(filename,'r')
   lines = file_.readlines()
   print(lines)
   data = lines[1].split(" ")
   print(data)
   normal[0] = float(data[4])
   normal[1] = float(data[5])
   normal[2] = float(data[6])

   return normal

def ExtractSectionBarycenter(section):
   barycenter = [0.0,0.0,0.0]

   for i in range(section.GetNumberOfPoints()):
      point = section.GetPoint(i)
      barycenter[0] = barycenter[0] + point[0]
      barycenter[1] = barycenter[1] + point[1]
      barycenter[2] = barycenter[2] + point[2]

   barycenter[0] = barycenter[0]/section.GetNumberOfPoints()
   barycenter[1] = barycenter[1]/section.GetNumberOfPoints()
   barycenter[2] = barycenter[2]/section.GetNumberOfPoints()
   return barycenter

def ComputeTargetPointOnCore(core,sourceId):
   voronoiSourceSeeds = vtk.vtkIdList()
   voronoiSourceSeeds.InsertNextId(sourceId) 

   voronoiFastMarching = vtkvmtk.vtkvmtkNonManifoldFastMarching()
   voronoiFastMarching.SetInputData(core)
   voronoiFastMarching.UnitSpeedOn()
   voronoiFastMarching.SetSolutionArrayName(eikonalSolutionArrayName)
   voronoiFastMarching.SeedsBoundaryConditionsOn()
   voronoiFastMarching.SetSeeds(voronoiSourceSeeds)
   voronoiFastMarching.Update()

   fastmarchingcore = voronoiFastMarching.GetOutput()

   eikonalSolutionArray = fastmarchingcore.GetPointData().GetArray(eikonalSolutionArrayName)
   coreRadiusArray = fastmarchingcore.GetPointData().GetArray(radiusArrayName)

   maxEikonalSolutionId = -1
   maxEikonalSolutionValue = 0.0
   for i in range(core.GetNumberOfPoints()):
       if eikonalSolutionArray.GetValue(i) + coreRadiusArray.GetValue(i) > maxEikonalSolutionValue:
           maxEikonalSolutionValue = eikonalSolutionArray.GetValue(i) + coreRadiusArray.GetValue(i)
           maxEikonalSolutionId = i

   voronoiTargetSeeds = vtk.vtkIdList()
   voronoiTargetSeeds.InsertNextId(maxEikonalSolutionId)
 
   return maxEikonalSolutionId


def ComputeEikonalCenterline(voronoi,sourceId,targetId):
   sourceSeeds = vtk.vtkIdList()
   sourceSeeds.InsertNextId(sourceId) 

   targetSeeds = vtk.vtkIdList()
   targetSeeds.InsertNextId(targetId)

   voronoiFastMarching = vtkvmtk.vtkvmtkNonManifoldFastMarching()
   voronoiFastMarching.SetInputData(voronoi)
   voronoiFastMarching.UnitSpeedOn()
   voronoiFastMarching.SetSolutionArrayName(eikonalSolutionArrayName)
   voronoiFastMarching.SeedsBoundaryConditionsOn()
   voronoiFastMarching.SetSeeds(sourceSeeds)
   voronoiFastMarching.Update()
  
   solutionVoronoi = voronoiFastMarching.GetOutput()

   radiusArray = solutionVoronoi.GetPointData().GetArray(radiusArrayName)

   costFunctionArray = vtk.vtkDoubleArray()
   costFunctionArray.SetName(costFunctionArrayName)
   costFunctionArray.SetNumberOfTuples(solutionVoronoi.GetNumberOfPoints())

   for i in range(solutionVoronoi.GetNumberOfPoints()):
      costFunctionArray.SetValue(i,1.0/radiusArray.GetValue(i))

   solutionVoronoi.GetPointData().AddArray(costFunctionArray)
    
   voronoiFastMarching = vtkvmtk.vtkvmtkNonManifoldFastMarching()
   voronoiFastMarching.SetInputData(solutionVoronoi)
   voronoiFastMarching.SetCostFunctionArrayName(costFunctionArrayName)
   voronoiFastMarching.SetSolutionArrayName(eikonalSolutionArrayName)
   voronoiFastMarching.SeedsBoundaryConditionsOn()
   voronoiFastMarching.SetSeeds(sourceSeeds)
   voronoiFastMarching.Update()

   voronoi = voronoiFastMarching.GetOutput()

   centerlineBacktracing = vtkvmtk.vtkvmtkSteepestDescentLineTracer()
   centerlineBacktracing.SetInputData(voronoi)
   centerlineBacktracing.SetDataArrayName(radiusArrayName)
   centerlineBacktracing.SetDescentArrayName(eikonalSolutionArrayName)
   centerlineBacktracing.SetEdgeArrayName(edgeArrayName)
   centerlineBacktracing.SetEdgePCoordArrayName(edgePCoordArrayName)
   centerlineBacktracing.SetSeeds(targetSeeds)
   centerlineBacktracing.MergePathsOff()
   centerlineBacktracing.StopOnTargetsOn()
   centerlineBacktracing.SetTargets(sourceSeeds)
   centerlineBacktracing.Update()

   centerlineSmoothing = vtkvmtk.vtkvmtkCenterlineSmoothing()
   centerlineSmoothing.SetInputData(centerlineBacktracing.GetOutput())
   centerlineSmoothing.SetNumberOfSmoothingIterations(20)
   centerlineSmoothing.SetSmoothingFactor(0.5)
   centerlineSmoothing.Update()

   centerline = centerlineSmoothing.GetOutput()
   centerlineLastPointRadius = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(0)

   return centerline, centerlineLastPointRadius

def ComputeCenterlineLength(centerline):
   length = 0.0 
   for i in range(centerline.GetNumberOfPoints()-1):
      point0 = centerline.GetPoint(i)
      point1 = centerline.GetPoint(i+1)
     
      d = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
      length = length+d

   return length  


# -----------------------------------------------------------------------------------------


## Program:	computeaneurysmsaccenterline.py	
## Language:	Python
## Date:	2012/02/27
## Version:	1.0
## Application: Cerebral Aneurysms - Aneurysm Sac Morphology 
## Author:	Marina Piccinelli

## Description: The aneurysm sac centerline is computed from the neck cross section 
##		barycenter on the Voronoi Diagram core previously extracted. The centerline  
##		vtp file is saved; centerline length is computed and saved in txt file. 


# -----------------------------------------------------------------------------------------

# VMTK COMMON DATA ARRAY
def computeaneurysmsaccenterline_run(inputfiledirectory, ID):
    #print('USAGE:')
    #print('      ./computeaneurysmsaccenterline.py inputfilesDirectory caseID')
    #print('')

    #print('Inputfiles Directory	', inputfiledirectory)
    #print('case ID			', ID)
    #print('')

    # input filenames
    voronoidiagramfilename     = inputfiledirectory +  '/' + ID + '_sacvoronoi.vtp'
    voronoidiagramcorefilename = inputfiledirectory +  '/' + ID + '_saccorevoronoi.vtp'
    firstclosedsectionfilename = inputfiledirectory +  '/' + ID + '_firstclosedsection.vtp'
    parametersfilename         = inputfiledirectory +  '/csvfile/' + ID + '_parameters.csv'

    # output filenames
    aneurysmsaccenterlinefilename = inputfiledirectory  + '/' + ID + '_saccl.vtp'

    voronoiDiagram = ReadPolyData(voronoidiagramfilename)
    voronoiDiagramCore = ReadPolyData(voronoidiagramcorefilename)
    aneurysmNeck = ReadPolyData(firstclosedsectionfilename)

    aneurysmNeckBarycenter = ExtractSectionBarycenter(aneurysmNeck)
    aneurysmNeckNormal = ExtractSectionNormal(parametersfilename)

    coreLocator = vtk.vtkPointLocator()
    coreLocator.SetDataSet(voronoiDiagramCore)
    coreLocator.BuildLocator()
    coreSourceId = coreLocator.FindClosestPoint(aneurysmNeckBarycenter)
    coreTargetId = ComputeTargetPointOnCore(voronoiDiagramCore,coreSourceId)

    voronoiLocator = vtk.vtkPointLocator()
    voronoiLocator.SetDataSet(voronoiDiagram)
    voronoiLocator.BuildLocator()
    voronoiSourceId = voronoiLocator.FindClosestPoint(aneurysmNeckBarycenter)
    voronoiTargetId = voronoiLocator.FindClosestPoint(voronoiDiagramCore.GetPoint(coreTargetId))

    print('Extracting centerline')
    sacCenterline,lastPointRadius = ComputeEikonalCenterline(voronoiDiagram,voronoiSourceId,voronoiTargetId)
    WritePolyData(sacCenterline,aneurysmsaccenterlinefilename)

    sacCenterlineLength = ComputeCenterlineLength(sacCenterline)

    lastSacCenterlinePoint = sacCenterline.GetPoint(0)

    print('Saving Data')
    csvdirectory = inputfiledirectory +  '/' + 'csvfile'
    if not (os.path.isdir(csvdirectory)):
        os.mkdir(csvdirectory)

    csvfilename = inputfiledirectory  + '/csvfile/' + ID + '_centerline.csv'
    file_ = open(csvfilename,'w')
    line = 'ID sacCenterlineLength sacLastMISR'+ '\n'
    file_.write(line)
    line = ID +' '+str(sacCenterlineLength)+' '+str(lastPointRadius)+'\n'
    file_.write(line)
    file_.close()
