#!/usr/bin/env python

import os,sys,vtk,math
from vmtk import vtkvmtk
from vmtk import pypes
from vmtk import vmtkscripts
from common import *


def GeneratePypeArgs(file, ofile):
   pypeargs = 'vmtkbranchextractor -radiusarray@ MaximumInscribedSphereRadius ' + \
              '-ifile %s ' % file + \
              '--pipe vmtkbifurcationreferencesystems ' + \
              '-ofile %s ' % ofile
   return pypeargs


def ComputeBifurcationReferenceSystem(filename, ofilename):
   pypeargs = GeneratePypeArgs(filename,ofilename)
   pipe = pypes.Pype()
   pipe.ExitOnError = 0
   pipe.Arguments = pypeargs.split()
   pipe.ParseArguments()
   pipe.Execute()

   return

def SaveParentArtery(centerlines):
   numberOfCells = centerlines.GetNumberOfCells()

   cell = centerlines.GetCell(1)
   numberOfArteryPoints = cell.GetNumberOfPoints()

   artery = vtk.vtkPolyData()
   arteryPoints = vtk.vtkPoints()
   arteryCellArray = vtk.vtkCellArray()

   radiusArray = vtk.vtkDoubleArray()
   radiusArray.SetName(radiusArrayName)
   radiusArray.SetNumberOfComponents(1)
   radiusArray.SetNumberOfTuples(numberOfArteryPoints)
   radiusArray.FillComponent(0,0.0)

   count = 0
   arteryCellArray.InsertNextCell(cell.GetNumberOfPoints())
   for j in range(cell.GetNumberOfPoints()):
      arteryPoints.InsertNextPoint(cell.GetPoints().GetPoint(j))
      radiusArray.SetTuple1(count,centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(j)))
      arteryCellArray.InsertCellPoint(count)
      count+=1

   artery.SetPoints(arteryPoints)
   artery.SetLines(arteryCellArray)
   artery.GetPointData().AddArray(radiusArray)
   return artery

def AddAttributesAndSmooth(centerlines):
   centerlineAttributesFilter = vtkvmtk.vtkvmtkCenterlineAttributesFilter()
   centerlineAttributesFilter.SetInputData(centerlines)
   centerlineAttributesFilter.SetAbscissasArrayName(abscissasArrayName)
   centerlineAttributesFilter.SetParallelTransportNormalsArrayName(parallelTransportNormalsArrayName)
   centerlineAttributesFilter.Update()

   centerlineGeometryFilter = vtkvmtk.vtkvmtkCenterlineGeometry()
   centerlineGeometryFilter.SetInputData(centerlineAttributesFilter.GetOutput())
   centerlineGeometryFilter.SetLineSmoothing(1)
   centerlineGeometryFilter.SetSmoothingFactor(1.0)
   centerlineGeometryFilter.SetLengthArrayName('Length')
   centerlineGeometryFilter.SetTortuosityArrayName('Tortuosity')
   centerlineGeometryFilter.SetTorsionArrayName('Torsion')
   centerlineGeometryFilter.SetCurvatureArrayName('Curvature')
   centerlineGeometryFilter.SetFrenetTangentArrayName('FrenetTangent')
   centerlineGeometryFilter.SetFrenetBinormalArrayName('FrenetBinormal')
   centerlineGeometryFilter.SetFrenetNormalArrayName('FrenetNormal')
   centerlineGeometryFilter.Update()

   return centerlineGeometryFilter.GetOutput()

def IdentifyClippingPointsIds(centerlines, points, filename, aneurysmType):
   pointLocator = vtk.vtkPointLocator()
   pointLocator.SetDataSet(centerlines)
   pointLocator.BuildLocator()

   if (aneurysmType == 'lateral'):
      upstreamPoint         = points.GetPoint(0)
      upstreamPointId       = pointLocator.FindClosestPoint(upstreamPoint)
      upstreamPointMISR     = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(upstreamPointId)
      upstreamPointAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(upstreamPointId)

      distance = upstreamPointAbscissa - 2.0 * upstreamPointMISR

      stopId = 0
      for i in range(upstreamPointId,0,-1):
         currentAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(i)
         if (currentAbscissa - distance <= 0.0):
            stopId = i
            break
      print('2 pointids ', upstreamPointId, stopId)

      upstreamMeanRadius = 0.0
      upstreamCenterlineTangent = [0.0,0.0,0.0]
      for j in range(upstreamPointId, stopId, -1):
         currentRadius = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
         currentTangent = centerlines.GetPointData().GetArray('FrenetTangent').GetTuple3(j)

         upstreamMeanRadius += currentRadius
         upstreamCenterlineTangent[0] += currentTangent[0]
         upstreamCenterlineTangent[1] += currentTangent[1]
         upstreamCenterlineTangent[2] += currentTangent[2]

      upstreamMeanRadius = upstreamMeanRadius/(upstreamPointId-stopId)
      upstreamCenterlineTangent[0] = upstreamCenterlineTangent[0]/(upstreamPointId-stopId)
      upstreamCenterlineTangent[1] = upstreamCenterlineTangent[1]/(upstreamPointId-stopId)
      upstreamCenterlineTangent[2] = upstreamCenterlineTangent[2]/(upstreamPointId-stopId)
      SaveDirectionAsPolyData(upstreamPoint,upstreamCenterlineTangent,filename)

      pvDiameter = 2.0 * upstreamMeanRadius
      print('Parent vessel diameter ', pvDiameter)

   else:
      commonPoint         = points.GetPoint(0)
      commonPointId       = pointLocator.FindClosestPoint(commonPoint)
      commonPointMISR     = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(commonPointId)
      commonPointAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(commonPointId)

      distance = commonPointAbscissa - 2.0 * commonPointMISR

      stopId = 0
      for i in range(commonPointId,0,-1):
         currentAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(i)
         if (currentAbscissa - distance <= 0.0):
            stopId = i
            break

      commonMeanRadius = 0.0
      commonCenterlineTangent = [0.0,0.0,0.0]
      for j in range(commonPointId,stopId,-1):
         currentRadius = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
         currentTangent = centerlines.GetPointData().GetArray('FrenetTangent').GetTuple3(j)

         commonMeanRadius += currentRadius
         commonCenterlineTangent[0] += currentTangent[0]
         commonCenterlineTangent[1] += currentTangent[1]
         commonCenterlineTangent[2] += currentTangent[2]

      commonMeanRadius = commonMeanRadius/(commonPointId-stopId)
      commonCenterlineTangent[0] = commonCenterlineTangent[0]/(commonPointId-stopId)
      commonCenterlineTangent[1] = commonCenterlineTangent[1]/(commonPointId-stopId)
      commonCenterlineTangent[2] = commonCenterlineTangent[2]/(commonPointId-stopId)
      SaveDirectionAsPolyData(commonPoint,commonCenterlineTangent,filename)

      # radius data from daughter artery 1
      dau1Point         = points.GetPoint(1)
      dau1PointId       = pointLocator.FindClosestPoint(dau1Point)
      dau1PointMISR     = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(dau1PointId)
      dau1PointAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(dau1PointId)

      distance1 = dau1PointAbscissa + 2.0 * dau1PointMISR

      stopId = 0
      for i in range(dau1PointId, dau1PointId+50):
         currentAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(i)
         if (distance1 - currentAbscissa <= 0.0):
            stopId = i
            break

      dau1MeanRadius = 0.0
      for j in range(dau1PointId,stopId):
          currentRadius = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
          dau1MeanRadius += currentRadius
      dau1MeanRadius = dau1MeanRadius/(stopId-dau1PointId)

      # radius data from daughter artery 2
      dau2Point         = points.GetPoint(2)
      dau2PointId       = pointLocator.FindClosestPoint(dau2Point)
      dau2PointMISR     = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(dau2PointId)
      dau2PointAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(dau2PointId)

      distance2 = dau2PointAbscissa + 2.0 * dau2PointMISR

      stopId = 0
      for i in range(dau2PointId, dau2PointId+50):
         currentAbscissa = centerlines.GetPointData().GetArray(abscissasArrayName).GetTuple1(i)
         if (distance2 - currentAbscissa <= 0.0):
            stopId = i
            break

      dau2MeanRadius = 0.0
      for j in range(dau2PointId,stopId):
          currentRadius = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
          dau2MeanRadius += currentRadius
      dau2MeanRadius = dau2MeanRadius/(stopId-dau2PointId)

      pvDiameter = 2.0 * ((commonMeanRadius + dau1MeanRadius + dau2MeanRadius)/3.0)
      print('Parent Vessel Diameter ', pvDiameter)

   return pvDiameter

def SaveDirectionAsPolyData(point,vector,filename):
   polydata = vtk.vtkPolyData()
   points = vtk.vtkPoints()
   cellArray = vtk.vtkCellArray()

   array = vtk.vtkDoubleArray()
   array.SetName('Direction')
   array.SetNumberOfComponents(3)
   array.SetNumberOfTuples(1)
   array.SetTuple3(0,vector[0],vector[1],vector[2])

   points.InsertNextPoint(point)

   cellArray.InsertNextCell(1)
   cellArray.InsertCellPoint(0)

   polydata.SetPoints(points)
   polydata.SetVerts(cellArray)
   polydata.GetPointData().AddArray(array)

   WritePolyData(polydata, filename)
   return


def ExtractSectionNormal(filename):
   normal = [0.0,0.0,0.0]
   file = open(filename,'r')
   lines = []
   for line in file.readlines():
       lines.append(line)
   data = lines[1].split()
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


def GeneratePypeArgs1(file):
   pypeargs = 'vmtkcenterlineresampling -length 0.1 ' + \
              '-ifile %s ' % file + \
              '-ofile %s ' % file
   return pypeargs
 

def ExtractSacCenterlineDirection(sacfilename):
   pypearg = GeneratePypeArgs1(sacfilename)
   pipe = pypes.Pype()
   pipe.ExitOnError = 0
   pipe.Arguments = pypearg.split()
   pipe.ParseArguments()
   pipe.Execute()

   print("before line")
   print(sacfilename)
   line = ReadPolyData(sacfilename)
   print(line.GetNumberOfPoints())
   centerlineGeometryFilter = vtkvmtk.vtkvmtkCenterlineGeometry()
   centerlineGeometryFilter.SetInputData(line)
   centerlineGeometryFilter.SetLineSmoothing(1)
   centerlineGeometryFilter.SetSmoothingFactor(1.0)
   centerlineGeometryFilter.SetLengthArrayName('Length')
   centerlineGeometryFilter.SetTortuosityArrayName('Tortuosity')
   centerlineGeometryFilter.SetTorsionArrayName('Torsion')
   centerlineGeometryFilter.SetCurvatureArrayName('Curvature')
   centerlineGeometryFilter.SetFrenetTangentArrayName('FrenetTangent')
   centerlineGeometryFilter.SetFrenetBinormalArrayName('FrenetBinormal')
   centerlineGeometryFilter.SetFrenetNormalArrayName('FrenetNormal')
   centerlineGeometryFilter.Update()

   print("After centerline")
   numberOfPoints = centerlineGeometryFilter.GetOutput().GetNumberOfPoints()
   print(numberOfPoints)
   averageDirection = [0.0,0.0,0.0]
   for i in range(numberOfPoints):
      tangent = centerlineGeometryFilter.GetOutput().GetPointData().GetArray('FrenetTangent').GetTuple3(i)
      averageDirection[0]+=-tangent[0]
      averageDirection[1]+=-tangent[1]
      averageDirection[2]+=-tangent[2]
   averageDirection[0] = averageDirection[0]/numberOfPoints
   averageDirection[1] = averageDirection[1]/numberOfPoints
   averageDirection[2] = averageDirection[2]/numberOfPoints

   return averageDirection

def ProjectVectorOnPlane(vector, planenormal):
   newvector = [0.0,0.0,0.0]

   dot = vtk.vtkMath().Dot(vector,planenormal)

   newvector[0] = vector[0] - dot*planenormal[0]
   newvector[1] = vector[1] - dot*planenormal[1]
   newvector[2] = vector[2] - dot*planenormal[2]
   return newvector

def ProjectPointOnPlane(point,center, normal):
  vector = [0.0,0.0,0.0]
  newpoint = [0.0,0.0,0.0]

  vector[0] = center[0] - point[0]
  vector[1] = center[1] - point[1]
  vector[2] = center[2] - point[2]

  dot = vtk.vtkMath().Dot(vector,normal)

  newpoint[0] = point[0]+dot*normal[0]
  newpoint[1] = point[1]+dot*normal[1]
  newpoint[2] = point[2]+dot*normal[2]
  return newpoint

def CreateAneurysmBifurcationPlane(point, normal, filename):
    hside1 = 5.5
    hside2 = 5.5

    plane = vtk.vtkPlaneSource()
    plane.SetOrigin(-hside1,-hside1,0)
    plane.SetPoint1(-hside1,hside2,0)
    plane.SetPoint2(hside1,-hside2,0)
    plane.SetCenter(point)
    plane.SetNormal(normal)
    plane.SetResolution(100, 100)
    plane.Update()

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(plane.GetOutput())
    writer.Write()

# -----------------------------------------------------------------------------------------


## Program:	extractparentvesselparameters.py	
## Language:	Python
## Date:	2012/02/27
## Version:	1.0
## Application: Cerebral Aneurysms - Aneurysm Sac Morphology 
## Author:	Marina Piccinelli

## Description:	Given the parent vessel centerlines, the clipping points, the sac centerline, 
##		the neck plane, the relative orientation between neck, sac and parent artery 
##		are computed. Angles are saved in a txt file. Parent vessel diameter is 
##		also computed and save in a txt file.


# -----------------------------------------------------------------------------------------

# VMTK COMMON DATA ARRAY
def extractparentvesselparameters_run(inputfiledirectory, ID, aneurysmType):
    radiusArrayName = 'MaximumInscribedSphereRadius'
    parallelTransportNormalsArrayName = 'ParallelTransportNormals'

    #print('USAGE:')
    #print('      ./extractparentvesselparameters.py inputfilesDirectory caseID aneurysmType')
    #print('')

    #print('Inputfiles Directory	', inputfiledirectory)
    #print('case ID			', ID)
    #print('Aneurysm Type 		', aneurysmType)
    #print('')

    if (aneurysmType == 'lateral'):
        parentvesselfilename   = inputfiledirectory  + '/' + ID + '_forwardcl.vtp'
        clippingpointsfilename = inputfiledirectory  + '/' + ID + '_clippingpoints.vtp'
        parametersfilename     = inputfiledirectory  + '/csvfile/' + ID + '_parameters.csv'
        necksectionfilename    = inputfiledirectory + '/' + ID + '_firstclosedsection.vtp'
        saccenterlinefilename  = inputfiledirectory + '/' + ID + '_saccl.vtp'

        aneurysmreferencesystemfilename  = inputfiledirectory  + '/' + ID + '_aneurysmreferencesystem.vtp'
        aneurysmbifurcationplanefilename = inputfiledirectory  + '/' + ID + '_aneurysmplane.vtp'
        upstreamdirectionfilename        = inputfiledirectory  + '/' + ID + '_upstreamdirection.vtp'
        neckdirectionfilename            = inputfiledirectory  + '/' + ID + '_neckdirection.vtp'
        saccenterlinedirectionfilename   = inputfiledirectory  + '/' + ID + '_saccldirection.vtp'
        upstreamdirectionpfilename       = inputfiledirectory  + '/' + ID + '_upstreamdirection-p.vtp'
        neckdirectionpfilename           = inputfiledirectory  + '/' + ID + '_neckdirection-p.vtp'
        saccenterlinedirectionpfilename  = inputfiledirectory  + '/' + ID + '_saccldirection-p.vtp'

        parentVesselCenterlines = ReadPolyData(parentvesselfilename)
        ComputeBifurcationReferenceSystem(parentvesselfilename, aneurysmreferencesystemfilename)
        parentVesselCenterlines = SaveParentArtery(parentVesselCenterlines)
        clippingPoints = ReadPolyData(clippingpointsfilename)
        neckSection = ReadPolyData(necksectionfilename)

        print('')
        print('Computing Parent Vessel Diameter')
        parentVesselCenterlines = AddAttributesAndSmooth(parentVesselCenterlines)
        parentVesselDiameter = IdentifyClippingPointsIds(parentVesselCenterlines,
                                                         clippingPoints, upstreamdirectionfilename,
                                                         aneurysmType)

    else:
        parentvesselfilename   = inputfiledirectory + '/' + ID + '_centerline_par.vtp'
        clippingpointsfilename = inputfiledirectory + '/' + ID + '_clippingpoints.vtp'
        parametersfilename     = inputfiledirectory + '/csvfile/' + ID + '_parameters.csv'
        necksectionfilename    = inputfiledirectory + '/' + ID + '_firstclosedsection.vtp'
        saccenterlinefilename  = inputfiledirectory + '/' + ID + '_saccl.vtp'

        aneurysmreferencesystemfilename  = inputfiledirectory + '/' + ID + '_aneurysmreferencesystem.vtp'
        aneurysmbifurcationplanefilename = inputfiledirectory + '/' + ID + '_aneurysmplane.vtp'
        commondirectionfilename          = inputfiledirectory + '/' + ID + '_commondirection.vtp'
        neckdirectionfilename            = inputfiledirectory + '/' + ID + '_neckdirection.vtp'
        saccenterlinedirectionfilename   = inputfiledirectory + '/' + ID + '_saccldirection.vtp'
        commondirectionpfilename         = inputfiledirectory + '/' + ID + '_commondirection-p.vtp'
        neckdirectionpfilename           = inputfiledirectory + '/' + ID + '_neckdirection-p.vtp'
        saccenterlinedirectionpfilename  = inputfiledirectory + '/' + ID + '_saccldirection-p.vtp'

        parentVesselCenterlines = ReadPolyData(parentvesselfilename)
        clippingPoints = ReadPolyData(clippingpointsfilename)
        neckSection = ReadPolyData(necksectionfilename)
        ComputeBifurcationReferenceSystem(parentvesselfilename, aneurysmreferencesystemfilename)

        print('')
        print('Computing Parent Vessel Diameter')
        parentVesselCenterlines = AddAttributesAndSmooth(parentVesselCenterlines)
        parentVesselDiameter = IdentifyClippingPointsIds(parentVesselCenterlines,
                                                         clippingPoints, commondirectionfilename,
                                                         aneurysmType)

    print('')
    print('Saving neck normal direction')
    neckBarycenter = ExtractSectionBarycenter(neckSection)
    neckNormal = ExtractSectionNormal(parametersfilename)
    SaveDirectionAsPolyData(neckBarycenter,neckNormal,neckdirectionfilename)

    print('')
    print('Saving sac centerline direction')
    sacCenterlineDirection = ExtractSacCenterlineDirection(saccenterlinefilename)
    SaveDirectionAsPolyData(neckBarycenter,sacCenterlineDirection,saccenterlinedirectionfilename)

    print('')
    print('Calculating direction projections and angles')
    parameters = get_parameters(inputfiledirectory)
    center = 1./3*(np.asarray(parameters["par"]["div_point"]) +
                   np.asarray(parameters["dau1"]["div_point"]) +
                   np.asarray(parameters["dau2"]["div_point"]))
    aneurysmBifurcationReferenceSystem = ReadPolyData(aneurysmreferencesystemfilename)
    numPoints = aneurysmBifurcationReferenceSystem.GetNumberOfPoints()
    tmp_bif_points = []
    for i in range(numPoints):
        tmp_bif_points.append(aneurysmBifurcationReferenceSystem.GetPoint(i))
    tmp_index = np.argsort(np.sum((np.asarray(tmp_bif_points) - center)**2, axis=1))[0]
    aneurysmBifurcationPoint = tmp_bif_points[int(tmp_index)]
    aneurysmBifurcationNormal = aneurysmBifurcationReferenceSystem.GetPointData().GetArray('Normal').GetTuple3(tmp_index)

    CreateAneurysmBifurcationPlane(aneurysmBifurcationPoint, aneurysmBifurcationNormal,
                                   aneurysmbifurcationplanefilename)

    neckBarycenterP = ProjectPointOnPlane(neckBarycenter, aneurysmBifurcationPoint, aneurysmBifurcationNormal)
    sacCenterlineDirectionP = ProjectVectorOnPlane(sacCenterlineDirection,aneurysmBifurcationNormal)
    vtk.vtkMath().Normalize(sacCenterlineDirectionP)
    SaveDirectionAsPolyData(neckBarycenterP,sacCenterlineDirectionP,saccenterlinedirectionpfilename)
    neckNormalP = ProjectVectorOnPlane(neckNormal,aneurysmBifurcationNormal)
    vtk.vtkMath().Normalize(neckNormalP)
    SaveDirectionAsPolyData(neckBarycenterP,neckNormalP,neckdirectionpfilename)

    if (aneurysmType == 'lateral'):
        upstreamDirectionPD = ReadPolyData(upstreamdirectionfilename)
        upstreamPoint = upstreamDirectionPD.GetPoint(0)
        upstreamDirection = upstreamDirectionPD.GetPointData().GetArray('Direction').GetTuple3(0)
        parentVesselDirectionP = ProjectVectorOnPlane(upstreamDirection,aneurysmBifurcationNormal)
        vtk.vtkMath().Normalize(parentVesselDirectionP)
        upstreamPointP = ProjectPointOnPlane(upstreamPoint,aneurysmBifurcationPoint,aneurysmBifurcationNormal)
        SaveDirectionAsPolyData(upstreamPointP,parentVesselDirectionP,upstreamdirectionpfilename)
    else:
        commonDirectionPD  = ReadPolyData(commondirectionfilename)
        commonPoint = commonDirectionPD.GetPoint(0)
        commonDirection = commonDirectionPD.GetPointData().GetArray('Direction').GetTuple3(0)
        parentVesselDirectionP = ProjectVectorOnPlane(commonDirection,aneurysmBifurcationNormal)
        vtk.vtkMath().Normalize(parentVesselDirectionP)
        commonPointP = ProjectPointOnPlane(commonPoint,aneurysmBifurcationPoint,aneurysmBifurcationNormal)
        SaveDirectionAsPolyData(commonPointP,parentVesselDirectionP,commondirectionpfilename)

    theta_sac_vessel = math.degrees(vtkvmtk.vtkvmtkMath.AngleBetweenNormals\
                                      (sacCenterlineDirectionP,parentVesselDirectionP))
    theta_neck_vessel = math.degrees(vtkvmtk.vtkvmtkMath.AngleBetweenNormals\
                                      (neckNormalP,parentVesselDirectionP))

    print('')
    print('Saving Data')
    csvdirectory = inputfiledirectory  + '/' + 'csvfile'
    if not (os.path.isdir(csvdirectory)):
        os.mkdir(csvdirectory)

    parentvesseldatafilename = inputfiledirectory +'/csvfile/'+ID+'_parentvessel.csv'
    anglesdatafilename       = inputfiledirectory +'/csvfile/'+ID+'_angles.csv'

    file1 = open(parentvesseldatafilename,'w')
    line1 = 'ID parentVesselDiameter'+'\n'
    file1.write(line1)
    line1 = str(ID)+' '+str(parentVesselDiameter)+'\n'
    file1.write(line1)
    file1.close()

    file2 = open(anglesdatafilename,'w')
    line2 = 'ID thetaSacVessel thetaNeckVessel'+'\n'
    file2.write(line2)
    line2 = str(ID)+' '+str(theta_sac_vessel)+' '+str(theta_neck_vessel)+'\n'
    file2.write(line2)
    file2.close()
