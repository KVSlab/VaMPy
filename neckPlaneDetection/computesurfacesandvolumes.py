#!/usr/bin/env python

import sys, os, math
import vtk
from vmtk import vtkvmtk
from vmtk import pypes
from common import *

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
   pointid = int(data[1])

   return normal,pointid

def ComputePolyBallEnvelope(voronoi):
   print('Computing Envelope')
   modeller = vtkvmtk.vtkvmtkPolyBallModeller()
   modeller.SetInputData(voronoi)
   modeller.SetRadiusArrayName(radiusArrayName)
   modeller.UsePolyBallLineOff()
   modeller.SetSampleDimensions([120, 120, 120])
   modeller.Update()

   marchingcubes = vtk.vtkMarchingCubes()
   marchingcubes.SetInputData(modeller.GetOutput())
   marchingcubes.SetValue(0,0.0)
   marchingcubes.Update()

   envelope =  marchingcubes.GetOutput()
   return envelope


def CutSurfaceWithPlane(surface, point, normal, capp=True):
   cutPlane = vtk.vtkPlane()
   cutPlane.SetOrigin(point)
   cutPlane.SetNormal(normal)

   clipper = vtk.vtkClipPolyData()
   clipper.SetInputData(surface)
   clipper.SetClipFunction(cutPlane)
   clipper.Update()

   clippedSurface = clipper.GetOutput()

   trianglefilter = vtk.vtkTriangleFilter()
   trianglefilter.SetInputData(clippedSurface)
   trianglefilter.PassLinesOff()
   trianglefilter.PassVertsOff()
   trianglefilter.Update()

   surfaceArea = ComputeSurfaceArea(trianglefilter.GetOutput())
   if capp:
       cappedSurface = CapSurface(trianglefilter.GetOutput())
   else:
       cappedSurface = trianglefilter.GetOutput()

   return surfaceArea, cappedSurface


def CapSurface(input):
   print('Capping Surface')
   entityIdsArrayName = 'CellEntityIds'
   cellEntityIdOffset = 1

   capper = vtkvmtk.vtkvmtkSmoothCapPolyData()
   capper.SetInputData(input)
   capper.SetConstraintFactor(0.0)
   capper.SetNumberOfRings(4)
   capper.Update()

   trianglefilter = vtk.vtkTriangleFilter()
   trianglefilter.SetInputData(capper.GetOutput())
   trianglefilter.Update()

   cappedInput = trianglefilter.GetOutput()

   return cappedInput


def GeneratePypeArgs(surfacefilename,gridfilename):
   pypeargs = 'vmtkdelaunayvoronoi ' + \
              '-ifile %s ' % surfacefilename + \
              '-delaunaytessellationfile %s ' % gridfilename
   return pypeargs

def ComputeVolume(grid):
    totvolume = 0.0
    for i in range(grid.GetNumberOfCells()):
        cell = vtk.vtkTetra()
        cell = grid.GetCell(i)
        volume = cell.ComputeVolume(cell.GetPoints().GetPoint(0),cell.GetPoints().GetPoint(1),cell.GetPoints().GetPoint(2),cell.GetPoints().GetPoint(3))
        totvolume +=volume
    return totvolume

def ComputeSurfaceArea(surface):
   cleaner = vtk.vtkCleanPolyData()
   cleaner.SetInputData(surface)
   cleaner.Update()

   trianglefilter = vtk.vtkTriangleFilter()
   trianglefilter.SetInputData(cleaner.GetOutput())
   trianglefilter.Update()

   envelope = trianglefilter.GetOutput()

   numberOfTriangles = envelope.GetNumberOfCells()
   totarea = 0
   for i in range(numberOfTriangles):
      cell = vtk.vtkTriangle()
      cell = envelope.GetCell(i)
      a = cell.TriangleArea(cell.GetPoints().GetPoint(0),cell.GetPoints().GetPoint(1),cell.GetPoints().GetPoint(2))
      totarea +=a
   return totarea


# -----------------------------------------------------------------------------------------


## Program:	computesurfacesandvolumes.py	
## Language:	Python
## Date:	2012/02/27
## Version:	1.0
## Application: Cerebral Aneurysms - Aneurysm Sac Morphology 
## Author:	Marina Piccinelli

## Description: The aneurysm sac volume and surface, the Voronoi Diagram core envelope 
##		volume and surface are computed. Computed values are saved in a txt file.


# -----------------------------------------------------------------------------------------

# VMTK COMMON DATA ARRAY
def computesurfacesandvolumes_run(inputfiledirectory, ID):
    radiusArrayName = 'MaximumInscribedSphereRadius'

    #print('USAGE:')
    #print('      ./computesurfacesandvolumes inputfilesDirectory caseID')
    #print('')

    #print('Inputfiles Directory   ', inputfiledirectory)
    #print('case ID		      ', ID)

    # input filenames
    aneurysmsacsurfacefilename = inputfiledirectory  + '/' + ID + '_aneurysmsurface.vtp'
    corevoronoidiagramfilename = inputfiledirectory  + '/' + ID + '_saccorevoronoi.vtp'
    sacsplinefilename          = inputfiledirectory  + '/' + ID + '_spline.vtp'
    parametersdatafilename     = inputfiledirectory  + '/csvfile/' + ID + '_parameters.csv'

    # output filenames
    sacsurfacefilename         = inputfiledirectory  + '/' + ID + '_sac.vtp'
    voronoicoresurfacefilename = inputfiledirectory  + '/' + ID + '_coreenvelope.vtp'
    sacgridfilename            = inputfiledirectory  + '/' + ID + '_sac.vtu'
    coreenvelopegridfilename   = inputfiledirectory   + '/' + ID + '_coreenvelope.vtu'

    aneurysmSac = ReadPolyData(aneurysmsacsurfacefilename)
    voronoiDiagramCore = ReadPolyData(corevoronoidiagramfilename)
    sacSpline   = ReadPolyData(sacsplinefilename)

    aneurysmNeckNormal, barycenterId = ExtractSectionNormal(parametersdatafilename)
    aneurysmNeckBarycenter = sacSpline.GetPoint(barycenterId)

    # if I already have the SAC FILE
    #openSacSurface = ReadPolyData(sacsurfacefilename)
    #sacSurfaceArea = ComputeSurfaceArea(openSacSurface)
    #sacSurface = CapSurface(openSacSurface)
    #WritePolyData(sacSurface,sacsurfacefilename)

    #original pipeline
    sacSurfaceArea, sacSurface = CutSurfaceWithPlane(aneurysmSac,aneurysmNeckBarycenter,aneurysmNeckNormal)
    WritePolyData(sacSurface,sacsurfacefilename)

    pypeargs = GeneratePypeArgs(sacsurfacefilename,sacgridfilename)
    pipe = pypes.Pype()
    pipe.ExitOnError = 0
    pipe.Arguments = pypeargs.split()
    pipe.ParseArguments()
    pipe.Execute()

    # if I don't want the cut
    #voronoiDiagramCoreEnvelope = ComputePolyBallEnvelope(voronoiDiagramCore)
    #openVoronoiDiagramCoreEnvelope = ReadPolyData(voronoicoresurfacefilename)
    #coreEnvelopeSurfaceArea    = ComputeSurfaceArea(openVoronoiDiagramCoreEnvelope)
    #voronoiDiagramCoreEnvelope = CapSurface(openVoronoiDiagramCoreEnvelope)
    #WritePolyData(voronoiDiagramCoreEnvelope,voronoicoresurfacefilename)

    #original pipeline
    voronoiDiagramCoreEnvelope = ComputePolyBallEnvelope(voronoiDiagramCore)
    coreEnvelopeSurfaceArea, voronoiDiagramCoreEnvelope = CutSurfaceWithPlane(voronoiDiagramCoreEnvelope,
                                                                              aneurysmNeckBarycenter,
                                                                              aneurysmNeckNormal)
    WritePolyData(voronoiDiagramCoreEnvelope,voronoicoresurfacefilename)

    pypeargs = GeneratePypeArgs(voronoicoresurfacefilename,coreenvelopegridfilename)
    pipe = pypes.Pype()
    pipe.ExitOnError = 0
    pipe.Arguments = pypeargs.split()
    pipe.ParseArguments()
    pipe.Execute()

    sacGrid = ReadPolyData(sacgridfilename)
    sacEnvelopeVolume = ComputeVolume(sacGrid)

    coreGrid = ReadPolyData(coreenvelopegridfilename)
    coreEnvelopeVolume = ComputeVolume(coreGrid)

    print('Saving Data')
    csvdirectory = inputfiledirectory  + '/' + 'csvfile'
    if not (os.path.isdir(csvdirectory)):
        os.mkdir(csvdirectory)

    datafile = inputfiledirectory  +'/csvfile/' + ID + '_surfacesandvolumes.csv'
    file_ = open(datafile, 'w')
    line = 'ID sacSurface sacVolume coreSurface coreVolume \n'
    file_.write(line)

    line = str(ID) + ' ' + str(sacSurfaceArea) + ' ' + str(sacEnvelopeVolume) + ' '\
            + str(coreEnvelopeSurfaceArea) + ' ' + str(coreEnvelopeVolume) + '\n'
    file_.write(line)
    file_.close()
