import sys
import math
import vtk

from vmtk import vtkvmtk

def FindMinimumAreaSection(areaArray,closeArray):
   numberOfElements = len(closeArray)
   min_ = 10e100
   index = 'nan'
   for i in range(numberOfElements):
      value = float(areaArray[i])
      close = closeArray[i]
      notanumber = math.isnan(value)
      if (notanumber == False) and (close == 1):
         if (value < min_):
            min_ = value
            index = i

   return min_, index


def ExtractFirstClosedSectionIndeces(areas,closed,phis,thetas):
   firstArea = -1
   firstAreaPoint = -1
   firstAreaPhi = -1
   firstAreaTheta = -1
   numberOfElements = len(areas)
   for i in range(numberOfElements):
      if (closed[i]==1):
         firstArea = areas[i]
         firstAreaPoint = i
         firstAreaPhi = phis[i]
         firstAreaTheta = thetas[i]
         break
   return firstArea, firstAreaPoint, firstAreaPhi, firstAreaTheta


def FindLocalMinimumClosedSectionIndeces(areas,closed,phis,thetas,firstId):
    minLocal = -1
    minLocalArea = areas[firstId]
    minLocalPoint = -1
    minLocalPhi = -1
    minLocalTheta = -1

    for i in range(firstId+1, min(len(closed)-firstId, firstId+10)):
        if (closed[i]==1):
            if areas[i]< minLocalArea:
                minLocal = 1
                minLocalArea = areas[i]
                minLocalPoint = i
                minLocalPhi = phis[i]
                minLocalTheta = thetas[i]

    return minLocal, minLocalArea, minLocalPoint, minLocalPhi, minLocalTheta


def FindMaxMinimumAreas(areas, closed, phis, thetas, firstId):
   max = -1
   maxArea = areas[firstId]
   maxPoint = -1
   maxPhi = -1
   maxTheta = -1
   numberOfElements = len(areas)
   for i in range(firstId+1,numberOfElements):
      if (closed[i]==1):
         if (areas[i]>maxArea):
            max = 1
            maxArea = areas[i]
            maxPoint = i
            maxPhi = phis[i]
            maxTheta = thetas[i]

   return max, maxArea, maxPoint, maxPhi, maxTheta


def ExtractMinimaAreasProfiles(profiledatafilename, numberOfPoints, numberOfPhis, thetaStep):
   profiledatafile = open(profiledatafilename,'r')

   results = []
   alllines = []
   for line in profiledatafile.readlines():
      alllines.append(line)

   numberOfTilts = 360.0/thetaStep

   minimaAreas = []
   minimaAreasPhis = []
   minimaAreasThetas = []
   minimaAreasClosed = []

   numberoflines = numberOfTilts * numberOfPhis + 1

   for point in range(numberOfPoints):
      initialLinesId = int(point * numberoflines +1)
      finalLinesId = int((point+1) * numberoflines +1)

      theta = []
      phi = []
      closed = []
      areas = []
      length = []
      minSize = []
      maxSize = []
      shape = []

      for id in range(initialLinesId,finalLinesId):
         line = alllines[id].split()
         if (line[5]=='NA'):
            phi.append('nan')
            theta.append('nan')
            areas.append('nan')
            length.append('nan')
            closed.append('nan')
            minSize.append('nan')
            maxSize.append('nan')
            shape.append('nan')
         else:
            phi.append(float(line[2]))
            theta.append(float(line[3]))
            closed.append(float(line[4]))
            areas.append(float(line[5]))
            length.append(float(line[6]))
            minSize.append(float(line[7]))
            maxSize.append(float(line[8]))
            shape.append(float(line[9]))

      minArea,index = FindMinimumAreaSection(areas,closed)
#      print(minArea, index

      if index == 'nan':
         minimaAreas.append('nan')
         minimaAreasPhis.append('nan')
         minimaAreasThetas.append('nan')
         minimaAreasClosed.append('nan')
      else:
         minimaAreas.append(minArea)
         minimaAreasPhis.append(phi[index])
         minimaAreasThetas.append(theta[index])
         minimaAreasClosed.append(1)

   firstArea,firstPoint,firstPhi,firstTheta = ExtractFirstClosedSectionIndeces(minimaAreas,minimaAreasClosed,minimaAreasPhis,minimaAreasThetas)

   results.append(firstPoint)
   results.append(firstArea)
   results.append(firstPhi)
   results.append(firstTheta)

   localMin,localMinArea,localMinPoint,localMinPhi,localMinTheta = FindLocalMinimumClosedSectionIndeces(minimaAreas,minimaAreasClosed,minimaAreasPhis,minimaAreasThetas,firstPoint)

   results.append(localMin)
   results.append(localMinPoint)
   results.append(localMinArea)
   results.append(localMinPhi)
   results.append(localMinTheta)

   max,maxArea,maxPoint,maxPhi,maxTheta = FindMaxMinimumAreas(minimaAreas,minimaAreasClosed,minimaAreasPhis,minimaAreasThetas,firstPoint)

   results.append(max)
   results.append(maxPoint)
   results.append(maxArea)
   results.append(maxPhi)
   results.append(maxTheta)

   return results
