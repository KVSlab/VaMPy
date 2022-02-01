import argparse

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def convert_to_polydata(model, nr, index):
    # create a new 'Xdmf3ReaderS'
    rRTxdmf = Xdmf3ReaderS(registrationName='{}.xdmf'.format(index), FileName=[
        '/Users/henriakj/PhD/OasisMove/results_moving_atrium/{}/data/{}/Solutions/{}.xdmf'.format(model, nr, index)])
    rRTxdmf.PointArrays = ['{}'.format(index)]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    rRTxdmfDisplay = Show(rRTxdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for '{}'.format(index)
    rRTLUT = GetColorTransferFunction('{}'.format(index))

    # get opacity transfer function/opacity map for '{}'.format(index)
    rRTPWF = GetOpacityTransferFunction('{}'.format(index))

    # trace defaults for the display properties.
    rRTxdmfDisplay.Representation = 'Surface'
    rRTxdmfDisplay.ColorArrayName = ['POINTS', '{}'.format(index)]
    rRTxdmfDisplay.LookupTable = rRTLUT
    rRTxdmfDisplay.SelectTCoordArray = 'None'
    rRTxdmfDisplay.SelectNormalArray = 'None'
    rRTxdmfDisplay.SelectTangentArray = 'None'
    rRTxdmfDisplay.OSPRayScaleArray = '{}'.format(index)
    rRTxdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    rRTxdmfDisplay.SelectOrientationVectors = 'None'
    rRTxdmfDisplay.ScaleFactor = 13.36917953491211
    rRTxdmfDisplay.SelectScaleArray = '{}'.format(index)
    rRTxdmfDisplay.GlyphType = 'Arrow'
    rRTxdmfDisplay.GlyphTableIndexArray = '{}'.format(index)
    rRTxdmfDisplay.GaussianRadius = 0.6684589767456055
    rRTxdmfDisplay.SetScaleArray = ['POINTS', '{}'.format(index)]
    rRTxdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    rRTxdmfDisplay.OpacityArray = ['POINTS', '{}'.format(index)]
    rRTxdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    rRTxdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    rRTxdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    rRTxdmfDisplay.ScalarOpacityFunction = rRTPWF
    rRTxdmfDisplay.ScalarOpacityUnitDistance = 7.29903693454099
    rRTxdmfDisplay.OpacityArrayName = ['POINTS', '{}'.format(index)]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    rRTxdmfDisplay.ScaleTransferFunction.Points = [1.165113091468811, 0.0, 0.5, 0.0, 135221.09375, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    rRTxdmfDisplay.OpacityTransferFunction.Points = [1.165113091468811, 0.0, 0.5, 0.0, 135221.09375, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    rRTxdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Extract Surface'
    extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=rRTxdmf)

    # show data in view
    extractSurface1Display = Show(extractSurface1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    extractSurface1Display.Representation = 'Surface'
    extractSurface1Display.ColorArrayName = ['POINTS', '{}'.format(index)]
    extractSurface1Display.LookupTable = rRTLUT
    extractSurface1Display.SelectTCoordArray = 'None'
    extractSurface1Display.SelectNormalArray = 'None'
    extractSurface1Display.SelectTangentArray = 'None'
    extractSurface1Display.OSPRayScaleArray = '{}'.format(index)
    extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSurface1Display.SelectOrientationVectors = 'None'
    extractSurface1Display.ScaleFactor = 13.36917953491211
    extractSurface1Display.SelectScaleArray = '{}'.format(index)
    extractSurface1Display.GlyphType = 'Arrow'
    extractSurface1Display.GlyphTableIndexArray = '{}'.format(index)
    extractSurface1Display.GaussianRadius = 0.6684589767456055
    extractSurface1Display.SetScaleArray = ['POINTS', '{}'.format(index)]
    extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.OpacityArray = ['POINTS', '{}'.format(index)]
    extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSurface1Display.ScaleTransferFunction.Points = [1.165113091468811, 0.0, 0.5, 0.0, 135221.09375, 1.0, 0.5,
                                                           0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSurface1Display.OpacityTransferFunction.Points = [1.165113091468811, 0.0, 0.5, 0.0, 135221.09375, 1.0, 0.5,
                                                             0.0]

    # hide data in view
    Hide(rRTxdmf, renderView1)

    # show color bar/color legend
    extractSurface1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # save data
    SaveData('/Users/henriakj/PhD/OasisMove/results_moving_atrium/{}/data/{}/Solutions/{}.vtp'.format(model, nr, index),
             proxy=extractSurface1, PointDataArrays=['{}'.format(index)])

    # ================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    # ================================================================

    # get layout
    layout1 = GetLayout()

    # --------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(2490, 1688)

    # -----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [14.070598602294922, 88.01505279541016, 214.50521989097683]
    renderView1.CameraFocalPoint = [14.070598602294922, 88.01505279541016, -186.3010025024414]
    renderView1.CameraParallelScale = 103.73628375101308

    # --------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--case")
    parser.add_argument("--nr")
    parser.add_argument("--id")
    args = parser.parse_args()
    case = args.case
    number = args.nr
    index = args.id

    convert_to_polydata(case, number, index)
