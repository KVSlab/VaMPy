from os import path

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, cycle):
    # Local
    solution_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics/'

    # On Oracle
    #solution_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/Hemodynamics"

    eCAP_cycle_05xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_05.xdmf', FileName=[
        path.join(solution_path, f'ECAP_cycle_{cycle:02d}.xdmf')])
    eCAP_cycle_05xdmf.PointArrays = ['ECAP']
    oSI_cycle_05xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_05.xdmf', FileName=[
        path.join(solution_path, f'OSI_cycle_{cycle:02d}.xdmf')])
    oSI_cycle_05xdmf.PointArrays = ['OSI']
    rRT_cycle_05xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_05.xdmf', FileName=[
        path.join(solution_path, f'RRT_cycle_{cycle:02d}.xdmf')])
    rRT_cycle_05xdmf.PointArrays = ['RRT']
    tAWSS_cycle_05xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_05.xdmf', FileName=[
        path.join(solution_path, f'TAWSS_cycle_{cycle:02d}.xdmf')])
    tAWSS_cycle_05xdmf.PointArrays = ['TAWSS']
    tWSSG_cycle_05xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_05.xdmf', FileName=[
        path.join(solution_path, f'TWSSG_cycle_{cycle:02d}.xdmf')])
    tWSSG_cycle_05xdmf.PointArrays = ['TWSSG']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    oSI_cycle_05xdmfDisplay = Show(oSI_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'OSI'
    oSILUT = GetColorTransferFunction('OSI')

    # get opacity transfer function/opacity map for 'OSI'
    oSIPWF = GetOpacityTransferFunction('OSI')

    # trace defaults for the display properties.
    oSI_cycle_05xdmfDisplay.Representation = 'Surface'
    oSI_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'OSI']
    oSI_cycle_05xdmfDisplay.LookupTable = oSILUT
    oSI_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    oSI_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    oSI_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    oSI_cycle_05xdmfDisplay.OSPRayScaleArray = 'OSI'
    oSI_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    oSI_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    oSI_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    oSI_cycle_05xdmfDisplay.SelectScaleArray = 'OSI'
    oSI_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    oSI_cycle_05xdmfDisplay.GlyphTableIndexArray = 'OSI'
    oSI_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    oSI_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'OSI']
    oSI_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    oSI_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'OSI']
    oSI_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    oSI_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    oSI_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    oSI_cycle_05xdmfDisplay.ScalarOpacityFunction = oSIPWF
    oSI_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 2.440939567077204
    oSI_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'OSI']
    oSI_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    oSI_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    oSI_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    oSI_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [-0.012716573663055897, 0.0, 0.5, 0.0, 0.5080192685127258,
                                                            1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    oSI_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [-0.012716573663055897, 0.0, 0.5, 0.0, 0.5080192685127258,
                                                              1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    oSI_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_05xdmfDisplay = Show(eCAP_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'ECAP'
    eCAPLUT = GetColorTransferFunction('ECAP')

    # get opacity transfer function/opacity map for 'ECAP'
    eCAPPWF = GetOpacityTransferFunction('ECAP')

    # trace defaults for the display properties.
    eCAP_cycle_05xdmfDisplay.Representation = 'Surface'
    eCAP_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'ECAP']
    eCAP_cycle_05xdmfDisplay.LookupTable = eCAPLUT
    eCAP_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    eCAP_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    eCAP_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    eCAP_cycle_05xdmfDisplay.OSPRayScaleArray = 'ECAP'
    eCAP_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    eCAP_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    eCAP_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    eCAP_cycle_05xdmfDisplay.SelectScaleArray = 'ECAP'
    eCAP_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    eCAP_cycle_05xdmfDisplay.GlyphTableIndexArray = 'ECAP'
    eCAP_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    eCAP_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'ECAP']
    eCAP_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    eCAP_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'ECAP']
    eCAP_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    eCAP_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    eCAP_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    eCAP_cycle_05xdmfDisplay.ScalarOpacityFunction = eCAPPWF
    eCAP_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 2.440939567077204
    eCAP_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'ECAP']
    eCAP_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    eCAP_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    eCAP_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    eCAP_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [-0.35945209860801697, 0.0, 0.5, 0.0, 159.9031219482422,
                                                             1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    eCAP_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [-0.35945209860801697, 0.0, 0.5, 0.0, 159.9031219482422,
                                                               1.0, 0.5, 0.0]

    # show color bar/color legend
    eCAP_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_05xdmfDisplay = Show(rRT_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'RRT'
    rRTLUT = GetColorTransferFunction('RRT')

    # get opacity transfer function/opacity map for 'RRT'
    rRTPWF = GetOpacityTransferFunction('RRT')

    # trace defaults for the display properties.
    rRT_cycle_05xdmfDisplay.Representation = 'Surface'
    rRT_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'RRT']
    rRT_cycle_05xdmfDisplay.LookupTable = rRTLUT
    rRT_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    rRT_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    rRT_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    rRT_cycle_05xdmfDisplay.OSPRayScaleArray = 'RRT'
    rRT_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    rRT_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    rRT_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    rRT_cycle_05xdmfDisplay.SelectScaleArray = 'RRT'
    rRT_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    rRT_cycle_05xdmfDisplay.GlyphTableIndexArray = 'RRT'
    rRT_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    rRT_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'RRT']
    rRT_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    rRT_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'RRT']
    rRT_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    rRT_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    rRT_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    rRT_cycle_05xdmfDisplay.ScalarOpacityFunction = rRTPWF
    rRT_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 2.440939567077204
    rRT_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'RRT']
    rRT_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    rRT_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    rRT_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    rRT_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [-7274674.0, 0.0, 0.5, 0.0, 3628691.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    rRT_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [-7274674.0, 0.0, 0.5, 0.0, 3628691.0, 1.0, 0.5, 0.0]

    # show color bar/color legend
    rRT_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_05xdmfDisplay = Show(tAWSS_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'TAWSS'
    tAWSSLUT = GetColorTransferFunction('TAWSS')

    # get opacity transfer function/opacity map for 'TAWSS'
    tAWSSPWF = GetOpacityTransferFunction('TAWSS')

    # trace defaults for the display properties.
    tAWSS_cycle_05xdmfDisplay.Representation = 'Surface'
    tAWSS_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'TAWSS']
    tAWSS_cycle_05xdmfDisplay.LookupTable = tAWSSLUT
    tAWSS_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    tAWSS_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    tAWSS_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    tAWSS_cycle_05xdmfDisplay.OSPRayScaleArray = 'TAWSS'
    tAWSS_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    tAWSS_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    tAWSS_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    tAWSS_cycle_05xdmfDisplay.SelectScaleArray = 'TAWSS'
    tAWSS_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    tAWSS_cycle_05xdmfDisplay.GlyphTableIndexArray = 'TAWSS'
    tAWSS_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    tAWSS_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'TAWSS']
    tAWSS_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    tAWSS_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'TAWSS']
    tAWSS_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    tAWSS_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    tAWSS_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    tAWSS_cycle_05xdmfDisplay.ScalarOpacityFunction = tAWSSPWF
    tAWSS_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 2.440939567077204
    tAWSS_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'TAWSS']
    tAWSS_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    tAWSS_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    tAWSS_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tAWSS_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [0.003126617521047592, 0.0, 0.5, 0.0, 3.4618287086486816,
                                                              1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tAWSS_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [0.003126617521047592, 0.0, 0.5, 0.0, 3.4618287086486816,
                                                                1.0, 0.5, 0.0]

    # show color bar/color legend
    tAWSS_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_05xdmfDisplay = Show(tWSSG_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'TWSSG'
    tWSSGLUT = GetColorTransferFunction('TWSSG')

    # get opacity transfer function/opacity map for 'TWSSG'
    tWSSGPWF = GetOpacityTransferFunction('TWSSG')

    # trace defaults for the display properties.
    tWSSG_cycle_05xdmfDisplay.Representation = 'Surface'
    tWSSG_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'TWSSG']
    tWSSG_cycle_05xdmfDisplay.LookupTable = tWSSGLUT
    tWSSG_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    tWSSG_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    tWSSG_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    tWSSG_cycle_05xdmfDisplay.OSPRayScaleArray = 'TWSSG'
    tWSSG_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    tWSSG_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    tWSSG_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    tWSSG_cycle_05xdmfDisplay.SelectScaleArray = 'TWSSG'
    tWSSG_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    tWSSG_cycle_05xdmfDisplay.GlyphTableIndexArray = 'TWSSG'
    tWSSG_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    tWSSG_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'TWSSG']
    tWSSG_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    tWSSG_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'TWSSG']
    tWSSG_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    tWSSG_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    tWSSG_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    tWSSG_cycle_05xdmfDisplay.ScalarOpacityFunction = tWSSGPWF
    tWSSG_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 2.440939567077204
    tWSSG_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'TWSSG']
    tWSSG_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    tWSSG_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    tWSSG_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tWSSG_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [0.002579913940280676, 0.0, 0.5, 0.0, 3.2623963356018066,
                                                              1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tWSSG_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [0.002579913940280676, 0.0, 0.5, 0.0, 3.2623963356018066,
                                                                1.0, 0.5, 0.0]

    # show color bar/color legend
    tWSSG_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get 2D transfer function for 'TWSSG'
    tWSSGTF2D = GetTransferFunction2D('TWSSG')

    # set active source
    SetActiveSource(eCAP_cycle_05xdmf)

    # get 2D transfer function for 'ECAP'
    eCAPTF2D = GetTransferFunction2D('ECAP')

    # set active source
    SetActiveSource(tWSSG_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[eCAP_cycle_05xdmf, oSI_cycle_05xdmf, rRT_cycle_05xdmf,
                                                tAWSS_cycle_05xdmf, tWSSG_cycle_05xdmf])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'
    appendAttributes1Display.ColorArrayName = ['POINTS', 'ECAP']
    appendAttributes1Display.LookupTable = eCAPLUT
    appendAttributes1Display.SelectTCoordArray = 'None'
    appendAttributes1Display.SelectNormalArray = 'None'
    appendAttributes1Display.SelectTangentArray = 'None'
    appendAttributes1Display.OSPRayScaleArray = 'ECAP'
    appendAttributes1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    appendAttributes1Display.SelectOrientationVectors = 'None'
    appendAttributes1Display.ScaleFactor = 8.22499008178711
    appendAttributes1Display.SelectScaleArray = 'ECAP'
    appendAttributes1Display.GlyphType = 'Arrow'
    appendAttributes1Display.GlyphTableIndexArray = 'ECAP'
    appendAttributes1Display.GaussianRadius = 0.41124950408935546
    appendAttributes1Display.SetScaleArray = ['POINTS', 'ECAP']
    appendAttributes1Display.ScaleTransferFunction = 'PiecewiseFunction'
    appendAttributes1Display.OpacityArray = ['POINTS', 'ECAP']
    appendAttributes1Display.OpacityTransferFunction = 'PiecewiseFunction'
    appendAttributes1Display.DataAxesGrid = 'GridAxesRepresentation'
    appendAttributes1Display.PolarAxes = 'PolarAxesRepresentation'
    appendAttributes1Display.ScalarOpacityFunction = eCAPPWF
    appendAttributes1Display.ScalarOpacityUnitDistance = 2.440939567077204
    appendAttributes1Display.OpacityArrayName = ['POINTS', 'ECAP']
    appendAttributes1Display.SelectInputVectors = [None, '']
    appendAttributes1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    appendAttributes1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    appendAttributes1Display.ScaleTransferFunction.Points = [-0.35945209860801697, 0.0, 0.5, 0.0, 159.9031219482422,
                                                             1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    appendAttributes1Display.OpacityTransferFunction.Points = [-0.35945209860801697, 0.0, 0.5, 0.0, 159.9031219482422,
                                                               1.0, 0.5, 0.0]

    # hide data in view
    Hide(oSI_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_05xdmf, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Extract Surface'
    extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=appendAttributes1)

    # show data in view
    extractSurface1Display = Show(extractSurface1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    extractSurface1Display.Representation = 'Surface'
    extractSurface1Display.ColorArrayName = ['POINTS', 'ECAP']
    extractSurface1Display.LookupTable = eCAPLUT
    extractSurface1Display.PointSize = 0.1
    extractSurface1Display.LineWidth = 2.0
    extractSurface1Display.RenderLinesAsTubes = 1
    extractSurface1Display.SelectTCoordArray = 'None'
    extractSurface1Display.SelectNormalArray = 'None'
    extractSurface1Display.SelectTangentArray = 'None'
    extractSurface1Display.OSPRayScaleArray = 'ECAP'
    extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSurface1Display.SelectOrientationVectors = 'None'
    extractSurface1Display.ScaleFactor = 8.22499008178711
    extractSurface1Display.SelectScaleArray = 'ECAP'
    extractSurface1Display.GlyphType = 'Arrow'
    extractSurface1Display.GlyphTableIndexArray = 'ECAP'
    extractSurface1Display.GaussianRadius = 0.41124950408935546
    extractSurface1Display.SetScaleArray = ['POINTS', 'ECAP']
    extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.OpacityArray = ['POINTS', 'ECAP']
    extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'
    extractSurface1Display.SelectInputVectors = [None, '']
    extractSurface1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    extractSurface1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSurface1Display.ScaleTransferFunction.Points = [-0.35945209860801697, 0.0, 0.5, 0.0, 159.9031219482422, 1.0,
                                                           0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSurface1Display.OpacityTransferFunction.Points = [-0.35945209860801697, 0.0, 0.5, 0.0, 159.9031219482422,
                                                             1.0, 0.5, 0.0]

    # hide data in view
    Hide(appendAttributes1, renderView1)

    # show color bar/color legend
    extractSurface1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # save data
    save_path = path.join(solution_path, f"hemodynamics_cycle_{cycle:02d}.vtp")
    SaveData(save_path, proxy=extractSurface1, PointDataArrays=['ECAP', 'OSI', 'RRT', 'TAWSS', 'TWSSG'])

    # ================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    # ================================================================

    # get layout
    layout1 = GetLayout()

    # --------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(1612, 1270)

    # -----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [-161.9212053970614, -9.620525906665149, -15.050279272305126]
    renderView1.CameraFocalPoint = [-7.099048614501934, 106.90110015869136, -184.287498474121]
    renderView1.CameraViewUp = [-0.4019679064885209, 0.8834884918920158, 0.24056160717705413]
    renderView1.CameraParallelScale = 66.58666397548878
    ResetSession()


def ResetSession():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()


if __name__ == '__main__':
    # TODO: UPDATE CASE PATH
    # TODO: UPDATE FILE PATH
    # TODO: UPDATE SAVE PATH
    # cases = sorted(listdir("/Users/henriakj/PhD/Code/VaMPy/models/models_inria"))
    # conditions = ["SR", "AF"]

    conditions = ['SR']
    cases = ['1029']
    cycle = 5
    for condition in conditions:
        for case in cases:
            print(f"Combining and converting HEMODYNAMICS from xdmf to vtu for {case} for condition {condition}")
            try:
                main(case, condition, cycle)
            except Exception as e:
                print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
