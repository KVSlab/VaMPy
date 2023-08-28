from paraview.simple import *
from os import getcwd, makedirs, path
from pathlib import Path
import os

#from dolfin import *
from mpi4py import MPI


paraview.simple._DisableFirstRenderCameraReset()


def convert_xdmf_to_vtp(index, path_xdmf):
    # create a new 'Xdmf3ReaderS'
    tWSSGxdmf = Xdmf3ReaderS(registrationName='{}.xdmf'.format(index), FileName=[
        '{}/{}.xdmf'.format(str(path_xdmf), str(index))])        
    tWSSGxdmf.PointArrays = ['{}'.format(index)]

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    tWSSGxdmfDisplay = Show(tWSSGxdmf, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for '{}'.format(index)
    tWSSGLUT = GetColorTransferFunction('{}'.format(index))

    # get opacity transfer function/opacity map for '{}'.format(index)
    tWSSGPWF = GetOpacityTransferFunction('{}'.format(index))

    # trace defaults for the display properties.
    tWSSGxdmfDisplay.Representation = 'Surface'
    tWSSGxdmfDisplay.ColorArrayName = ['POINTS', '{}'.format(index)]
    tWSSGxdmfDisplay.LookupTable = tWSSGLUT
    tWSSGxdmfDisplay.SelectTCoordArray = 'None'
    tWSSGxdmfDisplay.SelectNormalArray = 'None'
    tWSSGxdmfDisplay.SelectTangentArray = 'None'
    tWSSGxdmfDisplay.OSPRayScaleArray = '{}'.format(index)
    tWSSGxdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    tWSSGxdmfDisplay.SelectOrientationVectors = 'None'
    tWSSGxdmfDisplay.ScaleFactor = 10.267200088500978
    tWSSGxdmfDisplay.SelectScaleArray = '{}'.format(index)
    tWSSGxdmfDisplay.GlyphType = 'Arrow'
    tWSSGxdmfDisplay.GlyphTableIndexArray = '{}'.format(index)
    tWSSGxdmfDisplay.GaussianRadius = 0.5133600044250488
    tWSSGxdmfDisplay.SetScaleArray = ['POINTS', '{}'.format(index)]
    tWSSGxdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    tWSSGxdmfDisplay.OpacityArray = ['POINTS', '{}'.format(index)]
    tWSSGxdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    tWSSGxdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    tWSSGxdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    tWSSGxdmfDisplay.ScalarOpacityFunction = tWSSGPWF
    tWSSGxdmfDisplay.ScalarOpacityUnitDistance = 7.34230806856266
    tWSSGxdmfDisplay.OpacityArrayName = ['POINTS', '{}'.format(index)]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tWSSGxdmfDisplay.ScaleTransferFunction.Points = [7.836367876734585e-05, 0.0, 0.5, 0.0, 0.15637780725955963, 1.0,
                                                     0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tWSSGxdmfDisplay.OpacityTransferFunction.Points = [7.836367876734585e-05, 0.0, 0.5, 0.0, 0.15637780725955963, 1.0,
                                                       0.5, 0.0]

    # show color bar/color legend
    tWSSGxdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Extract Surface'
    extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=tWSSGxdmf)

    # show data in view
    extractSurface1Display = Show(extractSurface1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    extractSurface1Display.Representation = 'Surface'
    extractSurface1Display.ColorArrayName = ['POINTS', '{}'.format(index)]
    extractSurface1Display.LookupTable = tWSSGLUT
    extractSurface1Display.SelectTCoordArray = 'None'
    extractSurface1Display.SelectNormalArray = 'None'
    extractSurface1Display.SelectTangentArray = 'None'
    extractSurface1Display.OSPRayScaleArray = '{}'.format(index)
    extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSurface1Display.SelectOrientationVectors = 'None'
    extractSurface1Display.ScaleFactor = 0.2
    extractSurface1Display.SelectScaleArray = '{}'.format(index)
    extractSurface1Display.GlyphType = 'Arrow'
    extractSurface1Display.GlyphTableIndexArray = '{}'.format(index)
    extractSurface1Display.GaussianRadius = 0.01
    extractSurface1Display.SetScaleArray = ['POINTS', '{}'.format(index)]
    extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.OpacityArray = ['POINTS', '{}'.format(index)]
    extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSurface1Display.ScaleTransferFunction.Points = [0.7343099117279053, 0.0, 0.5, 0.0, 42.69994354248047, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSurface1Display.OpacityTransferFunction.Points = [0.7343099117279053, 0.0, 0.5, 0.0, 42.69994354248047, 1.0, 0.5, 0.0]

    #------------------------------------------------------------------------------------------------------    
    # Create folder to store the """.vtp""" files
    common_path_vtp = path.join(path_xdmf, 'vtp_files')
    #if MPI.rank(MPI.comm_world)==0:
    if not path.exists(common_path_vtp):
        makedirs(common_path_vtp)
    originial_vtp_path= common_path_vtp
    path_xdmf = Path(common_path_vtp)
    #------------------------------------------------------------------------------------------------------
    # Saveing data
    common_path_vtp = path.join(originial_vtp_path, '{}'.format(index))
    if not path.exists(common_path_vtp):
        makedirs(common_path_vtp)
    path_xdmf = Path(common_path_vtp)
    SaveData('{}/{}.vtp'.format(str(path_xdmf), str(index)),
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
    layout1.SetSize(1280, 804)

    # -----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [-7.729784319564094, 110.9843864440918, 73.12084965223914]
    renderView1.CameraFocalPoint = [-6.470470428466797, 110.9843864440918, -183.4224090576172]
    renderView1.CameraParallelScale = 66.3990812117987

    # --------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).

    ResetSession()

def ResetSession():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()


if __name__ == '__main__':

    indices = ['{}'.format(index) for index in ['TAWSS', 'OSI', 'RRT', 'ECAP', 'TWSSG']]
    path_xdmf = "/Users/gadursn/Oasis/Newtonian_LA_c167_100K_dt_0dot5AND1ms/C_167/data/3/PostProc/Diff_h5_files_NvsnonN_CY-37"
    for index in indices:
        print("-- Converting: {} index from xdmf to vtp format".format(index))
        convert_xdmf_to_vtp(index, path_xdmf)
