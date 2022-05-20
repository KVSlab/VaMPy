# trace generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Xdmf3ReaderS'
input_path_nu = '../nu.xdmf'
array_name = 'nu'
nuxdmf = Xdmf3ReaderS(registrationName='nu.xdmf', FileName=[input_path_nu])
nuxdmf.PointArrays = [array_name]

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

### Visualization part
# get active view
"""renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
nuxdmfDisplay = Show(nuxdmf, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'nu'
nuLUT = GetColorTransferFunction(array_name)
nuLUT.RGBPoints = [0.004047868773341179, 0.231373, 0.298039, 0.752941, 0.02986969333142042, 0.865003, 0.865003, 0.865003, 0.055691517889499664, 0.705882, 0.0156863, 0.14902]
nuLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'nu'
nuPWF = GetOpacityTransferFunction(array_name)
nuPWF.Points = [0.004047868773341179, 0.0, 0.5, 0.0, 0.055691517889499664, 1.0, 0.5, 0.0]
nuPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
nuxdmfDisplay.Representation = 'Surface'
nuxdmfDisplay.ColorArrayName = ['POINTS', array_name]
nuxdmfDisplay.LookupTable = nuLUT
nuxdmfDisplay.SelectTCoordArray = 'None'
nuxdmfDisplay.SelectNormalArray = 'None'
nuxdmfDisplay.SelectTangentArray = 'None'
nuxdmfDisplay.OSPRayScaleArray = array_name
nuxdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
nuxdmfDisplay.SelectOrientationVectors = 'None'
nuxdmfDisplay.ScaleFactor = 9.86716995239258
nuxdmfDisplay.SelectScaleArray = array_name
nuxdmfDisplay.GlyphType = 'Arrow'
nuxdmfDisplay.GlyphTableIndexArray = array_name
nuxdmfDisplay.GaussianRadius = 0.4933584976196289
nuxdmfDisplay.SetScaleArray = ['POINTS', array_name]
nuxdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
nuxdmfDisplay.OpacityArray = ['POINTS', array_name]
nuxdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
nuxdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
nuxdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
nuxdmfDisplay.ScalarOpacityFunction = nuPWF
nuxdmfDisplay.ScalarOpacityUnitDistance = 0.793608185079851
nuxdmfDisplay.OpacityArrayName = ['POINTS', array_name]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
nuxdmfDisplay.ScaleTransferFunction.Points = [0.004047868773341179, 0.0, 0.5, 0.0, 0.055691517889499664, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
nuxdmfDisplay.OpacityTransferFunction.Points = [0.004047868773341179, 0.0, 0.5, 0.0, 0.055691517889499664, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
nuxdmfDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color legend/bar for nuLUT in view renderView1
nuLUTColorBar = GetScalarBar(nuLUT, renderView1)
nuLUTColorBar.Title = array_name
nuLUTColorBar.ComponentTitle = ''

# change scalar bar placement
nuLUTColorBar.WindowLocation = 'AnyLocation'
nuLUTColorBar.Position = [0.5923076923076923, 0.27100840336134446]
nuLUTColorBar.ScalarBarLength = 0.32999999999999996"""

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=nuxdmf)
# Properties modified on calculator1
calculator1.ResultArrayName = 'nu_relative'
calculator1.Function = '{}/0.004'.format(str(array_name))

# show data in view
"""calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'nu_relative'
nu_relativeLUT = GetColorTransferFunction('nu_relative')
nu_relativeLUT.RGBPoints = [0.10119671933352947, 0.231373, 0.298039, 0.752941, 0.7467423332855105, 0.865003, 0.865003, 0.865003, 1.3922879472374916, 0.705882, 0.0156863, 0.14902]
nu_relativeLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'nu_relative'
nu_relativePWF = GetOpacityTransferFunction('nu_relative')
nu_relativePWF.Points = [0.10119671933352947, 0.0, 0.5, 0.0, 1.3922879472374916, 1.0, 0.5, 0.0]
nu_relativePWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'nu_relative']
calculator1Display.LookupTable = nu_relativeLUT
calculator1Display.SelectTCoordArray = 'None'
calculator1Display.SelectNormalArray = 'None'
calculator1Display.SelectTangentArray = 'None'
calculator1Display.OSPRayScaleArray = 'nu_relative'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 9.86716995239258
calculator1Display.SelectScaleArray = 'nu_relative'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'nu_relative'
calculator1Display.GaussianRadius = 0.4933584976196289
calculator1Display.SetScaleArray = ['POINTS', 'nu_relative']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'nu_relative']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = nu_relativePWF
calculator1Display.ScalarOpacityUnitDistance = 0.793608185079851
calculator1Display.OpacityArrayName = ['POINTS', 'nu_relative']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.10119671933352947, 0.0, 0.5, 0.0, 1.3922879472374916, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.10119671933352947, 0.0, 0.5, 0.0, 1.3922879472374916, 1.0, 0.5, 0.0]

# hide data in view
Hide(nuxdmf, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color legend/bar for nu_relativeLUT in view renderView1
nu_relativeLUTColorBar = GetScalarBar(nu_relativeLUT, renderView1)
nu_relativeLUTColorBar.Title = 'nu_relative'
nu_relativeLUTColorBar.ComponentTitle = ''

# change scalar bar placement
nu_relativeLUTColorBar.WindowLocation = 'AnyLocation'
nu_relativeLUTColorBar.Position = [0.6092307692307692, 0.2331932773109244]
nu_relativeLUTColorBar.ScalarBarLength = 0.32999999999999985"""

# save data
output_path_nu = './nu_relative.vtu'
SaveData(output_path_nu, proxy=calculator1, ChooseArraysToWrite=1,
    PointDataArrays=['nu_relative'],
    Writetimestepsasfileseries=1)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
"""layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1300, 476)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [-204.85947527687696, -66.36741796824411, 8.92781098754074]
renderView1.CameraFocalPoint = [51.29368105908454, 49.52752694382273, -40.06994759012861]
renderView1.CameraViewUp = [-0.1994959404234519, 0.02680570093873169, -0.9795319413637035]
renderView1.CameraParallelScale = 73.86409793442597"""

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).