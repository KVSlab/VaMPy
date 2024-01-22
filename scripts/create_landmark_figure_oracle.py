from os import path

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition):
    # create a new 'XML PolyData Reader'
    root = f'/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/Hemodynamics'
    hemodynamics_lavtp = XMLPolyDataReader(registrationName='hemodynamics_la.vtp', FileName=[
        path.join(root, 'hemodynamics_la.vtp')])

    # create a new 'XML PolyData Reader'
    hemodynamics_laavtp = XMLPolyDataReader(registrationName='hemodynamics_laa.vtp', FileName=[
        path.join(root, 'hemodynamics_laa.vtp')])

    # Properties modified on hemodynamics_laavtp
    hemodynamics_laavtp.TimeArray = 'None'

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    hemodynamics_laavtpDisplay = Show(hemodynamics_laavtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    hemodynamics_laavtpDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    hemodynamics_laavtpDisplay.SetScalarBarVisibility(renderView1, True)

    # Properties modified on hemodynamics_lavtp
    hemodynamics_lavtp.TimeArray = 'None'

    # show data in view
    hemodynamics_lavtpDisplay = Show(hemodynamics_lavtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    hemodynamics_lavtpDisplay.Representation = 'Surface'

    # show color bar/color legend
    hemodynamics_lavtpDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'RegionId'
    regionIdLUT = GetColorTransferFunction('RegionId')

    # get opacity transfer function/opacity map for 'RegionId'
    regionIdPWF = GetOpacityTransferFunction('RegionId')

    # get 2D transfer function for 'RegionId'
    regionIdTF2D = GetTransferFunction2D('RegionId')

    # turn off scalar coloring
    ColorBy(hemodynamics_laavtpDisplay, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(regionIdLUT, renderView1)

    # set active source
    SetActiveSource(hemodynamics_lavtp)

    # turn off scalar coloring
    ColorBy(hemodynamics_lavtpDisplay, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(regionIdLUT, renderView1)

    renderView1.ApplyIsometricView()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # reset view to fit data
    renderView1.ResetCamera(True, 0.9)

    # set active source
    SetActiveSource(hemodynamics_laavtp)

    # change solid color
    hemodynamics_laavtpDisplay.AmbientColor = [0.8392156862745098, 0.15294117647058825, 0.1568627450980392]
    hemodynamics_laavtpDisplay.DiffuseColor = [0.8392156862745098, 0.15294117647058825, 0.1568627450980392]

    # Properties modified on hemodynamics_laavtpDisplay
    hemodynamics_laavtpDisplay.Specular = 0.02

    # Properties modified on hemodynamics_laavtpDisplay
    hemodynamics_laavtpDisplay.Specular = 1.0

    # set active source
    SetActiveSource(hemodynamics_lavtp)

    # Properties modified on hemodynamics_lavtpDisplay
    hemodynamics_lavtpDisplay.Specular = 1.0

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1796, 1684)

    # current camera placement for renderView1
    renderView1.CameraPosition = [131.95300317629005, 247.0727198205771, -44.69630025998924]
    renderView1.CameraFocalPoint = [-3.725383758544922, 111.39433288574219, -180.37468719482422]
    renderView1.CameraViewUp = [-0.40824829046386296, 0.8164965809277263, -0.40824829046386285]
    renderView1.CameraViewAngle = 32.086302454473476
    renderView1.CameraParallelScale = 60.822956914549756

    # save screenshot
    SaveScreenshot(f'/app/VaMPy/scripts/landmark_figures/{condition}/{case}.png', renderView1, 16,
                   ImageResolution=[1796, 1684])

    ResetSession()


def ResetSession():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()


if __name__ == '__main__':
    cases = [f'{case:04d}' for case in
             [3, 4, 5, 6, 7, 8, 9, 19, 20, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 74, 76, 77, 78, 80,
              81, 1029, 1030, 1031, 1032, 1033, 1035, 1037, 1038, 1039, 2022]]
    conditions = ['AF', 'SR']
    for condition in conditions:
        for case in cases:
            print(f"Figure for {case} for condition {condition}")
            try:
                main(case, condition)
            except Exception as e:
                print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
