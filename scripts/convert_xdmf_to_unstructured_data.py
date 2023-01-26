from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(metric, cycle, conv, case, step=None):
    # find view
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
    if step is not None:
        steps = "_{:03d}".format(step)
    else:
        steps = ""
    # set active view
    SetActiveView(renderView1)
    if conv == "Mesh" or conv == "Time":
        lapath = "{}/".format(case)
    else:
        lapath = ""
    # create a new 'Xdmf3ReaderS'
    kinetic_energy_cycle_01xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_01.xdmf', FileName=[
        '/Users/henriakj/PhD/Code/OasisMove/results_moving_atrium/RigidVsMoving/{}/FlowMetrics/{}_cycle_{:02d}{}.xdmf'.format(case,
            #'/Users/henriakj/PhD/Code/OasisMove/results_moving_atrium/{}Convergence/FlowMetrics/{}{}_cycle_{:02d}{}.xdmf'.format(
                 metric, cycle, steps)])

    # show data in view
    kinetic_energy_cycle_01xdmfDisplay = Show(kinetic_energy_cycle_01xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # get color transfer function/color map for metric
    kinetic_energyLUT = GetColorTransferFunction(metric)

    # get opacity transfer function/opacity map for metric
    kinetic_energyPWF = GetOpacityTransferFunction(metric)

    # trace defaults for the display properties.
    kinetic_energy_cycle_01xdmfDisplay.Representation = 'Surface'
    kinetic_energy_cycle_01xdmfDisplay.ColorArrayName = ['POINTS', metric]
    kinetic_energy_cycle_01xdmfDisplay.LookupTable = kinetic_energyLUT
    kinetic_energy_cycle_01xdmfDisplay.SelectTCoordArray = 'None'
    kinetic_energy_cycle_01xdmfDisplay.SelectNormalArray = 'None'
    kinetic_energy_cycle_01xdmfDisplay.SelectTangentArray = 'None'
    kinetic_energy_cycle_01xdmfDisplay.OSPRayScaleArray = metric
    kinetic_energy_cycle_01xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    kinetic_energy_cycle_01xdmfDisplay.SelectOrientationVectors = 'None'
    kinetic_energy_cycle_01xdmfDisplay.ScaleFactor = 9.4589599609375
    kinetic_energy_cycle_01xdmfDisplay.SelectScaleArray = metric
    kinetic_energy_cycle_01xdmfDisplay.GlyphType = 'Arrow'
    kinetic_energy_cycle_01xdmfDisplay.GlyphTableIndexArray = metric
    kinetic_energy_cycle_01xdmfDisplay.GaussianRadius = 0.472947998046875
    kinetic_energy_cycle_01xdmfDisplay.SetScaleArray = ['POINTS', metric]
    kinetic_energy_cycle_01xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    kinetic_energy_cycle_01xdmfDisplay.OpacityArray = ['POINTS', metric]
    kinetic_energy_cycle_01xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    kinetic_energy_cycle_01xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    kinetic_energy_cycle_01xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    kinetic_energy_cycle_01xdmfDisplay.ScalarOpacityFunction = kinetic_energyPWF
    kinetic_energy_cycle_01xdmfDisplay.ScalarOpacityUnitDistance = 1.0462431131982175
    kinetic_energy_cycle_01xdmfDisplay.OpacityArrayName = ['POINTS', metric]

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    kinetic_energy_cycle_01xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1522868126630783, 1.0, 0.5,
                                                                     0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    kinetic_energy_cycle_01xdmfDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.08832293748855591, 1.0,
                                                                       0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    kinetic_energy_cycle_01xdmfDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.08832293748855591, 1.0,
                                                                         0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    kinetic_energy_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()
    if conv == "Mesh" or conv == "Time":
        laname = "{}_".format(case)
    else:
        laname = ""
    # save data
    SaveData(
        '/Users/henriakj/PhD/Code/OasisMove/results_moving_atrium/RigidVsMoving/{}/FlowMetrics/VTU/{}_cycle_{:02d}{}.vtu'.format(
            #'/Users/henriakj/PhD/Code/OasisMove/results_moving_atrium/{}Convergence/FlowMetrics/VTU/{}{}_cycle_{:02d}{}.vtu'.format(
              #  conv, laname,
            case, metric, cycle, steps),
        proxy=kinetic_energy_cycle_01xdmf, PointDataArrays=[metric])

    ResetSession()


def ResetSession():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()


if __name__ == '__main__':
    conv = "Save"  # Cycle, Time

    cycles = list(range(1, 21))
    cases = ["LA{:03d}".format(i) for i in [53]]  # , 67, 86, 105, 135, 175]]
    cases = ["DT{:03d}".format(i) for i in [10, 20, 50, 100, 200]]  # , 67, 86, 105, 135, 175]]
    metrics = ["kinetic_energy", "dissipation", "turbulent_kinetic_energy", "turbulent_dissipation"]
    cycle = 5
    steps = [4, 5, 7, 10, 16, 25, 50, 100]
    steps = [2]
    step = None
    cases =["Rigid"]
    cases = ["Moving", "Generic"]
    for metric in metrics:
        #for step in steps:
        #    case = "LA"
        for case in cases:
            # for cycle in cycles[1:]:
            print("Converting XDMF to VTU for {} for cycle {} for case {}".format(metric, cycle, case))
            main(metric, cycle, conv, case, step)