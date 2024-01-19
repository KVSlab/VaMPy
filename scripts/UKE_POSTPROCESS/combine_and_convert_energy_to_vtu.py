#### import the simple module from the paraview
import argparse
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, cycle, is_local=False):
    # create a new 'Xdmf3ReaderS'
    if is_local:
        # LOCAL
        filename_ke = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/kinetic_energy_cycle_{cycle:02d}.xdmf'
        filename_tke = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/turbulent_kinetic_energy_cycle_{cycle:02d}.xdmf'
    else:
        # ORACLE
        filename_ke = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/kinetic_energy_cycle_{cycle:02d}.xdmf"
        filename_tke = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/turbulent_kinetic_energy_cycle_{cycle:02d}.xdmf"

    kinetic_energy_cycle_05xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_05.xdmf', FileName=[filename_ke])
    turbulent_kinetic_energy_cycle_05xdmf = Xdmf3ReaderS(registrationName='turbulent_kinetic_energy_cycle_05.xdmf',
                                                         FileName=[filename_tke])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    turbulent_kinetic_energy_cycle_05xdmfDisplay = Show(turbulent_kinetic_energy_cycle_05xdmf, renderView1,
                                                        'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'turbulent_kinetic_energy'
    turbulent_kinetic_energyLUT = GetColorTransferFunction('turbulent_kinetic_energy')

    # get opacity transfer function/opacity map for 'turbulent_kinetic_energy'
    turbulent_kinetic_energyPWF = GetOpacityTransferFunction('turbulent_kinetic_energy')

    # trace defaults for the display properties.
    turbulent_kinetic_energy_cycle_05xdmfDisplay.Representation = 'Surface'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'turbulent_kinetic_energy']
    turbulent_kinetic_energy_cycle_05xdmfDisplay.LookupTable = turbulent_kinetic_energyLUT
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OSPRayScaleArray = 'turbulent_kinetic_energy'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SelectScaleArray = 'turbulent_kinetic_energy'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.GlyphTableIndexArray = 'turbulent_kinetic_energy'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'turbulent_kinetic_energy']
    turbulent_kinetic_energy_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'turbulent_kinetic_energy']
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.ScalarOpacityFunction = turbulent_kinetic_energyPWF
    turbulent_kinetic_energy_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 0.632171218531951
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'turbulent_kinetic_energy']
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    turbulent_kinetic_energy_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875,
                                                                               1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0,
                                                                                 0.008310319855809212, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    turbulent_kinetic_energy_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0,
                                                                                   0.008310319855809212, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    kinetic_energy_cycle_05xdmfDisplay = Show(kinetic_energy_cycle_05xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'kinetic_energy'
    kinetic_energyLUT = GetColorTransferFunction('kinetic_energy')

    # get opacity transfer function/opacity map for 'kinetic_energy'
    kinetic_energyPWF = GetOpacityTransferFunction('kinetic_energy')

    # trace defaults for the display properties.
    kinetic_energy_cycle_05xdmfDisplay.Representation = 'Surface'
    kinetic_energy_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'kinetic_energy']
    kinetic_energy_cycle_05xdmfDisplay.LookupTable = kinetic_energyLUT
    kinetic_energy_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    kinetic_energy_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    kinetic_energy_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    kinetic_energy_cycle_05xdmfDisplay.OSPRayScaleArray = 'kinetic_energy'
    kinetic_energy_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    kinetic_energy_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    kinetic_energy_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    kinetic_energy_cycle_05xdmfDisplay.SelectScaleArray = 'kinetic_energy'
    kinetic_energy_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    kinetic_energy_cycle_05xdmfDisplay.GlyphTableIndexArray = 'kinetic_energy'
    kinetic_energy_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    kinetic_energy_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'kinetic_energy']
    kinetic_energy_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    kinetic_energy_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'kinetic_energy']
    kinetic_energy_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    kinetic_energy_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    kinetic_energy_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    kinetic_energy_cycle_05xdmfDisplay.ScalarOpacityFunction = kinetic_energyPWF
    kinetic_energy_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 0.632171218531951
    kinetic_energy_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'kinetic_energy']
    kinetic_energy_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    kinetic_energy_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    kinetic_energy_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5,
                                                                     0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    kinetic_energy_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.06549811363220215, 1.0,
                                                                       0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    kinetic_energy_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.06549811363220215, 1.0,
                                                                         0.5, 0.0]

    # show color bar/color legend
    kinetic_energy_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get 2D transfer function for 'turbulent_kinetic_energy'
    turbulent_kinetic_energyTF2D = GetTransferFunction2D('turbulent_kinetic_energy')

    # set active source
    SetActiveSource(kinetic_energy_cycle_05xdmf)

    # get 2D transfer function for 'kinetic_energy'
    kinetic_energyTF2D = GetTransferFunction2D('kinetic_energy')

    # set active source
    SetActiveSource(turbulent_kinetic_energy_cycle_05xdmf)

    # set active source
    SetActiveSource(kinetic_energy_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[kinetic_energy_cycle_05xdmf, turbulent_kinetic_energy_cycle_05xdmf])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'
    appendAttributes1Display.ColorArrayName = ['POINTS', 'kinetic_energy']
    appendAttributes1Display.LookupTable = kinetic_energyLUT
    appendAttributes1Display.SelectTCoordArray = 'None'
    appendAttributes1Display.SelectNormalArray = 'None'
    appendAttributes1Display.SelectTangentArray = 'None'
    appendAttributes1Display.OSPRayScaleArray = 'kinetic_energy'
    appendAttributes1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    appendAttributes1Display.SelectOrientationVectors = 'None'
    appendAttributes1Display.ScaleFactor = 8.22499008178711
    appendAttributes1Display.SelectScaleArray = 'kinetic_energy'
    appendAttributes1Display.GlyphType = 'Arrow'
    appendAttributes1Display.GlyphTableIndexArray = 'kinetic_energy'
    appendAttributes1Display.GaussianRadius = 0.41124950408935546
    appendAttributes1Display.SetScaleArray = ['POINTS', 'kinetic_energy']
    appendAttributes1Display.ScaleTransferFunction = 'PiecewiseFunction'
    appendAttributes1Display.OpacityArray = ['POINTS', 'kinetic_energy']
    appendAttributes1Display.OpacityTransferFunction = 'PiecewiseFunction'
    appendAttributes1Display.DataAxesGrid = 'GridAxesRepresentation'
    appendAttributes1Display.PolarAxes = 'PolarAxesRepresentation'
    appendAttributes1Display.ScalarOpacityFunction = kinetic_energyPWF
    appendAttributes1Display.ScalarOpacityUnitDistance = 0.632171218531951
    appendAttributes1Display.OpacityArrayName = ['POINTS', 'kinetic_energy']
    appendAttributes1Display.SelectInputVectors = [None, '']
    appendAttributes1Display.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    appendAttributes1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    appendAttributes1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.06549811363220215, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    appendAttributes1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.06549811363220215, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_05xdmf, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # save data
    if is_local:
        save_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/energy_cycle_{cycle:02d}.vtu'
    else:
        save_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/energy_cycle_{cycle:02d}.vtu"

    SaveData(save_path,
             proxy=appendAttributes1, PointDataArrays=['kinetic_energy', 'turbulent_kinetic_energy'])

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

    parser = argparse.ArgumentParser()
    parser.add_argument('--case', help='Description for foo argument', required=True)
    parser.add_argument('--condition', help='Description for bar argument', required=True)
    args = parser.parse_args()

    conditions = [args.condition]
    cases = [args.case]
    cycle = [1,2,3,4,5]
    for condition in conditions:
        for case in cases:
            print(f"Combining and converting KE & TKE from xdmf to vtu for {case} for condition {condition}")
            try:
                main(case, condition, cycle)
            except Exception as e:
                print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
