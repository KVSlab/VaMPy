import argparse

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, cycle, is_local=False):
    # create a new 'Xdmf3ReaderS'
    if is_local:
        # LOCAL
        filename = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/blood_residence_time_cycle_{cycle:02d}.xdmf'
    else:
        # ORACLE
        filename = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/blood_residence_time_cycle_{cycle:02d}.xdmf"

    blood_residence_time_cycle_05xdmf = Xdmf3ReaderS(registrationName='blood_residence_time_cycle_05.xdmf',
                                                     FileName=[filename])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    blood_residence_time_cycle_05xdmfDisplay = Show(blood_residence_time_cycle_05xdmf, renderView1,
                                                    'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'blood_residence_time'
    blood_residence_timeLUT = GetColorTransferFunction('blood_residence_time')

    # get opacity transfer function/opacity map for 'blood_residence_time'
    blood_residence_timePWF = GetOpacityTransferFunction('blood_residence_time')

    # trace defaults for the display properties.
    blood_residence_time_cycle_05xdmfDisplay.Representation = 'Surface'
    blood_residence_time_cycle_05xdmfDisplay.ColorArrayName = ['POINTS', 'blood_residence_time']
    blood_residence_time_cycle_05xdmfDisplay.LookupTable = blood_residence_timeLUT
    blood_residence_time_cycle_05xdmfDisplay.SelectTCoordArray = 'None'
    blood_residence_time_cycle_05xdmfDisplay.SelectNormalArray = 'None'
    blood_residence_time_cycle_05xdmfDisplay.SelectTangentArray = 'None'
    blood_residence_time_cycle_05xdmfDisplay.OSPRayScaleArray = 'blood_residence_time'
    blood_residence_time_cycle_05xdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    blood_residence_time_cycle_05xdmfDisplay.SelectOrientationVectors = 'None'
    blood_residence_time_cycle_05xdmfDisplay.ScaleFactor = 8.22499008178711
    blood_residence_time_cycle_05xdmfDisplay.SelectScaleArray = 'blood_residence_time'
    blood_residence_time_cycle_05xdmfDisplay.GlyphType = 'Arrow'
    blood_residence_time_cycle_05xdmfDisplay.GlyphTableIndexArray = 'blood_residence_time'
    blood_residence_time_cycle_05xdmfDisplay.GaussianRadius = 0.41124950408935546
    blood_residence_time_cycle_05xdmfDisplay.SetScaleArray = ['POINTS', 'blood_residence_time']
    blood_residence_time_cycle_05xdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    blood_residence_time_cycle_05xdmfDisplay.OpacityArray = ['POINTS', 'blood_residence_time']
    blood_residence_time_cycle_05xdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    blood_residence_time_cycle_05xdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
    blood_residence_time_cycle_05xdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
    blood_residence_time_cycle_05xdmfDisplay.ScalarOpacityFunction = blood_residence_timePWF
    blood_residence_time_cycle_05xdmfDisplay.ScalarOpacityUnitDistance = 0.632171218531951
    blood_residence_time_cycle_05xdmfDisplay.OpacityArrayName = ['POINTS', 'blood_residence_time']
    blood_residence_time_cycle_05xdmfDisplay.SelectInputVectors = [None, '']
    blood_residence_time_cycle_05xdmfDisplay.WriteLog = ''

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    blood_residence_time_cycle_05xdmfDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1123.3660888671875, 1.0,
                                                                           0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    blood_residence_time_cycle_05xdmfDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 16736728.0, 1.0, 0.5,
                                                                             0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    blood_residence_time_cycle_05xdmfDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 16736728.0, 1.0, 0.5,
                                                                               0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    blood_residence_time_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get 2D transfer function for 'blood_residence_time'
    blood_residence_timeTF2D = GetTransferFunction2D('blood_residence_time')

    # save data
    if is_local:
        save_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/blood_residence_time_cycle_{cycle:02d}.vtu'
    else:
        save_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/blood_residence_time_cycle_{cycle:02d}.vtu"

    SaveData(save_path,
             proxy=blood_residence_time_cycle_05xdmf, PointDataArrays=['blood_residence_time'])

    layout1 = GetLayout()
    layout1.SetSize(2856, 1270)
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
    cycle = [1, 2, 3, 4, 5]

    for condition in conditions:
        for case in cases:
            print(f"Converting BRT xdmf to vtu for {case} for condition {condition}")
            try:
                main(case, condition, cycle)
            except Exception as e:
                print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
