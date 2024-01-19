import argparse
from os import path

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, is_local=False):
    # create a new 'Xdmf3ReaderS'
    if is_local:
        # LOCAL
        solution_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/blood_residence_time_cycle_{cycle:02d}.xdmf'
    else:
        # ORACLE
        solution_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics"
    print(f"-- Reading solutions from {solution_path} --")

    # create a new 'Xdmf3 Reader S'
    blood_residence_time_cycle_01xdmf = Xdmf3ReaderS(registrationName='blood_residence_time_cycle_01.xdmf', FileName=[
        path.join(solution_path, 'blood_residence_time_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    blood_residence_time_cycle_02xdmf = Xdmf3ReaderS(registrationName='blood_residence_time_cycle_02.xdmf', FileName=[
        path.join(solution_path, 'blood_residence_time_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    blood_residence_time_cycle_03xdmf = Xdmf3ReaderS(registrationName='blood_residence_time_cycle_03.xdmf', FileName=[
        path.join(solution_path, 'blood_residence_time_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    blood_residence_time_cycle_04xdmf = Xdmf3ReaderS(registrationName='blood_residence_time_cycle_04.xdmf', FileName=[
        path.join(solution_path, 'blood_residence_time_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    blood_residence_time_cycle_05xdmf = Xdmf3ReaderS(registrationName='blood_residence_time_cycle_05.xdmf', FileName=[
        path.join(solution_path, 'blood_residence_time_cycle_05.xdmf')])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    blood_residence_time_cycle_01xdmfDisplay = Show(blood_residence_time_cycle_01xdmf, renderView1,
                                                    'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    blood_residence_time_cycle_01xdmfDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    blood_residence_time_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    blood_residence_time_cycle_04xdmfDisplay = Show(blood_residence_time_cycle_04xdmf, renderView1,
                                                    'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    blood_residence_time_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    blood_residence_time_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    blood_residence_time_cycle_05xdmfDisplay = Show(blood_residence_time_cycle_05xdmf, renderView1,
                                                    'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    blood_residence_time_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    blood_residence_time_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    blood_residence_time_cycle_02xdmfDisplay = Show(blood_residence_time_cycle_02xdmf, renderView1,
                                                    'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    blood_residence_time_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    blood_residence_time_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    blood_residence_time_cycle_03xdmfDisplay = Show(blood_residence_time_cycle_03xdmf, renderView1,
                                                    'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    blood_residence_time_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    blood_residence_time_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'blood_residence_time'
    blood_residence_timeLUT = GetColorTransferFunction('blood_residence_time')

    # get opacity transfer function/opacity map for 'blood_residence_time'
    blood_residence_timePWF = GetOpacityTransferFunction('blood_residence_time')

    # get 2D transfer function for 'blood_residence_time'
    blood_residence_timeTF2D = GetTransferFunction2D('blood_residence_time')

    # set active source
    SetActiveSource(blood_residence_time_cycle_01xdmf)

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # set active source
    SetActiveSource(blood_residence_time_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[blood_residence_time_cycle_01xdmf, blood_residence_time_cycle_02xdmf,
                                                blood_residence_time_cycle_03xdmf, blood_residence_time_cycle_04xdmf,
                                                blood_residence_time_cycle_05xdmf])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'

    # hide data in view
    Hide(blood_residence_time_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(blood_residence_time_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(blood_residence_time_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(blood_residence_time_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(blood_residence_time_cycle_04xdmf, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    if is_local:
        save_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/blood_residence_time.vtu'
    else:
        save_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/blood_residence_time.vtu"
    # save data
    SaveData(save_path,
             proxy=appendAttributes1,
             PointDataArrays=['blood_residence_time', 'blood_residence_time_input_1', 'blood_residence_time_input_2',
                              'blood_residence_time_input_3', 'blood_residence_time_input_4'])

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

    for condition in conditions:
        for case in cases:
            print(f"Converting BRT xdmf to vtu for {case} for condition {condition}")
            try:
                main(case, condition)
            except Exception as e:
                print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
