#### import the simple module from the paraview
import argparse
from os import path

from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, cycle, is_local=False):
    # create a new 'Xdmf3ReaderS'
    if is_local:
        # LOCAL
        solution_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/'
    else:
        # ORACLE
        solution_path = f'/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/'

        # create a new 'Xdmf3 Reader S'
    kinetic_energy_cycle_01xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_01.xdmf', FileName=[
        path.join(solution_path, 'kinetic_energy_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    kinetic_energy_cycle_02xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_02.xdmf', FileName=[
        path.join(solution_path, 'kinetic_energy_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    kinetic_energy_cycle_03xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_03.xdmf', FileName=[
        path.join(solution_path, 'kinetic_energy_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    kinetic_energy_cycle_04xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_04.xdmf', FileName=[
        path.join(solution_path, 'kinetic_energy_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    kinetic_energy_cycle_05xdmf = Xdmf3ReaderS(registrationName='kinetic_energy_cycle_05.xdmf', FileName=[
        path.join(solution_path, 'kinetic_energy_cycle_05.xdmf')])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    kinetic_energy_cycle_03xdmfDisplay = Show(kinetic_energy_cycle_03xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    kinetic_energy_cycle_03xdmfDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    kinetic_energy_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    kinetic_energy_cycle_05xdmfDisplay = Show(kinetic_energy_cycle_05xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    kinetic_energy_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    kinetic_energy_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    kinetic_energy_cycle_04xdmfDisplay = Show(kinetic_energy_cycle_04xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    kinetic_energy_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    kinetic_energy_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    kinetic_energy_cycle_02xdmfDisplay = Show(kinetic_energy_cycle_02xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    kinetic_energy_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    kinetic_energy_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    kinetic_energy_cycle_01xdmfDisplay = Show(kinetic_energy_cycle_01xdmf, renderView1,
                                              'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    kinetic_energy_cycle_01xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    kinetic_energy_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'kinetic_energy'
    kinetic_energyLUT = GetColorTransferFunction('kinetic_energy')

    # get opacity transfer function/opacity map for 'kinetic_energy'
    kinetic_energyPWF = GetOpacityTransferFunction('kinetic_energy')

    # get 2D transfer function for 'kinetic_energy'
    kinetic_energyTF2D = GetTransferFunction2D('kinetic_energy')

    # hide data in view
    Hide(kinetic_energy_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_04xdmf, renderView1)

    # set active source
    SetActiveSource(kinetic_energy_cycle_01xdmf)

    # set active source
    SetActiveSource(kinetic_energy_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[kinetic_energy_cycle_01xdmf, kinetic_energy_cycle_02xdmf,
                                                kinetic_energy_cycle_03xdmf, kinetic_energy_cycle_04xdmf,
                                                kinetic_energy_cycle_05xdmf])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # hide data in view
    Hide(kinetic_energy_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(kinetic_energy_cycle_05xdmf, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # hide data in view
    Hide(appendAttributes1, renderView1)

    # create a new 'Xdmf3 Reader S'
    turbulent_kinetic_energy_cycle_01xdmf = Xdmf3ReaderS(registrationName='turbulent_kinetic_energy_cycle_01.xdmf',
                                                         FileName=[
                                                             path.join(solution_path,
                                                                       'turbulent_kinetic_energy_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    turbulent_kinetic_energy_cycle_02xdmf = Xdmf3ReaderS(registrationName='turbulent_kinetic_energy_cycle_02.xdmf',
                                                         FileName=[
                                                             path.join(solution_path,
                                                                       'turbulent_kinetic_energy_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    turbulent_kinetic_energy_cycle_03xdmf = Xdmf3ReaderS(registrationName='turbulent_kinetic_energy_cycle_03.xdmf',
                                                         FileName=[
                                                             path.join(solution_path,
                                                                       'turbulent_kinetic_energy_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    turbulent_kinetic_energy_cycle_04xdmf = Xdmf3ReaderS(registrationName='turbulent_kinetic_energy_cycle_04.xdmf',
                                                         FileName=[
                                                             path.join(solution_path,
                                                                       'turbulent_kinetic_energy_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    turbulent_kinetic_energy_cycle_05xdmf = Xdmf3ReaderS(registrationName='turbulent_kinetic_energy_cycle_05.xdmf',
                                                         FileName=[
                                                             path.join(solution_path,
                                                                       'turbulent_kinetic_energy_cycle_05.xdmf')])

    # show data in view
    turbulent_kinetic_energy_cycle_05xdmfDisplay = Show(turbulent_kinetic_energy_cycle_05xdmf, renderView1,
                                                        'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    turbulent_kinetic_energy_cycle_05xdmfDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # show color bar/color legend
    turbulent_kinetic_energy_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    turbulent_kinetic_energy_cycle_04xdmfDisplay = Show(turbulent_kinetic_energy_cycle_04xdmf, renderView1,
                                                        'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    turbulent_kinetic_energy_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    turbulent_kinetic_energy_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    turbulent_kinetic_energy_cycle_03xdmfDisplay = Show(turbulent_kinetic_energy_cycle_03xdmf, renderView1,
                                                        'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    turbulent_kinetic_energy_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    turbulent_kinetic_energy_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    turbulent_kinetic_energy_cycle_01xdmfDisplay = Show(turbulent_kinetic_energy_cycle_01xdmf, renderView1,
                                                        'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    turbulent_kinetic_energy_cycle_01xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    turbulent_kinetic_energy_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    turbulent_kinetic_energy_cycle_02xdmfDisplay = Show(turbulent_kinetic_energy_cycle_02xdmf, renderView1,
                                                        'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    turbulent_kinetic_energy_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    turbulent_kinetic_energy_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(kinetic_energy_cycle_04xdmf)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_05xdmf, renderView1)

    # set active source
    SetActiveSource(turbulent_kinetic_energy_cycle_01xdmf)

    # get color transfer function/color map for 'turbulent_kinetic_energy'
    turbulent_kinetic_energyLUT = GetColorTransferFunction('turbulent_kinetic_energy')

    # get opacity transfer function/opacity map for 'turbulent_kinetic_energy'
    turbulent_kinetic_energyPWF = GetOpacityTransferFunction('turbulent_kinetic_energy')

    # get 2D transfer function for 'turbulent_kinetic_energy'
    turbulent_kinetic_energyTF2D = GetTransferFunction2D('turbulent_kinetic_energy')

    # set active source
    SetActiveSource(turbulent_kinetic_energy_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes2 = AppendAttributes(registrationName='AppendAttributes2',
                                         Input=[turbulent_kinetic_energy_cycle_01xdmf,
                                                turbulent_kinetic_energy_cycle_02xdmf,
                                                turbulent_kinetic_energy_cycle_03xdmf,
                                                turbulent_kinetic_energy_cycle_04xdmf,
                                                turbulent_kinetic_energy_cycle_05xdmf])

    # show data in view
    appendAttributes2Display = Show(appendAttributes2, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes2Display.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(turbulent_kinetic_energy_cycle_02xdmf, renderView1)

    # show color bar/color legend
    appendAttributes2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # hide data in view
    Hide(appendAttributes2, renderView1)

    # set active source
    SetActiveSource(appendAttributes1)

    # set active source
    SetActiveSource(appendAttributes2)

    # set active source
    SetActiveSource(appendAttributes1)

    # set active source
    SetActiveSource(appendAttributes2)

    # set active source
    SetActiveSource(appendAttributes1)

    # set active source
    SetActiveSource(appendAttributes2)

    # create a new 'Append Attributes'
    appendAttributes3 = AppendAttributes(registrationName='AppendAttributes3',
                                         Input=[appendAttributes1, appendAttributes2])

    # show data in view
    appendAttributes3Display = Show(appendAttributes3, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes3Display.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # hide data in view
    Hide(appendAttributes1, renderView1)

    # hide data in view
    Hide(appendAttributes2, renderView1)

    # show color bar/color legend
    appendAttributes3Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # hide data in view
    Hide(appendAttributes3, renderView1)

    # save data
    if is_local:
        save_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/energy.vtu'
    else:
        save_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/energy.vtu"

    # save data
    SaveData(save_path,
             proxy=appendAttributes3,
             PointDataArrays=['kinetic_energy', 'kinetic_energy_input_1', 'kinetic_energy_input_2',
                              'kinetic_energy_input_3', 'kinetic_energy_input_4', 'turbulent_kinetic_energy',
                              'turbulent_kinetic_energy_input_1', 'turbulent_kinetic_energy_input_2',
                              'turbulent_kinetic_energy_input_3', 'turbulent_kinetic_energy_input_4'])

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
    cycles = [1, 2, 3, 4, 5]
    for condition in conditions:
        for case in cases:
            for cycle in cycles:
                print(f"Combining and converting KE & TKE from xdmf to vtu for {case} for condition {condition}")
                try:
                    main(case, condition, cycle)
                except Exception as e:
                    print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
