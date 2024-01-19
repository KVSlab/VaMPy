import argparse
from os import path

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, is_local=False):
    # create a new 'XML Unstructured Grid Reader'
    if is_local:
        # LOCAL
        solution_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/'
    else:
        # ORACLE
        solution_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/"

    energy_cycle_01vtu = XMLUnstructuredGridReader(registrationName='energy_cycle_01.vtu', FileName=[
        path.join(solution_path, 'energy_cycle_01.vtu')])

    # create a new 'XML Unstructured Grid Reader'
    energy_cycle_02vtu = XMLUnstructuredGridReader(registrationName='energy_cycle_02.vtu', FileName=[
        path.join(solution_path, 'energy_cycle_02.vtu')])

    # create a new 'XML Unstructured Grid Reader'
    energy_cycle_03vtu = XMLUnstructuredGridReader(registrationName='energy_cycle_03.vtu', FileName=[
        path.join(solution_path, 'energy_cycle_03.vtu')])

    # create a new 'XML Unstructured Grid Reader'
    energy_cycle_04vtu = XMLUnstructuredGridReader(registrationName='energy_cycle_04.vtu', FileName=[
        path.join(solution_path, 'energy_cycle_04.vtu')])

    # create a new 'XML Unstructured Grid Reader'
    energy_cycle_05vtu = XMLUnstructuredGridReader(registrationName='energy_cycle_05.vtu', FileName=[
        path.join(solution_path, 'energy_cycle_05.vtu')])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    energy_cycle_02vtuDisplay = Show(energy_cycle_02vtu, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    energy_cycle_02vtuDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    energy_cycle_02vtuDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    energy_cycle_04vtuDisplay = Show(energy_cycle_04vtu, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    energy_cycle_04vtuDisplay.Representation = 'Surface'

    # show color bar/color legend
    energy_cycle_04vtuDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    energy_cycle_05vtuDisplay = Show(energy_cycle_05vtu, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    energy_cycle_05vtuDisplay.Representation = 'Surface'

    # show color bar/color legend
    energy_cycle_05vtuDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    energy_cycle_03vtuDisplay = Show(energy_cycle_03vtu, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    energy_cycle_03vtuDisplay.Representation = 'Surface'

    # show color bar/color legend
    energy_cycle_03vtuDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    energy_cycle_01vtuDisplay = Show(energy_cycle_01vtu, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    energy_cycle_01vtuDisplay.Representation = 'Surface'

    # show color bar/color legend
    energy_cycle_01vtuDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'kinetic_energy'
    kinetic_energyLUT = GetColorTransferFunction('kinetic_energy')

    # get opacity transfer function/opacity map for 'kinetic_energy'
    kinetic_energyPWF = GetOpacityTransferFunction('kinetic_energy')

    # get 2D transfer function for 'kinetic_energy'
    kinetic_energyTF2D = GetTransferFunction2D('kinetic_energy')

    # hide data in view
    Hide(energy_cycle_01vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_05vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_04vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_03vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_02vtu, renderView1)

    # set active source
    SetActiveSource(energy_cycle_01vtu)

    # set active source
    SetActiveSource(energy_cycle_05vtu)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[energy_cycle_01vtu, energy_cycle_02vtu, energy_cycle_03vtu,
                                                energy_cycle_04vtu, energy_cycle_05vtu])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # hide data in view
    Hide(energy_cycle_02vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_03vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_01vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_04vtu, renderView1)

    # hide data in view
    Hide(energy_cycle_05vtu, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # save data
    # save data
    if is_local:
        save_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/FlowMetrics/energy.vtu'
    else:
        save_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/FlowMetrics/energy.vtu"

    SaveData(save_path,
             proxy=appendAttributes1,
             PointDataArrays=['kinetic_energy', 'kinetic_energy_input_1', 'kinetic_energy_input_2',
                              'kinetic_energy_input_3', 'kinetic_energy_input_4', 'turbulent_kinetic_energy',
                              'turbulent_kinetic_energy_input_1', 'turbulent_kinetic_energy_input_2',
                              'turbulent_kinetic_energy_input_3', 'turbulent_kinetic_energy_input_4'],
             FieldDataArrays=['TimeValue', 'TimeValue_input_1', 'TimeValue_input_2', 'TimeValue_input_3',
                              'TimeValue_input_4'])

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
