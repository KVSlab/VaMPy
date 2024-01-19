import argparse
from os import path

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition):
    # Local
    # solution_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics/'

    # On Oracle
    solution_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/Hemodynamics"
    print(f"-- Reading solutions from {solution_path} --")

    # Local
    # solution_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics/'

    # On Oracle
    solution_path = f"/home/opc/Simulation40/{condition}/{case}/results_moving_atrium/data/1/Hemodynamics"
    print(f"-- Reading solutions from {solution_path} --")

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_01xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_01.xdmf', FileName=[
        path.join(solution_path, '/ECAP_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_02xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_02.xdmf', FileName=[
        path.join(solution_path, '/ECAP_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_03xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_03.xdmf', FileName=[
        path.join(solution_path, '/ECAP_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_04xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_04.xdmf', FileName=[
        path.join(solution_path, '/ECAP_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_05xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_05.xdmf', FileName=[
        path.join(solution_path, '/ECAP_cycle_05.xdmf')])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    eCAP_cycle_01xdmfDisplay = Show(eCAP_cycle_01xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_01xdmfDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    eCAP_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_05xdmfDisplay = Show(eCAP_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_04xdmfDisplay = Show(eCAP_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_03xdmfDisplay = Show(eCAP_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_02xdmfDisplay = Show(eCAP_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'ECAP'
    eCAPLUT = GetColorTransferFunction('ECAP')

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # get opacity transfer function/opacity map for 'ECAP'
    eCAPPWF = GetOpacityTransferFunction('ECAP')

    # get 2D transfer function for 'ECAP'
    eCAPTF2D = GetTransferFunction2D('ECAP')

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_01xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_01.xdmf', FileName=[
        path.join(solution_path, '/OSI_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_02xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_02.xdmf', FileName=[
        path.join(solution_path, '/OSI_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_03xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_03.xdmf', FileName=[
        path.join(solution_path, '/OSI_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_04xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_04.xdmf', FileName=[
        path.join(solution_path, '/OSI_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_05xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_05.xdmf', FileName=[
        path.join(solution_path, '/OSI_cycle_05.xdmf')])

    # show data in view
    oSI_cycle_01xdmfDisplay = Show(oSI_cycle_01xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_01xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_03xdmfDisplay = Show(oSI_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_05xdmfDisplay = Show(oSI_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_04xdmfDisplay = Show(oSI_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_02xdmfDisplay = Show(oSI_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # get color transfer function/color map for 'OSI'
    oSILUT = GetColorTransferFunction('OSI')

    # get opacity transfer function/opacity map for 'OSI'
    oSIPWF = GetOpacityTransferFunction('OSI')

    # get 2D transfer function for 'OSI'
    oSITF2D = GetTransferFunction2D('OSI')

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_01xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_01.xdmf', FileName=[
        path.join(solution_path, '/RRT_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_02xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_02.xdmf', FileName=[
        path.join(solution_path, '/RRT_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_03xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_03.xdmf', FileName=[
        path.join(solution_path, '/RRT_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_04xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_04.xdmf', FileName=[
        path.join(solution_path, '/RRT_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_05xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_05.xdmf', FileName=[
        path.join(solution_path, '/RRT_cycle_05.xdmf')])

    # show data in view
    rRT_cycle_01xdmfDisplay = Show(rRT_cycle_01xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_01xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_02xdmfDisplay = Show(rRT_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_04xdmfDisplay = Show(rRT_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_03xdmfDisplay = Show(rRT_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_05xdmfDisplay = Show(rRT_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # get color transfer function/color map for 'RRT'
    rRTLUT = GetColorTransferFunction('RRT')

    # get opacity transfer function/opacity map for 'RRT'
    rRTPWF = GetOpacityTransferFunction('RRT')

    # get 2D transfer function for 'RRT'
    rRTTF2D = GetTransferFunction2D('RRT')

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_01xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_01.xdmf', FileName=[
        path.join(solution_path, '/TAWSS_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_02xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_02.xdmf', FileName=[
        path.join(solution_path, '/TAWSS_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_03xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_03.xdmf', FileName=[
        path.join(solution_path, '/TAWSS_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_04xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_04.xdmf', FileName=[
        path.join(solution_path, '/TAWSS_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_05xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_05.xdmf', FileName=[
        path.join(solution_path, '/TAWSS_cycle_05.xdmf')])

    # show data in view
    tAWSS_cycle_01xdmfDisplay = Show(tAWSS_cycle_01xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_01xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_03xdmfDisplay = Show(tAWSS_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_02xdmfDisplay = Show(tAWSS_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_05xdmfDisplay = Show(tAWSS_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_04xdmfDisplay = Show(tAWSS_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # get color transfer function/color map for 'TAWSS'
    tAWSSLUT = GetColorTransferFunction('TAWSS')

    # get opacity transfer function/opacity map for 'TAWSS'
    tAWSSPWF = GetOpacityTransferFunction('TAWSS')

    # get 2D transfer function for 'TAWSS'
    tAWSSTF2D = GetTransferFunction2D('TAWSS')

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_01xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_01.xdmf', FileName=[
        path.join(solution_path, '/TWSSG_cycle_01.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_02xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_02.xdmf', FileName=[
        path.join(solution_path, '/TWSSG_cycle_02.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_03xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_03.xdmf', FileName=[
        path.join(solution_path, '/TWSSG_cycle_03.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_04xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_04.xdmf', FileName=[
        path.join(solution_path, '/TWSSG_cycle_04.xdmf')])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_05xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_05.xdmf', FileName=[
        path.join(solution_path, '/TWSSG_cycle_05.xdmf')])

    # show data in view
    tWSSG_cycle_01xdmfDisplay = Show(tWSSG_cycle_01xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_01xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_01xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_04xdmfDisplay = Show(tWSSG_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_05xdmfDisplay = Show(tWSSG_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_02xdmfDisplay = Show(tWSSG_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_03xdmfDisplay = Show(tWSSG_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # get color transfer function/color map for 'TWSSG'
    tWSSGLUT = GetColorTransferFunction('TWSSG')

    # get opacity transfer function/opacity map for 'TWSSG'
    tWSSGPWF = GetOpacityTransferFunction('TWSSG')

    # get 2D transfer function for 'TWSSG'
    tWSSGTF2D = GetTransferFunction2D('TWSSG')

    # set active source
    SetActiveSource(eCAP_cycle_01xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # set active source
    SetActiveSource(eCAP_cycle_05xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.6033542156219482, 159.9031219482422)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[eCAP_cycle_01xdmf, eCAP_cycle_02xdmf, eCAP_cycle_03xdmf,
                                                eCAP_cycle_04xdmf, eCAP_cycle_05xdmf])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'

    # hide data in view
    Hide(eCAP_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_04xdmf, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(oSI_cycle_01xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(oSI_cycle_05xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # create a new 'Append Attributes'
    appendAttributes2 = AppendAttributes(registrationName='AppendAttributes2',
                                         Input=[oSI_cycle_01xdmf, oSI_cycle_02xdmf, oSI_cycle_03xdmf, oSI_cycle_04xdmf,
                                                oSI_cycle_05xdmf])

    # show data in view
    appendAttributes2Display = Show(appendAttributes2, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes2Display.Representation = 'Surface'

    # hide data in view
    Hide(oSI_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_03xdmf, renderView1)

    # show color bar/color legend
    appendAttributes2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(rRT_cycle_01xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(rRT_cycle_05xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # create a new 'Append Attributes'
    appendAttributes3 = AppendAttributes(registrationName='AppendAttributes3',
                                         Input=[rRT_cycle_01xdmf, rRT_cycle_02xdmf, rRT_cycle_03xdmf, rRT_cycle_04xdmf,
                                                rRT_cycle_05xdmf])

    # show data in view
    appendAttributes3Display = Show(appendAttributes3, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes3Display.Representation = 'Surface'

    # hide data in view
    Hide(rRT_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_03xdmf, renderView1)

    # show color bar/color legend
    appendAttributes3Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(tAWSS_cycle_01xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(tAWSS_cycle_05xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # create a new 'Append Attributes'
    appendAttributes4 = AppendAttributes(registrationName='AppendAttributes4',
                                         Input=[tAWSS_cycle_01xdmf, tAWSS_cycle_02xdmf, tAWSS_cycle_03xdmf,
                                                tAWSS_cycle_04xdmf, tAWSS_cycle_05xdmf])

    # show data in view
    appendAttributes4Display = Show(appendAttributes4, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes4Display.Representation = 'Surface'

    # hide data in view
    Hide(tAWSS_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_04xdmf, renderView1)

    # show color bar/color legend
    appendAttributes4Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(tWSSG_cycle_01xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(tWSSG_cycle_05xdmf)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # create a new 'Append Attributes'
    appendAttributes5 = AppendAttributes(registrationName='AppendAttributes5',
                                         Input=[tWSSG_cycle_01xdmf, tWSSG_cycle_02xdmf, tWSSG_cycle_03xdmf,
                                                tWSSG_cycle_04xdmf, tWSSG_cycle_05xdmf])

    # show data in view
    appendAttributes5Display = Show(appendAttributes5, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes5Display.Representation = 'Surface'

    # hide data in view
    Hide(tWSSG_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_01xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_02xdmf, renderView1)

    # show color bar/color legend
    appendAttributes5Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(appendAttributes1)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # set active source
    SetActiveSource(appendAttributes5)

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # create a new 'Append Attributes'
    appendAttributes6 = AppendAttributes(registrationName='AppendAttributes6',
                                         Input=[appendAttributes1, appendAttributes2, appendAttributes3,
                                                appendAttributes4, appendAttributes5])

    # show data in view
    appendAttributes6Display = Show(appendAttributes6, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes6Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes3, renderView1)

    # hide data in view
    Hide(appendAttributes5, renderView1)

    # hide data in view
    Hide(appendAttributes2, renderView1)

    # hide data in view
    Hide(appendAttributes1, renderView1)

    # hide data in view
    Hide(appendAttributes4, renderView1)

    # show color bar/color legend
    appendAttributes6Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # create a new 'Extract Surface'
    extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=appendAttributes6)

    # show data in view
    extractSurface1Display = Show(extractSurface1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    extractSurface1Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes6, renderView1)

    # show color bar/color legend
    extractSurface1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    eCAPLUT.RescaleTransferFunction(-0.3172625005245209, 159.9031219482421)

    # save data
    save_path = path.join(solution_path, f"hemodynamics.vtp")
    SaveData(save_path,
             proxy=extractSurface1,
             PointDataArrays=['ECAP', 'ECAP_input_1', 'ECAP_input_2', 'ECAP_input_3', 'ECAP_input_4', 'OSI',
                              'OSI_input_1', 'OSI_input_2', 'OSI_input_3', 'OSI_input_4', 'RRT', 'RRT_input_1',
                              'RRT_input_2', 'RRT_input_3', 'RRT_input_4', 'TAWSS', 'TAWSS_input_1', 'TAWSS_input_2',
                              'TAWSS_input_3', 'TAWSS_input_4', 'TWSSG', 'TWSSG_input_1', 'TWSSG_input_2',
                              'TWSSG_input_3', 'TWSSG_input_4'])

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
            print(
                f"Combining and converting HEMODYNAMICS from xdmf to vtu for {case} for condition {condition}")
            try:
                main(case, condition)
            except Exception as e:
                print(f"-- FAILED for case {case}, condition {condition}), error: {e}")
