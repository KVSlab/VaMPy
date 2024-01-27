# trace generated using paraview version 5.12.0-RC2
# import paraview
# paraview.compatibility.major = 5
# paraview.compatibility.minor = 12

#### import the simple module from the paraview
import argparse

from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# create a new 'Xdmf3 Reader S'
def main(case):
    eCAP_cycle_02xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_02.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/ECAP_cycle_02.xdmf'])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_03xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_03.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/ECAP_cycle_03.xdmf'])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_04xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_04.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/ECAP_cycle_04.xdmf'])

    # create a new 'Xdmf3 Reader S'
    eCAP_cycle_05xdmf = Xdmf3ReaderS(registrationName='ECAP_cycle_05.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/ECAP_cycle_05.xdmf'])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_02xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_02.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/OSI_cycle_02.xdmf'])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_03xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_03.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/OSI_cycle_03.xdmf'])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_04xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_04.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/OSI_cycle_04.xdmf'])

    # create a new 'Xdmf3 Reader S'
    oSI_cycle_05xdmf = Xdmf3ReaderS(registrationName='OSI_cycle_05.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/OSI_cycle_05.xdmf'])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_02xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_02.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/RRT_cycle_02.xdmf'])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_03xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_03.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/RRT_cycle_03.xdmf'])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_04xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_04.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/RRT_cycle_04.xdmf'])

    # create a new 'Xdmf3 Reader S'
    rRT_cycle_05xdmf = Xdmf3ReaderS(registrationName='RRT_cycle_05.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/RRT_cycle_05.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_02xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_02.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TAWSS_cycle_02.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_03xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_03.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TAWSS_cycle_03.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_04xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_04.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TAWSS_cycle_04.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tAWSS_cycle_05xdmf = Xdmf3ReaderS(registrationName='TAWSS_cycle_05.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TAWSS_cycle_05.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_02xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_02.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TWSSG_cycle_02.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_03xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_03.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TWSSG_cycle_03.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_04xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_04.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TWSSG_cycle_04.xdmf'])

    # create a new 'Xdmf3 Reader S'
    tWSSG_cycle_05xdmf = Xdmf3ReaderS(registrationName='TWSSG_cycle_05.xdmf', FileName=[
        f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/TWSSG_cycle_05.xdmf'])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    eCAP_cycle_03xdmfDisplay = Show(eCAP_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_03xdmfDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # show color bar/color legend
    eCAP_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_03xdmfDisplay = Show(oSI_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_02xdmfDisplay = Show(rRT_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_02xdmfDisplay = Show(oSI_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_05xdmfDisplay = Show(oSI_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_02xdmfDisplay = Show(eCAP_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_05xdmfDisplay = Show(tAWSS_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_03xdmfDisplay = Show(tAWSS_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tAWSS_cycle_04xdmfDisplay = Show(tAWSS_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_04xdmfDisplay = Show(tWSSG_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    tWSSG_cycle_03xdmfDisplay = Show(tWSSG_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tWSSG_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tWSSG_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_03xdmfDisplay = Show(rRT_cycle_03xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_03xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_03xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    oSI_cycle_04xdmfDisplay = Show(oSI_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    oSI_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    oSI_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_05xdmfDisplay = Show(eCAP_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_04xdmfDisplay = Show(rRT_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    rRT_cycle_05xdmfDisplay = Show(rRT_cycle_05xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    rRT_cycle_05xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    rRT_cycle_05xdmfDisplay.SetScalarBarVisibility(renderView1, True)

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
    tAWSS_cycle_02xdmfDisplay = Show(tAWSS_cycle_02xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    tAWSS_cycle_02xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    tAWSS_cycle_02xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # show data in view
    eCAP_cycle_04xdmfDisplay = Show(eCAP_cycle_04xdmf, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    eCAP_cycle_04xdmfDisplay.Representation = 'Surface'

    # show color bar/color legend
    eCAP_cycle_04xdmfDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'TWSSG'
    tWSSGLUT = GetColorTransferFunction('TWSSG')

    # get opacity transfer function/opacity map for 'TWSSG'
    tWSSGPWF = GetOpacityTransferFunction('TWSSG')

    # get 2D transfer function for 'TWSSG'
    tWSSGTF2D = GetTransferFunction2D('TWSSG')

    # set active source
    SetActiveSource(eCAP_cycle_02xdmf)

    # get color transfer function/color map for 'ECAP'
    eCAPLUT = GetColorTransferFunction('ECAP')

    # get opacity transfer function/opacity map for 'ECAP'
    eCAPPWF = GetOpacityTransferFunction('ECAP')

    # get 2D transfer function for 'ECAP'
    eCAPTF2D = GetTransferFunction2D('ECAP')

    # set active source
    SetActiveSource(eCAP_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1',
                                         Input=[eCAP_cycle_02xdmf, eCAP_cycle_03xdmf, eCAP_cycle_04xdmf,
                                                eCAP_cycle_05xdmf])

    # set active source
    SetActiveSource(oSI_cycle_02xdmf)

    # get color transfer function/color map for 'OSI'
    oSILUT = GetColorTransferFunction('OSI')

    # get opacity transfer function/opacity map for 'OSI'
    oSIPWF = GetOpacityTransferFunction('OSI')

    # get 2D transfer function for 'OSI'
    oSITF2D = GetTransferFunction2D('OSI')

    # set active source
    SetActiveSource(oSI_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes2 = AppendAttributes(registrationName='AppendAttributes2',
                                         Input=[oSI_cycle_02xdmf, oSI_cycle_03xdmf, oSI_cycle_04xdmf, oSI_cycle_05xdmf])

    # set active source
    SetActiveSource(rRT_cycle_02xdmf)

    # get color transfer function/color map for 'RRT'
    rRTLUT = GetColorTransferFunction('RRT')

    # get opacity transfer function/opacity map for 'RRT'
    rRTPWF = GetOpacityTransferFunction('RRT')

    # get 2D transfer function for 'RRT'
    rRTTF2D = GetTransferFunction2D('RRT')

    # set active source
    SetActiveSource(rRT_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes3 = AppendAttributes(registrationName='AppendAttributes3',
                                         Input=[rRT_cycle_02xdmf, rRT_cycle_03xdmf, rRT_cycle_04xdmf, rRT_cycle_05xdmf])

    # set active source
    SetActiveSource(tAWSS_cycle_02xdmf)

    # get color transfer function/color map for 'TAWSS'
    tAWSSLUT = GetColorTransferFunction('TAWSS')

    # get opacity transfer function/opacity map for 'TAWSS'
    tAWSSPWF = GetOpacityTransferFunction('TAWSS')

    # get 2D transfer function for 'TAWSS'
    tAWSSTF2D = GetTransferFunction2D('TAWSS')

    # set active source
    SetActiveSource(tAWSS_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes4 = AppendAttributes(registrationName='AppendAttributes4',
                                         Input=[tAWSS_cycle_02xdmf, tAWSS_cycle_03xdmf, tAWSS_cycle_04xdmf,
                                                tAWSS_cycle_05xdmf])

    # set active source
    SetActiveSource(tWSSG_cycle_02xdmf)

    # set active source
    SetActiveSource(tWSSG_cycle_05xdmf)

    # create a new 'Append Attributes'
    appendAttributes5 = AppendAttributes(registrationName='AppendAttributes5',
                                         Input=[tWSSG_cycle_02xdmf, tWSSG_cycle_03xdmf, tWSSG_cycle_04xdmf,
                                                tWSSG_cycle_05xdmf])

    # show data in view
    appendAttributes1Display = Show(appendAttributes1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes1Display.Representation = 'Surface'

    # hide data in view
    Hide(eCAP_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(eCAP_cycle_02xdmf, renderView1)

    # show color bar/color legend
    appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

    # show data in view
    appendAttributes3Display = Show(appendAttributes3, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes3Display.Representation = 'Surface'

    # hide data in view
    Hide(rRT_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(rRT_cycle_04xdmf, renderView1)

    # show color bar/color legend
    appendAttributes3Display.SetScalarBarVisibility(renderView1, True)

    # show data in view
    appendAttributes2Display = Show(appendAttributes2, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes2Display.Representation = 'Surface'

    # hide data in view
    Hide(oSI_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(oSI_cycle_03xdmf, renderView1)

    # show color bar/color legend
    appendAttributes2Display.SetScalarBarVisibility(renderView1, True)

    # show data in view
    appendAttributes4Display = Show(appendAttributes4, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes4Display.Representation = 'Surface'

    # hide data in view
    Hide(tAWSS_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_05xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(tAWSS_cycle_03xdmf, renderView1)

    # show color bar/color legend
    appendAttributes4Display.SetScalarBarVisibility(renderView1, True)

    # show data in view
    appendAttributes5Display = Show(appendAttributes5, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes5Display.Representation = 'Surface'

    # hide data in view
    Hide(tWSSG_cycle_03xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_04xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_02xdmf, renderView1)

    # hide data in view
    Hide(tWSSG_cycle_05xdmf, renderView1)

    # show color bar/color legend
    appendAttributes5Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(appendAttributes1)

    # create a new 'Calculator'
    calculator1 = Calculator(registrationName='Calculator1', Input=appendAttributes1)

    # Properties modified on calculator1
    calculator1.ResultArrayName = 'ecap'
    calculator1.Function = '0.25*(ECAP+ECAP_input_1+ECAP_input_2+ECAP_input_3)'

    # show data in view
    calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes1, renderView1)

    # show color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'ecap'
    ecapLUT = GetColorTransferFunction('ecap')

    # Rescale transfer function
    ecapLUT.RescaleTransferFunction(-0.17174672428518534, 272209252.0)

    # get opacity transfer function/opacity map for 'ecap'
    ecapPWF = GetOpacityTransferFunction('ecap')

    # get 2D transfer function for 'ecap'
    ecapTF2D = GetTransferFunction2D('ecap')

    # set active source
    SetActiveSource(appendAttributes2)

    # create a new 'Calculator'
    calculator2 = Calculator(registrationName='Calculator2', Input=appendAttributes2)

    # Properties modified on calculator2
    calculator2.ResultArrayName = 'osi'
    calculator2.Function = '0.25*(OSI+OSI_input_1+OSI_input_2+OSI_input_3)'

    # show data in view
    calculator2Display = Show(calculator2, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    calculator2Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes2, renderView1)

    # show color bar/color legend
    calculator2Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    ecapLUT.RescaleTransferFunction(-0.17174672428518534, 272209252.0)

    # get color transfer function/color map for 'osi'
    osiLUT = GetColorTransferFunction('osi')

    # get opacity transfer function/opacity map for 'osi'
    osiPWF = GetOpacityTransferFunction('osi')

    # get 2D transfer function for 'osi'
    osiTF2D = GetTransferFunction2D('osi')

    # set active source
    SetActiveSource(appendAttributes3)

    # create a new 'Calculator'
    calculator3 = Calculator(registrationName='Calculator3', Input=appendAttributes3)

    # Properties modified on calculator3
    calculator3.ResultArrayName = 'rrt'
    calculator3.Function = '0.25*(RRT+RRT_input_1+RRT_input_2+RRT_input_3)'

    # show data in view
    calculator3Display = Show(calculator3, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    calculator3Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes3, renderView1)

    # show color bar/color legend
    calculator3Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    ecapLUT.RescaleTransferFunction(-0.17174672428518534, 272209252.0)

    # get color transfer function/color map for 'rrt'
    rrtLUT = GetColorTransferFunction('rrt')

    # get opacity transfer function/opacity map for 'rrt'
    rrtPWF = GetOpacityTransferFunction('rrt')

    # get 2D transfer function for 'rrt'
    rrtTF2D = GetTransferFunction2D('rrt')

    # set active source
    SetActiveSource(appendAttributes4)

    # create a new 'Calculator'
    calculator4 = Calculator(registrationName='Calculator4', Input=appendAttributes4)

    # Properties modified on calculator4
    calculator4.ResultArrayName = 'tawss'
    calculator4.Function = '0.25*(TAWSS+TAWSS_input_1+TAWSS_input_2+TAWSS_input_3)'

    # show data in view
    calculator4Display = Show(calculator4, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    calculator4Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes4, renderView1)

    # show color bar/color legend
    calculator4Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    ecapLUT.RescaleTransferFunction(-0.17174672428518534, 272209252.0)

    # get color transfer function/color map for 'tawss'
    tawssLUT = GetColorTransferFunction('tawss')

    # get opacity transfer function/opacity map for 'tawss'
    tawssPWF = GetOpacityTransferFunction('tawss')

    # get 2D transfer function for 'tawss'
    tawssTF2D = GetTransferFunction2D('tawss')

    # set active source
    SetActiveSource(appendAttributes5)

    # create a new 'Calculator'
    calculator5 = Calculator(registrationName='Calculator5', Input=appendAttributes5)

    # Properties modified on calculator5
    calculator5.ResultArrayName = 'twssg'
    calculator5.Function = '0.25*(TWSSG+TWSSG_input_1+TWSSG_input_2+TWSSG_input_3)'

    # show data in view
    calculator5Display = Show(calculator5, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    calculator5Display.Representation = 'Surface'

    # hide data in view
    Hide(appendAttributes5, renderView1)

    # show color bar/color legend
    calculator5Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    ecapLUT.RescaleTransferFunction(-0.17174672428518534, 272209252.0)

    # get color transfer function/color map for 'twssg'
    twssgLUT = GetColorTransferFunction('twssg')

    # get opacity transfer function/opacity map for 'twssg'
    twssgPWF = GetOpacityTransferFunction('twssg')

    # get 2D transfer function for 'twssg'
    twssgTF2D = GetTransferFunction2D('twssg')

    # set active source
    SetActiveSource(calculator1)

    # set active source
    SetActiveSource(calculator2)

    # set active source
    SetActiveSource(calculator1)

    # set active source
    SetActiveSource(calculator2)

    # set active source
    SetActiveSource(calculator1)

    # set active source
    SetActiveSource(calculator3)

    # set active source
    SetActiveSource(calculator4)

    # set active source
    SetActiveSource(calculator5)

    # create a new 'Append Attributes'
    appendAttributes6 = AppendAttributes(registrationName='AppendAttributes6',
                                         Input=[calculator2, calculator1, calculator3, calculator4, calculator5])

    # show data in view
    appendAttributes6Display = Show(appendAttributes6, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    appendAttributes6Display.Representation = 'Surface'

    # hide data in view
    Hide(calculator3, renderView1)

    # hide data in view
    Hide(calculator2, renderView1)

    # hide data in view
    Hide(calculator1, renderView1)

    # hide data in view
    Hide(calculator4, renderView1)

    # hide data in view
    Hide(calculator5, renderView1)

    # show color bar/color legend
    appendAttributes6Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Extract Surface'
    extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=appendAttributes6)

    # set active source
    SetActiveSource(appendAttributes6)

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

    # set active source
    SetActiveSource(extractSurface1)

    # save data
    save_path = f'/home/opc/Simulation1/{case}/results_moving_atrium/data/1/Hemodynamics/avghemodynamics.vtp'
    SaveData(save_path, proxy=extractSurface1, ChooseArraysToWrite=1,
             PointDataArrays=['ECAP', 'ECAP_input_1', 'ECAP_input_2', 'ECAP_input_3', 'OSI', 'OSI_input_1',
                              'OSI_input_2', 'OSI_input_3', 'RRT', 'RRT_input_1', 'RRT_input_2', 'RRT_input_3', 'TAWSS',
                              'TAWSS_input_1', 'TAWSS_input_2', 'TAWSS_input_3', 'TWSSG', 'TWSSG_input_1',
                              'TWSSG_input_2', 'TWSSG_input_3', 'ecap', 'osi', 'rrt', 'tawss', 'twssg'])

    ResetSession()


def ResetSession():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', help='Description for foo argument', required=True)
    parser.add_argument('--condition', help='Description for bar argument', required=True)
    args = parser.parse_args()

    conditions = [args.condition]
    cases = [args.case]
    for case in cases:
        main(case)
