# Read in probe ids
import csv
import json
from collections import defaultdict
from os import path

import numpy as np
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()


def main(case, condition, probe,probes):
    # create a new 'XML PolyData Reader'
    hemodynamics_laavtp = XMLPolyDataReader(registrationName='hemodynamics_laa.vtp', FileName=[
        f'/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics/hemodynamics_laa.vtp'])

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

    # update the view to ensure updated data information
    renderView1.Update()

    # get color transfer function/color map for 'RegionId'
    regionIdLUT = GetColorTransferFunction('RegionId')

    # get opacity transfer function/opacity map for 'RegionId'
    regionIdPWF = GetOpacityTransferFunction('RegionId')

    # get 2D transfer function for 'RegionId'
    regionIdTF2D = GetTransferFunction2D('RegionId')

    renderView1.ApplyIsometricView()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToNegativeX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveY()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToNegativeY()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # turn off scalar coloring
    ColorBy(hemodynamics_laavtpDisplay, None)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(regionIdLUT, renderView1)

    # Properties modified on hemodynamics_laavtpDisplay
    hemodynamics_laavtpDisplay.Opacity = 0.2

    # Properties modified on hemodynamics_laavtpDisplay
    hemodynamics_laavtpDisplay.Opacity = 0.4

    renderView1.ResetActiveCameraToPositiveX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToNegativeX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveY()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # create a new 'XML PolyData Reader'
    model_centerlinesvtp = XMLPolyDataReader(registrationName='model_centerlines.vtp', FileName=[
        f'/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/model_centerlines.vtp'])

    # create a new 'XML PolyData Reader'
    model_region_centerline_0vtp = XMLPolyDataReader(registrationName='model_region_centerline_0.vtp', FileName=[
        f'/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/model_region_centerline_0.vtp'])

    # Properties modified on model_region_centerline_0vtp
    model_region_centerline_0vtp.TimeArray = 'None'

    # show data in view
    model_region_centerline_0vtpDisplay = Show(model_region_centerline_0vtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    model_region_centerline_0vtpDisplay.Representation = 'Surface'

    # Properties modified on model_centerlinesvtp
    model_centerlinesvtp.TimeArray = 'None'

    # show data in view
    model_centerlinesvtpDisplay = Show(model_centerlinesvtp, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    model_centerlinesvtpDisplay.Representation = 'Surface'

    # update the view to ensure updated data information
    renderView1.Update()

    # reset view to fit data bounds
    renderView1.ResetCamera(13.773507118225098, 25.916683197021484, 94.66801452636719, 116.58678436279297,
                            -193.3843994140625, -136.79444885253906, False, 0.9)

    # reset view to fit data bounds
    renderView1.ResetCamera(13.773507118225098, 25.916683197021484, 94.66801452636719, 116.58678436279297,
                            -193.3843994140625, -136.79444885253906, False, 0.9)

    # set active source
    SetActiveSource(model_centerlinesvtp)

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ApplyIsometricView()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToNegativeX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveY()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToNegativeY()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveZ()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToNegativeZ()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    renderView1.ResetActiveCameraToPositiveX()

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # Properties modified on model_centerlinesvtpDisplay
    model_centerlinesvtpDisplay.LineWidth = 10.0

    # set active source
    SetActiveSource(model_region_centerline_0vtp)

    # Properties modified on model_region_centerline_0vtpDisplay
    model_region_centerline_0vtpDisplay.LineWidth = 10.0

    # change solid color
    model_region_centerline_0vtpDisplay.AmbientColor = [0.6666666666666666, 0.0, 0.0]
    model_region_centerline_0vtpDisplay.DiffuseColor = [0.6666666666666666, 0.0, 0.0]

    # create a new 'Sphere'
    for i,probe_tmp in emumerate(probes):
        sphere_tmp = Sphere(registrationName=f'Sphere{i+2}')
        sphere1Display_tmp = Show(sphere_tmp, renderView1, 'GeometryRepresentation')
        sphere1Display_tmp.Representation = 'Surface'
        SetActiveSource(sphere_tmp)
        sphere_tmp.Center = probe_tmp
        sphere_tmp.Radius = 1.0
        sphere1Display_tmp.AmbientColor = [0.0, 0., 0.0]
        sphere1Display_tmp.DiffuseColor = [0.0, 0., 0.0]
        renderView1.Update()

    sphere1 = Sphere(registrationName='Sphere1')
    sphere1Display = Show(sphere1, renderView1, 'GeometryRepresentation')
    sphere1Display.Representation = 'Surface'
    SetActiveSource(sphere1)
    sphere1.Center = probe
    sphere1.Radius = 2.0
    sphere1Display.AmbientColor = [0.0, 0.6666666666666666, 0.0]
    sphere1Display.DiffuseColor = [0.0, 0.6666666666666666, 0.0]

    # Properties modified on sphere1
    sphere1.ThetaResolution = 25

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(model_centerlinesvtp)

    # reset view to fit data bounds
    renderView1.ResetCamera(-50.78473663330078, 33.083797454833984, 80.89457702636719, 137.28143310546875,
                            -211.2445526123047, -173.57327270507812, False, 0.9)

    # reset view to fit data
    renderView1.ResetCamera(False, 0.9)

    # create a new 'Text'
    text1 = Text(registrationName='Text1')

    # show data in view
    text1Display = Show(text1, renderView1, 'TextSourceRepresentation')

    # update the view to ensure updated data information
    renderView1.Update()

    # Properties modified on text1Display
    text1Display.WindowLocation = 'Upper Center'

    # Properties modified on text1
    # Properties modified on text1
    text1.Text = f"""


    CASE {case} """

    text1Display.FontSize = 200

    # update the view to ensure updated data information
    renderView1.Update()

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1592, 1684)

    # current camera placement for renderView1
    renderView1.CameraPosition = [-207.78465279022632, 21.682427487704473, -43.37114608783702]
    renderView1.CameraFocalPoint = [-3.7945690155029297, 109.08800506591797, -173.43577575683594]
    renderView1.CameraViewUp = [0.2645033170921009, 0.5556310404803709, 0.7882335580917473]
    renderView1.CameraParallelScale = 66.57663888089196

    # save screenshot
    save_path = path.join("/app/VaMPy/scripts/ProbeViz", condition.upper(), f"probe_{case}.png")
    SaveScreenshot(save_path, renderView1, 16, ImageResolution=[1592, 1684])

    ResetSession()


def ResetSession():
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()


if __name__ == '__main__':

    cases = ["0003", "0004", "0005", "0006", "0007", "0008", "0009", "0019",
             "0020", "0021", "0023", "0024", "0025", "0026", "0027", "0028",
             "0029", "0030", "0031", "0032", "0033", "0034", "0035", "0074",
             "0076", "0077", "0078", "0080", "0081", "1029", "1030", "1031",
             "1032", "1033", "1035", "1037", "1038", "1039", "2022"]

    conditions = ['af', 'sr']
    local = False
    for condition in conditions:
        if local:
            probe_id_path = f"/Users/henriakj/PhD/Code/OasisMove/results_34case/Morphology/probe_ids_{condition}.csv"
        else:
            probe_id_path = f"/home/opc/probe_ids_{condition}.csv"

        probe_df = defaultdict(list)  # each value in each column is appended to a list

        with open(probe_id_path) as f:
            reader = csv.DictReader(f)  # read rows into a dictionary format
            for row in reader:  # read a row as {column1: value1, column2: value2,...}
                for (k, v) in row.items():  # go over each column name and value
                    probe_df[k].append(v)  # append the value into the appropriate list
                    # based on column name k
        ids = probe_df['laa_probe_id']
        # With pandas
        # probe_df = pd.read_csv(probe_id_path)
        # ids = probe_df.laa_probe_id

        for id, case in zip(ids, cases):
            # Read in probes (points)
            print(f"-- Viz probe for case {case} condition {condition.upper()}")
            if local:
                probe_point_path = f"/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_{condition}/{case}/model_probe_point.json"
            else:
                probe_point_path = f"/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/model_probe_point.json"

            with open(probe_point_path, 'r') as infile:
                probe_points = np.array(json.load(infile))
            try:
                probe = probe_points[int(id)]

                main(case, condition, probe, probe_points)
            except:
                print(f"Failed for case {case} condtion {condition}")
            # Match probe id with point
