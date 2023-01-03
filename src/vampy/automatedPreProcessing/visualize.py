import vtk

from DisplayData import DisplayModel, VtkPointCloud

version = vtk.vtkVersion().GetVTKMajorVersion()


def visualize(network_elements, probe_points, output_surface, mean_inflow_rate):
    """
    Visualize surface model / mesh with distributed flow rates at each inlet/outlet
    in an interactive VTK window.

    Args:
        network_elements (list): List of inlet/outlet elements within a surface model network
        probe_points (ndarray): Array of probe points
        output_surface (vtkPolyData): Surface model to visualize
        mean_inflow_rate (float): Mean inflow rate of input model
    """
    points = vtk.vtkPoints()
    scalar = vtk.vtkDoubleArray()
    scalar.SetNumberOfComponents(1)
    for element in network_elements:
        if not (element.IsAnOutlet()):
            continue
        points.InsertNextPoint(element.GetOutPointsx1()[0])
        scalar.InsertNextValue(100.0 * element.GetBeta())
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(scalar)

    pointInlet = vtk.vtkPoints()
    scalarInlet = vtk.vtkDoubleArray()
    scalarInlet.SetNumberOfComponents(1)
    for element in network_elements:
        if element.IsAnInlet():
            pointInlet.InsertNextPoint(element.GetInPointsx0()[0])
            scalarInlet.InsertNextValue(mean_inflow_rate)
            break

    polydataInlet = vtk.vtkPolyData()
    polydataInlet.SetPoints(pointInlet)
    polydataInlet.GetPointData().SetScalars(scalarInlet)

    labelMapperInlet = vtk.vtkLabeledDataMapper()
    labelMapper = vtk.vtkLabeledDataMapper()
    if version < 6:
        labelMapperInlet.SetInputConnection(polydataInlet.GetProducerPort())
        labelMapper.SetInputConnection(polydata.GetProducerPort())
    else:
        labelMapperInlet.SetInputData(polydataInlet)
        labelMapper.SetInputData(polydata)
    labelMapperInlet.SetLabelModeToLabelScalars()
    labelMapper.SetLabelModeToLabelScalars()

    labelProperties = labelMapperInlet.GetLabelTextProperty()
    labelProperties.SetFontFamilyToArial()
    labelProperties.SetFontSize(30)

    labelProperties = labelMapper.GetLabelTextProperty()
    labelProperties.SetFontFamilyToArial()
    labelProperties.SetFontSize(30)

    labelsInlet = vtk.vtkActor2D()
    labelsInlet.SetMapper(labelMapperInlet)
    labelMapperInlet.SetLabelFormat("%2.2f ml/s")
    labels = vtk.vtkActor2D()
    labels.SetMapper(labelMapper)
    labelMapper.SetLabelFormat("%2.2f %%")

    # Create a cloud of points from the list of probe points.
    pointCloud = VtkPointCloud()
    for i in range(probe_points.shape[0]):
        pointCloud.addPoint(probe_points[i])

    # Create the renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(labels)
    renderer.AddActor(labelsInlet)
    renderer.AddActor(pointCloud.vtkActor)

    opacity = 0.3
    renderer.AddActor(DisplayModel().polyDataToActor(output_surface, opacity))
    renderer.SetBackground(0.1, 0.1, 0.2)

    # Set the lights of the renderer
    DisplayModel().setLight(renderer)

    # Create the RenderWindow and RenderWindowInteractor
    windowTitle = "Mean inlet flow rate and percents of it at the inflows"
    DisplayModel().renderWindow(renderer, windowTitle)
