import numpy as np
import vtk
from morphman import read_polydata
from vtk.numpy_interface import dataset_adapter as dsa


def convert_to_cell_data(model):
    pointToCellFilter = vtk.vtkPointDataToCellData()
    pointToCellFilter.SetInputData(model)
    pointToCellFilter.PassPointDataOn()
    pointToCellFilter.Update()
    outputData = pointToCellFilter.GetOutput()
    return outputData


def get_volume(model, volume=True, area=False):
    cell_size_filter = vtk.vtkCellSizeFilter()
    cell_size_filter.SetInputData(model)
    cell_size_filter.SetComputeVolume(volume)
    cell_size_filter.SetComputeArea(area)
    cell_size_filter.SetComputeLength(False)
    cell_size_filter.SetComputeVertexCount(False)
    cell_size_filter.Update()
    return cell_size_filter.GetOutput()


def filter_metric(metric, area, m):
    # M = name
    if m == 'TAWSS':
        ids = np.where(metric >= 0)[0]
    elif m == "TWSSG":
        ids = np.where(metric >= 0)[0]
    elif m == "OSI":
        ids = np.where((0 <= metric) & (metric <= 0.5))[0]
    elif m == "RRT":
        ids = np.where((0 <= metric) & (metric <= 1E6))[0]
    elif m == "ECAP":
        ids = np.where((0 <= metric) & (metric <= 1E5))[0]
    else:
        pass

    return metric[ids], area[ids]


def main_surface(case, condition, cycle, is_local):
    # Compute mean, median, min, max, and volume average
    print(f"-- Loading case {case} condition {condition} cycle {cycle}")
    if is_local:
        mesh_path_laa = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics/hemodynamics_cycle_{cycle:02d}_laa.vtp'
        mesh_path_la = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics/hemodynamics_cycle_{cycle:02d}_la.vtp'
    else:
        mesh_path_laa = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics/hemodynamics_cycle_{cycle:02d}_laa.vtp"
        mesh_path_la = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics/hemodynamics_cycle_{cycle:02d}_la.vtp"

    for name, mesh_path in zip(['LA', 'LAA'], [mesh_path_la, mesh_path_laa]):
        data = read_polydata(mesh_path)
        point_to_cell_data = convert_to_cell_data(data)
        vtk_model = get_volume(point_to_cell_data, volume=False, area=True)
        model = dsa.WrapDataObject(vtk_model)
        print('PART, METRIC, Volume avg, avg, median, min, max')
        for m in ['TAWSS', 'TWSSG', 'OSI', 'RRT', 'ECAP']:
            area = model.CellData.GetArray('Area')
            metric = model.CellData.GetArray(m)

            # Filter negatives
            metric, area = filter_metric(metric, area, m)
            positive_area = np.abs(area)
            total_area = np.sum(positive_area)

            area_avg = (np.dot(positive_area, metric)) / total_area
            avg = np.mean(metric)
            median = np.median(metric)
            minimum = np.min(metric)
            maximum = np.max(metric)

            # TODO: dimension?
            values = [name, m, area_avg, avg, median, minimum, maximum]
            print(values)


def main_volume(case, condition, cycle, metric, is_local):
    # Compute mean, median, min, max, and volume average

    print(f"-- Loading case {case} condition {condition} cycle {cycle} metric {metric}")
    metric_name = 'blood_residence_time' if metric == 'blood_residence_time' else 'energy'
    if is_local:
        mesh_path_laa = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_1029_SR/FlowMetrics/{metric_name}_cycle_{cycle:02d}_laa.vtu'
        mesh_path_la = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_1029_SR/FlowMetrics/{metric_name}_cycle_{cycle:02d}_la.vtu'
    else:
        mesh_path_laa = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/FlowMetrics/hemodynamics_cycle_{cycle:02d}_laa.vtu"
        mesh_path_la = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/FlowMetrics/hemodynamics_cycle_{cycle:02d}_la.vtu"

    for name, mesh_path in zip(['LA', 'LAA'], [mesh_path_la, mesh_path_laa]):
        data = read_polydata(mesh_path)
        point_to_cell_data = convert_to_cell_data(data)
        vtk_model = get_volume(point_to_cell_data)
        model = dsa.WrapDataObject(vtk_model)
        volume = model.CellData.GetArray('Volume')
        quantity = model.CellData.GetArray(metric)
        positive_volume = np.abs(volume)
        total_volume = np.sum(positive_volume)

        # TODO: dimension?
        volume_avg = (np.dot(positive_volume, quantity)) / total_volume
        avg = np.mean(quantity)
        median = np.median(quantity)
        minimum = np.min(quantity)
        maximum = np.max(quantity)
        values = [name, volume_avg, avg, median, minimum, maximum]
        print('PART, Volume avg, avg, median, min, max')
        print(values)


if __name__ == '__main__':
    conditions = ["sr"]
    cases = ['1029']
    cycle = 3
    is_local = True
    metric = "blood_residence_time"
    metric = "kinetic_energy"
    metric = 'hemodynamics'
    metrics = ['kinetic_energy', 'turbulent_kinetic_energy', 'blood_residence_time', 'hemodynamics']
    for metric in metrics:
        for case in cases:
            for condition in conditions:
                if metric in ['kinetic_energy', 'turublent_kinetic_energy', 'blood_residence_time']:
                    main_volume(case, condition, cycle, metric, is_local)
                elif metric == 'hemodynamics':
                    main_surface(case, condition, cycle, is_local)
                else:
                    print(f'Not valie metric {metric}')
