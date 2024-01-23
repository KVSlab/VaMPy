import numpy as np
import pandas as pd
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


def main_surface(case, condition, cycle, region, metric_name, is_local):
    # Compute mean, median, min, max, and volume average
    print(f"-- Loading case {case} condition {condition} cycle {cycle} metric {metric_name}")
    if is_local:
        mesh_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_{case}_{condition}/Hemodynamics/hemodynamics_{region}.vtp'
    else:
        mesh_path = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics/hemodynamics_{region}.vtp"

    if cycle != 1:
        metric_array_name = f"{metric_name}_input_{cycle - 1}"
    else:
        metric_array_name = metric_name

    data = read_polydata(mesh_path)
    point_to_cell_data = convert_to_cell_data(data)
    vtk_model = get_volume(point_to_cell_data, volume=False, area=True)
    model = dsa.WrapDataObject(vtk_model)

    area = model.CellData.GetArray('Area')
    metric = model.CellData.GetArray(metric_array_name)

    # Filter negatives
    metric, area = filter_metric(metric, area, metric_name)
    positive_area = np.abs(area)
    total_area = np.sum(positive_area)

    area_avg = (np.dot(positive_area, metric)) / total_area
    avg = np.mean(metric)
    median = np.median(metric)
    minimum = np.min(metric)
    maximum = np.max(metric)

    return area_avg, avg, median, minimum, maximum


def main_volume(case, condition, cycle, metric, region, is_local):
    # Compute mean, median, min, max, and volume average

    print(f"-- Loading case {case} condition {condition} cycle {cycle} metric {metric}")
    metric_name = 'blood_residence_time' if metric == 'blood_residence_time' else 'energy'
    if is_local:
        mesh_path = f'/Users/henriakj/PhD/Code/OasisMove/results_34case/results_1029_SR/FlowMetrics/{metric_name}_{region}.vtu'
    else:
        mesh_path = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/FlowMetrics/{metric_name}_{region}.vtu"
    if cycle != 1:
        metric_name = f"{metric}_input_{cycle - 1}"
    else:
        metric_name = metric
    data = read_polydata(mesh_path)
    point_to_cell_data = convert_to_cell_data(data)
    vtk_model = get_volume(point_to_cell_data)
    model = dsa.WrapDataObject(vtk_model)
    volume = model.CellData.GetArray('Volume')
    quantity = model.CellData.GetArray(metric_name)
    positive_volume = np.abs(volume)
    total_volume = np.sum(positive_volume)

    # TODO: dimension?
    volume_avg = (np.dot(positive_volume, quantity)) / total_volume
    avg = np.mean(quantity)
    median = np.median(quantity)
    minimum = np.min(quantity)
    maximum = np.max(quantity)

    return volume_avg, avg, median, minimum, maximum


if __name__ == '__main__':
    conditions = ["sr", 'af']
    cases = ["0003", "0004", "0005", "0006", "0007", "0008", "0009", "0019",
             "0020", "0021", "0023", "0024", "0025", "0026", "0027", "0028",
             "0029", "0030", "0031", "0032", "0033", "0034", "0035", "0074",
             "0076", "0077", "0078", "0080", "0081", "1029", "1030", "1031",
             "1032", "1033", "1035", "1037", "1038", "1039", "2022"][:4]
    cycles = [1, 2, 3, 4, 5]
    is_local = False
    metrics = ['kinetic_energy', 'turbulent_kinetic_energy', 'blood_residence_time', 'hemodynamics']
    regions = ['laa', 'la']
    for condition in conditions:
        for region in regions:
            save_path = f'metrics_{condition}_{region}.csv'
            data = {
                'case_id': [],
                'metric': [],
                'cycle': [],
                'volume_avg': [],
                'avg': [],
                'median': [],
                'minimum': [],
                'maximum':[]
            }
            for case in cases:
                for metric in metrics:
                    for cycle in cycles:
                        if metric in ['kinetic_energy', 'turbulent_kinetic_energy', 'blood_residence_time']:
                            volume_avg, avg, median,minimum,maximum = main_volume(case, condition, cycle, metric, region, is_local)
                            # Append the current case, metric, and cycle information to the data dictionary
                            data['case_id'].append(case)
                            data['metric'].append(metric)
                            data['cycle'].append(cycle)
                            data['volume_avg'].append(volume_avg)
                            data['avg'].append(avg)
                            data['median'].append(median)
                            data['minimum'].append(minimum)
                            data['maximum'].append(maximum)
                        elif metric == 'hemodynamics':
                            for m in ['TAWSS', 'TWSSG', 'OSI', 'RRT', 'ECAP']:
                                volume_avg, avg, median,minimum,maximum = main_surface(case, condition, cycle, region, m, is_local)
                                data['case_id'].append(case)
                                data['metric'].append(metric)
                                data['cycle'].append(cycle)
                                data['volume_avg'].append(volume_avg)
                                data['avg'].append(avg)
                                data['median'].append(median)
                                data['minimum'].append(minimum)
                                data['maximum'].append(maximum)
                        else:
                            print(f'Not valid metric {metric}')
            df = pd.DataFrame(data)
            df.to_csv(save_path)
