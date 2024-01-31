import pandas as pd
from morphman import *
from vtk.numpy_interface import dataset_adapter as dsa



def main():
    alt_ids_1_9 = [f"cardct{case}" for case in ["1001", "1003", "1004", "1005", "1006", "1007", "1008", "1009"]]
    alt_ids_19_31 = [str(f) for f in
                     [988210734, 988914897, 982060010, 973883493, 980899385, 988848277, 988860504, 988714701, 988659161,
                      988831992, 968252017, 974406807, 975361340]]
    alt_ids_32_35 = ["2 B", "3 M", "4 T", "1 A"]
    alt_ids_74_81 = [f"cardct10{int(case)}" for case in
                     ["0074", "0075", "0076", "0077", "0078", "0079", "0080", "0081"]]
    alt_ids_1029_1039 = [
        f"cardct{case}" for case in
        ["1029", "1030", "1031", "1032", "1033", "1034", "1035", "1036", "1037", "1038", "1039"]]
    alt_ids_2022 = ["none"]
    sources = ["henrik (2023)"] * len(alt_ids_1_9) + ["henrik (2023)"] * len(alt_ids_19_31) + ["meraj (2022)"] * len(
        alt_ids_32_35) + ["carla (2023)"] * len(alt_ids_74_81) + ["carla (2023)"] * len(alt_ids_1029_1039) + [
                  "henrik (2022)"] * len(alt_ids_2022)
    alt_ids = alt_ids_1_9 + alt_ids_19_31 + alt_ids_32_35 + alt_ids_74_81 + alt_ids_1029_1039 + alt_ids_2022
    cases = ["0001", "0003", "0004", "0005", "0006", "0007", "0008", "0009", "0019", "0020", "0021", "0022", "0023",
             "0024", "0025", "0026", "0027", "0028", "0029", "0030", "0031", "0032", "0033", "0034", "0035", "0074",
             "0075", "0076", "0077", "0078", "0079", "0080", "0081", "1029", "1030", "1031", "1032", "1033", "1034",
             "1035", "1036", "1037", "1038", "1039", "2022"]

    conditions = ["af", 'sr']  # , "sr"]
    for condition in conditions:
        data = {
            'case_id': [],
            'laa_probe_id': []
        }
        for j, case in enumerate(cases):
            # print(f"-- Loading case {case} condition {condition}")
            try:
                local = False
                if local:
                    model_path_laa = f"/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_landmarked/{condition}/{case}/model_laa.vtp"
                    laa_cl_path = f"/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_{condition}/{case}/model_laa_centerline.vtp"
                    if condition == "sr":
                        laa_cl_path = f"/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_{condition}/{case}/model_region_centerline_0.vtp"
                    point_path = f"/Users/henriakj/PhD/Code/VaMPy/models/models_inria/models_{condition}/{case}/model_probe_point.json"
                else:
                    model_path_laa = f"/home/opc/Simulation40/{condition.upper()}/{case}/results_moving_atrium/data/1/Hemodynamics/hemodynamics_laa.vtp"
                    laa_cl_path = f"/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/model_region_centerline_0.vtp"
                    point_path = f"/app/OasisMove/src/oasismove/mesh/UKE_{condition.upper()}/{case}/model_probe_point.json"

                laa = read_polydata(model_path_laa)

                # LAA SURFACE AREA
                surface_area = vtk_compute_mass_properties(laa)

                # (1) FIND ORIFICE
                boundary, laa_orifice_area = compute_laa_orifice_area(laa)

                min_id = compute_laa_centerline_stats(boundary, case, laa_cl_path, point_path)

                data['case_id'].append(case)
                data['laa_probe_id'].append(min_id)
            except Exception as e:
                print(f"-- Failed for case {case}: {e}")

        df = pd.DataFrame(data)
        df.to_csv(f"probe_ids_{condition}.csv")


def compute_laa_centerline_stats(boundary, case, laa_cl_path, point_path):
    # (") Get orifice point
    laa_centerline = read_polydata(laa_cl_path)
    boundary_points = dsa.WrapDataObject(boundary).Points
    p_orifice = np.mean(boundary_points, axis=0)
    locator = get_vtk_point_locator(laa_centerline)
    p_orifice_id = locator.FindClosestPoint(p_orifice)
    p_orifice_on_cl = laa_centerline.GetPoint(p_orifice_id)
    id_tmp = p_orifice_id
    p_orifice_on_cl = np.array(p_orifice_on_cl)
    while id_tmp < laa_centerline.GetNumberOfPoints() - 1:
        id_tmp += 1
        p_cl = laa_centerline.GetPoint(id_tmp)
        p_cl = np.array(p_cl)

        if np.linalg.norm(p_cl - p_orifice_on_cl) >= 10:
            break

    with open(point_path, 'r') as file:
        probe_points = json.load(file)
    print(f"case {case}, number of points: {len(probe_points)}")
    min_id = np.argmin(np.linalg.norm(probe_points - p_cl, axis=1))
    return min_id


def compute_laa_orifice_area(laa):
    boundary = vtk_extract_feature_edges(laa)
    delaunay = vtk.vtkDelaunay2D()
    delaunay.SetInputData(boundary)
    delaunay.SetTolerance(0.05)
    delaunay.Update()
    laa_orifice = delaunay.GetOutput()
    orifice_area = vtk_compute_mass_properties(laa_orifice)
    return boundary, orifice_area


if __name__ == '__main__':
    main()
