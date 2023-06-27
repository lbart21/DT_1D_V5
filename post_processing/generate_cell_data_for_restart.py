from Algorithms.DT_1D_V5.post_processing.process_spatial_cell_data import ProcessSpatialCellData

def generate_cell_data_for_component(component_data_file):
    detailed_data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
    detailed_component_data = detailed_data.component_data
    t_final = detailed_data.t_final
    cell_data = [{} for _ in range(detailed_data.n_cells)]
    multi_species = False

    spcs_names = []

    for column_name in detailed_component_data.columns:
        if "massf" in column_name:
            multi_species = True
            spcs_names.append(column_name)

    for i in range(detailed_data.n_cells):
        cell_data[i]["p"] = detailed_component_data["p"][i]
        cell_data[i]["T"] = detailed_component_data["T"][i]
        cell_data[i]["vel_x"] = detailed_component_data["vel_x"][i]

        if multi_species:
            massf = [None] * len(spcs_names)
            for ind, spcs in enumerate(spcs_names):
                massf[ind] = detailed_component_data[spcs][i]
            cell_data[i]["massf"] = massf

    return cell_data, t_final

