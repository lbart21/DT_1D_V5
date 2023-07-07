import os
def print_DT_1D_V5_spatial_cell_data_file_names():
    """
    spatial cell data file name convention: "Sim" + int(sim_number) + "SpatialCellDataAt" + float(time) + "ForComponent" + str(component_name) + ".txt"
    - Additional data is sim_number, time and component_name.
    """
    current_dir = os.getcwd()
    for file in os.listdir(current_dir):
        print(type(file))

print(print_DT_1D_V5_spatial_cell_data_file_names())