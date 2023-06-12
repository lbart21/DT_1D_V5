"""
Function:
Author: Luke Bartholomew
Edits:
"""
import pandas as pd
from Algorithms.DT_1D_V5.post_processing.data_file_to_structured_data \
            import GenerateDataObject

def calculate_net_thrust_from_pressures_notinterpolated_areas(component_data_files, spatial_interface_data):
    data = None

    for component_data_file in component_data_files:
        component_data = GenerateDataObject(data_file_name = component_data_file)
        if data is None:
            data = component_data.component_data
            sim_number = component_data.sim_number
            t_final = component_data.t_final
        else:
            data = pd.concat([data, component_data.component_data])
            data = data.sort_values(by = ["pos_x"])
    
def calculate_net_thrust_from_pressures_interpolated_areas(component_data_files, spatial_interface_data):
    pass