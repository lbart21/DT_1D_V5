import os
import pandas as pd
from numpy import pi
from Algorithms.DT_1D_V5.post_processing.process_spatial_cell_data import ProcessSpatialCellData
### Read 1D files in to read out pos_x values
def generate_eilmer_extract_line_files(spatial_cell_data_files, n_points, job_name, x_scale):
    base_eilmer_extract_str = 'e4shared --post --job='+ job_name + ' --vtk-xml --add-vars="mach,pitot,total-p,total-T,total-h" '
    data_object = None
    for component_data_file in spatial_cell_data_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
        if data_object is None:
            data_object = data.component_data
        else:
            data_object = pd.concat([data_object, data.component_data], \
                                        axis = 0, ignore_index = True)
    data_object = data_object.sort_values(by = ["pos_x"])

    for ind, x in enumerate(data_object["pos_x"]):
        print(x)
        y_max = (data_object["A_c"][ind] / pi) ** 0.5
        start_pos = str(x_scale*x) + ",0.0,0.0,"
        end_pos = str(x_scale*x) + "," + str(y_max) + ",0.0,"
        extract_line_points_str = '--extract-line=' + start_pos + end_pos + str(n_points)
        output_file_name_str=' --output-file=./data/' + "extract_line_at_x" + str(x) + ".txt > /dev/null"
        final_str = base_eilmer_extract_str + extract_line_points_str + output_file_name_str
        os.system(final_str)
    