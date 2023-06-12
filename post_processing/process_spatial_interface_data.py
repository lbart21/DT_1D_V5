"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
import pandas

class ProcessSpatialInterfaceData():
    def __init__(self, spatial_interface_data_file) -> None:
        cwd = os.getcwd()
        file_length = sum(1 for line in open(cwd + "/data/" + spatial_interface_data_file))
        file = open(cwd + "/data/" + spatial_interface_data_file)
        variable_flag = -1
        data_start = False
        for row_ind, row in enumerate(file):
            split_row = row.split(": ")
            if split_row[0] == "Time":
                self.t_final = float(split_row[1])
            elif split_row[0] == "Sim":
                self.sim_number = int(split_row[1])
            elif split_row[0] == "Variables":
                variable_flag = row_ind + 1
            if row_ind == variable_flag:
                variable_names = list(row[:-1].split(' '))
                data_start_ind = row_ind + 1
                self.n_interfaces = file_length - data_start_ind
                self.interface_data = pandas.DataFrame(columns = variable_names, index = range(self.n_interfaces))
                data_start = True
                continue
            if data_start:
                for ind, name in enumerate(variable_names):
                    self.interface_data[name][row_ind - data_start_ind] = float(row.split(' ')[ind])
        file.close()
        self.interface_data = self.interface_data.sort_values(by = ["pos_x"])
