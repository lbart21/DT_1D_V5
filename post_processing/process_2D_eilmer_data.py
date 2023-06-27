import os
import pandas as pd
class Process2DEilmerData():
    __slots__ = ["cell_data", "t_final"]
    def __init__(self, eilmer_data_file) -> None:
        cwd = os.getcwd()

        file_length = sum(1 for line in open(cwd + "/data/" + eilmer_data_file))
        with open(cwd + "/data/" + eilmer_data_file) as f:
            data_start = False
            variable_flag = -1
            for row_ind, row in enumerate(f):
                split_row = row.split(": ")
                if split_row[0] == "sim_time":
                    self.t_final = split_row[1]
                elif split_row[0] == "variables":
                    variable_flag = row_ind + 1
                elif split_row[0] == "nkcell":
                    data_start = True
                    data_start_ind = row_ind + 1
                    n_cells = file_length - data_start_ind
                    print("n_cells:", n_cells)
                    self.cell_data = pd.DataFrame(index=range(n_cells), \
                                        columns = variable_names)
                    continue
                elif row_ind == variable_flag:
                    variable_names = list(row[1:-1].replace('"', "").split(' '))
                    for var_ind, var in enumerate(variable_names):
                        if var == "pos.x":
                            variable_names[var_ind] = "pos_x"
                        elif var == "pos.y":
                            variable_names[var_ind] = "pos_y"
                        elif var == "vel.x":
                            variable_names[var_ind] = "vel_x"
                        elif var == "vel.y":
                            variable_names[var_ind] = "vel_y"
                        elif len(var) > 5 and var[:5] == "massf":
                            variable_names[var_ind] = "massf_" + var.split("-")[1]
                    print(variable_names)
                    
                elif data_start:
                    data = row[1:].split(' ')
                    for var_ind, var in enumerate(variable_names):
                        self.cell_data[var][row_ind - data_start_ind] = float(data[var_ind])

