import os
import pandas as pd
class ProcessExtractLineFiles():
    __slots__ = ["cell_data", "n_cells"]
    def __init__(self, extract_line_file) -> None:
        cwd = os.getcwd()
        file_length = sum(1 for line in open(cwd + "/data/" + extract_line_file))
        self.n_cells = file_length - 1
        with open(cwd + "/data/" + extract_line_file) as f:
            for row_ind, row in enumerate(f):
                if row_ind == 0:
                    messy_property_names = row[:-1].split(' ')
                    property_names = [prop.split(":")[1] for prop in messy_property_names]
                    for ind, prop in enumerate(property_names):
                        if prop == "pos.x":
                            property_names[ind] = "pos_x"
                        elif prop == "pos.y":
                            property_names[ind] = "pos_y"
                        elif prop == "vel.x":
                            property_names[ind] = "vel_x"
                        elif prop == "vel.y":
                            property_names[ind] = "vel_y"
                        elif prop == "total_p":
                            property_names[ind] = "p_t"
                        elif prop == "total_h":
                            property_names[ind] = "h_t"
                        elif prop == "total_T":
                            property_names[ind] = "T_t"
                        elif len(prop) > 5 and prop[:5] == "massf":
                            property_names[ind] = "massf_" + property_names[ind].split("-")[1]
                    self.cell_data = pd.DataFrame(index=range(self.n_cells), \
                                        columns = property_names)
                else:
                    data = row[:-1].split(" ")
                    cell_idx = row_ind - 1
                    for ind, prop in enumerate(property_names):
                        self.cell_data[prop][cell_idx] = float(data[ind])
