from Algorithms.DT_1D_V5.post_processing.process_2D_eilmer_data import Process2DEilmerData
import pandas as pd

class Averaged2DEilmerDataInto1D():
    __slots__ = ["eilmer_cell_data", "averaged_data"]
    def __init__(self, eilmer_cell_data_files, properties, pos_x_array, x_tol) -> None:
        """
        - eilmer_cell_data_files: list(str)
        - properties: list(str)
        - x_pos_array: list(float)
        - x_tol: float
        inputs: 
            - eilmer cell data files 
            - properties that are wanted to be averaged
            - list of x positions of cells to be mapped onto
            - tolerance on x values to get average around desired x value
        output: 
            - pandas data frame with pos_x + properties as column names
            - entries into data frame are averaged values over "columns" of cells 
                about a desired x value
        * Average is determined by volume average (val = sum(vol_i * prop_i) / sum(vol_i))
        """
        self.eilmer_cell_data = None

        for eilmer_cell_data_file in eilmer_cell_data_files:
            eilmer_data_total = Process2DEilmerData(eilmer_data_file = eilmer_cell_data_file)
            if self.eilmer_cell_data is None:
                self.eilmer_cell_data = eilmer_data_total.cell_data
            else:
                self.eilmer_cell_data = pd.concat([self.eilmer_cell_data, eilmer_data_total.cell_data], \
                                                axis = 0, ignore_index = True)
        self.eilmer_cell_data = self.eilmer_cell_data.sort_values(by = ["pos_x"])

        self.averaged_data = pd.DataFrame(index = range(len(pos_x_array)), \
                                        columns = ["pos_x"] + properties)
        for x_ind, x in enumerate(pos_x_array):
            neighbour_cell_ids = []
            for cell_ind in range(len(self.eilmer_cell_data["pos_x"])):
                if abs(self.eilmer_cell_data["pos_x"][cell_ind] - x) < x_tol:
                    neighbour_cell_ids.append(cell_ind)

            for prop in properties:
                average_prop = 0.0
                vol = 0.0

                for cell_id in neighbour_cell_ids:
                    if prop == "vel":
                        pass
                    else:
                        average_prop += self.eilmer_cell_data[prop][cell_id] \
                                        * self.eilmer_cell_data["volume"][cell_id]
                        vol += self.eilmer_cell_data["volume"][cell_id]

                try:
                    self.averaged_data[prop][x_ind] = average_prop / vol
                except:
                    print("No cells found within bounds of x = ", x)
                    self.averaged_data[prop][x_ind] = 0.0
            self.averaged_data["pos_x"][x_ind] = x   