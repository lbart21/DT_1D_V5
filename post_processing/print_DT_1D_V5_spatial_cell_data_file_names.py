import os
import re
def print_DT_1D_V5_spatial_cell_data_file_names_for_multiple_times(sim_number, labels, times):
    """
    spatial cell data file name convention: "Sim" + int(sim_number) + "SpatialCellDataAt" + float(time) + "ForComponent" + str(component_name) + ".txt"
    - Additional data is sim_number, time and component_name.
    """
    current_dir = os.getcwd()
    file_name_list = [[] * len(times)]
    time_list = []
    for file_name in os.listdir(current_dir + "/data"):
        if "Sim" + str(sim_number) == file_name[:len("Sim" + str(sim_number))]:
            # We've got the right simulation 
            if file_name[len("Sim" + str(sim_number)):len("Sim" + str(sim_number)) + 17] == "SpatialCellDataAt":
                # We've got the right data file type for the simulation
                file_time = float(re.findall("\d+\.\d+", file_name)[0])
                if len(times) == 0: # Want to  find all files regardless of time value
                    if file_time not in time_list:
                        time_list.append(file_time)
                    time_idx = time_list.index(file_time)
                    comp_name_location = file_name.find("Component") + 9
                    for label in labels:
                        # This is a file we want to extract
                        if file_name[comp_name_location:-4] == label:
                            file_name_list[time_idx].append(file_name)
                else:
                    for time_idx, time in enumerate(times):
                        if file_time == time:
                            # We've got the right time snapshot for the right data file type for the right simulation
                            comp_name_location = file_name.find("Component") + 9
                            for label in labels:
                                # This is a file we want to extract
                                if file_name[comp_name_location:-4] == label:
                                    file_name_list[time_idx].append(file_name)
    print(file_name_list)
