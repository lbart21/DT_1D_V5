import os
import re
def print_eilmer_extract_line_file_names():
    """
    Arguments
    --------
    None

    Intermediate Outputs
    -----------

    Returns
    ------------
    None
    """
    current_dir = os.getcwd()
    file_name_list = []
    position_list = []
    for file_name in os.listdir(current_dir + "/data"):
        if len(file_name) > 17 and file_name[:17] == "extract_line_at_x":
            position = float(re.findall("\d+\.\d+", file_name)[0])
            position_list.append(position)
            file_name_list.append(file_name)
    sorted_lists = [name for (pos, name) in sorted(zip(position_list, file_name_list), key=lambda pair: pair[0])]
    print(sorted_lists)
