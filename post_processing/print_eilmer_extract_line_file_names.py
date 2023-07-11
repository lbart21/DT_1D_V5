import os
import re
def print_eilmer_extract_line_file_names():
    """
    Want to print Eilmer extract file names in order of ascending pos_x
    file name convention is: "extract_line_at_x" + float(pos_x) + ".txt"
    - Need to look for files who start with "extract_line_at_x"
        - This will require the file names to already be long enough for this to be true
    - Output needs to be in the form '"file_with_pos_x_1", "file_with_pos_x_2", ..., "file_with_pos_x_N"'
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
