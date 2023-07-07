import os
def print_eilmer_extract_line_file_names():
    """
    Want to print Eilmer extract file names in order of ascending pos_x
    file name convention is: "extract_line_at_x" + float(pos_x) + ".txt"
    - Need to look for files who start with "extract_line_at_x"
        - This will require the file names to already be long enough for this to be true
    - Output needs to be in the form '"file_with_pos_x_1", "file_with_pos_x_2", ..., "file_with_pos_x_N"'
    """
