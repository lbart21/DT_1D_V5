from numpy import linspace, pi
from Algorithms.DT_1D_V5.post_processing.process_extract_line_files import ProcessExtractLineFiles
def generate_average_from_eilmer_extract(line_extract_data_file, x, ymin, ymax, n, extract_vars):
    """
    Pass in data file that Eilmer --extract-line post processing generates.
    Produce data object from this file
    Pass in x position, ymin, ymax and number of sample points.
    Generate line of sample points between ymin and ymax at the given x value
    At the same time, generate array of positions for "cell" boundaries used for area calculation. -> [ymin, y1, y2, ..., ymax] yi = 0.5 * ave(consectutive sample points)
    Iterate over sample points and find the closest cell
    average = sum(phi * pi * (y2^2 - y1^2))/sum(pi*(y2^2 - y1^2)))
    """
    # TODO include flag to do average average or mass flux average
    # mass flux average = sum(phi * rho * vel_x * pi * (y2^2 - y1^2)) / sum(rho * vel_x * pi*(y2^2 - y1^2))
    y_sample = linspace(ymin, ymax, n+2)
    y_sample = y_sample[1:-1]
    y_boundaries = 0.5 * (y_sample[:-1] + y_sample[1:])
    y_sample = y_sample.tolist()
    y_boundaries = [ymin] + y_boundaries.tolist() + [ymax]

    data = ProcessExtractLineFiles(extract_line_file = line_extract_data_file)
    cell_data = data.cell_data
    n_cells = data.n_cells
    nearest_cell_idxs = [None] * n
    for i in range(n):
        nearest_cell_idxs[i] = find_nearest_cell_idx(sample_location = [x, y_sample[i]], cell_array = cell_data)
    
    if "vel" in extract_vars:
        cell_data["vel"] = (cell_data["vel_x"] ** 2.0 + cell_data["vel_y"] ** 2.0) ** 0.5
    if "Ma" in extract_vars:
        cell_data["Ma"] = (cell_data["vel_x"] ** 2.0 + cell_data["vel_y"] ** 2.0) ** 0.5 / cell_data["a"]
    
    average_props = {}
    max_props = {}
    min_props = {}
    max_vals = cell_data.max()
    min_vals = cell_data.min()
    for var in extract_vars:
        max_props[var] = max_vals[var]
        min_props[var] = min_vals[var]
        average = 0.0
        area = 0.0
        for i in range(n):
            area_i = pi * (y_boundaries[i + 1] ** 2.0 - y_boundaries[i] ** 2.0)
            area += area_i
            if var == "Ma":
                vel_x = cell_data["vel_x"][nearest_cell_idxs[i]]
                vel_y = cell_data["vel_y"][nearest_cell_idxs[i]]
                a = cell_data["a"][nearest_cell_idxs[i]]
                val = (vel_x ** 2.0 + vel_y ** 2.0) ** 0.5 / a
                average += val * area_i
                
            elif var == "vel":
                vel_x = cell_data["vel_x"][nearest_cell_idxs[i]]
                vel_y = cell_data["vel_y"][nearest_cell_idxs[i]]
                val = (vel_x ** 2.0 + vel_y ** 2.0) ** 0.5
                average += val * area_i

            else:
                val = cell_data[var][nearest_cell_idxs[i]]
                average += val * area_i

        average_props[var] = average / area
    return average_props, min_props, max_props

def find_nearest_cell_idx(sample_location, cell_array):
    pos_x, pos_y = sample_location
    nearest_cell_idx = 0
    min_distance = 1e6 #placeholder value
    for i in range(cell_array.shape[0]):
        cell_pos_x = cell_array["pos_x"][i]
        cell_pos_y = cell_array["pos_y"][i]
        dist = ((pos_x - cell_pos_x) ** 2.0 + (pos_y - cell_pos_y) ** 2.0)
        if dist < min_distance:
            nearest_cell_idx = i
            min_distance = dist
    return nearest_cell_idx