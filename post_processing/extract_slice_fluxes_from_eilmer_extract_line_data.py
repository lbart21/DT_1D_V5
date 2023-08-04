from numpy import linspace, pi
from Algorithms.DT_1D_V5.post_processing.process_extract_line_files import ProcessExtractLineFiles
from Algorithms.DT_1D_V5.post_processing.generate_average_from_eilmer_extract_line import find_nearest_cell_idx
def extract_slice_fluxes_from_eilmer_extract_line_data(line_extract_data_file, x, \
                                                       ymin, ymax, n, extract_vars):
    """_summary_

    Args:
        line_extract_data_file (_type_): _description_
        x (_type_): _description_
        ymax (_type_): _description_
        n (_type_): _description_
        extract_vars (_type_): _description_

    Returns:
        _type_: _description_
    """
    y_sample = linspace(ymin, ymax, n)
    y_boundaries = 0.5 * (y_sample[:-1] + y_sample[1:])
    y_sample = y_sample.tolist()
    y_boundaries = [ymin] + y_boundaries.tolist() + [ymax]
    data = ProcessExtractLineFiles(extract_line_file = line_extract_data_file)
    cell_data = data.cell_data
    n_cells = data.n_cells
    nearest_cell_idxs = [None] * n

    for i in range(n):
        nearest_cell_idxs[i] = find_nearest_cell_idx(sample_location = [x, y_sample[i]], cell_array = cell_data)
    
    fluxes = {key:0.0 for key in extract_vars}

    for i in range(n):
        area_i = pi * (y_boundaries[i + 1] ** 2.0 - y_boundaries[i] ** 2.0)
        rho_nc = cell_data["rho"][nearest_cell_idxs[i]]
        vel_x_nc = cell_data["vel_x"][nearest_cell_idxs[i]]
        if "mass_flux" in extract_vars:
            fluxes["mass_flux"] += rho_nc * vel_x_nc * area_i
            
        if "xMom_flux" in extract_vars:
            vel_y_nc = cell_data["vel_y"][nearest_cell_idxs[i]]
            vel_nc = (vel_x_nc ** 2.0 + vel_y_nc ** 2.0 ) ** 0.5
            fluxes["xMom_flux"] += rho_nc * vel_nc * area_i * vel_x_nc

        if "energy_flux" in extract_vars:
            u_nc = cell_data["u"][nearest_cell_idxs[i]]
            p_nc = cell_data["p"][nearest_cell_idxs[i]]
            vel_y_nc = cell_data["vel_y"][nearest_cell_idxs[i]]
            vel_nc = (vel_x_nc ** 2.0 + vel_y_nc ** 2.0 ) ** 0.5
            fluxes["energy_flux"] += rho_nc * (u_nc + p_nc / rho_nc + 0.5 * vel_nc ** 2.0) * area_i * vel_x_nc
    
    return fluxes
