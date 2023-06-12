"""
Function:
Author: Luke Bartholomew
Edits:
"""
from numpy import pi
from Algorithms.DT_1D_V5.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject
from Algorithms.DT_1D_V5.models.interface_models.single_phase_single_species_interface \
            import SinglePhaseSingleSpeciesInterface
from Algorithms.DT_1D_V5.models.cell_models.single_phase_single_species_nonreactive_cell \
            import SinglePhaseSingleSpeciesNonreactiveCell

class SinglePhaseSingleSpeciesNonreactivePipe(SingleInlet1DMeshObject):
    __slots__ = []
    def __init__(self, restart_flag, init_flow_state, comp_label, geometry, n_cells, limiter, \
                        recon_scheme, recon_props, update_from, flux_scheme, \
                        reverse_direction_for_ghost_cells = False) -> None:
        super().__init__(n_cells, reverse_direction_for_ghost_cells)
        # if restart_flag[0] is True, init_fs will be an empty flow_state object with the desired gas model
        #                               restart_flag[1] will then be a list of dictionaries with 
        #                               p, T and vel_x values for each cell
        # geometry = [D, L]

        self.component_labels = [comp_label]

        self.initialise_cells(geometry = geometry, n_cells = n_cells, \
                                init_flow_state = init_flow_state, \
                                restart_flag = restart_flag, comp_label = comp_label)

        self.initialise_interfaces(geometry = geometry, n_cells = n_cells, limiter = limiter, \
                                    recon_scheme = recon_scheme, recon_props = recon_props, \
                                    update_from = update_from, flux_scheme = flux_scheme, \
                                    init_flow_state = init_flow_state)
    
    def initialise_cells(self, geometry, n_cells, init_flow_state, restart_flag, comp_label):
        [D, L] = geometry

        for cell in range(n_cells):
            flow_state_object = init_flow_state.__class__
            fluid_state_object = init_flow_state.fluid_state.__class__
            fluid_model_filename = init_flow_state.fluid_state.gmodel.file_name
            fluid_model_object = init_flow_state.fluid_state.gmodel.__class__

            gm = fluid_model_object(fluid_model_filename)
            gs = fluid_state_object(gm)
            fs = flow_state_object(gs)

            if restart_flag[0]:
                fs.fluid_state.p = restart_flag[1][cell]["p"]
                fs.fluid_state.T = restart_flag[1][cell]["T"]
                fs.fluid_state.update_thermo_from_pT()
                fs.vel_x = restart_flag[1][cell]["vel_x"]

            else:
                fs.fluid_state.copy_values(init_flow_state.fluid_state)
                fs.vel_x = init_flow_state.vel_x
            
            cell_object = SinglePhaseSingleSpeciesNonreactiveCell(cell_id = cell, \
                                                                      label = comp_label)
            geo = {
                "dx"    :   L / n_cells,
                "dV"    :   0.25 * pi * D ** 2 * L / n_cells,
                "A_c"   :   0.25 * pi * D ** 2,
                "A_s"   :   pi * D * L / n_cells,
                "pos_x" :   (0.5 + cell) * L / n_cells
            }
            cell_object.fill_geometry(geometry = geo)
            cell_object.flow_state = fs
            cell_object.initialise_conserved_quantities() 
            self.cell_array[cell] = cell_object

    def initialise_interfaces(self, geometry, n_cells, limiter, recon_scheme, recon_props, \
                                    update_from, flux_scheme, init_flow_state):
        [D, L] = geometry
        for interface in range(n_cells + 1):
            flow_state_object = init_flow_state.__class__
            fluid_state_object = init_flow_state.fluid_state.__class__
            fluid_model_filename = init_flow_state.fluid_state.gmodel.file_name
            fluid_model_object = init_flow_state.fluid_state.gmodel.__class__
            gm_lft = fluid_model_object(fluid_model_filename)
            gm__rght = fluid_model_object(fluid_model_filename)
            gs_lft = fluid_state_object(gm_lft)
            gs_rght = fluid_state_object(gm__rght)
            fs_lft = flow_state_object(gs_lft)
            fs_rght = flow_state_object(gs_rght)

            interface_object = SinglePhaseSingleSpeciesInterface(\
                                    interface_id = interface, \
                                    flux_scheme = flux_scheme, \
                                    recon_scheme = recon_scheme, \
                                    limiter = limiter, \
                                    recon_props = recon_props, \
                                    update_from = update_from)
            pos_x = interface * L / n_cells
            geo = {"A"      :   0.25 * pi * D ** 2,
                    "pos_x" :   pos_x}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            self.interface_array[interface] = interface_object

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)