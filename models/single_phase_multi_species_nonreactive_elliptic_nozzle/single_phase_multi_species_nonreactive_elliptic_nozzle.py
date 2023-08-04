"""
Function:
Author: Luke Bartholomew
Edits:
"""
from numpy import array, pi
from Algorithms.DT_1D_V5.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject
from Algorithms.DT_1D_V5.models.cell_models.single_phase_multi_species_nonreactive_cell \
        import SinglePhaseMultiSpeciesNonreactiveCell
from Algorithms.DT_1D_V5.models.interface_models.single_phase_multi_species_interface \
        import SinglePhaseMultiSpeciesInterface

class SinglePhaseMultiSpeciesNonreactiveEllipticNozzle(SingleInlet1DMeshObject):
    __slots__ = []
    def __init__(self, n_cells, geometry, init_flow_state, restart_flag, flux_scheme, limiter, \
                        recon_scheme, recon_props, update_from, comp_label, reverse_direction_for_ghost_cells = False) -> None:
        super().__init__(n_cells, reverse_direction_for_ghost_cells)
        # geometry = [D_L, D_R, L]
        self.component_labels = [comp_label]

        self.initialise_cells(geometry = geometry, n_cells = n_cells, \
                                init_flow_state = init_flow_state, comp_label = comp_label, \
                                restart_flag = restart_flag)

        self.initialise_interfaces(geometry = geometry, n_cells = n_cells, \
                                    init_flow_state = init_flow_state, flux_scheme = flux_scheme, \
                                    recon_scheme = recon_scheme, limiter = limiter, \
                                    recon_props = recon_props, update_from = update_from)

    def initialise_cells(self, geometry, n_cells, init_flow_state, comp_label, restart_flag):
        [_, _,  L] = geometry
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
                massf = restart_flag[1][cell]["massf"]
                if abs(sum(massf) - 1.0) > 1e-12:
                    massf = (array(massf) / sum(massf)).tolist()
                fs.fluid_state.massf = massf
                fs.fluid_state.update_thermo_from_pT()
                fs.vel_x = restart_flag[1][cell]["vel_x"]

            else:
                fs.fluid_state.copy_values(init_flow_state.fluid_state)
                fs.vel_x = init_flow_state.vel_x

            pos_x_w = cell * L / n_cells
            pos_x_c = (cell + 0.5) * L / n_cells
            pos_x_e = (cell + 1.0) * L / n_cells

            D_w = self.find_diameter_at_x(geometry = geometry, x = pos_x_w)
            #D_c = self.find_diameter_at_x(geometry = geometry, x = pos_x_c)
            D_e = self.find_diameter_at_x(geometry = geometry, x = pos_x_e)
            D_c = 0.5 * (D_w + D_e)

            dx = L / n_cells
            geo = {
                "dx"    :   dx,
                "dV"    :   pi * (D_w ** 2.0 + D_w * D_e + D_e ** 2.0) * dx / 12.0,
                "A_c"   :   0.25 * pi * D_c ** 2.0,
                "A_s"   :   0.25 * pi * (D_w + D_e) * (4.0 * dx ** 2.0 + (D_w - D_e) ** 2.0) ** 0.5,
                "pos_x" :   pos_x_c
            }

            cell_object = SinglePhaseMultiSpeciesNonreactiveCell(cell_id = cell, \
                                                                                label = comp_label)
            
            cell_object.fill_geometry(geometry = geo)
            cell_object.flow_state = fs
            cell_object.initialise_conserved_quantities() 
            self.cell_array[cell] = cell_object

    def initialise_interfaces(self, geometry, n_cells, init_flow_state, flux_scheme, \
                                    recon_scheme, limiter, recon_props, update_from):
        [_, _, L] = geometry
        for interface in range(n_cells + 1):
            pos_x = interface * L / n_cells
            D = self.find_diameter_at_x(geometry = geometry, x = pos_x)
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

            interface_object = SinglePhaseMultiSpeciesInterface(interface_id = interface, \
                                        flux_scheme = flux_scheme, recon_scheme = recon_scheme, \
                                        limiter = limiter, recon_props = recon_props, \
                                        update_from = update_from)
            geo = {"A"      :   0.25 * pi * D ** 2,
                    "pos_x" :   pos_x}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            self.interface_array[interface] = interface_object

    def find_diameter_at_x(self, geometry, x):
        [D_L, D_R, L] = geometry
        return ((D_L ** 2.0 - D_R ** 2.0) * (1.0 - (x / L) ** 2.0) + D_R ** 2.0) ** 0.5
    
    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)