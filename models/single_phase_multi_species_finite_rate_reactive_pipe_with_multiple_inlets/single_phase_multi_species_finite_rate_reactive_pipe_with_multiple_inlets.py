"""
Function:
Author: Luke Bartholomew
Edits:
"""
from numpy import pi, array
from Algorithms.DT_1D_V5.models.prefilled_multi_inlet_mesh_object import MultiInlet1DMeshObject
from Algorithms.DT_1D_V5.models.cell_models.single_phase_multi_species_finite_rate_reactive_cell \
            import SinglePhaseMultiSpeciesFiniteRateReactiveCell
from Algorithms.DT_1D_V5.models.interface_models.single_phase_multi_species_interface \
            import SinglePhaseMultiSpeciesInterface

class SinglePhaseMultiSpeciesFiniteRateReactivePipeWithMultipleInlets(MultiInlet1DMeshObject):
    __slots__ = []
    def __init__(self, n_cells, n_inlets, bulk_geometry, inlet_areas, init_flow_state, \
                        comp_label, flux_scheme, recon_scheme, limiter, recon_props, \
                        update_from, restart_flag, \
                        reverse_direction_for_mirrored_flow = False) -> None:
        super().__init__(n_cells, n_inlets, reverse_direction_for_mirrored_flow)
        # bulk_geometry = [D, L]
        # inlet_areas = [A_1, A_2, ..., A_n]

        self.component_labels = [comp_label]

        self.initialise_cells(bulk_geometry = bulk_geometry, n_cells = n_cells, \
                                init_flow_state = init_flow_state, comp_label = comp_label, \
                                restart_flag = restart_flag)

        self.initialise_interfaces(bulk_geometry = bulk_geometry, inlet_areas = inlet_areas, \
                                    n_cells = n_cells, init_flow_state = init_flow_state, \
                                    flux_scheme = flux_scheme, recon_scheme = recon_scheme, \
                                    limiter = limiter, recon_props = recon_props, \
                                    update_from = update_from, \
                                    reverse_direction = reverse_direction_for_mirrored_flow)

    def initialise_cells(self, bulk_geometry, n_cells, init_flow_state, comp_label, restart_flag):
        [D, L] = bulk_geometry

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

            cell_object = SinglePhaseMultiSpeciesFiniteRateReactiveCell(\
                                            cell_id = cell, label = comp_label)
            
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

    def initialise_interfaces(self, bulk_geometry, inlet_areas, n_cells, init_flow_state, \
                                    flux_scheme, recon_scheme, limiter, recon_props, \
                                    update_from, reverse_direction):
        [D, L] = bulk_geometry
        inlet_interface_array = [None] * len(inlet_areas)
        bulk_interface_array = [None] * n_cells

        for inlet_ind, area in enumerate(inlet_areas):
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

            if reverse_direction: #Inlets are on the east boundary
                interface_id = n_cells + inlet_ind
                pos_x = L
            else:
                interface_id = inlet_ind
                pos_x = 0.0

            interface_object = SinglePhaseMultiSpeciesInterface(interface_id = interface_id, \
                                        flux_scheme = flux_scheme, recon_scheme = recon_scheme, \
                                        limiter = limiter, recon_props = recon_props, \
                                        update_from = update_from)

            geo = {"A"      :   area,
                    "pos_x" :   pos_x}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            inlet_interface_array[inlet_ind] = interface_object

            gm_lft = fluid_model_object(fluid_model_filename)
            gm__rght = fluid_model_object(fluid_model_filename)
            gs_lft = fluid_state_object(gm_lft)
            gs_rght = fluid_state_object(gm__rght)
            fs_lft = flow_state_object(gs_lft)
            fs_rght = flow_state_object(gs_rght)
        
        for interface in range(n_cells):
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

            if reverse_direction: #Inlets are on the east boundary
                interface_id = interface
                pos_x = interface * L / n_cells
            else:
                interface_id = len(inlet_areas) + interface
                pos_x = (interface + 1) * L / n_cells

            interface_object = SinglePhaseMultiSpeciesInterface(interface_id = interface_id, \
                                    flux_scheme = flux_scheme, recon_scheme = recon_scheme, \
                                    limiter = limiter, recon_props = recon_props, \
                                    update_from = update_from)

            geo = {"A"      :   0.25 * pi * D ** 2,
                    "pos_x" :   pos_x}
            interface_object.fill_geometry(geometry = geo)
            interface_object.lft_state = fs_lft
            interface_object.rght_state = fs_rght
            bulk_interface_array[interface] = interface_object
        
        if reverse_direction:
            self.interface_array = bulk_interface_array + inlet_interface_array
        else:
            self.interface_array = inlet_interface_array + bulk_interface_array

    def add_boundary_conditions(self, BC):
        self.boundary_conditions.append(BC)