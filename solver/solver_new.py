"""
Function:
Author: Luke Bartholomew
Edits:
"""
from copy import deepcopy
from time import perf_counter
from Algorithms.DT_1D_V5.solver.write_simulation_description import write_simulation_description
from Algorithms.DT_1D_V5.solver.write_spatial_cell_data_to_file import write_spatial_cell_data_to_file
from Algorithms.DT_1D_V5.solver.write_spatial_interface_data_to_file import write_spatial_interface_data_to_file
from Algorithms.DT_1D_V5.solver.write_transient_cell_data_to_file import write_transient_cell_data_to_file
from Algorithms.DT_1D_V5.solver.write_transient_interface_data_to_file import write_transient_interface_data_to_file
from Algorithms.DT_1D_V5.solver.write_interface_geometry_data import write_interface_geometry_data
from Algorithms.DT_1D_V5.integrate.integrate_new import Integrate
from Algorithms.DT_1D_V5.extras.join_blocks import JointBlock
from Algorithms.DT_1D_V5.models.prefilled_single_inlet_mesh_object import SingleInlet1DMeshObject

from Algorithms.DT_1D_V5.reconstruction.locate_neighbouring_cell_indices import find_idx_of_cell_recursively
from Algorithms.DT_1D_V5.reconstruction.locate_neighbouring_interface_indices import find_idx_of_interface_recursively

class Solver():
    __slots__ = ["sim_number", "simulation_description", "restart_flag", "n_ghost_cells", "mesh", \
                    "t_current", \
                    "transient_cell_flow_property_variables_to_write", \
                    "transient_interface_flow_property_variables_to_write", \
                    "spatial_cell_flow_property_variables_to_write", \
                    "spatial_interface_flow_property_varaible_to_write"]
    def __init__(self, mesh_object, transient_cell_flow_property_variables_to_write, \
                        restart_flag, transient_interface_flow_property_variables_to_write, \
                        sim_number, spatial_cell_flow_property_variables_to_write, \
                        spatial_interface_flow_property_varaible_to_write, \
                        simulation_description) -> None:
        """
        meshObject = object with attributes: cell_array, interface_array, 
                                                map_cell_id_to_west_interface_idx, 
                                                map_cell_id_to_east_interface_idx,
                                                map_interface_id_to_west_cell_idx,
                                                map_interface_id_to_east_cell_idx,
                                                boundary_conditions,
                                                boundary_interface_ids,
                                                component_labels,
                                                cell_idx_to_track,
                                                interface_idx_to_track
        cfl_flag = [Bool, float]
        t_final = 
        data_save_dt = 
        transient_cell_flow_property_variables_to_write = 
        transient_interface_flow_property_variables_to_write = 
        spatial_cell_flow_property_variables_to_write = 
        sim_number = 
        restart_flag = 
        rapid_data_save_steps = 
        simulation_description = 
        """
        self.sim_number = sim_number
        self.simulation_description = simulation_description

        self.restart_flag = restart_flag
        self.n_ghost_cells = 0

        self.mesh = mesh_object

        self.t_current = 0.0
        if self.restart_flag[0]:
            self.t_current = self.restart_flag[1]

        self.transient_cell_flow_property_variables_to_write = transient_cell_flow_property_variables_to_write
        self.transient_interface_flow_property_variables_to_write = transient_interface_flow_property_variables_to_write
        self.spatial_cell_flow_property_variables_to_write = spatial_cell_flow_property_variables_to_write
        self.spatial_interface_flow_property_varaible_to_write = spatial_interface_flow_property_varaible_to_write

        ### Write data to files
        write_spatial_cell_data_to_file(cell_array = mesh_object.cell_array, \
                        time = self.t_current, labels = mesh_object.component_labels, \
                        flow_property_variables = self.spatial_cell_flow_property_variables_to_write, \
                        sim_number = sim_number, \
                        boundary_conditions = mesh_object.boundary_conditions, \
                        simulation_description = self.simulation_description)
        
        for cell_id in mesh_object.cell_idx_to_track:
            write_transient_cell_data_to_file(cell = mesh_object.cell_array[cell_id], \
                        time = self.t_current, \
                        flow_property_variables = self.transient_cell_flow_property_variables_to_write, \
                        sim_number = sim_number, simulation_description = self.simulation_description)
        
        write_interface_geometry_data(interface_array = mesh_object.interface_array, \
                                        sim_number = sim_number)


        ### Add ghost cells to mesh, but don't put the flow properties into the cells
        self.add_ghost_cells()


    
    def start_simulation(self, cfl_flag, t_final, spatial_data_save_dt, temporal_data_save_dt, n_steps_print = 100):
        write_simulation_description(sim_number = self.sim_number, \
                                        simulation_description = self.simulation_description, \
                                        boundary_conditions = self.mesh.boundary_conditions, \
                                        cfl_flag = cfl_flag, t_final = t_final)

        time_tol = 1e-12
        t_write_tol = 1e-12

        current_step = 0
        spatial_data_written = False
        temporal_data_written = False
        spatial_t_write = self.t_current + spatial_data_save_dt
        temporal_t_write = self.t_current + temporal_data_save_dt

        while self.t_current < t_final and abs(self.t_current - t_final) > time_tol:
            spatial_data_written = False
            temporal_data_written = False
            if (current_step + 1) % n_steps_print == 0:
                time_counter_start = perf_counter()
            new_data = Integrate(mesh = self.mesh, cfl_flag = cfl_flag, t_current = self.t_current, current_step = current_step)
            if (current_step + 1) % n_steps_print == 0:
                time_counter_end = perf_counter()
            dt = new_data.dt_total
            self.t_current += dt
            current_step += 1
            if current_step % n_steps_print == 0:
                print("t: ", self.t_current, "WCTFT:", (t_final - self.t_current) / dt * (time_counter_end - time_counter_start), "CFL:", new_data.cfl_used)
            self.mesh = new_data.mesh

            # Check if we've just done the first time step so we should write out the spatial 
            # interface data
            if current_step == 1:
                write_spatial_interface_data_to_file(\
                    interface_array = self.mesh.interface_array[:-1 * self.n_ghost_cells], \
                    sim_number = self.sim_number, time = self.t_current - dt, \
                    simulation_description = self.simulation_description, \
                    interface_property_variables = self.spatial_interface_flow_property_varaible_to_write)
                
            # Check if we're due to write out transient interface and cell data
            if self.t_current > temporal_t_write - t_write_tol:
                for cell_id in self.mesh.cell_idx_to_track:
                    write_transient_cell_data_to_file(cell = self.mesh.cell_array[cell_id], time = self.t_current, sim_number = self.sim_number, \
                                        flow_property_variables = self.transient_cell_flow_property_variables_to_write, \
                                        simulation_description = self.simulation_description)
                
                for interface_id in self.mesh.interface_idx_to_track:
                    write_transient_interface_data_to_file(interface = self.mesh.interface_array[interface_id], time = self.t_current - dt, \
                                                flow_property_variables = self.transient_interface_flow_property_variables_to_write, sim_number = self.sim_number, \
                                                simulation_description = self.simulation_description)
                temporal_data_written = True
                temporal_t_write += temporal_data_save_dt
            
            # Check if we're due to write out snapshot
            if self.t_current > spatial_t_write - t_write_tol:
                print("Writing data, t = ", self.t_current)
                write_spatial_cell_data_to_file(cell_array = self.mesh.cell_array[:-1 * self.n_ghost_cells], \
                    time = self.t_current, labels = self.mesh.component_labels, \
                    flow_property_variables = self.spatial_cell_flow_property_variables_to_write, \
                    sim_number = self.sim_number, \
                    boundary_conditions = self.mesh.boundary_conditions, \
                    simulation_description = self.simulation_description)
                write_spatial_interface_data_to_file(\
                    interface_array = self.mesh.interface_array[:-1 * self.n_ghost_cells], \
                    sim_number = self.sim_number, time = self.t_current - dt, \
                    simulation_description = self.simulation_description, \
                    interface_property_variables = self.spatial_interface_flow_property_varaible_to_write)
                spatial_data_written = True
                spatial_t_write += spatial_data_save_dt
            
        if not spatial_data_written:
            write_spatial_cell_data_to_file(cell_array = self.mesh.cell_array[:-1 * self.n_ghost_cells], time = self.t_current, \
                            labels = self.mesh.component_labels, sim_number = self.sim_number, \
                            flow_property_variables = self.spatial_cell_flow_property_variables_to_write, \
                            boundary_conditions = self.mesh.boundary_conditions, \
                            simulation_description = self.simulation_description)
            
        if not temporal_data_written:
            for cell_id in new_data.mesh.cell_idx_to_track:
                write_transient_cell_data_to_file(cell = self.mesh.cell_array[cell_id], time = self.t_current, sim_number = self.sim_number, \
                                    flow_property_variables = self.transient_cell_flow_property_variables_to_write, \
                                    simulation_description = self.simulation_description)
            for interface_id in new_data.mesh.interface_idx_to_track:
                write_transient_interface_data_to_file(interface = new_data.mesh.interface_array[interface_id], \
                                    time = self.t_current - new_data.dt_total, sim_number = self.sim_number, \
                                    flow_property_variables = self.transient_interface_flow_property_variables_to_write, \
                                    simulation_description = self.simulation_description)

    def add_ghost_cells(self):

        for bc_ind, bc in enumerate(self.mesh.boundary_conditions): # bc = [int(interface_id), int(n_layers), str(east or west), [str(condition type), list(additional info)]]
            interface_id = bc[0]
            n_gc_layers = bc[1]
            boundary_edge = bc[2]

            cell_pairs = [] #Going to be a list of tuple pairs, first element will be ghost cell id, 
                            # second will be interior cell id that the ghost cell mirrors
            
            n_cells_current = len(self.mesh.cell_array)
            
            self.n_ghost_cells += n_gc_layers
            if boundary_edge == "East":
                gc_layer = SingleInlet1DMeshObject(n_cells = n_gc_layers, reverse_direction_for_ghost_cells = False)
                for cell in range(n_gc_layers):
                    inner_cell_id, _ = find_idx_of_cell_recursively(interface_id = interface_id, \
                                            recursion_depth = cell + 1, direction = "West", \
                                            map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                            map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                            map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                            map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx)
                    inner_cell = self.mesh.cell_array[inner_cell_id]
                    ghost_cell = deepcopy(inner_cell)

                    ### Create clone of cell
                    gm_class = inner_cell.flow_state.fluid_state.gmodel.__class__
                    gm_filename = inner_cell.flow_state.fluid_state.gmodel.file_name
                    gs_class = inner_cell.flow_state.fluid_state.__class__
                    fs_class = inner_cell.flow_state.__class__

                    gm = gm_class(gm_filename)
                    gs = gs_class(gm)
                    gs.copy_values(inner_cell.flow_state.fluid_state)
                    fs = fs_class(gs, inner_cell.flow_state.vel_x)

                    ghost_cell.flow_state = fs
                    ghost_cell.interior_cell_flag = False

                    gc_layer.cell_array[cell] = ghost_cell

                    cell_pairs.append((n_cells_current + cell, inner_cell_id))
                
                for interface in range(n_gc_layers + 1):
                    inner_interface_id = find_idx_of_interface_recursively(interface_id = interface_id, \
                                            recursion_depth = interface, direction = "West", \
                                            map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                            map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                            map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                            map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx)
                    
                    inner_interface = self.mesh.interface_array[inner_interface_id]
                    ghost_interface = deepcopy(inner_interface)

                    #lft_gm_class = inner_interface.lft_state.fluid_state.gmodel.__class__
                    #lft_gm_filename = inner_interface.lft_state.fluid_state.gmodel.file_name
                    #lft_gs_class = inner_interface.lft_state.fluid_state.__class__
                    #lft_fs_class = inner_interface.lft_state.__class__

                    #lft_gm = lft_gm_class(lft_gm_filename)
                    #lft_gs = lft_gs_class(lft_gm)
                    #lft_gs.copy_values(inner_interface.lft_state.fluid_state)
                    #lft_fs = lft_fs_class(lft_gs, inner_interface.lft_state.vel_x)


                    #rght_gm_class = inner_interface.rght_state.fluid_state.gmodel.__class__
                    #rght_gm_filename = inner_interface.rght_state.fluid_state.gmodel.file_name
                    #rght_gs_class = inner_interface.rght_state.fluid_state.__class__
                    #rght_fs_class = inner_interface.rght_state.__class__
                    
                    #rght_gm = rght_gm_class(rght_gm_filename)
                    #rght_gs = rght_gs_class(rght_gm)
                    #rght_gs.copy_values(inner_interface.rght_state.fluid_state)
                    #rght_fs = rght_fs_class(rght_gs, inner_interface.rght_state.vel_x)

                    #ghost_interface.lft_state = lft_fs
                    #ghost_interface.rght_state = rght_fs

                    ghost_interface.flux_flag = False

                    gc_layer.interface_array[interface] = ghost_interface

            elif boundary_edge == "West":
                gc_layer = SingleInlet1DMeshObject(n_cells = n_gc_layers, reverse_direction_for_ghost_cells = True)
                for cell in range(n_gc_layers):
                    inner_cell_id, _ = find_idx_of_cell_recursively(interface_id = interface_id, \
                                            recursion_depth = cell + 1, direction = "East", \
                                            map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                            map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                            map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                            map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx)
                    
                    inner_cell = self.mesh.cell_array[inner_cell_id]
                    ghost_cell = deepcopy(inner_cell)

                    ### Create clone of cell
                    gm_class = inner_cell.flow_state.fluid_state.gmodel.__class__
                    gm_filename = inner_cell.flow_state.fluid_state.gmodel.file_name
                    gs_class = inner_cell.flow_state.fluid_state.__class__
                    fs_class = inner_cell.flow_state.__class__

                    gm = gm_class(gm_filename)
                    gs = gs_class(gm)
                    gs.copy_values(inner_cell.flow_state.fluid_state)
                    fs = fs_class(gs, inner_cell.flow_state.vel_x)

                    ghost_cell.flow_state = fs
                    ghost_cell.interior_cell_flag = False

                    gc_layer.cell_array[cell] = ghost_cell

                    cell_pairs.append((n_cells_current + cell, inner_cell_id))
                
                for interface in range(n_gc_layers + 1):
                    inner_interface_id = find_idx_of_interface_recursively(interface_id = interface_id, \
                                            recursion_depth = interface, direction = "East", \
                                            map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                            map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                            map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                            map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx)
                    
                    inner_interface = self.mesh.interface_array[inner_interface_id]
                    ghost_interface = deepcopy(inner_interface)

                    #lft_gm_class = inner_interface.lft_state.fluid_state.gmodel.__class__
                    #lft_gm_filename = inner_interface.lft_state.fluid_state.gmodel.file_name
                    #lft_gs_class = inner_interface.lft_state.fluid_state.__class__
                    #lft_fs_class = inner_interface.lft_state.__class__

                    #lft_gm = lft_gm_class(lft_gm_filename)
                    #lft_gs = lft_gs_class(lft_gm)
                    #lft_gs.copy_values(inner_interface.lft_state.fluid_state)
                    #lft_fs = lft_fs_class(lft_gs, inner_interface.lft_state.vel_x)


                    #rght_gm_class = inner_interface.rght_state.fluid_state.gmodel.__class__
                    #rght_gm_filename = inner_interface.rght_state.fluid_state.gmodel.file_name
                    #rght_gs_class = inner_interface.rght_state.fluid_state.__class__
                    #rght_fs_class = inner_interface.rght_state.__class__
                    
                    #rght_gm = rght_gm_class(rght_gm_filename)
                    #rght_gs = rght_gs_class(rght_gm)
                    #rght_gs.copy_values(inner_interface.rght_state.fluid_state)
                    #rght_fs = rght_fs_class(rght_gs, inner_interface.rght_state.vel_x)

                    #ghost_interface.lft_state = lft_fs
                    #ghost_interface.rght_state = rght_fs

                    ghost_interface.flux_flag = False

                    gc_layer.interface_array[interface] = ghost_interface
                
            self.mesh = JointBlock(mesh_object_1 = self.mesh, mesh_object_2 = gc_layer, \
                                block_1_interface_id_being_replaced = interface_id, \
                                block_2_interface_id_being_replaced = 0, \
                                new_interface = deepcopy(self.mesh.interface_array[interface_id]), \
                                adding_ghost_cells_bool = True)
            self.mesh.boundary_conditions[bc_ind].append(cell_pairs)
