"""
Function:
Author: Luke Bartholomew
Edits:
"""
import math as m

from Algorithms.DT_1D_V5.boundary_conditions.form_boundary_condition_information_new import apply_boundary_conditions

class Integrate():
    def __init__(self, mesh, cfl_flag, t_current, current_step) -> None:
        self.dt_total = 1_000_000
        
        ### Add boundary conditions
        self.mesh = apply_boundary_conditions(mesh = mesh)

        ### Find out inviscid flux dt from CFL
        if not cfl_flag[0]: #Using fixed dt, cfl_flag has form [False, float dt_fixed]
            self.dt_total = cfl_flag[1]
        else:
            if cfl_flag[1] == "no_ramping":
                # Using fixed cfl
                # cfl_flag has form [True, "no_ramping", float cfl_fixed]
                min_dx_over_u = 1e6
                for cell in self.mesh.cell_array:
                    min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                self.dt_total = cfl_flag[2] * min_dx_over_u
            elif cfl_flag[1] == "cfl_ramp_from_tCurrent": 
                # Ramping cfl based on tCurrent
                # cfl_flag has form [True, "cfl_ramp_from_tCurrent", [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
                # If we are outside of t_n, use cfl_n. Else, find first instance tCurrent < t_i, then use cfl_i
                if t_current > cfl_flag[2][-1][1]:
                    cfl_max = cfl_flag[2][-1][0]
                    #print("cfl value used: ", cfl_max)
                    min_dx_over_u = 1e6
                    for cell in self.mesh.cell_array:
                        min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                    self.dt_total = cfl_max * min_dx_over_u
                else:
                    for cfl_pair in range(len(cfl_flag[2])):
                        if t_current > cfl_flag[2][cfl_pair][1]:
                            pass
                        else:
                            cfl_current = cfl_flag[2][cfl_pair][0]
                            #print("cfl value used: ", cfl_current)
                            min_dx_over_u = 1e6
                            for cell in self.mesh.cell_array:
                                min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                            self.dt_total = cfl_current * min_dx_over_u
                            break
                
            elif cfl_flag[1] == "cfl_ramp_from_step":
                # Ramping cfl based on currentStep
                # cfl_flag has form [True, "cfl_ramp_from_step", [[cfl_1, step_1], [cfl_2, step_2], ..., [cfl_n, step_n]]]
                # If we are outside of step_n, use cfl_n. Else, find first instance currentStep < step_i, then use cfl_i
                if current_step > cfl_flag[2][-1][1]: 
                    cfl_max = cfl_flag[2][-1][0]
                    #print("step: ", currentStep, " cfl used: ", cfl_max)
                    min_dx_over_u = 1e6
                    for cell in self.mesh.cell_array:
                        min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                    self.dt_total = cfl_max * min_dx_over_u
                else:
                    for cfl_pair in range(len(cfl_flag[2])):
                        if current_step > cfl_flag[2][cfl_pair][1]:
                            pass
                        
                        else:
                            cfl_current = cfl_flag[2][cfl_pair][0]
                            #print("step: ", currentStep, " cfl used: ", cfl_current)
                            min_dx_over_u = 1e6
                            for cell in self.mesh.cell_array:
                                min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                            self.dt_total = cfl_current * min_dx_over_u
                            break
                
            elif cfl_flag[1] == "dt_ramp_then_cfl_ramp_from_step":
                # Ramping dt initially, then use cfl based on current step number
                # cfl_flag has the form [True, "dt_ramp_then_cfl_ramp_from_step", [dt_init, scale], [[cfl_1, step_1], [cfl_2, step_2], ..., [cfl_n, step_n]]]
                # Starting at dt = dt_init, ramp this every time step by scale until the time step is greater than what comes from
                # cfl_1. Then ramp cfl similarly to "cfl_ramp_from_step" method.
                [dt_init, scale] = cfl_flag[2]
                min_dx_over_u = 1e6
                for cell in self.mesh.cell_array:
                    min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                min_dt_from_cfl = cfl_flag[3][0][0] * min_dx_over_u
                start_step = m.ceil( m.log(min_dt_from_cfl / dt_init, scale) )
                if current_step < start_step:
                    dt_ramped = dt_init * scale ** current_step
                    self.dt_total = dt_ramped
                
                else:
                    if current_step > cfl_flag[3][-1][1]:
                        min_dx_over_u = 1e6
                        for cell in self.mesh.cell_array:
                            min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                        self.dt_total = cfl_flag[3][-1][0] * min_dx_over_u
                        #print("step = ", current_step, " final cfl criteria used, cfl = ", cfl_flag[3][-1][0])
                    else:
                        for cfl_pair in range(len(cfl_flag[3])):
                            if current_step > cfl_flag[3][cfl_pair][1]:
                                pass
                            else:
                                cfl_current = cfl_flag[3][cfl_pair][0]
                                #print("step = ", current_step, " Intermediate cfl value used = ", cfl_current)
                                min_dx_over_u = 1e6
                                for cell in self.mesh.cell_array:
                                    min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                                self.dt_total = cfl_current * min_dx_over_u
                                break

            elif cfl_flag[1] == "dt_ramp_then_cfl_ramp_from_tCurrent":
                # Ramping dt initially, then use cfl based on current time.
                # cfl_flag has the form [True, "dt_ramp_then_cfl_ramp_from_tCurrent", [dt_init, scale], [[cfl_1, t_1], [cfl_2, t_2], ..., [cfl_n, t_n]]]
                # Starting at dt = dt_init, ramp this every time step by scale until the time step is greater than what comes from
                # cfl_1. Then ramp cfl similarly to "cfl_ramp_from_tCurrent".
                [dt_init, scale] = cfl_flag[2]
                min_dx_over_u = 1e6
                for cell in self.mesh.cell_array:
                    min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                min_dt_from_cfl = cfl_flag[3][0][0] * min_dx_over_u
                start_step = m.ceil( m.log(min_dt_from_cfl / dt_init, scale) )
                if current_step < start_step:
                    dt_ramped = dt_init * scale ** current_step
                    self.dt_total = dt_ramped
                
                else:
                    if t_current > cfl_flag[3][-1][1]:
                        min_dx_over_u = 1e6
                        for cell in self.mesh.cell_array:
                            min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                        self.dt_total = cfl_flag[3][-1][0] * min_dx_over_u
                    else:
                        for cfl_pair in range(len(cfl_flag[3])):
                            if t_current > cfl_flag[3][cfl_pair][1]:
                                pass
                            else:
                                cfl_current = cfl_flag[3][cfl_pair][0]
                                min_dx_over_u = 1e6 
                                for cell in self.mesh.cell_array:
                                    min_dx_over_u = min(min_dx_over_u, cell.max_allowable_dt())
                                self.dt_total = cfl_current * min_dx_over_u
                                break
                            
            else:
                print("Invalid CFL criteria given.")

        for interface in self.mesh.interface_array:
            ### Do all the interface specific methods that are defined within each interface (reconstruction, fluxes, perform conserved property updates etc etc)
            interface.complete_interface_methods(cell_array = self.mesh.cell_array, \
                                                interface_array = self.mesh.interface_array, \
                                                map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx, \
                                                map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                                map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                                map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                                dt_inv = self.dt_total)

        for ind, cell in enumerate(self.mesh.cell_array):
            cell.decode_to_primative_properties_after_interface_methods()
        
        ### Perform iteration over cells and update conserved properties
        for cell in self.mesh.cell_array:
            cell.complete_cell_methods(cell_array = self.mesh.cell_array, \
                                        interface_array = self.mesh.interface_array, \
                                        map_cell_id_to_west_interface_idx = self.mesh.map_cell_id_to_west_interface_idx, \
                                        map_cell_id_to_east_interface_idx = self.mesh.map_cell_id_to_east_interface_idx, \
                                        map_interface_id_to_west_cell_idx = self.mesh.map_interface_id_to_west_cell_idx, \
                                        map_interface_id_to_east_cell_idx = self.mesh.map_interface_id_to_east_cell_idx, \
                                        dt_inv = self.dt_total)
        
        for cell in self.mesh.cell_array:
            cell.decode_to_primative_properties_after_cell_methods()

            