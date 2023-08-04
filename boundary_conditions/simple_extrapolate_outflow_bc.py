"""
Function:
Author: Luke Bartholomew
Edits:
"""
def simple_extrapolate_outflow_bc(mesh, bc_ind):
    if len(mesh.boundary_conditions[bc_ind][4]) == 1:
        print("This boundary condition is not suitable for only 1 layer of ghost cells, resorting to zero-gradient copy boundary condition.")
        gc_id, intc_id = mesh.boundary_conditions[bc_ind][4][0]
        mesh.cell_array[gc_id].flow_state.fluid_state.copy_values(mesh.cell_array[intc_id].flow_state.fluid_state)
        mesh.cell_array[gc_id].flow_state.vel_x = mesh.cell_array[intc_id].flow_state.vel_x
    else:
        int_cell_1 = mesh.cell_array[mesh.boundary_conditions[bc_ind][4][0][1]]
        int_cell_2 = mesh.cell_array[mesh.boundary_conditions[bc_ind][4][1][1]]

        rho_1 = int_cell_1.flow_state.fluid_state.rho
        u_1 = int_cell_1.flow_state.fluid_state.u
        vel_x_1 = int_cell_1.flow_state.vel_x
        massf_1 = int_cell_1.flow_state.fluid_state.massf

        rho_2 = int_cell_2.flow_state.fluid_state.rho
        u_2 = int_cell_2.flow_state.fluid_state.u
        vel_x_2 = int_cell_2.flow_state.vel_x
        massf_2 = int_cell_2.flow_state.fluid_state.massf

        dx_base = 0.5 * (int_cell_1.geo["dx"] + int_cell_2.geo["dx"])

        for gc_pair_ind in range(len(mesh.boundary_conditions[bc_ind][4])):
            gc_id, intc_id = mesh.boundary_conditions[bc_ind][4][gc_pair_ind]
            gc_pos_x = 0.5 * int_cell_1.geo["dx"] + int_cell_2.geo["dx"] + 0.5 * mesh.cell_array[gc_id].geo["dx"]
            for prev_gc_cell_ind in range(gc_pair_ind):
                gc_id_prev, _ = mesh.boundary_conditions[bc_ind][4][prev_gc_cell_ind]
                gc_pos_x += mesh.cell_array[gc_id_prev].geo["dx"]
            rho_extrapolated = (rho_1 - rho_2) / dx_base * gc_pos_x + rho_2
            u_extrapolated = (u_1 - u_2) / dx_base * gc_pos_x + u_2
            vel_x_extrapolated = (vel_x_1 - vel_x_2) / dx_base * gc_pos_x + vel_x_2
            massf_extrapolated = massf_2.copy()
            for spcs_ind, _ in enumerate(massf_extrapolated):
                massf_extrapolated[spcs_ind] = (massf_1[spcs_ind] - massf_2[spcs_ind]) / dx_base * gc_pos_x + massf_2[spcs_ind]
            mesh.cell_array[gc_id].flow_state.fluid_state.rho = rho_extrapolated
            mesh.cell_array[gc_id].flow_state.fluid_state.u = u_extrapolated
            mesh.cell_array[gc_id].flow_state.fluid_state.massf = massf_extrapolated
            mesh.cell_array[gc_id].flow_state.fluid_state.update_thermo_from_rhou()
            mesh.cell_array[gc_id].flow_state.vel_x = vel_x_extrapolated
    return mesh