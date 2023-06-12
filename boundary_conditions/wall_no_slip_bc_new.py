"""
Function:
Author: Luke Bartholomew
Edits:
"""
def wall_no_slip_bc(mesh, bc_ind):
    for gc_id, intc_id in mesh.boundary_conditions[bc_ind][4]:
        mesh.cell_array[gc_id].flow_state.fluid_state.copy_values(mesh.cell_array[intc_id].flow_state.fluid_state)
        mesh.cell_array[gc_id].flow_state.vel_x = -1.0 * mesh.cell_array[intc_id].flow_state.vel_x
    return mesh