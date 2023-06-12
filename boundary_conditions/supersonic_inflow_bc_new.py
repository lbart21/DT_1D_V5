"""
Function:
Author: Luke Bartholomew
Edits:
"""
def supersonic_inflow_bc(mesh, bc_ind):
    # bc = [interface_id, n_layers, str(West or East), ["SupersonicInFlow_BC", fs], [(gc_id1, intc_id1), (gc_id2, intc_id2), ...]]
    fs = mesh.boundary_conditions[bc_ind][3][1]
    for gc_id, _ in mesh.boundary_conditions[bc_ind][4]:
        mesh.cell_array[gc_id].flow_state.fluid_state.copy_values(fs.fluid_state)
        mesh.cell_array[gc_id].flow_state.vel_x = fs.vel_x

    return mesh