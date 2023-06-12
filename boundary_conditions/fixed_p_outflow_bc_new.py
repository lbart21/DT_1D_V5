"""
Function:
Author: Luke Bartholomew
Edits:
"""
def fixed_p_outflow_bc(mesh, bc_ind):
    # bc = [interface_id, n_layers, str(West or East), ["FixedPOutFlow_BC", gs], [(gc_id1, intc_id1), (gc_id2, intc_id2), ...]]
    gs = mesh.boundary_conditions[bc_ind][3][1]
    for gc_id, intc_id in mesh.boundary_conditions[bc_ind][4]:
        mesh.cell_array[gc_id].flow_state.fluid_state.copy_values(mesh.cell_array[intc_id].flow_state.fluid_state)
        mesh.cell_array[gc_id].flow_state.fluid_state.p = gs.p
        mesh.cell_array[gc_id].flow_state.fluid_state.update_thermo_from_pT()
        mesh.cell_array[gc_id].flow_state.vel_x = mesh.cell_array[intc_id].flow_state.vel_x
    return mesh