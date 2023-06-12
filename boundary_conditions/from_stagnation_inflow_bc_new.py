"""
Function:
Author: Luke Bartholomew
Edits:
"""
def from_stagnation_inflow_bc(mesh, bc_ind):
    # bc = [interface_id, n_layers, str(West or East), ["FromStagnationInFlow_BC", gs], [(gc_id1, intc_id1), (gc_id2, intc_id2), ...]]
    boundary_edge = mesh.boundary_conditions[bc_ind][2]
    gs = mesh.boundary_conditions[bc_ind][3][1]
    interior_cell_id = mesh.boundary_conditions[bc_ind][4][0][1]
    if boundary_edge == "West":
        bulk_speed = max(0.0, mesh.cell_array[interior_cell_id].flow_state.vel_x)
        
    elif boundary_edge == "East":
        bulk_speed = min(0.0, mesh.cell_array[interior_cell_id].flow_state.vel_x)
        
    stag_enthalpy = gs.enthalpy
    stag_entropy = gs.entropy
    static_enthalpy = stag_enthalpy - 0.5 * bulk_speed ** 2.0
    gs.update_thermo_from_hs(h = static_enthalpy, s = stag_entropy)

    for gc_id, _ in mesh.boundary_conditions[bc_ind][4]:
        mesh.cell_array[gc_id].flow_state.fluid_state.copy_values(gs)
        mesh.cell_array[gc_id].flow_state.vel_x = bulk_speed
        
    gs.update_thermo_from_hs(h = stag_enthalpy, s = stag_entropy)
    return mesh