"""
Function:
Author: Luke Bartholomew
Edits:
"""

def mdot_from_stagnation_inflow_bc(mesh, bc_ind, t_current):
    relaxation_factor = 0.1
    # bc = [interface_id, n_layers, str(West or East), 
    #           ["FromStagnationWithMassFlowRateInFlow_BC", gs, mass_flux, [refresh_period, refresh_time, corrector_flag]], 
    #           [(gc_id1, intc_id1), (gc_id2, intc_id2), ...]]
    boundary_edge = mesh.boundary_conditions[bc_ind][2]
    gs = mesh.boundary_conditions[bc_ind][3][1]
    mass_flux = mesh.boundary_conditions[bc_ind][3][2]
    [refresh_period, refresh_time, corrector_flag] = mesh.boundary_conditions[bc_ind][3][3]

    interior_cell_id = mesh.boundary_conditions[bc_ind][4][0][1]

    p_stag = gs.p
    T_stag = gs.T

    if mass_flux > 0.0:
        if (not corrector_flag) and (t_current > refresh_time):
            print("Switching on corrector step for stagnation boundary condition.")
            corrector_flag = True
            mesh.boundary_conditions[bc_ind][3][3][2] = True

        if corrector_flag:
            p0_min = 0.1 * p_stag
            p0_max = 10.0 * p_stag

            A_boundary = mesh.interface_array[mesh.boundary_conditions[bc_ind][0]].geo["A"]

            rho_boundary = mesh.cell_array[interior_cell_id].flow_state.fluid_state.rho
            rho_v_boundary = rho_boundary * mesh.cell_array[interior_cell_id].flow_state.vel_x
            p_boundary = mesh.cell_array[interior_cell_id].flow_state.fluid_state.p

            dp_over_p = 0.5 * relaxation_factor / rho_boundary * ( (mass_flux / A_boundary) ** 2.0 \
                                                                    - rho_v_boundary * abs(rho_v_boundary) ) / p_boundary
            
            p_stag = (1.0 + dp_over_p) * p_stag
            p_stag = min(max(p_stag, p0_min), p0_max)

            if abs(dp_over_p) < 0.00001: # Allow min 0.1% change
                print("Switching off corrector step in stagnation boundary condition.")
                mesh.boundary_conditions[bc_ind][3][3][2] = False # Turn off corrector step for now
                mesh.boundary_conditions[bc_ind][3][3][1] = t_current + refresh_period # Set new refresh time
                #mesh.boundary_conditions[bc_ind][3][2] = 0.0

            gs.p = p_stag
            gs.update_thermo_from_pT()
    
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
    
    gs.p = p_stag
    gs.T = T_stag
    gs.update_thermo_from_pT()

    return mesh
