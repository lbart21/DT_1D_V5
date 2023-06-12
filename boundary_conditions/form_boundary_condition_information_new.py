"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DT_1D_V5.boundary_conditions.simple_outflow_bc_new import simple_outflow_bc
from Algorithms.DT_1D_V5.boundary_conditions.fixed_p_outflow_bc_new import fixed_p_outflow_bc
from Algorithms.DT_1D_V5.boundary_conditions.fixed_pt_outflow_bc_new import fixed_pt_outflow_bc
from Algorithms.DT_1D_V5.boundary_conditions.from_stagnation_inflow_bc_new import from_stagnation_inflow_bc
from Algorithms.DT_1D_V5.boundary_conditions.mdot_from_stagnation_inflow_bc_new import mdot_from_stagnation_inflow_bc
from Algorithms.DT_1D_V5.boundary_conditions.supersonic_inflow_bc_new import supersonic_inflow_bc
from Algorithms.DT_1D_V5.boundary_conditions.wall_no_slip_bc_new import wall_no_slip_bc
from Algorithms.DT_1D_V5.boundary_conditions.wall_with_slip_bc_new import wall_with_slip_bc

def apply_boundary_conditions(mesh):
    for bc_ind, bc in enumerate(mesh.boundary_conditions):
        if bc[3][0] == "SimpleOutFlow_BC":
            mesh = simple_outflow_bc(mesh = mesh, bc_ind = bc_ind)

        elif bc[3][0] == "SimpleExtrapolateOutFlow_BC":
            print("SimpleExtrapolateOutFlow_BC has not been implemented yet.")
            
        elif bc[3][0] == "WallNoSlip_BC":
            mesh = wall_no_slip_bc(mesh = mesh, bc_ind = bc_ind)

        elif bc[3][0] == "WallWithSlip_BC":
            mesh = wall_with_slip_bc(mesh = mesh, bc_ind = bc_ind)
            
        elif bc[3][0] == "SupersonicInFlow_BC":
            mesh = supersonic_inflow_bc(mesh = mesh, bc_ind = bc_ind)
            
        elif bc[3][0] == "FromStagnationInFlow_BC":
            mesh = from_stagnation_inflow_bc(mesh = mesh, bc_ind = bc_ind)
            
        elif bc[3][0] == "FromStagnationWithMassFlowRateInFlow_BC":
            mesh = mdot_from_stagnation_inflow_bc(mesh = mesh, bc_ind = bc_ind)
            
        elif bc[3][0] == "FixedPOutFlow_BC":
            mesh = fixed_p_outflow_bc(mesh = mesh, bc_ind = bc_ind)
            
        elif bc[3][0] == "FixedPTOutFlow_BC":
            mesh = fixed_pt_outflow_bc(mesh = mesh, bc_ind = bc_ind)

    return mesh
