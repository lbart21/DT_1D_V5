"""
Function:
Author: Luke Bartholomew
Edits:
"""
import sys
from numpy import array
from Algorithms.DT_1D_V5.models.cell_methods.single_phase_multi_species_decoding \
        import multi_species_decode_to_primative_properties
from Algorithms.DT_1D_V5.models.cell_methods.single_phase_multi_species_encode_cqs \
        import encode_multi_species_cqs
from Algorithms.DT_1D_V5.models.cell_methods.bulk_reaction_fast_chemistry_method \
        import bulk_reaction_fast_chemsitry
from Algorithms.DT_1D_V5.models.cell_methods.goal_massf_solver \
        import goal_massf_solver, form_goal_massf_solver_matrices

class SinglePhaseMultiSpeciesNonuniformBulkReactionFastChemistryCubicNozzleCell():
    __slots__ = ["cell_id", "label", "geo", "flow_state", "phase", "interior_cell_flag", \
                    "cqs", "bulk_reaction_parameters", "massf_solver_matrices_been_formed", \
                    "LHS_A", "RHS_A", "source"]
    def __init__(self, cell_id, label) -> None:
        self.cell_id = cell_id
        self.label = label
        self.geo = {}
        self.flow_state = None
        self.phase = "Single"
        self.interior_cell_flag = True

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0,
            "spcs_mass" : []
        }

        self.bulk_reaction_parameters = [] #Will be [x, goal_massf_dict, reaction]
        self.massf_solver_matrices_been_formed = False
        self.LHS_A = None
        self.RHS_A = None
        self.source = None
    
    def fill_geometry(self, geometry):
        self.geo = geometry
    
    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho, u, massf, vel_x = multi_species_decode_to_primative_properties(cqs = self.cqs)
            if rho < 0.0:
                print("Failed to decode to rho, u, massf, vel_x due to negative mass.")
                print("ID:", self.cell_id)
                print("cqs:", self.cqs)
                print("rho, u, massf:", rho, u, massf)
                sys.exit(1)

            self.flow_state.fluid_state.rho = rho
            self.flow_state.fluid_state.u = u
            self.flow_state.fluid_state.massf = massf
            self.flow_state.vel_x = vel_x
            try:
                self.flow_state.fluid_state.update_thermo_from_rhou()
            except:
                print("Failed to update from rho, u. ID:", self.cell_id)
                print("cqs:", self.cqs)
                print("rho, u, massf:", rho, u, massf)
                sys.exit(1)
    
    def initialise_conserved_quantities(self):
        self.cqs = encode_multi_species_cqs(flow_state = self.flow_state)

    def max_allowable_dt(self):
        return self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)
    
    def complete_cell_methods(self, **kwargs):
        if not self.massf_solver_matrices_been_formed:
            self.LHS_A, self.RHS_A, self.source = form_goal_massf_solver_matrices(\
                                    bulk_reaction_parameters = self.bulk_reaction_parameters, \
                                    molar_masses = self.flow_state.fluid_state.gmodel.mol_masses, \
                                    species_names = self.flow_state.fluid_state.gmodel.species_names)
            self.massf_solver_matrices_been_formed = True
            
        goal_massf = goal_massf_solver(massf_current = self.flow_state.fluid_state.massf, \
                                       LHS_A = self.LHS_A, RHS_A = self.RHS_A, source = self.source)
        new_massf = bulk_reaction_fast_chemsitry(massf_current = self.flow_state.fluid_state.massf, \
                                                x = self.bulk_reaction_parameters[0], \
                                                massf_goal = goal_massf)
        self.flow_state.fluid_state.massf = new_massf
        self.reinitialise_massf_conserved_quantity()
    
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()

    def reinitialise_massf_conserved_quantity(self):
        self.cqs["spcs_mass"] = (self.flow_state.fluid_state.rho * array(self.flow_state.fluid_state.massf)).tolist()