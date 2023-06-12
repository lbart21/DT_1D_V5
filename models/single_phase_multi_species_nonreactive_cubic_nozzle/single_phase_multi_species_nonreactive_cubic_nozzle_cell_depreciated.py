"""
Function:
Author: Luke Bartholomew
Edits:
"""
import sys
from Algorithms.DT_1D_V5.models.cell_methods.single_phase_multi_species_decoding \
        import multi_species_decode_to_primative_properties
from Algorithms.DT_1D_V5.models.cell_methods.single_phase_multi_species_encode_cqs \
        import encode_multi_species_cqs

class SinglePhaseMultiSpeciesNonreactiveCubicNozzleCell():
    __slots__ = ["cell_id", "label", "geo", "flow_state", "phase", "interior_cell_flag", "cqs"]
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
        pass
    
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        pass