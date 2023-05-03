"""
Function:
Author: Luke Bartholomew
Edits:
"""
from Algorithms.DT_1D_V4.models.cell_methods.single_phase_single_species_decoding \
        import single_species_decode_to_primative_properties

from Algorithms.DT_1D_V4.models.cell_methods.single_phase_single_species_encode_cqs \
        import encode_single_species_cqs

class SinglePhaseMultiSpeciesNonReactivePipeWithHeatFluxCell():
    def __init__(self, cell_id, label, heat_flux) -> None:
        self.geo = {}
        self.flow_state = None
        self.cell_id = cell_id
        self.label = label
        self.phase = "Single"
        self.interior_cell_flag = True

        self.heat_flux = heat_flux

        self.cqs = {
            "mass"      : 0.0,
            "xMom"      : 0.0,
            "energy"    : 0.0
        }

        self.flow_state = None
    
    def fill_geometry(self, geometry):
        self.geo = geometry
    
    def decode_to_primative_properties(self):
        if self.interior_cell_flag:
            rho, u, vel_x = single_species_decode_to_primative_properties(cqs = self.cqs)

            self.flow_state.fluid_state.rho = rho
            self.flow_state.fluid_state.u = u

            self.flow_state.vel_x = vel_x

            self.flow_state.fluid_state.update_thermo_from_rhou()
    
    def initialise_conserved_quantities(self):
        self.cqs = encode_single_species_cqs(flow_state = self.flow_state)

    def max_allowable_dt(self, cfl):
        return cfl * self.geo["dx"] / (abs(self.flow_state.vel_x) + self.flow_state.fluid_state.a)

    def complete_cell_methods(self, **kwargs):
        dt_inv = kwargs["dt_inv"]
        self.cqs["energy"] += dt_inv * self.geo["A_s"] * self.heat_flux / self.geo["dV"]
    
    def decode_to_primative_properties_after_interface_methods(self):
        self.decode_to_primative_properties()

    def decode_to_primative_properties_after_cell_methods(self):
        self.decode_to_primative_properties()