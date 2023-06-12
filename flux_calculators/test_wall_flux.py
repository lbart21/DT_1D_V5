from Algorithms.DT_1D_V5.flux_calculators.fluid_fluxes import AUSMPlusupORIGINAL, ausm_plus_up_original
from gdtk.gas import GasModel, GasState
from Algorithms.DT_1D_V5.fluid_models.flow_state import FlowState
import time


class call_class():
    def __init__(self, fs_L, fs_R) -> None:
        self.fs_L = fs_L
        self.fs_R = fs_R

    def do_thing(self):
        fluxes = AUSMPlusupORIGINAL(self.fs_L, self.fs_R, False, False)
        
class call_function():
    def __init__(self, fs_L, fs_R) -> None:
        self.fs_L = fs_L
        self.fs_R = fs_R
    
    def do_thing(self):
        fluxes = ausm_plus_up_original(self.fs_L, self.fs_R, False, False)
"""
gm_L = GasModel("ideal-air-gas-model.lua")
gm_R = GasModel("ideal-air-gas-model.lua")
#print(numba.typeof(gm_L))
gs_L = GasState(gm_L)
gs_R = GasState(gm_R)

p = 1e5
T = 300.0

gs_L.p = p
gs_L.T = T
gs_L.update_thermo_from_pT()

gs_R.p = p
gs_R.T = T
gs_R.update_thermo_from_pT()

vel_x = 100.0
vel_x_L = vel_x
vel_x_R = -vel_x

fluxes = ausm_plus_up_original(p_L = gs_L.p, p_R = gs_R.p, rho_L = gs_L.rho, rho_R = gs_R.rho, \
                                a_L = gs_L.a, a_L_for_Ma = gs_L.a, a_R = gs_R.a, a_R_for_Ma = gs_R.a, \
                                u_L = gs_L.u, u_R = gs_R.u, vel_x_L = vel_x_L, vel_x_R = vel_x_R, \
                                multi_species_flux = False, multi_phase = False)
"""
n_tests = 10000
class_time = 0.0
function_time = 0.0

class_data = [None] * n_tests
function_data = [None] * n_tests

p = 1e5
T = 300.0

for i in range(n_tests):
    gm_L = GasModel("ideal-air-gas-model.lua")
    gm_R = GasModel("ideal-air-gas-model.lua")

    gs_L = GasState(gm_L)
    gs_R = GasState(gm_R)

    p += 1e3
    T += 10.0

    gs_L.p = p
    gs_L.T = T
    gs_L.update_thermo_from_pT()

    gs_R.p = p
    gs_R.T = T
    gs_R.update_thermo_from_pT()

    vel_x = 100.0
    fs_L = FlowState(state = gs_L, vel_x = vel_x)
    fs_R = FlowState(state = gs_R, vel_x = -vel_x)

    class_data[i] = call_class(fs_L, fs_R)
    function_data[i] = call_function(fs_L, fs_R)

for test in range(n_tests):

    class_time_start = time.perf_counter()
    class_data[i].do_thing()
    class_time_end = time.perf_counter()
    class_time += (class_time_end - class_time_start)

    function_time_start = time.perf_counter()
    function_data[i].do_thing()
    function_time_end = time.perf_counter()
    function_time += (function_time_end - function_time_start)

print("Average class time:", class_time / n_tests)
print("Average function time:", function_time / n_tests)


