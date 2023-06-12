"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
def write_simulation_description(sim_number, simulation_description, boundary_conditions, cfl_flag, t_final):
    cwd = os.getcwd()
    file_name = "Sim" + str(sim_number) + "SimulationDescription.txt"
    file = open(cwd + "/data/" + file_name, "w")
    file.write("Sim: " + str(sim_number) + "\n")
    file.write("Description: " + simulation_description + "\n")
    file.write("Boundary conditions: {" + "\n")
    for bc in boundary_conditions:
        file.write("[" + str(bc[0]) + ", " + str(bc[1]) + ", " + str(bc[2]) + ", [")

        if bc[3][0] == "WallNoSlip_BC":
            file.write("'WallNoSlip_BC']]")

        elif bc[3][0] == "SupersonicInFlow_BC":
            fs = bc[3][1]
            p = fs.fluid_state.p
            T = fs.fluid_state.T
            massf = fs.fluid_state.massf
            vel_x = fs.vel_x

            file.write("'SupersonicInFlow_BC', p: %f, T: %f, vel_x: %.3f, massf: [" % (p, T, vel_x))
            for spcs in massf:
                file.write("%.9f, " % spcs)
            file.write("]]]")

        elif bc[3][0] == "FromStagnationInFlow_BC":
            gs = fs = bc[3][1]
            p = gs.p
            T = gs.T
            massf = gs.massf

            file.write("'FromStagnationInFlow_BC', p:, %f, T: %f, massf: [" % (p, T))
            for spcs in massf:
                file.write("%.9f, " % spcs)
            file.write("]]]")
            
        elif bc[3][0] == "FromStagnationWithMassFlowRateInFlow_BC":
            gs = fs = bc[3][1]
            mass_flux = bc[3][2]
            p = gs.p
            T = gs.T
            massf = gs.massf

            file.write("'FromStagnationInFlow_BC', p:, %f, T: %f, mass_flux: %f, massf: [" % (p, T, mass_flux))
            for spcs in massf:
                file.write("%.9f, " % spcs)
            file.write("]]]")
            
        elif bc[3][0] == "SimpleOutFlow_BC":
            file.write("'SimpleOutFlow_BC']]")
            
        elif bc[3][0] == "SimpleExtrapolateOutFlow_BC":
            file.write("'SimpleExtrapolateOutFlow_BC']]")
            
        elif bc[3][0] == "FixedPOutFlow_BC":
            gs = fs = bc[3][1]
            p = gs.p
            file.write("'FixedPOutFlow_BC', p: %f]]" % p)
            
        elif bc[3][0] == "FixedPTOutFlow_BC":
            gs = fs = bc[3][1]
            p = gs.p
            T = gs.T
            file.write("'FixedPTOutFlow_BC', p: %f, T: %f]]" % (p, T))
        file.write("}\n")

    file.write("cfl_flag: ")
    cfl_flag_string = convert_list_to_string(cfl_flag)
    file.write(cfl_flag_string + "\n")
    file.write("t_final: %f" % t_final)
    
    file.close()

def convert_list_to_string(list_object):

    string_list = [None] * len(list_object)
    for ind, element in enumerate(list_object):
        if type(element) == type([]):
            element = convert_list_to_string(element)
        string_list[ind] = str(element)
    
    string_list = "["+ ', '.join(string_list) + "]"

    return string_list
