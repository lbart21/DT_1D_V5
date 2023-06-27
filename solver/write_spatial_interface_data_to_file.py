"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
def write_spatial_interface_data_to_file(interface_array, sim_number, time, simulation_description, \
                                         interface_property_variables):
    cwd = os.getcwd()
    file_name = "Sim" + str(sim_number) + "SpatialInterfaceDataAt" + str(format(time, ".9f")) + ".txt"
    with open(cwd + "/data/" + file_name, "w") as file:
        file.write("Sim: " + str(sim_number) + "\n")
        file.write("Time: " + str(time) + "\n")
        file.write("Description: " + simulation_description + "\n")
        n_interfaces = len(interface_array)
        first_interface_in_component = True
        for interface_idx in range(n_interfaces):
            interface_data = {}
            if "mass_flux" in interface_property_variables:
                interface_data["mass_flux"] = interface_array[interface_idx].boundary_fluxes["mass"]
            if "xMom_flux" in interface_property_variables:
                interface_data["xMom_flux"] = interface_array[interface_idx].boundary_fluxes["xMom"]
            if "energy_flux" in interface_property_variables:
                interface_data["energy_flux"] = interface_array[interface_idx].boundary_fluxes["energy"]
            if "Ma" in interface_property_variables:
                interface_data["Ma"] = interface_array[interface_idx].boundary_fluxes["Ma"]
            if "vel_x" in interface_property_variables:
                interface_data["vel_x"] = interface_array[interface_idx].boundary_fluxes["vel_x"]
            if "p" in interface_property_variables:
                interface_data["p"] = interface_array[interface_idx].boundary_fluxes["p"]
            if "spcs_mass" in interface_property_variables:
                spcs_names = interface_array[interface_idx].lft_state.fluid_state.gmodel.species_names
                spcs_names = ["massf_" + name for name in spcs_names]
                for ind, name in enumerate(spcs_names):
                    interface_data[name] = interface_array[interface_idx].boundary_fluxes["massf"][ind]
            
            joint_interface_data = {**interface_array[interface_idx].geo, **interface_data}
            variable_names = list(joint_interface_data.keys())
            if  first_interface_in_component:
                file.write("Variables: " + str(len(variable_names)) + "\n")
                file.write(" ".join(variable_names) + "\n")
                first_interface_in_component = False
            for name in variable_names:
                file.write(str(format(joint_interface_data[name], ".9f")) + " ")
            file.write("\n")
