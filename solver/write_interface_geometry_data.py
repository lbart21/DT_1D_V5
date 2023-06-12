import os
def write_interface_geometry_data(interface_array, sim_number):
    cwd = os.getcwd()
    filename = "Sim" + str(sim_number) + "InterfaceGeometryData.txt"
    with open(cwd + "/data/" + filename, "w") as file:
        file.write("Sim: " + str(sim_number) + "\n")
        variable_names = ["pos_x", "A"]
        file.write("Variables: " + str(len(variable_names)) + "\n")
        file.write(" ".join(variable_names) + "\n")
        for interface in interface_array:
            for variable in variable_names:
                file.write(str(format(interface.geo[variable], ".9f")) + " ")
            file.write("\n")

    file.close()            