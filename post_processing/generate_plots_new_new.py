"""
Function:
Author: Luke Bartholomew
Edits:
"""
import os
import re
import time
import numpy as np
import pandas as pd

from Algorithms.DT_1D_V5.post_processing.process_spatial_cell_data import ProcessSpatialCellData
from Algorithms.DT_1D_V5.post_processing.process_spatial_interface_data \
        import ProcessSpatialInterfaceData
from Algorithms.DT_1D_V5.post_processing.process_transient_cell_data import ProcessTransientCellData
from Algorithms.DT_1D_V5.post_processing.process_transient_interface_data \
            import ProcessTransientInterfaceData
from Algorithms.DT_1D_V5.post_processing.process_eilmer_data import ProcessEilmerData
from Algorithms.DT_1D_V5.post_processing.SI_units_dictionary import SI_UNITS
from Algorithms.DT_1D_V5.post_processing.symbols import SYMBOLS

from Algorithms.DT_0D_V2.post_processing.interface_data_file_to_object \
            import FormInterfaceDataFromFile as interface_data_object_0d

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation
####################################################################################################
def generate_single_spatial_cell_data_plots(data_files, plot_vars, label):
    """
    - Plots with a single y-axis
    - Options for plot_vars elements:
        - list of mole or mass fraction names
        - Individual properties (p, T, vel_x, stagnation condition etc, as well as specific mass 
            or mole fractions)
        - "massf", "molef" and "conc" will plot all available mass fractions, mole fractions 
            and concentrations.
    - data_files: list(str)
    - plot_vars: list(str, list(str))
    - label: str
    """
    data_object = None
    for component_data_file in data_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
        if data_object is None:
            data_object = data.component_data
            t_final = data.t_final
            sim_number = data.sim_number
        else:
            data_object = pd.concat([data_object, data.component_data], axis = 0)
    data_object = data_object.sort_values(by = ["pos_x"])

    for var in plot_vars:
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)
        formatted_title_time = f"{t_final / 1e-6:.3f}"
        formatted_file_name_time = f"{t_final:.9f}"

        #formatted_title_time = '{:.3f}'.format(t_final / 1e-6)
        #formatted_file_name_time = '{:.9f}'.format(t_final)
        
        if var == "massf": # Plot all mass fractions
            column_names = list(data_object.columns)
            massf_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "massf"]
            species_names = [species.split("_")[1] for species in massf_species_names] 
            for index, species in enumerate(species_names):
                split_name = re.findall('[a-zA-Z][^A-Z]*', species)
                #split_name = re.findall('(\d+|[A-Za-z]+)', species)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                ax.scatter(data_object["pos_x"], \
                            data_object[massf_species_names[index]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
            ax.legend()
            ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
            filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
        
        elif var == "molef":
            column_names = list(data_object.columns)
            molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
            species_names = [species.split("_")[1] for species in molef_species_names]
            for index, species in enumerate(species_names):
                #split_name = re.findall('(\d+|[A-Za-z]+)', species)
                split_name = re.findall('[a-zA-Z][^A-Z]*', species)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                ax.scatter(data_object["pos_x"], \
                            data_object[molef_species_names[index]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
            ax.legend()
            ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
            filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
        
        elif var == "conc":
            column_names = list(data_object.columns)
            conc_species_names = [species for species in column_names \
                                            if len(species) > 4 and species[:4] == "conc"]
            species_names = [species.split("_")[1] for species in conc_species_names]
            for index, species in enumerate(species_names):
                if val.isnumeric():
                    split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                ax.scatter(data_object["pos_x"], \
                            data_object[conc_species_names[index]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
            ax.legend()
            ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
            filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
        
        elif len(var) > 5 and var[:5] == "massf": # Plot specific mass fractions
            species_name = var[6:] # Trim off starting massf_
            #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
            for ind, val in enumerate(split_name):
                if val.isnumeric():
                    split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'
            ax.scatter(data_object["pos_x"], data_object[var], marker = '.')
            ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
            filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

        elif len(var) > 5 and var[:5] == "molef": # Plot specific mass fractions
            species_name = var[6:] # Trim off starting molef_
            #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
            for ind, val in enumerate(split_name):
                if val.isnumeric():
                    split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                    + ' ' + r'$Fraction$'
            ax.scatter(data_object["pos_x"], data_object[var], marker = '.')
            ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
            ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
            filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

        elif isinstance(var, list):
            if var[0][:5] == "massf":
                for species in var:
                    species_name = species[6:] # Trim off starting massf_
                    #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[species] = r'$-$'
                    SYMBOLS[species] = r"$" + formatted_species_name + "$" + \
                                    ' ' + r'$Mass$' + ' ' + r'$Fraction$'
                    ax.scatter(data_object["pos_x"], data_object[species], marker = '.', \
                                label = r"$" + formatted_species_name + "$")
                ax.legend()
                ax.set_ylabel(SYMBOLS["massf"] + " (" + SI_UNITS["massf"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Multiple Mass Fractions at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
                filename = "Sim " + str(sim_number) + ' ' + label + ' ' + ' '.join(var) \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

            elif var[0][:5] == "molef":
                for species in var:
                    species_name = species[6:] # Trim off starting molef_
                    #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
                    split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[species] = r'$-$'
                    SYMBOLS[species] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'
                    ax.scatter(data_object["pos_x"], data_object[species], marker = '.', \
                                label = r"$" + formatted_species_name + "$")
                ax.legend()
                ax.set_ylabel(SYMBOLS["molef"] + " (" + SI_UNITS["molef"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Multiple Mole Fractions at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
                filename = "Sim " + str(sim_number) + ' ' + label + ' ' + ' '.join(var) \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

            elif var[0][:4] == "conc":
                pass

        
        else: # All other properties
            ax.scatter(data_object["pos_x"], data_object[var], marker = '.')
            ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
            ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                + formatted_title_time + r'$\mu$' + "s")

            filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
        
        ax.set_xlabel("Position (m)")
        #ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    #rotation = "horizontal", ha = "right")
        #ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            #+ formatted_title_time + r'$\mu$' + "s")

        #filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    #+ " distribution at t = " + formatted_file_name_time + ".jpg"

        ax.grid()
        current_dir = os.getcwd()
        fig.savefig(current_dir + "/plots/" + filename, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_spatial_cell_data_waterfall_plots(data_files, plot_vars, label):
    """
    - data_files: list(list(str))
    - plot_vars: 
    - label: str
    """
    data_objects = [None] * len(data_files)
    time_list = [None] * len(data_files)

    for ind, snapshot in enumerate(data_files):
        data_made = False
        for file in snapshot:
            data = ProcessSpatialCellData(data_file_name = file)
            if not data_made:
                data_objects[ind] = data.component_data
                time_list[ind] = data.t_final
                sim_number = data.sim_number
                data_made = True
            else:
                data_objects[ind] = pd.concat([data_objects[ind], data.component_data], \
                                                axis = 0)
        
        data_objects[ind] = data_objects[ind].sort_values(by = ["pos_x"])
    
    for var in plot_vars:
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)

        if len(var) > 5 and var[:5] == "massf":
            species_name = var[6:] # Trim off starting massf_
            #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'

        elif len(var) > 5 and var[:5] == "molef":
            species_name = var[6:] # Trim off starting molef_
            #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'
        
        for ind, t_snapshot in enumerate(time_list):
            #formatted_time = '{:.3f}'.format(t_snapshot / 1e-6)
            formatted_time = f'{t_snapshot / 1e-6:.3f}'
            """
            if var == "massf":
                column_names = list(data_objects[ind].columns)
                massf_species_names = [species for species in column_names \
                                            if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names] 
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    ax.scatter(data_objects[ind]["pos_x"], \
                                data_objects[ind][massf_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$" \
                                            + " at t = " + formatted_time + r'$\mu$' + "s")
                ax.legend()
            
            elif var == "molef":
                column_names = list(data_objects[ind].columns)
                molef_species_names = [species for species in column_names \
                                            if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names] 
                for index, species in enumerate(species_names):
                    split_name = re.findall('(\d+|[A-Za-z]+)', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    ax.scatter(data_objects[ind]["pos_x"], \
                                data_objects[ind][molef_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$" \
                                            + " at t = " + formatted_time + r'$\mu$' + "s")
                ax.legend()
            """
            #else:
            ax.scatter(data_objects[ind]["pos_x"], data_objects[ind][var], \
                    label = "Distribution at t = " + formatted_time + r'$\mu$' + "s", \
                    marker = ".")
        
        ax.set_title("Distribution of " + SYMBOLS[var] + " at Various Times")
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
        ax.set_xlabel("Position (m)")
        ax.legend()
        ax.grid()
        filename = "Sim " + str(sim_number) + ' ' + label \
                        + ' ' + var + " distribution at multiple times.jpg"
        current_dir = os.getcwd()
        fig.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
####################################################################################################
def generate_transient_gross_thrust_plot(interface_data_file, max_line = None):
    """
    - interface_data_file: str
    - max_line: int
    """
    interface_data_total = ProcessTransientInterfaceData(data_file_name = interface_data_file, \
                                                    max_line = max_line)
    interface_data = interface_data_total.interface_data
    sim_number = interface_data_total.sim_number
    m_dot = interface_data["mass_flux"]
    p_exit = interface_data["p"]
    vel_x_e = interface_data["vel_x"]
    A_e = interface_data["A"]
    interface_data["Thrust"] = m_dot * vel_x_e * A_e + p_exit * A_e

    fig, ax = plt.subplots()
    fig.set_size_inches(15, 5)

    ax.set_title("Transient Gross Thrust Profile")
    ax.set_ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
    ax.set_xlabel("Time (ms)")
    ax.scatter(interface_data["time"] * 1e3, interface_data["Thrust"], marker = '.')
    ax.grid()

    filename = "Sim"+ str(sim_number) + "GrossThrustProfile.jpg"
    current_dir = os.getcwd()
    fig.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
    plt.close()
####################################################################################################
def generate_net_thrust_plot(inlet_interface_data_file, outlet_interface_data_file, \
                                max_line = None):
    """
    - inlet_interface_data_file: str
    - outlet_interface_data_file: str
    - max_line: int
    """
    inlet_interface_data_total = ProcessTransientInterfaceData(\
                                    data_file_name = inlet_interface_data_file, \
                                    max_line = max_line)
    outlet_interface_data_total = ProcessTransientInterfaceData(\
                                    data_file_name = outlet_interface_data_file, \
                                    max_line = max_line)
    
    inlet_interface_data = inlet_interface_data_total.interface_data
    outlet_interface_data = outlet_interface_data_total.interface_data

    sim_number = outlet_interface_data_total.sim_number

    m_dot_exit = outlet_interface_data["mass_flux"]
    m_dot_inlet = inlet_interface_data["mass_flux"]

    vel_x_exit = outlet_interface_data["vel_x"]
    vel_x_inlet = inlet_interface_data["vel_x"]

    p_exit = outlet_interface_data["p"]
    p_inlet = inlet_interface_data["p"]

    A_exit = outlet_interface_data["A"]
    A_inlet = inlet_interface_data["A"]

    outlet_interface_data["Thrust"] = m_dot_exit * vel_x_exit * A_exit \
                                        - m_dot_inlet * vel_x_inlet * A_inlet \
                                        + p_exit * A_exit - p_inlet * A_inlet
    
    fig, ax = plt.subplots()
    fig.set_size_inches(15, 5)

    ax.set_title("Transient Net Thrust Profile")
    ax.set_ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
    ax.set_xlabel("Time (ms)")
    ax.scatter(outlet_interface_data["time"] * 1e3, outlet_interface_data["Thrust"], marker = '.')
    ax.grid()

    filename = "Sim"+ str(sim_number) + "NetThrustProfile.jpg"
    current_dir = os.getcwd()
    fig.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
    plt.close()
####################################################################################################
def generate_transient_cell_property_plots(cell_data_file, plot_vars, max_line = None):
    """
    - cell_data_file: str
    - plot_vars: 
    - max_line: int
    """
    cell_data_total = ProcessTransientCellData(data_file_name = cell_data_file, \
                                           max_line = max_line)
    cell_data = cell_data_total.cell_data
    cell_id = cell_data_total.cell_id
    sim_number = cell_data_total.sim_number

    for var in plot_vars:
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)

        ax.set_title("Transient Development of " + SYMBOLS[var] + " at Cell " \
                            + str(cell_id))
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                            rotation = "horizontal", ha = "right")
        ax.set_xlabel("Time (ms)")

        if len(var) > 5 and var[:5] == "massf":
            species_name = var[6:] # Trim off starting massf_
            #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'
    
        elif len(var) > 5 and var[:5] == "molef":
            species_name = var[6:] # Trim off starting molef_
            #split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            split_name = re.findall('[a-zA-Z][^A-Z]*', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'

        if var == "massf":
            massf_names = [name for name in cell_data.columns if "massf" in name]
            for species in massf_names:
                name = species[6:]
                split_name = re.findall('[a-zA-Z][^A-Z]*', name)
                #split_name = re.findall('(\d+|[A-Za-z]+)', name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_name = ''.join(split_name)
                ax.scatter(cell_data["time"] * 1e3, cell_data[species], marker = '.', \
                           label = r"$" + formatted_name + "$")
            ax.legend()
        
        elif var == "molef":
            molef_names = [name for name in cell_data.columns if "molef" in name]
            for species in molef_names:
                name = species[6:]
                #split_name = re.findall('(\d+|[A-Za-z]+)', name)
                split_name = re.findall('[a-zA-Z][^A-Z]*', name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_name = ''.join(split_name)
                ax.scatter(cell_data["time"] * 1e3, cell_data[species], marker = '.', \
                            label = r"$" + formatted_name + "$")
            ax.legend()
        
        else:
            ax.scatter(cell_data["time"] * 1e3, cell_data[var], marker = '.')
        
        ax.grid()
        current_dir = os.getcwd()
        filename = "Sim " + str(sim_number) + " Transient Development of " + \
                            var + " at Cell " + str(cell_id) + ".jpg"
        fig.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
####################################################################################################
def generate_transient_interface_property_plots(interface_data_file, plot_vars, \
                                                    max_line = None):
    """
    - interface_data_file: str
    - plot_vars:
    - max_line: int
    """
    interface_data_total = ProcessTransientInterfaceData(data_file_name = interface_data_file, \
                                                            max_line = max_line)
    interface_data = interface_data_total.interface_data
    sim_number = interface_data_total.sim_number
    interface_id = interface_data_total.interface_id

    if "mass_flux" in plot_vars:
        interface_data["mass_flux"] *= interface_data["A"]
    if "energy_flux" in plot_vars:
        interface_data["energy_flux"] *= interface_data["A"]

    for var in plot_vars:
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)
        
        ax.scatter(interface_data["time"] * 1e3, interface_data[var], marker = '.')
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
        ax.set_title("Transient Development of " + SYMBOLS[var]+ \
                                " at Interface " + str(interface_id))
        ax.grid()
        filename = "Sim " + str(sim_number) + "Transient Development of " + var + " at Interface " \
                        + str(interface_id) + ".jpg"
        current_dir = os.getcwd()
        fig.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
####################################################################################################
def generate_thrust_contribution_plot_without_interpolation(spatial_cell_data_file, \
                                                            interface_geometry_data_file):
    """
    - spatial_cell_data_file: list(str)
    - interface_geometry_data_file: str
    """
    pass
