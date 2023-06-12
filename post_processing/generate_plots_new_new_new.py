"""
Function:
Author: Luke Bartholomew
Edits:
"""

import os
import re
import numpy as np
import pandas as pd

from Algorithms.DT_1D_V5.post_processing.process_spatial_cell_data \
        import ProcessSpatialCellData
from Algorithms.DT_1D_V5.post_processing.process_spatial_interface_data \
        import ProcessSpatialInterfaceData
from Algorithms.DT_1D_V5.post_processing.process_transient_cell_data \
        import ProcessTransientCellData
from Algorithms.DT_1D_V5.post_processing.process_transient_interface_data \
        import ProcessTransientInterfaceData

from Algorithms.DT_1D_V5.post_processing.symbols import SYMBOLS
from Algorithms.DT_1D_V5.post_processing.SI_units_dictionary import SI_UNITS

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import animation

####################################################################################################
def generate_single_spatial_cell_data_plots(data_files, plot_vars, label):
    """
    - data_files: list(str)
    - plot_vars: list(list(str), str)
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

        if isinstance(var, list): #either list("massf_" or "molef_")
            ismassf = False
            ismolef = True

            if var[0][:6] == "massf_":
                ismassf = True
                ismolef = False

            for prop in var:
                spcs_name = prop[6:]
                split_name = re.findall('[0-9]+|[A-Za-z]+', spcs_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[prop] = r'$-$'
                if ismassf:
                    SYMBOLS[prop] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'
                elif ismolef:
                    SYMBOLS[prop] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'
                ax.scatter(data_object["pos_x"], data_object[prop], marker = '.', \
                                label = r"$" + formatted_species_name + "$")
            ax.legend()

            if ismassf:
                ax.set_ylabel(SYMBOLS["massf"] + " (" + SI_UNITS["massf"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Multiple Mass Fractions at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + ' '.join(var) \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
            elif ismolef:
                ax.set_ylabel(SYMBOLS["molef"] + " (" + SI_UNITS["molef"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Multiple Mole Fractions at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + ' '.join(var) \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
                
        else: #either regular property, specific massf or molef for a species or "massf"/"molef"
            if var == "massf": #plot all mass fractions
                column_names = list(data_object.columns)
                massf_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names]
                for idx, species in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    ax.scatter(data_object["pos_x"], \
                            data_object[massf_species_names[idx]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
                
            elif var == "molef": #plot all mole fractions
                column_names = list(data_object.columns)
                molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names]
                for idx, species in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    ax.scatter(data_object["pos_x"], \
                            data_object[molef_species_names[idx]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                            + formatted_title_time + r'$\mu$' + "s")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

            elif len(var) > 5 and var[:5] == "massf": #plot specific mass fraction
                species_name = var[6:] # Trim off starting massf_
                split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
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
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

            elif len(var) > 5 and var[:5] == "molef": #plot specific mole fraction
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
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
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

            else: #plot regular property
                ax.scatter(data_object["pos_x"], data_object[var], marker = '.')
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of " + SYMBOLS[var] + " at t = " \
                                + formatted_title_time + r'$\mu$' + "s")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"

        ax.set_xlabel("Position (m)")
        ax.grid()
        current_dir = os.getcwd()
        fig.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_spatial_cell_data_waterfall_plots(data_files, plot_vars, label):
    """
    - data_files: list(list(str))
    - plot_vars: list(list(str), str)
    - label: str
    """
    data_objects = [None] * len(data_files)
    time_list = [None] * len(data_files)

    for ind, snapshot in enumerate(data_files):
        data_made = False
        for file in snapshot:
            data = ProcessSpatialCellData(spatial_cell_data_file = file)
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

        if isinstance(var, list):
            ismassf = False
            ismolef = True

            if var[0][:6] == "massf_":
                ismassf = True
                ismolef = False
            
            for prop in var:
                spcs_name = prop[6:]
                split_name = re.findall('[0-9]+|[A-Za-z]+', spcs_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[prop] = r'$-$'
                if ismassf:
                    SYMBOLS[prop] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'
                elif ismolef:
                    SYMBOLS[prop] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'
                for ind, time in enumerate(time_list):
                    formatted_time = f'{time / 1e-6:.3f}'
                    ax.scatter(data_objects[ind]["pos_x"], data_objects[ind][prop], \
                            marker = ".", label = SYMBOLS[prop] + "Distribution at " + \
                                formatted_time + r'$\mu$' + "s")
            ax.legend()

            if ismassf:
                ax.set_ylabel(SYMBOLS["massf"] + " (" + SI_UNITS["massf"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Mass Fractions at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + ' '.join(var) \
                    + " distribution at various times.jpg"

            elif ismolef:
                ax.set_ylabel(SYMBOLS["molef"] + " (" + SI_UNITS["molef"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Mole Fractions at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + ' '.join(var) \
                    + " distribution at various times.jpg"

        else:
            if var == "massf":
                column_names = list(data_objects[0].columns)
                massf_species_names = [species for species in column_names \
                                            if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names]
                for idx, species in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    for ind, time in enumerate(time_list):
                        formatted_time = f'{time / 1e-6:.3f}'
                        ax.scatter(data_objects[ind]["pos_x"], \
                                    data_objects[ind][massf_species_names[idx]], \
                                    marker = ".", label = r"$" + formatted_species_name + "$" \
                                                            + "at " + formatted_time \
                                                            + r'$\mu$' + "s")
                ax.legend()
                ax.set_ylabel(SYMBOLS["massf"] + " (" + SI_UNITS["massf"] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Mass Fractions at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                                + " distribution at various times.jpg"

            elif var == "molef":
                column_names = list(data_objects[0].columns)
                molef_species_names = [species for species in column_names \
                                            if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names]
                for idx, species in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    for ind, time in enumerate(time_list):
                        formatted_time = f'{time / 1e-6:.3f}'
                        ax.scatter(data_objects[ind]["pos_x"], \
                                    data_objects[ind][molef_species_names[idx]], \
                                    marker = ".", label = r"$" + formatted_species_name + "$" \
                                                            + "at " + formatted_time \
                                                            + r'$\mu$' + "s")
                ax.legend()
                ax.set_ylabel(SYMBOLS["molef"] + " (" + SI_UNITS["molef"] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of Mole Fractions at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                                + " distribution at various times.jpg"

            elif len(var) > 5 and var[:5] == "massf":
                species_name = var[6:] # Trim off starting massf_
                split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'
                for ind, time in enumerate(time_list):
                    formatted_time = f'{time / 1e-6:.3f}'
                    ax.scatter(data_objects[ind]["pos_x"], data_objects[ind][var], \
                                    marker = ".", label = r"$" + formatted_species_name + "$" \
                                                            + "at " + formatted_time \
                                                            + r'$\mu$' + "s")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of " + SYMBOLS[var] + " at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                                + " distribution at various times.jpg"

            elif len(var) > 5 and var[:] == "molef":
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'
                for ind, time in enumerate(time_list):
                    formatted_time = f'{time / 1e-6:.3f}'
                    ax.scatter(data_objects[ind]["pos_x"], data_objects[ind][var], \
                                    marker = ".", label = r"$" + formatted_species_name + "$" \
                                                            + "at " + formatted_time \
                                                            + r'$\mu$' + "s")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of " + SYMBOLS[var] + " at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                                + " distribution at various times.jpg"

            else:
                for ind, time in enumerate(time_list):
                    formatted_time = f'{time / 1e-6:.3f}'
                    ax.scatter(data_objects[ind]["pos_x"], data_objects[ind][var], \
                                    marker = ".", label = SYMBOLS[var] + " at " \
                                                            + formatted_time \
                                                            + r'$\mu$' + "s")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Distribution of " + SYMBOLS[var] + " at Various Times")
                file_name = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                                + " distribution at various times.jpg"

        ax.set_xlabel("Position (m)")
        ax.grid()
        current_dir = os.getcwd()
        fig.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_transient_cell_property_plots(cell_data_file, plot_vars, max_line = None):
    """
    - cell_data_file: str
    - plot_vars: list(list(str), str)
    - max_line: int if not None
    """
    cell_data_object = ProcessTransientCellData(data_file_name = cell_data_file, \
                                                    max_line = max_line)
    cell_data = cell_data_object.cell_data
    cell_data = cell_data.sort_values(by = ["time"])
    cell_id = cell_data_object.cell_id
    sim_number = cell_data_object.sim_number

    for var in plot_vars:
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)

        if isinstance(var, list):
            ismolef = True
            ismassf = False

            if var[0][:6] == "massf_":
                ismassf = True
                ismolef = False
            
            for prop in var:
                spcs_name = prop[6:]
                split_name = re.findall('[0-9]+|[A-Za-z]+', spcs_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[prop] = r'$-$'
                if ismassf:
                    SYMBOLS[prop] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'
                elif ismolef:
                    SYMBOLS[prop] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'
                ax.scatter(cell_data["time"], cell_data[prop], marker = ".", \
                                label = r"$" + formatted_species_name + "$")
            ax.legend()

            if ismassf:
                ax.set_ylabel(SYMBOLS["massf"] + " (" + SI_UNITS["massf"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of Multiple Mass Fractions at Cell ID " + \
                                str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + ' '.join(var) \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"

            elif ismolef:
                ax.set_ylabel(SYMBOLS["molef"] + " (" + SI_UNITS["molef"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of Multiple Mole Fractions at Cell ID " + \
                                str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + ' '.join(var) \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"
        
        else:
            if var == "massf":
                column_names = list(cell_data.columns)
                massf_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names]
                for idx, species in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    ax.scatter(cell_data["time"], cell_data[massf_species_names[idx]], \
                                marker = ".", label = r"$" + formatted_species_name + "$")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of Mass Fractions for Cell ID " \
                                    + str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + var \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"

            elif var == "molef":
                column_names = list(cell_data.columns)
                molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names]
                for idx, species in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    ax.scatter(cell_data["time"], cell_data[molef_species_names[idx]], \
                                marker = ".", label = r"$" + formatted_species_name + "$")
                ax.legend()
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                                rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of Mole Fractions for Cell ID " \
                                    + str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + var \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"

            elif len(var) > 5 and var[:5] == "massf":
                species_name = var[6:] # Trim off starting massf_
                split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$ $Fraction$'
                ax.scatter(cell_data["time"], cell_data[var], marker = '.')
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of " + SYMBOLS[var] + " for Cell ID " \
                                    + str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + var \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"

            elif len(var) > 5 and var[:5] == "molef":
                species_name = var[6:] # Trim off starting molef_
                split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                SI_UNITS[var] = r'$-$'
                SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$ $Fraction'
                ax.scatter(cell_data["time"], cell_data[var], marker = '.')
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of " + SYMBOLS[var] + " for Cell ID " \
                                    + str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + var \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"

            else:
                ax.scatter(cell_data["time"], cell_data[var], marker = '.')
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of " + SYMBOLS[var] + " for Cell ID " \
                                    + str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + var \
                                +  " Transient Evolultion for Cell ID " + str(cell_id) + ".jpg"

        ax.set_xlabel("Time (s)")
        ax.grid()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_transient_interface_data_plots(interface_data_file, plot_vars, max_line = None):
    """
    - interface_data_file: str
    - plot_vars: list(list(str), str)
    - max_line: int if not None
    """
    interface_data_object = ProcessTransientInterfaceData(\
                                data_file_name = interface_data_file, \
                                max_line = max_line)
    interface_data = interface_data_object.interface_data
    interface_data = interface_data.sort_values(by = ["time"])
    interface_id = interface_data_object.interface_id
    sim_number = interface_data_object.sim_number

    for var in plot_vars:
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)

        if var == "mass_flux" or var == "xMom_flux" or var == "energy_flux":
            ax.scatter(interface_data["time"], interface_data[var] * interface_data["A"], \
                        marker = '.')
        else:
            ax.scatter(interface_data["time"], interface_data[var], marker = '.')
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
        ax.set_title("Transient Evolution of " + SYMBOLS[var] + " for Interface ID " \
                                    + str(interface_id))
        file_name = "Sim " + str(sim_number) + " " + var \
                                + " Transient Evolultion for Interface ID " \
                                + str(interface_id) + ".jpg"
        ax.set_xlabel("Time (s)")
        ax.grid()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_spatial_cell_data_plots_with_multiple_y_axes(spatial_data_files, plot_vars, \
                                                          visible_axes):
    """
    - spatial_data_files: list(str)
    - plot_vars: list(list(list(str)))
    - visible_axes: list(list(int))
    """

    data_object = None
    for component_data_file in spatial_data_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
        if data_object is None:
            data_object = data.component_data
            t_final = data.t_final
            sim_number = data.sim_number
        else:
            data_object = pd.concat([data_object, data.component_data], axis = 0)
    data_object = data_object.sort_values(by = ["pos_x"])

    for plot in plot_vars: # plot_vars is a list of lists. Each sub-list describes what properties
                            # are wanted on each plot. Within these, there are lists defining what
                            # properties get plotted on the same axis
                            