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
from sympy.abc import M

from Algorithms.DT_1D_V5.post_processing.process_spatial_cell_data \
        import ProcessSpatialCellData
from Algorithms.DT_1D_V5.post_processing.process_spatial_interface_data \
        import ProcessSpatialInterfaceData
from Algorithms.DT_1D_V5.post_processing.process_transient_cell_data \
        import ProcessTransientCellData
from Algorithms.DT_1D_V5.post_processing.process_transient_interface_data \
        import ProcessTransientInterfaceData
from Algorithms.DT_1D_V5.post_processing.generate_1D_averaged_data_from_2D_eilmer_data \
        import Averaged2DEilmerDataInto1D
from Algorithms.DT_1D_V5.post_processing.generate_average_from_eilmer_extract_line import generate_average_from_eilmer_extract
from Algorithms.MultiDimensionNewton.MultiDimensionNS import NewtonSolver


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
            data_object = pd.concat([data_object, data.component_data], \
                                        axis = 0, ignore_index = True)
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
                                                axis = 0, ignore_index = True)
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
                                                            + " at " + formatted_time \
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
                                                            + " at " + formatted_time \
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
                                                            + " at " + formatted_time \
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
                                                            + " at " + formatted_time \
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
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"

            elif ismolef:
                ax.set_ylabel(SYMBOLS["molef"] + " (" + SI_UNITS["molef"] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of Multiple Mole Fractions at Cell ID " + \
                                str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + ' '.join(var) \
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"
        
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
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"

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
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"

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
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"

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
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"

            else:
                ax.scatter(cell_data["time"], cell_data[var], marker = '.')
                ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
                ax.set_title("Transient Evolution of " + SYMBOLS[var] + " for Cell ID " \
                                    + str(cell_id))
                file_name = "Sim " + str(sim_number) + " " + var \
                                +  " Transient Evolultion for Cell ID " \
                                + str(cell_id) + ".jpg"

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
def generate_spatial_cell_data_plots_with_multiple_y_axes(spatial_cell_data_files, \
                                                            plot_vars, visible_axes, label):
    """
    - spatial_cell_data_files: list(str)
    - plot_vars: list(list(list(str)))
    - visible_axes: list(list(int))
    """

    data_object = None
    for component_data_file in spatial_cell_data_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
        if data_object is None:
            data_object = data.component_data
            t_final = data.t_final
            sim_number = data.sim_number
        else:
            data_object = pd.concat([data_object, data.component_data], \
                                        axis = 0, ignore_index = True)
    data_object = data_object.sort_values(by = ["pos_x"])

    formatted_title_time = '{:.3f}'.format(t_final / 1e-6)
    formatted_file_name_time = '{:.9f}'.format(t_final)

    marker_list = [".", "x", "+", "*", "|", "_", "^", "s", "o"]

    for plot_num in range(len(plot_vars)): # plot_vars is a list of lists. 
                            # Each sub-list describes what properties are wanted on each plot. 
                            # Within these, there are lists defining what properties 
                            # get plotted on the same axis
        num_active_axes = sum(visible_axes[plot_num])

        fig, ax = plt.subplots()
        axes = [ax] + [ax.twinx() for _ in range(num_active_axes - 1)]
        fig.set_size_inches(15, 5)

        fig.subplots_adjust(right=0.9 - 0.05 * (num_active_axes - 1))

        if ["D"] in plot_vars[plot_num]:
            A_c = data_object["A_c"]
            data_object["D"] = (4.0 * A_c / np.pi) ** 0.5

        new_axes_location = 1.0
        right_spine_taken = False

        for axis_plot_num in range(1, len(plot_vars[plot_num])):
            if visible_axes[plot_num][axis_plot_num] == 1:
                if not right_spine_taken:
                    right_spine_taken = True
                else:
                    new_axes_location += 0.1
                    axes[axis_plot_num].spines["right"].\
                                        set_position(('axes', new_axes_location))
                    axes[axis_plot_num].set_frame_on(True)
                    axes[axis_plot_num].patch.set_visible(False)
            
            else:
                axes[axis_plot_num].spines['right'].set_visible(False)
                axes[axis_plot_num].get_yaxis().set_visible(False)

        plot_list = []
        marker_location = 0

        for axis_plot_ind, axis_plot_vars in enumerate(plot_vars[plot_num]):
            if axis_plot_vars == ["D"]:
                axes[axis_plot_ind].set_ylim((0.0, 3.0 * np.max(data_object["D"])))
                p_ind, = axes[axis_plot_ind].plot(data_object["pos_x"], data_object["D"], \
                                                    label = SYMBOLS["D"], color = "black")
                if visible_axes[plot_num][axis_plot_ind] == 1:
                    plot_list.append(p_ind)
                    axes[axis_plot_ind].set_ylabel(SYMBOLS["D"] + " (" + SI_UNITS["D"] +")", \
                            ha = "right")
            elif axis_plot_vars == ["massf"]: #plot all mass fractions on the 1 axis
                column_names = list(data_object.columns)
                massf_species_names = [species for species in column_names \
                                            if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names]
                for sp_idx, sp in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', sp)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[axis_plot_ind].scatter(data_object["pos_x"], \
                                            data_object[massf_species_names[sp_idx]], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[axis_plot_ind].set_ylabel(SYMBOLS["massf"] + " (" \
                                    + SI_UNITS["massf"] +")", ha = "right")
            elif axis_plot_vars == ["molef"]: #plot all mole fractions on the 1 axis
                column_names = list(data_object.columns)
                molef_species_names = [species for species in column_names \
                                            if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names]
                for sp_idx, sp in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', sp)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[axis_plot_ind].scatter(data_object["pos_x"], \
                                            data_object[molef_species_names[sp_idx]], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[axis_plot_ind].set_ylabel(SYMBOLS["molef"] + " (" \
                                    + SI_UNITS["molef"] +")", ha = "right")
            
            elif len(axis_plot_vars[0]) > 5 and axis_plot_vars[0][:5] == "massf":
                massf_label_list = []
                for massf_species in axis_plot_vars:
                    species_name = massf_species[6:] # Trim off starting massf_
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[massf_species] = r'$-$'
                    SYMBOLS[massf_species] = r"$f_{" + formatted_species_name + r'}$'
                    p_ind = axes[axis_plot_ind].scatter(data_object["pos_x"], \
                                            data_object[massf_species], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    massf_label_list.append(SYMBOLS[massf_species])
                    marker_location += 1
                axes[axis_plot_ind].set_ylabel(', '.join(massf_label_list) + " (" + \
                                    SI_UNITS[axis_plot_vars[0]] +")", ha = "right")

            elif len(axis_plot_vars[0]) > 5 and axis_plot_vars[0][:5] == "molef":
                molef_label_list = []
                for molef_species in axis_plot_vars:
                    species_name = molef_species[6:] # Trim off starting molef_
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[molef_species] = r'$-$'
                    SYMBOLS[molef_species] = r"$f_{" + formatted_species_name + r'}$'
                    p_ind = axes[axis_plot_ind].scatter(data_object["pos_x"], \
                                            data_object[molef_species], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    molef_label_list.append(SYMBOLS[molef_species])
                    marker_location += 1
                axes[axis_plot_ind].set_ylabel(', '.join(molef_label_list) + " (" + \
                                    SI_UNITS[axis_plot_vars[0]] +")", ha = "right")

            else:
                for axis_plot_var in axis_plot_vars:
                    p_ind = axes[axis_plot_ind].scatter(data_object["pos_x"], \
                                                    data_object[axis_plot_var], \
                                                    marker = marker_list[marker_location], \
                                                    label = SYMBOLS[axis_plot_var])
                    marker_location += 1
                    if visible_axes[plot_num][axis_plot_ind] == 1:
                        plot_list.append(p_ind)
                if visible_axes[plot_num][axis_plot_ind] == 1:
                    axes[axis_plot_ind].set_ylabel(', '.join([SYMBOLS[var] \
                                                for var in axis_plot_vars]) \
                                                + " (" + SI_UNITS[axis_plot_vars[0]] +")", \
                                                ha = "right")
        ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), \
                    handles = plot_list, ncol = len(plot_list))
        axes[0].set_xlabel("Position (m)")
        axes[0].grid(visible = True, axis = 'y', which = "both")
        axes[0].minorticks_on()
        ax.set_title("Distribution of Multiple Properties at t = " \
                                    + formatted_title_time + r'$\mu$' + "s")
        file_name = "Sim " + str(sim_number) + " plot" + str(plot_num) + label + \
            " multiple y axes spatial cell distributions at t=" \
                + formatted_file_name_time + ".jpg"
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()

####################################################################################################
def generate_thrust_contribution_plot(spatial_cell_data_files, spatial_interface_data_file):
    """
    spatial_cell_data_files: list(str)
    spatial_interface_data_file: str
    """
    cell_data_object = None
    for component_cell_data_file in spatial_cell_data_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_cell_data_file)
        if cell_data_object is None:
            cell_data_object = data.component_data
            t_final = data.t_final
            sim_number = data.sim_number
        else:
            cell_data_object = pd.concat([cell_data_object, data.component_data], \
                                            axis = 0, ignore_index = True)
    cell_data_object = cell_data_object.sort_values(by = ["pos_x"])
    
    cell_data_object["net_thrust"] = None
    cell_data_object["D"] = (4.0 * cell_data_object["A_c"] / np.pi) ** 0.5
    interface_data_object = ProcessSpatialInterfaceData( \
                                spatial_interface_data_file = spatial_interface_data_file)
    interface_data = interface_data_object.interface_data
    interface_data = interface_data.sort_values(by = ["pos_x"])

    formatted_title_time = '{:.3f}'.format(t_final / 1e-6)
    formatted_file_name_time = '{:.9f}'.format(t_final)
    
    for ind, cell in enumerate(cell_data_object["pos_x"]):
        net_thrust_cell = cell_data_object["p"][ind] * \
                        (interface_data["A"][ind + 1] - interface_data["A"][ind])
        cell_data_object["net_thrust"][ind] = net_thrust_cell
    
    cell_data_object["cumulative_thrust"] = cell_data_object["net_thrust"].cumsum()

    fig, ax = plt.subplots()
    fig.set_size_inches(15, 5)
    axes = [ax] + [ax.twinx()]

    plot_list = []

    axes[1].spines["right"].set_visible(False) #Turn off diameter spine
    axes[1].get_yaxis().set_visible(False)

    p_ind0 = axes[0].scatter(cell_data_object["pos_x"], \
                        cell_data_object["net_thrust"], \
                        label = r"$Net$ $Thrust$ $(N)$", marker = ".")
    p_ind1 = axes[0].scatter(cell_data_object["pos_x"], \
                        cell_data_object["cumulative_thrust"], \
                        label = r"$Cumulative$ $Thrust$ $(N)$", marker = "x")
    plot_list.append(p_ind0)
    plot_list.append(p_ind1)
    axes[0].set_ylabel(r"$Thrust$ $(N)$")
    axes[1].plot(cell_data_object["pos_x"], cell_data_object["D"])
    axes[1].set_ylim((0.0, 3.0 * np.max(cell_data_object["D"])))
    ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), \
                handles = plot_list, ncol = 2)
        
    axes[0].set_xlabel("Position (m)")
    axes[0].grid(visible = True, axis = 'y', which = "both")
    axes[0].minorticks_on()
    ax.set_title("Net and Cumulative Thrust Contribution at " + formatted_title_time \
                        + r"$\mu$" + "s")
    file_name = "Sim" + str(sim_number) + "Net and Cumulative Thrust Contribution at " \
        + formatted_file_name_time +".jpg"
    
    current_dir = os.getcwd()
    plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
    plt.close()
    last_cell_index = len(cell_data_object["cumulative_thrust"]) - 1
    total_cumulative_thrust = cell_data_object["cumulative_thrust"][last_cell_index]
    print("Total cumulative thrust:", total_cumulative_thrust, "(N)")

####################################################################################################
def generate_transient_gross_thrust_plot(transient_interface_data_file, max_line = None):
    """
    - transient_interface_data_file: str
    """
    transient_interface_data_object = ProcessTransientInterfaceData(\
                    data_file_name = transient_interface_data_file, max_line = max_line)
    transient_interface_data = transient_interface_data_object.interface_data
    transient_interface_data = transient_interface_data.sort_values(by = ["time"])
    sim_number = transient_interface_data_object.sim_number

    A_e = transient_interface_data["A"]
    vel_x_e = transient_interface_data["vel_x"]
    m_dot = transient_interface_data["mass_flux"]
    p_e = transient_interface_data["p"]

    transient_interface_data["Thrust"] = A_e * (m_dot * vel_x_e + p_e)

    fig, ax = plt.subplots()
    fig.set_size_inches(15, 5)

    ax.scatter(transient_interface_data["time"] * 1e3, transient_interface_data["Thrust"], \
                marker = ".")
    ax.set_ylabel(r"$Thrust$ $(N)$", rotation = "horizontal", ha = "right")
    ax.set_title("Transient Gross Thrust Profile")
    ax.set_xlabel("Time (ms)")
    plt.grid()
    file_name = "Sim" + str(sim_number) + "TransientGrossThrustProfile.jpg"
    current_dir = os.getcwd()
    plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
    plt.close()
####################################################################################################
def generate_spatial_interface_data_plots_with_multiple_y_axes(spatial_interface_data_file, \
                                                            plot_vars, visible_axes, label):
    """
    - spatial_interface_data_file: str
    - plot_vars: list(list(list(str)))
    - visible_axes: list(list(int))
    - label: str
    """

    interface_data_object = ProcessSpatialInterfaceData(spatial_interface_data_file = \
                                                            spatial_interface_data_file)
    interface_data = interface_data_object.interface_data
    interface_data = interface_data.sort_values(by = ["pos_x"])

    t_final = interface_data_object.t_final
    sim_number = interface_data_object.sim_number

    formatted_title_time = '{:.3f}'.format(t_final / 1e-6)
    formatted_file_name_time = '{:.9f}'.format(t_final)

    marker_list = [".", "x", "+", "*", "|", "_", "^", "s", "o"]

    for plot_num in range(len(plot_vars)):
        num_active_axes = sum(visible_axes[plot_num])
        fig, ax = plt.subplots()
        axes = [ax] + [ax.twinx() for _ in range(num_active_axes - 1)]
        fig.set_size_inches(15, 5)

        fig.subplots_adjust(right=0.9 - 0.05 * (num_active_axes - 1))

        if ["D"] in plot_vars[plot_num]:
            A_c = interface_data["A_c"]
            interface_data["D"] = (4.0 * A_c / np.pi) ** 0.5
        
        new_axes_location = 1.0
        right_spine_taken = False

        for axis_plot_num in range(1, len(plot_vars[plot_num])):
            if visible_axes[plot_num][axis_plot_num] == 1:
                if not right_spine_taken:
                    right_spine_taken = True
                else:
                    new_axes_location += 0.1
                    axes[axis_plot_num].spines["right"].\
                                        set_position(('axes', new_axes_location))
                    axes[axis_plot_num].set_frame_on(True)
                    axes[axis_plot_num].patch.set_visible(False)
            
            else:
                axes[axis_plot_num].spines['right'].set_visible(False)
                axes[axis_plot_num].get_yaxis().set_visible(False)

        plot_list = []
        marker_location = 0
        for axis_plot_ind, axis_plot_vars in enumerate(plot_vars[plot_num]):
            if axis_plot_vars == ["D"]:
                axes[axis_plot_ind].set_ylim((0.0, 3.0 * np.max(interface_data["D"])))
                p_ind, = axes[axis_plot_ind].plot(interface_data["pos_x"], \
                                                    interface_data["D"], \
                                                    label = SYMBOLS["D"], color = "black")
                if visible_axes[plot_num][axis_plot_ind] == 1:
                    plot_list.append(p_ind)
                    axes[axis_plot_ind].set_ylabel(SYMBOLS["D"] + " (" + SI_UNITS["D"] +")", \
                            ha = "right")
            
            else:
                for axis_plot_var in axis_plot_vars:
                    if axis_plot_var in set(("mass_flux", "xMom_flux", "energy_flux")):
                        interface_data[axis_plot_var] *= axis_plot_var["A"]
                    p_ind = axes[axis_plot_ind].scatter(interface_data["pos_x"], \
                                                    interface_data[axis_plot_var], \
                                                    marker = marker_list[marker_location], \
                                                    label = SYMBOLS[axis_plot_var])
                    marker_location += 1
                    if visible_axes[plot_num][axis_plot_ind] == 1:
                        plot_list.append(p_ind)
                if visible_axes[plot_num][axis_plot_ind] == 1:
                    axes[axis_plot_ind].set_ylabel(', '.join([SYMBOLS[var] \
                                                for var in axis_plot_vars]) \
                                                + " (" + SI_UNITS[axis_plot_vars[0]] +")", \
                                                ha = "right")
        ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), \
                    handles = plot_list, ncol = len(plot_list))
        axes[0].set_xlabel("Position (m)")
        axes[0].grid(visible = True, axis = 'y', which = "both")
        axes[0].minorticks_on()
        ax.set_title("Distribution of Multiple Properties at t = " \
                                    + formatted_title_time + r'$\mu$' + "s")
        file_name = "Sim " + str(sim_number) + " plot" + str(plot_num) + label + \
            " multiple y axes spatial interface distributions at t=" \
                + formatted_file_name_time + ".jpg"
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
class GenerateTransientAnimationOfSpatialCellData():
    def __init__(self, spatial_cell_data_files, plot_vars, slow_down_factor, label) -> None:
        self.data_objects = [None] * len(spatial_cell_data_files)
        self.time_list = [None] * len(spatial_cell_data_files)

        for snapshot_ind, snapshot in enumerate(spatial_cell_data_files):
            for component_data_file in snapshot:
                component_data_object = ProcessSpatialCellData(spatial_cell_data_file = \
                                                               component_data_file)
                if self.data_objects[snapshot_ind] is None:
                    self.data_objects[snapshot_ind] = component_data_object.component_data
                    self.time_list[snapshot_ind] = component_data_object.t_final
                    sim_number = component_data_object.sim_number
                else:
                    self.data_objects[snapshot_ind] = \
                                            pd.concat([self.data_objects[snapshot_ind], \
                                                    component_data_object.component_data], \
                                                    axis = 0, ignore_index = True)
            self.data_objects[snapshot_ind] = self.data_objects[snapshot_ind].sort_values(\
                                                                            by = ["pos_x"])

        self.dt_list = (np.array(self.time_list[1:]) - np.array(self.time_list[:-1])).tolist()

        ffwriter = animation.FFMpegWriter()
        current_dir = os.getcwd()
        for var in plot_vars:
            self.fig, self.ax = plt.subplots()
            self.fig.set_size_inches(15, 5)

            y_min = min([min(df[var]) for df in self.data_objects])
            y_max = max([max(df[var]) for df in self.data_objects])
            dy = y_max - y_min

            y_min, y_max = y_min - 0.1 * dy, y_max + 0.1 * dy

            anim = FuncAnimation(self.fig, self.update, frames = len(self.time_list), \
                                interval = 1, fargs = [var, slow_down_factor, y_min, y_max])
            file_name = "Sim" + str(sim_number) + label + "SpatialCellAnimationOf" + var \
                            + ".mp4"
            anim.save(current_dir + "/plots/" + file_name, writer = ffwriter)
            
    def update(self, frame, *fargs):
        self.ax.clear()
        var, slow_down_factor, y_min, y_max = fargs
        self.ax.set_ylim([y_min, y_max])
        self.ax.grid()
        self.ax.set_xlabel("Position (m)")
        self.ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                            rotation = "horizontal", ha = "right")
        self.ax.set_title("Transient Evolution of " + SYMBOLS[var] + " Distribution")
        x = self.data_objects[frame]["pos_x"]
        y = self.data_objects[frame][var]
        scat = self.ax.scatter(x, y, marker = ".")
        if frame < len(self.data_objects) - 1:
            time.sleep(self.dt_list[frame] * slow_down_factor)
        return (scat, )
####################################################################################################
def compare_to_1D_analytical_result(spatial_cell_data_files, A_throat, x_throat, gamma, R):
    data_object = None
    for component_data_file in spatial_cell_data_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
        if data_object is None:
            data_object = data.component_data
            sim_number = data.sim_number
        else:
            data_object = pd.concat([data_object, data.component_data], \
                                        axis = 0, ignore_index = True)
    data_object = data_object.sort_values(by = ["pos_x"])

    Ma_anal = [None] * len(data_object["pos_x"])

    def Ma_function(M, kwargs):
        gamma = kwargs["gamma"]
        A_star = kwargs["A_star"]
        A_x = kwargs["A_x"]
        return M ** -1 * (2.0 / (gamma + 1.0) * (1.0 + 0.5 * (gamma - 1.0) * M ** 2.0)) \
                            ** (0.5 * (gamma + 1.0) / (gamma - 1.0)) - A_x / A_star

    for ind, x_i in enumerate(data_object["pos_x"]):
        if x_i < x_throat:
            Ma_guess = 0.1
        else:
            Ma_guess = 1.5
        result = NewtonSolver(initial_guess = [Ma_guess], f = [Ma_function], \
                                jacobian_method = "Analytical", \
                                residual_smoothing = False, variables = [M], \
                                verbose = False, A_star = A_throat, \
                                A_x = data_object["A_c"][ind], gamma = gamma)
        Ma_anal[ind] = result.final_result[0]
    
    p_t = data_object["p_t"][0] 
    T_t = data_object["T_t"][0] 
    T_anal = (T_t / (1.0 + 0.5 * (gamma - 1.0) * np.array(Ma_anal) ** 2.0)).tolist()
    p_anal = (p_t / \
            (1.0 + 0.5 * (gamma - 1.0) * np.array(Ma_anal) ** 2.0) \
            ** (gamma / (gamma - 1.0))).tolist()
    vel_x_anal = (np.array(Ma_anal) * (gamma * R * np.array(T_anal)) ** 0.5).tolist()
    rho_anal = (np.array(p_anal) / (R * np.array(T_anal))).tolist()
    props = ["p", "T", "rho", "Ma", "vel_x"]
    analytical_props = [p_anal, T_anal, rho_anal, Ma_anal, vel_x_anal]
    current_dir = os.getcwd()
    for i in range(len(props)):
        fig, ax = plt.subplots()
        fig.set_size_inches(15, 5)

        ax.scatter(data_object["pos_x"], data_object[props[i]], marker = ".",
                        label = "Simulation")
        ax.plot(data_object["pos_x"], analytical_props[i], \
                        label = "Analytical")

        ax.set_ylabel(SYMBOLS[props[i]] + " (" + SI_UNITS[props[i]] +")", \
                    rotation = "horizontal", ha = "right")
        ax.set_xlabel("Position (m)")
        ax.set_title("Comparison to Analytical Solution for " + SYMBOLS[props[i]])
        ax.legend()
        ax.grid()

        file_name = "Sim" + str(sim_number) + "Comparison to Analytical Solution For " + \
                       props[i] + ".jpg"
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches="tight")
        plt.close()
####################################################################################################
def generate_transient_cell_property_plots_with_multiple_y_axes(transient_cell_data_file, \
                                                                plot_vars, visible_axes, \
                                                                max_line = None):
    cell_data_total = ProcessTransientCellData(data_file_name = transient_cell_data_file, \
                                            max_line = max_line)
    sim_number = cell_data_total.sim_number
    cell_id = cell_data_total.cell_id
    cell_data = cell_data_total.cell_data

    cell_data = cell_data.sort_values(by = ["time"])

    marker_list = [".", "x", "+", "*", "|", "_", "^", "s", "o"]

    for plot_num in range(len(plot_vars)):
        num_active_axes = sum(visible_axes[plot_num])

        fig, ax = plt.subplots(figsize=(15,5))
        axes = [ax] + [ax.twinx() for _ in range(num_active_axes - 1)]
        #fig.set_size_inches(15, 5)

        fig.subplots_adjust(right=0.9 - 0.05 * (num_active_axes - 1))

        new_axes_location = 1.0
        right_spline_taken = False

        for axis_plot_num in range(1, len(plot_vars[plot_num])):
            if visible_axes[plot_num][axis_plot_num] == 1:
                if not right_spline_taken:
                    right_spline_taken = True
                else:
                    new_axes_location += 0.1
                    axes[axis_plot_num].spines["right"].\
                                        set_position(('axes', new_axes_location))
                    axes[axis_plot_num].set_frame_on(True)
                    axes[axis_plot_num].patch.set_visible(False)
            
            else:
                axes[axis_plot_num].spines["right"].set_visible(False)
                axes[axis_plot_num].get_yaxis().set_visible(False)
        
        plot_list = []
        marker_location = 0

        for axis_plot_ind, axis_plot_vars in enumerate(plot_vars[plot_num]):
            if axis_plot_vars == ["massf"]:
                column_names = list(cell_data.columns)
                massf_species_names = [species for species in column_names \
                                        if len(species) > 5 and species[:5] == "massf"]
                species_names = [species.split("_")[1] for species in massf_species_names]
                for sp_idx, sp in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', sp)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[axis_plot_ind].scatter(cell_data["time"] * 1e3, \
                                            cell_data[massf_species_names[sp_idx]], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[axis_plot_ind].set_ylabel(SYMBOLS["massf"] + " (" \
                                                    + SI_UNITS["massf"] +")", ha = "right")
            
            elif axis_plot_vars == ["molef"]:
                column_names = list(cell_data.columns)
                molef_species_names = [species for species in column_names \
                                        if len(species) > 5 and species[:5] == "molef"]
                species_names = [species.split("_")[1] for species in molef_species_names]
                for sp_idx, sp in enumerate(species_names):
                    split_name = re.findall('[0-9]+|[A-Za-z]+', sp)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    p_ind = axes[axis_plot_ind].scatter(cell_data["time"] * 1e3, \
                                            cell_data[molef_species_names[sp_idx]], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    plot_list.append(p_ind)
                    marker_location += 1
                    axes[axis_plot_ind].set_ylabel(SYMBOLS["molef"] + " (" \
                                                    + SI_UNITS["molef"] +")", ha = "right")
                    
            elif len(axis_plot_vars[0]) > 5 and axis_plot_vars[0][:5] == "massf":
                massf_label_list = []
                for massf_species in axis_plot_vars:
                    species_name = massf_species[6:] # Trim off starting _massf
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[massf_species] = r'$-$'
                    SYMBOLS[massf_species] = r"$f_{" + formatted_species_name + r'}$'
                    p_ind = axes[axis_plot_ind].scatter(cell_data["time"] * 1e3, \
                                            cell_data[massf_species], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    massf_label_list.append(SYMBOLS[massf_species])
                    marker_location += 1
                axes[axis_plot_ind].set_ylabel(', '.join(massf_label_list) + " (" + \
                                    SI_UNITS[axis_plot_vars[0]] +")", ha = "right")
            
            elif len(axis_plot_vars[0]) > 5 and axis_plot_vars[0][:5] == "molef":
                molef_label_list = []
                for molef_species in axis_plot_vars:
                    species_name = molef_species[6:] # Trim off starting _molef
                    split_name = re.findall('[0-9]+|[A-Za-z]+', species_name)
                    for ind, val in enumerate(split_name):
                        if val.isnumeric():
                            split_name[ind] = "_{" + split_name[ind] + "}"
                    formatted_species_name = ''.join(split_name)
                    SI_UNITS[molef_species] = r'$-$'
                    SYMBOLS[molef_species] = r"$f_{" + formatted_species_name + r'}$'
                    p_ind = axes[axis_plot_ind].scatter(cell_data["time"] * 1e3, \
                                            cell_data[molef_species], \
                                            marker = marker_list[marker_location], \
                                            label = r"$f_{" + formatted_species_name + "}$")
                    molef_label_list.append(SYMBOLS[molef_species])
                    marker_location += 1
                axes[axis_plot_ind].set_ylabel(', '.join(molef_label_list) + " (" + \
                                    SI_UNITS[axis_plot_vars[0]] +")", ha = "right")
            
            else:
                for axis_plot_var in axis_plot_vars:
                    p_ind = axes[axis_plot_ind].scatter(cell_data["time"] * 1e3, \
                                                cell_data[axis_plot_var], \
                                                marker = marker_list[marker_location], \
                                                label = SYMBOLS[axis_plot_var])
                    marker_location += 1
                    if visible_axes[plot_num][axis_plot_ind] == 1:
                        plot_list.append(p_ind)
                if visible_axes[plot_num][axis_plot_ind] == 1:
                    axes[axis_plot_ind].set_ylabel(', '.join([SYMBOLS[var] \
                                                for var in axis_plot_vars]) \
                                                + " (" + SI_UNITS[axis_plot_vars[0]] +")", \
                                                ha = "right")
        ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), \
                    handles = plot_list, ncol = len(plot_list))
        axes[0].set_xlabel("Time (ms)")
        axes[0].grid(visible = True, axis = 'y', which = "both")
        axes[0].minorticks_on()    
        ax.set_title("Transient Evolution of Multiple Properties for Cell ID " + str(cell_id))
        file_name = "Sim " + str(sim_number) + " plot" + str(plot_num) + \
            " multiple y axes transient cell property evolution at cell id" \
                + str(cell_id) + ".jpg"
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_transient_interface_property_plots_with_multiple_y_axes(\
                                                            transient_interface_data_file, \
                                                            plot_vars, visible_axes, \
                                                            max_line = None):
    interface_data_total = ProcessTransientInterfaceData(\
                                data_file_name = transient_interface_data_file, \
                                max_line = max_line)
    interface_data = interface_data_total.interface_data
    sim_number = interface_data_total.sim_number
    interface_id = interface_data_total.interface_id

    interface_data = interface_data.sort_values(by = ["time"])

    marker_list = [".", "x", "+", "*", "|", "_", "^", "s", "o"]

    for plot_num in range(len(plot_vars)):
        num_active_axes = sum(visible_axes[plot_num])

        fig, ax = plt.subplots(figsize=(15,5))
        axes = [ax] + [ax.twinx() for _ in range(num_active_axes - 1)]
        fig.set_size_inches(15, 5)

        fig.subplots_adjust(right=0.9 - 0.05 * (num_active_axes - 1))

        new_axes_location = 1.0
        right_spline_taken = False

        for axis_plot_num in range(1, len(plot_vars[plot_num])):
            if visible_axes[plot_num][axis_plot_num] == 1:
                if not right_spline_taken:
                    right_spline_taken = True
                else:
                    new_axes_location += 0.1
                    axes[axis_plot_num].spines["right"].\
                                        set_position(('axes', new_axes_location))
                    axes[axis_plot_num].set_frame_on(True)
                    axes[axis_plot_num].patch.set_visible(False)
            
            else:
                axes[axis_plot_num].spines["right"].set_visible(False)
                axes[axis_plot_num].get_yaxis().set_visible(False)
        
        plot_list = []
        marker_location = 0

        for axis_plot_ind, axis_plot_vars in enumerate(plot_vars[plot_num]):
            for axis_plot_var in axis_plot_vars:
                if axis_plot_var in set(("mass_flux", "xMom_flux", "energy_flux")):
                    interface_data[axis_plot_var] *= interface_data["A"]
                p_ind = axes[axis_plot_ind].scatter(interface_data["time"] * 1e3, \
                                                    interface_data[axis_plot_var], \
                                                    marker = marker_list[marker_location], \
                                                    label = SYMBOLS[axis_plot_var])
                marker_location += 1
                if visible_axes[plot_num][axis_plot_ind] == 1:
                    plot_list.append(p_ind)
            if visible_axes[plot_num][axis_plot_ind] == 1:
                axes[axis_plot_ind].set_ylabel(', '.join([SYMBOLS[var] \
                                            for var in axis_plot_vars]) \
                                            + " (" + SI_UNITS[axis_plot_vars[0]] +")", \
                                            ha = "right")
    ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.12), \
                    handles = plot_list, ncol = len(plot_list))
    axes[0].set_xlabel("Time (ms)")
    axes[0].grid(visible = True, axis = 'y', which = "both")
    axes[0].minorticks_on()
    ax.set_title("Transient Evolution of Multiple Properties for Interface ID " \
                    + str(interface_id))
    file_name = "Sim " + str(sim_number) + " plot" + str(plot_num) + \
            " multiple y axes transient interface property evolution at interface id" \
                + str(interface_id) + ".jpg"
    current_dir = os.getcwd()
    plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
    plt.close()
####################################################################################################
def generate_1D_volume_averaged_plots_from_eilmer_data(eilmer_cell_data_files, properties, \
                                                pos_x_array, x_tol, plot_vars):
    eilmer_data = Averaged2DEilmerDataInto1D(eilmer_cell_data_files = eilmer_cell_data_files, \
                                             properties = properties, \
                                            pos_x_array = pos_x_array, \
                                            x_tol = x_tol)
    formatted_title_time = f"{eilmer_data.t_final / 1e-6:.3f}"
    formatted_file_name_time = f"{eilmer_data.t_final:.9f}"
    current_dir = os.getcwd()
    for var in plot_vars:
        fig, ax = plt.subplots()
        #fig.set_size_inches(15, 5)
        ax.scatter(eilmer_data.averaged_data["pos_x"], eilmer_data.averaged_data[var], marker = ".")
        ax.set_xlabel("Position (m)")
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
        ax.set_title("Averaged 1D Distribution of " + SYMBOLS[var] + " at t = " \
                        + formatted_title_time + r"$\mu$" + "s")
        ax.grid()
        file_name = "Averaged 1D plot of " + var + " from eilmer simulation at " + \
                                formatted_file_name_time  + ".jpg"
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def generate_1D_area_averaged_plots_from_eilmer_data(eilmer_extract_line_files, plot_vars, n_points, pos_x_array, ymax_array):
    """
    - eilmer_extract_line_files: list(str)
    - plot_vars: list(str)
    - n_points: int -> number of points to sample along for each x value
    - pos_x_array: list(float)
    - ymax_array: list(float)
    """
    n_cells = len(pos_x_array)
    data_df = pd.DataFrame(index=range(n_cells), columns=plot_vars)
    for ind, x_val in enumerate(pos_x_array):
        average_data, min_data, max_data = generate_average_from_eilmer_extract(line_extract_data_file = eilmer_extract_line_files[ind], x = x_val, ymin = 0.0, ymax = ymax_array[ind], n = n_points, extract_vars = plot_vars)
        for column in plot_vars:
            data_df[column][ind] = average_data[column]
    current_dir = os.getcwd()
    for var in plot_vars:
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.scatter(pos_x_array, data_df[var], marker = ".")
        ax.set_xlabel("Position (m)")
        ax.set_title("Area Averaged Spatial Distribution of " + SYMBOLS[var])
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
        
        file_name = "Area averaged spatial distribution of " + var + ".jpg"
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()
####################################################################################################
def compare_1D_simulation_to_area_averaged_eilmer_extract(eilmer_extract_line_files, simulation_files, plot_vars, n_points, add_bounding_lines):
    """
    - eilmer_extract_line_files: list(str)
    - simulation_files: list(str)
    - plot_vars: list(str)
    - n_points: int -> number of points to sample along for each x value
    """
    data_object = None
    for component_data_file in simulation_files:
        data = ProcessSpatialCellData(spatial_cell_data_file = component_data_file)
        if data_object is None:
            data_object = data.component_data
            t_final = data.t_final
            sim_number = data.sim_number
        else:
            data_object = pd.concat([data_object, data.component_data], \
                                        axis = 0, ignore_index = True)
    data_object = data_object.sort_values(by = ["pos_x"])
    pos_x_list = data_object["pos_x"].tolist()
    ymax_list = ((data_object["A_c"] / np.pi) ** 0.5).tolist()
    n_cells = data_object.shape[0]
    ave_data_df = pd.DataFrame(index=range(n_cells), columns=plot_vars)
    if add_bounding_lines:
        min_data_df = pd.DataFrame(index=range(n_cells), columns=plot_vars)
        max_data_df = pd.DataFrame(index=range(n_cells), columns=plot_vars)
    for ind, x_val in enumerate(pos_x_list):
        average_data, min_data, max_data = generate_average_from_eilmer_extract(line_extract_data_file = eilmer_extract_line_files[ind], x = x_val, ymin = 0.0, ymax = ymax_list[ind], n = n_points, extract_vars = plot_vars)
        for column in plot_vars:
            ave_data_df[column][ind] = average_data[column]
            if add_bounding_lines:
                min_data_df[column][ind] = min_data[column]
                max_data_df[column][ind] = max_data[column]
    current_dir = os.getcwd()
    for var in plot_vars:
        print(var)
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.scatter(pos_x_list, ave_data_df[var], marker = ".", label = "2D Average")
        if var == "vel":
            ax.scatter(pos_x_list, data_object["vel_x"], marker = ".", label = "1D Simulation")
        else:
            ax.scatter(pos_x_list, data_object[var], marker = ".", label = "1D Simulation")
        if add_bounding_lines:
            ax.fill_between(pos_x_list, min_data_df[var].tolist(), max_data_df[var].tolist(), alpha=0.2)
            ax.plot(pos_x_list, min_data_df[var].tolist(), c="k")
            ax.plot(pos_x_list, max_data_df[var].tolist(), c="k")
            #ax.scatter(pos_x_list, min_data_df[var], marker = ".", label = "2D Minimum")
            #ax.scatter(pos_x_list, max_data_df[var], marker = ".", label = "2D Maximum")
        ax.legend()
        ax.grid()
        ax.set_xlabel("Position (m)")
        ax.set_title("Comparison Between Area Averaged 2D Data and 1D Simulation Result of " + SYMBOLS[var])
        ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
        if add_bounding_lines:
            file_name = "Comparison of area averaged to 1D simulation spatial distribution of " + var + " with bounding lines.jpg"
        else:
            file_name = "Comparison of area averaged to 1D simulation spatial distribution of " + var + " without bounding lines.jpg"
        plt.savefig(current_dir + "/plots/" + file_name, bbox_inches = 'tight')
        plt.close()