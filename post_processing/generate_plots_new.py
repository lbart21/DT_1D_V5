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

from Algorithms.DT_1D_V5.post_processing.data_file_to_structured_data import GenerateDataObject
from Algorithms.DT_1D_V5.post_processing.cell_data_file_to_object import FormCellDataFromFile
from Algorithms.DT_1D_V5.post_processing.interface_data_file_to_object \
            import FormInterfaceDataFromFile
from Algorithms.DT_1D_V5.post_processing.process_eilmer_data import ProcessEilmerData
from Algorithms.DT_1D_V5.post_processing.SI_units_dictionary import SI_UNITS
from Algorithms.DT_1D_V5.post_processing.symbols import SYMBOLS

from Algorithms.DT_0D_V2.post_processing.interface_data_file_to_object \
            import FormInterfaceDataFromFile as interface_data_object_0d

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation


def generate_single_plots(data_files, plot_vars, label):
    data_object = None
    for component_data_file in data_files:
        data = GenerateDataObject(data_file_name = component_data_file)
        if data_object is None:
            data_object = data.component_data
            time = data.t_final
            sim_number = data.sim_number
        else:
            data_object = pd.concat([data_object, data.component_data], axis = 0)
    data_object = data_object.sort_values(by = ["pos_x"])

    for var in plot_vars:
        fig = plt.figure(figsize=(15, 5))
        formatted_title_time = '{:.3f}'.format(time / 1e-6)
        formatted_file_name_time = '{:.9f}'.format(time)
            
        if var == "massf": # Plot all mass fractions
            column_names = list(data_object.columns)
            massf_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "massf"]
            species_names = [species.split("_")[1] for species in massf_species_names] 
            for index, species in enumerate(species_names):
                split_name = re.findall('(\d+|[A-Za-z]+)', species)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                plt.scatter(data_object["pos_x"], \
                            data_object[massf_species_names[index]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
            plt.legend()

        elif len(var) > 5 and var[:5] == "massf": # Plot specific mass fractions
            species_name = var[6:] # Trim off starting massf_
            split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' + ' ' + r'$Fraction$'
            plt.scatter(data_object["pos_x"], \
                        data_object[var], marker = '.')
                
        elif var == "molef":
            column_names = list(data_object.columns)
            molef_species_names = [species for species in column_names \
                                                if len(species) > 5 and species[:5] == "molef"]
            species_names = [species.split("_")[1] for species in molef_species_names]
            for index, species in enumerate(species_names):
                split_name = re.findall('(\d+|[A-Za-z]+)', species)
                for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                plt.scatter(data_object["pos_x"], \
                            data_object[molef_species_names[index]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")
            plt.legend()

        elif len(var) > 5 and var[:5] == "molef": # Plot specific mass fractions
            species_name = var[6:] # Trim off starting molef_
            split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' + ' ' + r'$Fraction$'
            plt.scatter(data_object["pos_x"], \
                        data_object[var], marker = '.')
            
        elif var == "conc":
            column_names = list(data_object.columns)
            conc_species_names = [species for species in column_names \
                                            if len(species) > 4 and species[:4] == "conc"]
            species_names = [species.split("_")[1] for species in conc_species_names]
            for index, species in enumerate(species_names):
                if val.isnumeric():
                    split_name[ind] = "_{" + split_name[ind] + "}"
                formatted_species_name = ''.join(split_name)
                plt.scatter(data_object["pos_x"], \
                            data_object[conc_species_names[index]], \
                            marker = '.', label = r"$" + formatted_species_name + "$")

        else: # All other properties
            plt.scatter(data_object["pos_x"], data_object[var], marker = '.')
        plt.xlabel("Position (m)")
        plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
        plt.title("Distribution of " + SYMBOLS[var] + " at t = " \
                                                + formatted_title_time + r'$\mu$' + "s")
            
        filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var \
                    + " distribution at t = " + formatted_file_name_time + ".jpg"
        plt.grid()
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()

def generate_waterfall_plots(data_files, plot_vars, label):
    data_objects = [None] * len(data_files)
    time_list = [None] * len(data_files)

    for ind, snapshot in enumerate(data_files):
        data_made = False
        for file in snapshot:
            data = GenerateDataObject(data_file_name = file)
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
        fig = plt.figure(figsize=(15, 5))
        if len(var) > 5 and var[:5] == "massf":
            species_name = var[6:] # Trim off starting massf_
            split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mass$' \
                                + ' ' + r'$Fraction$'

        elif len(var) > 5 and var[:5] == "molef":
            species_name = var[6:] # Trim off starting molef_
            split_name = re.findall('(\d+|[A-Za-z]+)', species_name)
            for ind, val in enumerate(split_name):
                    if val.isnumeric():
                        split_name[ind] = "_{" + split_name[ind] + "}"
            formatted_species_name = ''.join(split_name)
            SI_UNITS[var] = r'$-$'
            SYMBOLS[var] = r"$" + formatted_species_name + "$" + ' ' + r'$Mole$' \
                                + ' ' + r'$Fraction$'

        for ind, time in enumerate(time_list):
            formatted_time = '{:.3f}'.format(time / 1e-6)
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
                    plt.scatter(data_objects[ind]["pos_x"], \
                                data_objects[ind][massf_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$" \
                                                + " at t = " + formatted_time + r'$\mu$' + "s")
                plt.legend()
                
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
                    plt.scatter(data_objects[ind]["pos_x"], \
                                data_objects[ind][molef_species_names[index]], \
                                marker = '.', label = r"$" + formatted_species_name + "$" \
                                                + " at t = " + formatted_time + r'$\mu$' + "s")
                plt.legend()

            else:
                plt.scatter(data_objects[ind]["pos_x"], data_objects[ind][var], \
                            label = "Distribution at t = " + formatted_time + r'$\mu$' + "s", \
                            marker = ".")
        plt.title("Distribution of " + SYMBOLS[var] + " at Various Times")
        plt.ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
        plt.xlabel("Position (m)")
        plt.legend()
        plt.grid()
        filename = "Sim " + str(sim_number) + ' ' + label + ' ' + var + " distribution at multiple times.jpg"
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()

def generate_thrust_plot(interface_data_file, max_line = None):
    interface_data_total = FormInterfaceDataFromFile(data_file_name = interface_data_file, max_line = max_line)
    interface_data = interface_data_total.interface_data
    sim_number = interface_data_total.sim_number
    m_dot = interface_data["mass_flux"]
    p_exit = interface_data["p"]
    vel_x_e = interface_data["vel_x"]
    A_e = interface_data["A"]
    interface_data["Thrust"] = m_dot * vel_x_e * A_e + p_exit * A_e

    fig = plt.figure(figsize=(15, 5))
    plt.title("Transient Thrust Profile")
    plt.ylabel("Thrust (N)", rotation = "horizontal", ha = "right")
    plt.xlabel("Time (ms)")
    plt.scatter(interface_data["time"] * 1e3, interface_data["Thrust"], marker = '.')
    plt.grid()
    filename = "Sim"+ str(sim_number) + "ThrustProfile.jpg"
    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    current_dir = os.getcwd()
    plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
    plt.close()

def generate_transient_cell_property_plots(cell_file_name, plot_vars, max_line = None):
    cell_data = FormCellDataFromFile(data_file_name = cell_file_name, max_line = max_line)
    sim_number = cell_data.sim_number
    for var in plot_vars:
        
        pass

def generate_transient_interface_property_plots():
    pass

def generate_spatial_plots_with_multiple_y_axes():
    pass

def generate_transient_cell_plots_with_multiple_y_axes():
    pass

def generate_transient_interface_plots_with_multiple_y_axes():
    pass

def generate_custom_spatial_plots_from_multiple_sims():
    pass

def generate_custom_transient_cell_plots_from_multiple_sims():
    pass

def generate_custom_transient_interface_plots_from_multiple_sims():
    pass

class AnimateData():
    def __init__(self, data_files, slow_down_factor, plot_vars, label) -> None:
        self.dfs = [None] * len(data_files)
        self.time_list = [None] * len(data_files)
        for ind, snapshot in enumerate(data_files):
            for component in snapshot:
                data = GenerateDataObject(data_file_name = component)
                if self.dfs[ind] is None:
                    self.dfs[ind] = data.component_data
                    self.time_list[ind] = data.t_final
                    self.sim_number = data.sim_number
                else:
                    self.dfs[ind] = pd.concat([self.dfs[ind], data.component_data], axis = 0)
            self.dfs[ind] = self.dfs[ind].sort_values(by = ["pos_x"])
        
        self.dt_list = (np.array(self.time_list[1:]) - np.array(self.time_list[:-1])).tolist()

        FFwriter = animation.FFMpegWriter()
        current_dir = os.getcwd()
        for var in plot_vars:
            self.fig, self.ax = plt.subplots(figsize = (15, 5))
            min_val = min([min(df[var]) for df in self.dfs])
            max_val = max([max(df[var]) for df in self.dfs])
            y_min = min_val - 0.1 * (max_val - min_val)
            y_max = max_val + 0.1 * (max_val - min_val)
            print(y_min, y_max)
            print("Inside __init__:", var)
            anim = FuncAnimation(plt.gcf(), self.update, frames = len(self.dfs), interval = 1, \
                                    fargs = [[var, slow_down_factor, y_min, y_max]])
            filename = "Sim" + str(self.sim_number) + var + "AnimationForComponent" + label + ".mp4"
            anim.save(current_dir + "/plots/" + filename, writer = FFwriter)
        
    def update(self, frame, fargs):
        self.ax.clear()
        [var, slow_down_factor, y_min, y_max] = fargs

        df = self.dfs[frame]
        self.ax.grid()
        self.ax.set_ylim([y_min, y_max])
        plt.scatter(df["pos_x"], df[var], marker = ".")
        self.ax.set_xlabel("Position (m)")
        self.ax.set_ylabel(SYMBOLS[var] + " (" + SI_UNITS[var] +")", \
                    rotation = "horizontal", ha = "right")
        self.ax.set_title("Transient Evolution of " + SYMBOLS[var] + " Distribution")
        if frame < len(self.dfs) - 1:
            time.sleep(self.dt_list[frame] * slow_down_factor)
        return plt

def generate_x_t_contour_plots(data_file_names, plot_vars, label):
    dfs = [None] * len(data_file_names)
    time_list = [None] * len(data_file_names)
    
    for ind, snapshot in enumerate(data_file_names):
        for component in snapshot:
            data = GenerateDataObject(data_file_name = component)
            if dfs[ind] is None:
                dfs[ind] = data.component_data
                time_list[ind] = data.t_final
                sim_number = data.sim_number
            else:
                dfs[ind] = pd.concat([dfs[ind], data.component_data], axis = 0)
        dfs[ind] = dfs[ind].sort_values(by = ["pos_x"])
    pos_x = dfs[0]["pos_x"]
    X, T = np.meshgrid(pos_x, time_list)

    for var in plot_vars:
        fig = plt.figure(figsize=(15, 5))
        var_mesh = np.array([df[var].tolist() for df in dfs])
        CS = plt.contourf(X, T, var_mesh, levels = 30)
        cbar = fig.colorbar(CS)
        cbar.ax.set_ylabel("(" + SI_UNITS[var] +")", \
                        rotation = "horizontal", ha = "right")
        plt.title("x-t contour plot of " + SYMBOLS[var])
        plt.xlabel("Position (m)")
        plt.ylabel("Time (s)")
        filename = "Sim" + str(sim_number) + "ContourPlotOf" + var + \
                            "ForComponent" + label + ".jpg"
        mng = plt.get_current_fig_manager()
        mng.full_screen_toggle()
        current_dir = os.getcwd()
        plt.savefig(current_dir + "/plots/" + filename, bbox_inches="tight")
        plt.close()
