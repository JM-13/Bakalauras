import os
import re
import pandas as pd
import numpy as np
import json

import constants as c
from tools.catalogue_data import Catalogue
from tools.extract_values import Extract
from tools.angstroms_and_degrees import Convert_measurement_type


Working_dir = os.getenv("PWD", os.getcwd())
Data_folder = os.path.join(Working_dir, 'Data')
os.makedirs(Data_folder, exist_ok=True)


D_endings = ['optS0', 'tdS0', 'optS1', 'optR1']
S_endings = ['fa', 'fb']

Catalogue(folders = [c.S0_folder, c.S1_folder],
          filename_endings = D_endings,
          save_as_json = True,
          save_location = Working_dir).Data_Files()

Catalogue(folders = [c.Scan_folder],
          filename_endings = S_endings,
          save_as_json = True,
          save_location = Working_dir).Scan_Files()


with open("Data_files.json", "r") as file:
    Data_files = json.load(file)


Save_coordinates = True
Save_optimized_coordinates = False
To_angstroms_and_degrees = False

for solute in Data_files:

    Save_folder_name = c.Filename_to_Solute[solute]
    Save_folder_path = os.path.join(Data_folder, Save_folder_name)
    os.makedirs(Save_folder_path, exist_ok=True)

    for solvent in Data_files[solute]:

        print(f"Getting data from files {solute}-{solvent}-...log")

        Save_file_name = f"{c.Filename_to_Solvent[solvent]}_{c.Filename_to_Solute[solute]}_DATA.txt"
        Save_file_path = os.path.join(Save_folder_path, Save_file_name)

        Energys_combined = {}
        Coordinates_combined = {}
        Optimized_coordinates_combined = {}

        for file_type in Data_files[solute][solvent]:

            Extractor = Extract(Data_files[solute][solvent][file_type])

            if file_type == 'optS0':
                Energys = Extractor.S0()
                Coordinates = Extractor.RAD()
                Optimized_Coordinates = Extractor.Optimized_RAD()

            elif file_type == 'tdS0':
                Energys = Extractor.S1()

            elif file_type == 'optS1':
                Energys = Extractor.S1()
                Coordinates = Extractor.RAD()
                Optimized_Coordinates = Extractor.Optimized_RAD()

            elif file_type == 'optR1':
                Energys = Extractor.S1()
                Coordinates = Extractor.RAD()
                Optimized_Coordinates = Extractor.Optimized_RAD()

            Energys_combined[file_type] = Energys[-1]

            if file_type != 'tdS0':

                if len(Coordinates_combined) == 0:
                    Coordinates_combined['Name']                 = list(Optimized_Coordinates[-1]['Name'])
                    Coordinates_combined['Definition']           = list(Optimized_Coordinates[-1]['Definition'])
                    Optimized_coordinates_combined['Name']       = list(Optimized_Coordinates[-1]['Name'])
                    Optimized_coordinates_combined['Definition'] = list(Optimized_Coordinates[-1]['Definition'])

                relavent_coordinates = Coordinates[-1]
                relavent_optimized_coordinates = Optimized_Coordinates[-1]

                if To_angstroms_and_degrees:
                    relavent_coordinates['New X'] = relavent_coordinates.apply(Convert_measurement_type, axis=1)

                Coordinates_combined[file_type]           = list(relavent_coordinates['New X'])
                Optimized_coordinates_combined[file_type] = list(relavent_optimized_coordinates['Value'])

        Energys_combined               = pd.DataFrame([Energys_combined])[['optS0', 'tdS0', 'optS1', 'optR1']]
        Coordinates_combined           = pd.DataFrame(Coordinates_combined)[['Name', 'Definition', 'optS0', 'optS1', 'optR1']]
        Optimized_coordinates_combined = pd.DataFrame(Optimized_coordinates_combined)[['Name', 'Definition', 'optS0', 'optS1', 'optR1']]

        with open(Save_file_path, 'w') as f:

            for file_type in Energys_combined:
                f.write(f"{file_type} = {Energys_combined[file_type].iloc[0]}\n")

            if Save_coordinates:
                f.write("\n")
                f.write(Coordinates_combined.to_string(index=False))
                f.write("\n")

            if Save_optimized_coordinates:
                f.write("\n")
                f.write(Optimized_coordinates_combined.to_string(index=False))
                f.write("\n")

        # print(f'Written data to file: {Save_file_name} in folder: {Save_folder_path}\n')





with open("Scan_files.json", "r") as file:
    Data_files = json.load(file)


for solute in Data_files:

    Save_folder_name = c.Filename_to_Solute[solute]
    Save_folder_path = os.path.join(Data_folder, Save_folder_name)
    os.makedirs(Save_folder_path, exist_ok=True)

    for solvent in Data_files[solute]:

        print(f"Getting data from Scan files {solute}-{solvent}-...log")

        Save_file_name = f"{c.Filename_to_Solvent[solvent]}_{c.Filename_to_Solute[solute]}_SCAN.txt"
        Save_file_path = os.path.join(Save_folder_path, Save_file_name)

        Data_combined = {}


        for direction in Data_files[solute][solvent]:

            Extractor = Extract(Data_files[solute][solvent][direction])

            Fixed_val = Extractor.Scan_fixed()
            Fixed_val = ', '.join(str(item) for item in Fixed_val)

            Energies = Extractor.Scan_Optimized_Energy()
            RAD_fixed_vals = Extractor.Scan_RAD_fixed_values()
            RAD_averages = Extractor.Scan_RAD_fixed_average(return_fixed_opposit=True)

            if direction == 'fa':

                Data_dataframe = pd.DataFrame({'S0':Energies[0],
                                               'S1':Energies[1],
                                               'D_fixed':RAD_fixed_vals,
                                               'D_average':RAD_averages[1]})

                Base_S0 = Data_dataframe['S0'].min()
                Data_dataframe['S0'] = Data_dataframe['S0'] - Base_S0
                Data_dataframe['S1'] = Data_dataframe['S1'] - Base_S0

                Data_combined['fa'] = {'Fixed':[Fixed_val, RAD_averages[0]], 'Data':Data_dataframe}

            elif direction == 'fb':

                Data_dataframe = pd.DataFrame({'S0':reversed(Energies[0]),
                                               'S1':reversed(Energies[1]),
                                               'D_fixed':reversed(RAD_fixed_vals),
                                               'D_average':reversed(RAD_averages[1])})

                Base_S0 = Data_dataframe['S0'].min()
                Data_dataframe['S0'] = Data_dataframe['S0'] - Base_S0
                Data_dataframe['S1'] = Data_dataframe['S1'] - Base_S0

                Data_combined['fb'] = {'Fixed':[Fixed_val, RAD_averages[0]], 'Data':Data_dataframe}


        with open(Save_file_path, 'w') as f:
            f.write("From fa scan:\n")
            f.write(f"{'Fixed'}   = {Data_combined['fa']['Fixed'][0]}\n")
            f.write(f"{'Opposit'} = {Data_combined['fa']['Fixed'][1]}\n\n")
            f.write(Data_combined['fa']['Data'].to_string(index=False))

            f.write("\n\n")

            f.write("From fb scan:\n")
            f.write(f"{'Fixed'}   = {Data_combined['fb']['Fixed'][0]}\n")
            f.write(f"{'Opposit'} = {Data_combined['fb']['Fixed'][1]}\n\n")
            f.write(Data_combined['fb']['Data'].to_string(index=False))

        # print(f'Written data to file: {Save_file_name} in folder: {Save_folder_path}\n')

processed_data_folder = []
for i in list(c.Filename_to_Solute.values()):
    folder_path = os.path.join(Data_folder, i)
    processed_data_folder.append(folder_path)

Catalogue(folders = processed_data_folder,
          save_as_json = True,
          save_location = Data_folder).Processed_Data_Files([1,2,0])

print(f'\nAll data saved to folder {Data_folder} sorted into folders based on used solute')


