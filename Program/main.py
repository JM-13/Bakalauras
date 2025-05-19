import os
import re
import pandas as pd
import numpy as np
import json

import constants as c
from tools.catalogue_data import Catalogue
from tools.extract_values import Extract


def Convert_measurement_type(row):

    rounding = [True, 6]

    if 'R' in row.iloc[0]:
        converted_value = row.iloc[-1] * 0.529177249 #To angstroms
    else:
        converted_value = np.degrees(row.iloc[-1])   #To degrees

    if rounding[0]:
        converted_value = round(converted_value, rounding[1])

    return converted_value


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


Save_coordinates = False
Save_optimized_coordinates = True
To_angstroms_and_degrees = False

for solute in list(Data_files.keys()):

    Save_folder_name = c.Filename_to_Solute[solute]
    Save_folder_path = os.path.join(c.Procesed_Data_Folder, Save_folder_name)
    os.makedirs(Save_folder_path, exist_ok=True)

    for solvent in list(Data_files[solute].keys()):

        print(f"Getting data from files {solute}-{solvent}-...log")

        Save_file_name = c.Filename_to_Solvent[solvent] + '_' + c.Filename_to_Solute[solute] + '.txt'
        Save_file_path = os.path.join(Save_folder_path, Save_file_name)

        Energys_combined = {}
        Coordinates_combined = {}
        Optimized_coordinates_combined = {}

        for file_type in list(Data_files[solute][solvent].keys()):

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
                    Coordinates_combined['Name'] = list(Optimized_Coordinates[-1]['Name'])
                    Coordinates_combined['Definition'] = list(Optimized_Coordinates[-1]['Definition'])
                    Optimized_coordinates_combined['Name'] = list(Optimized_Coordinates[-1]['Name'])
                    Optimized_coordinates_combined['Definition'] = list(Optimized_Coordinates[-1]['Definition'])

                relavent_coordinates = Coordinates[-1]
                relavent_optimized_coordinates = Optimized_Coordinates[-1]

                if To_angstroms_and_degrees:
                    relavent_coordinates['New X'] = relavent_coordinates.apply(Convert_measurement_type, axis=1)

                Coordinates_combined[file_type] = list(relavent_coordinates['New X'])
                Optimized_coordinates_combined[file_type] = list(relavent_optimized_coordinates['Value'])

        Energys_combined = pd.DataFrame([Energys_combined])[['optS0', 'tdS0', 'optS1', 'optR1']]
        Coordinates_combined = pd.DataFrame(Coordinates_combined)[['Name', 'Definition', 'optS0', 'optS1', 'optR1']]
        Optimized_coordinates_combined = pd.DataFrame(Optimized_coordinates_combined)[['Name', 'Definition', 'optS0', 'optS1', 'optR1']]

        with open(Save_file_path, 'w') as f:

            for file_type in list(Energys_combined.keys()):
                f.write(f"{file_type} = {Energys_combined[file_type].iloc[0]}\n")

            if Save_coordinates:
                f.write("\n")
                f.write(Coordinates_combined.to_string(index=False))
                f.write("\n")

            if Save_optimized_coordinates:
                f.write("\n")
                f.write(Optimized_coordinates_combined.to_string(index=False))
                f.write("\n")

        print(f'Written data to file: {Save_file_name} in folder: {Save_folder_path}\n')
