import os
import re
import pandas as pd
import numpy as np
import json

import constants as c
from tools.catalogue_data import Catalogue
from tools.extract_values import Extract
from tools.analyze import Analyze


Working_dir = os.getenv("PWD", os.getcwd())
Data_folder = os.path.join(Working_dir, 'Data')
os.makedirs(Data_folder, exist_ok=True)

D_endings = ['optS0', 'tdS0', 'optS1', 'optR1']
S_endings = ['fa', 'fb']


def user_input_check(input_text):
    print(f"")
    user_choice = None
    while not user_choice:
        user_choice = input(f"{input_text} yes/no\n")
        if user_choice == "yes":
            return True
        elif user_choice == "no":
            return False
        else:
            print(f"{user_choice} is not a yes/no answer, please try again")
            user_choice = None


to_extract_data = user_input_check('Do you want to extract data from the raw .log files?')
if to_extract_data:
    print(f'\nCataloging .log data in {c.S0_folder}, {c.S1_folder} folders')
    Catalogue(folders = [c.S0_folder, c.S1_folder],
            filename_endings = D_endings,
            save_as_json = True,
            save_location = Working_dir).Data_Files()

    with open("Data_files.json", "r") as file:
        Data_files = json.load(file)


    print(f'\nExtracting data from .log files:')
    for solute in Data_files:

        Save_folder_name = c.Filename_to_Solute[solute]
        Save_folder_path = os.path.join(Data_folder, Save_folder_name)
        os.makedirs(Save_folder_path, exist_ok=True)

        for solvent in Data_files[solute]:

            print(f"extracting {solute}-{solvent} file data")

            Save_file_name = f"{c.Filename_to_Solvent[solvent]}_{c.Filename_to_Solute[solute]}_DATA.txt"
            Save_file_path = os.path.join(Save_folder_path, Save_file_name)

            Energies_combined = {}
            Coordinates_combined = {}
            Optimized_coordinates_combined = {}

            for file_type in Data_files[solute][solvent]:

                Extractor = Extract(Data_files[solute][solvent][file_type])

                if file_type == 'optS0':
                    Energies = Extractor.S0()
                    Coordinates = Extractor.RAD()
                    Optimized_Coordinates = Extractor.Optimized_RAD()

                elif file_type == 'tdS0':
                    Energies = Extractor.S1()

                elif file_type == 'optS1':
                    Energies = Extractor.S1()
                    Coordinates = Extractor.RAD()
                    Optimized_Coordinates = Extractor.Optimized_RAD()

                elif file_type == 'optR1':
                    Energies = Extractor.S1()
                    Coordinates = Extractor.RAD()
                    Optimized_Coordinates = Extractor.Optimized_RAD()

                Energies_combined[file_type] = Energies[-1]

                if file_type != 'tdS0':

                    if len(Coordinates_combined) == 0:
                        Coordinates_combined['Name']                 = list(Optimized_Coordinates[-1]['Name'])
                        Coordinates_combined['Definition']           = list(Optimized_Coordinates[-1]['Definition'])
                        Optimized_coordinates_combined['Name']       = list(Optimized_Coordinates[-1]['Name'])
                        Optimized_coordinates_combined['Definition'] = list(Optimized_Coordinates[-1]['Definition'])

                    relavent_coordinates = Coordinates[-1]
                    relavent_optimized_coordinates = Optimized_Coordinates[-1]

                    Coordinates_combined[file_type]           = list(relavent_coordinates['New X'])
                    Optimized_coordinates_combined[file_type] = list(relavent_optimized_coordinates['Value'])

            Energies_combined              = pd.DataFrame([Energies_combined])[['optS0', 'tdS0', 'optS1', 'optR1']]
            Coordinates_combined           = pd.DataFrame(Coordinates_combined)[['Name', 'Definition', 'optS0', 'optS1', 'optR1']]
            Optimized_coordinates_combined = pd.DataFrame(Optimized_coordinates_combined)[['Name', 'Definition', 'optS0', 'optS1', 'optR1']]

            with open(Save_file_path, 'w') as f:

                for file_type in Energies_combined:
                    f.write(f"{file_type} = {Energies_combined[file_type].iloc[0]}\n")

                if c.Save_optimized_coordinates:
                    f.write("\n")
                    f.write(Optimized_coordinates_combined.to_string(index=False))
                    f.write("\n")

                else:
                    f.write("\n")
                    f.write(Coordinates_combined.to_string(index=False))
                    f.write("\n")

            print(f'written data to {Save_file_path}')


to_extract_scan = user_input_check('Do you want to extract data from the raw scan .log files?')
if to_extract_scan:
    print(f'\nCataloging scan .log data in {c.Scan_folder} folder')
    Catalogue(folders = [c.Scan_folder],
            filename_endings = S_endings,
            save_as_json = True,
            save_location = Working_dir).Scan_Files()

    with open("Scan_files.json", "r") as file:
        Data_files = json.load(file)

    print(f'\nExtracting data from scan .log files:')
    for solute in Data_files:

        Save_folder_name = c.Filename_to_Solute[solute]
        Save_folder_path = os.path.join(Data_folder, Save_folder_name)
        os.makedirs(Save_folder_path, exist_ok=True)

        for solvent in Data_files[solute]:

            Save_file_name = f"{c.Filename_to_Solvent[solvent]}_{c.Filename_to_Solute[solute]}_SCAN.txt"
            Save_file_path = os.path.join(Save_folder_path, Save_file_name)

            Data_combined = {}

            for direction in Data_files[solute][solvent]:

                print(f"extracting {solute}-{solvent}-S1-fscan-{direction} file data")

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

            print(f'written data to {Save_file_path}')


if to_extract_data or to_extract_scan:
    print(f'\nCataloging processed data in {Data_folder} folder')
    processed_data_folders = []
    for i in list(c.Filename_to_Solute.values()):
        folder_path = os.path.join(Data_folder, i)
        processed_data_folders.append(folder_path)

    Catalogue(folders = processed_data_folders,
            save_as_json = True,
            save_location = Data_folder).Processed_Data_Files([1,2,0])

    print(f'\nAll data has been extracted and saved to folder {Data_folder}')


to_analyze = user_input_check("Do you want to analyze the extracted data?")
if to_analyze:
    analasys = Analyze(Data_folder,
                       c.difference_function,
                       c.Solute_to_shorten.keys(),
                       c.Solvent_to_shorten.keys(),
                       c.Solute_to_shorten.values(),
                       c.Solvent_to_shorten.values())

    to_use_angstroms = user_input_check("Convert distances to angstroms for the analasys?")
    to_use_degrees   = user_input_check("Convert angles to degrees for the analasys?")

    to_calc_solvent_diff = user_input_check("Do you want to calculate differences between solvents?")
    if to_calc_solvent_diff:
        analasys.solvent_differences(use_angstroms=to_use_angstroms, use_degrees=to_use_degrees)
        print(f'\tDifferences calculated')

        to_display_solvent_diff = user_input_check("Display solvent differences?")
        if to_display_solvent_diff:
            analasys.display_solvent_differences()

        to_generate_latex = user_input_check("Do you want to save the differences to a pdf file?")
        if to_generate_latex:
            file_name = input("Please input a filename with no file type (if no name is given a default name will be used)")
            analasys.generate_latex_results_document(file_name=file_name, differences="Solvent", use_solvent_by_solute=False)

        to_display_solvent_diff_by_solute = user_input_check("Display solvent differences by solute?")
        if to_display_solvent_diff_by_solute:
            analasys.display_solvent_differences_by_solute()

        to_generate_latex = user_input_check("Do you want to save the differences to a pdf file?")
        if to_generate_latex:
            file_name = input("Please input a filename with no file type (if no name is given a default name will be used)")
            analasys.generate_latex_results_document(file_name=file_name, differences="Solvent", use_solvent_by_solute=True)

    to_calc_solute_diff = user_input_check("Do you want to calculate differences between solutes?")
    if to_calc_solute_diff:
        bdp_central = input(f"Pick which central atoms to compare:\n{"A"} - All 21 atoms\n{"H"} - Exclude replaced Hydrogen atoms\n{"NoH"} - No Hydrogen atoms\n")
        analasys.solute_differences(use_angstroms=to_use_angstroms, use_degrees=to_use_degrees, bdp_central=bdp_central)
        print(f'\tDifferences calculated')

        to_display_solute_diff = user_input_check("Display solute differences?")
        if to_display_solute_diff:
            analasys.display_solute_differences()

        to_generate_latex = user_input_check("Do you want to save the differences to a pdf file?")
        if to_generate_latex:
            file_name = input("Please input a filename with no file type (if no name is given a default name will be used)")
            analasys.generate_latex_results_document(file_name=file_name, differences="Solute")








