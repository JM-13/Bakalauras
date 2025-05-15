import os
import re
import pandas as pd
import numpy as np
import constants as c
from extract import Extract

Suffixes = [['-optS0.log', '-tdS0.log'], ['-optS1.log', '-optR1.log']]
Data = {}

for folder, sufixes in zip([c.S0_folder, c.S1_folder], Suffixes):
    for file in os.listdir(folder):
        for sufix in sufixes:
            if file.endswith(sufix):
                file_path = os.path.join(folder, file)
                rest_of_name = file.removesuffix(sufix)
                key_value = sufix.removesuffix('.log')

                Extractor = Extract(file_path)

                if rest_of_name not in Data:
                    Data[rest_of_name] = {}

                if sufix == Suffixes[0][0]:
                    Energys = Extractor.S0()
                else:
                    Energys = Extractor.S1()

                if sufix != Suffixes[0][1]:
                    Coordinates = Extractor.RAD()
                    Optimized_Coordinates = Extractor.Optimized_RAD()
                else:
                    Coordinates = [None, None]
                    Optimized_Coordinates = [None, None]

                Data[rest_of_name][key_value] = {'Energys':Energys, 'Coordinates':Coordinates, 'Optimized_Coordinates':Optimized_Coordinates}
                print(f'Gotten data from file: {file} in folder: {folder}')
                break

print("\n")

def Convert_measurement_type(row):

    rounding = [True, 6]

    if 'R' in row.iloc[0]:
        converted_value = row.iloc[-1] * 0.529177249 #To angstroms
    else:
        converted_value = np.degrees(row.iloc[-1])   #To degrees

    if rounding[0]:
        converted_value = round(converted_value, rounding[1])

    return converted_value


To_angstroms_and_degrees = False
Refined_Data = {}
for file_name in Data.keys():
    if file_name not in Refined_Data:
        Refined_Data[file_name] = {}

    Coordinates_combined = {}
    Optimized_coordinates_combined = {}

    for file_type in Data[file_name].keys():

        energy = Data[file_name][file_type]['Energys'][-1]
        coordinate = Data[file_name][file_type]['Coordinates'][-1]
        optimized_coordinate = Data[file_name][file_type]['Optimized_Coordinates'][-1]

        Refined_Data[file_name][file_type] = energy

        if file_type != '-tdS0':

            if 'Variable' not in Coordinates_combined:
                Coordinates_combined['Variable'] = coordinate['Variable']

            if To_angstroms_and_degrees:
                coordinate['New X'] = coordinate.apply(Convert_measurement_type, axis=1)

            Coordinates_combined[f'{file_type}_Value'] = coordinate['New X']


            if 'Name' or 'Definition' not in Optimized_coordinates_combined:
                Optimized_coordinates_combined['Name'] = optimized_coordinate['Name']
                Optimized_coordinates_combined['Definition'] = optimized_coordinate['Definition']

            Optimized_coordinates_combined[f'{file_type}_Value'] = optimized_coordinate['Value']


    Refined_Data[file_name]['Coordinates'] = pd.DataFrame(Coordinates_combined)[['Variable', '-optS0_Value', '-optS1_Value', '-optR1_Value']]
    Refined_Data[file_name]['Optimized_Coordinates'] = pd.DataFrame(Optimized_coordinates_combined)[['Name', 'Definition', '-optS0_Value', '-optS1_Value', '-optR1_Value']]


Working_dir = os.getenv("PWD", os.getcwd())
Data_folder = os.path.join(Working_dir, 'Data')
os.makedirs(Data_folder, exist_ok=True)


Save_coordinates = False
Save_optimized_coordinates = True
for material_name, material_abriviation in c.Material_to_fname.items():
    Material_folder = os.path.join(Data_folder, material_name)
    os.makedirs(Material_folder, exist_ok=True)

    if material_abriviation == 'ppp':
        carefull = ['ppp-ome-out', 'ppp-ome-in']
    else:
        carefull = ['99999999', '99999999']

    for material_solvent in Refined_Data.keys():

        if material_abriviation in material_solvent and not any(x in material_solvent for x in carefull):
            for solvent_name, solvent_abriviation in c.Solvent_to_fname.items():

                if solvent_abriviation in material_solvent:

                    save_file_name = f'{solvent_name}_{material_name}.txt'
                    save_file_path = os.path.join(Material_folder, save_file_name)

                    with open(save_file_path, 'w') as f:
                        f.write(f"{'optS0'} = {Refined_Data[material_solvent]['-optS0']}\n")
                        f.write(f"{'tdS0'}  = {Refined_Data[material_solvent]['-tdS0']}\n")
                        f.write(f"{'optS1'} = {Refined_Data[material_solvent]['-optS1']}\n")
                        f.write(f"{'optR1'} = {Refined_Data[material_solvent]['-optR1']}\n")

                        if Save_coordinates:
                            f.write("\n")
                            f.write(Refined_Data[material_solvent]['Coordinates'].to_string(index=False))
                            f.write("\n")

                        if Save_optimized_coordinates:
                            f.write("\n")
                            f.write(Refined_Data[material_solvent]['Optimized_Coordinates'].to_string(index=False))
                            f.write("\n")

                    print(f'Written data to file: {save_file_name} in folder: {Material_folder}')
                    break


Scan_Data = {}
Scan_direction = ['fscan-fa', 'fscan-fb']
track = 0
for file in os.listdir(c.Scan_folder):
    if file.endswith('.log'):
        print(f'Extracting data from file: {file}')
        file_path = os.path.join(c.Scan_folder, file)

        Extractor = Extract(file_path)

        Fixed_val = Extractor.Scan_fixed()
        Energies = Extractor.Scan_Optimized_Energy()
        RAD_fixed_vals = Extractor.Scan_RAD_fixed_values()
        RAD_averages = Extractor.Scan_RAD_fixed_average(return_fixed_opposit=True)



        for direction in Scan_direction:
            if direction in file:
                for material_name, material_abriviation in c.Material_to_fname.items():

                    if material_abriviation == 'ppp':
                        carefull = ['ppp-ome-out', 'ppp-ome-in']
                    else:
                        carefull = ['99999999', '99999999']

                    if material_abriviation in file and not any(x in file for x in carefull):
                        for solvent_name, solvent_abriviation in c.Solvent_to_fname.items():
                            if solvent_abriviation in file:
                                track+=1

                                if material_name not in Scan_Data:
                                    Scan_Data[material_name] = {}
                                if solvent_name not in Scan_Data[material_name]:
                                    Scan_Data[material_name][solvent_name] = {}
                                if direction not in Scan_Data[material_name][solvent_name]:
                                    Scan_Data[material_name][solvent_name][direction] = {}

                                if direction == Scan_direction[0]:
                                    Dat = {'S0':Energies[0], 'S1':Energies[1], 'D_fixed':RAD_fixed_vals, 'D_average':RAD_averages[1]}
                                else:
                                    Dat = {'S0':reversed(Energies[0]), 'S1':reversed(Energies[1]), 'D_fixed':reversed(RAD_fixed_vals), 'D_average':reversed(RAD_averages[1])}

                                df = pd.DataFrame(Dat)
                                min_S0 = df['S0'].min()
                                df['S0'] = df['S0'] - min_S0
                                df['S1'] = df['S1'] - min_S0

                                Scan_Data[material_name][solvent_name][direction] = {'Fixed':[Fixed_val, RAD_averages[0]], 'Data':df}

                                break
                        break
                break
print("\n")
Working_dir = os.getenv("PWD", os.getcwd())
Data_folder = os.path.join(Working_dir, 'Data')
os.makedirs(Data_folder, exist_ok=True)
for material_name in Scan_Data.keys():
    Material_folder = os.path.join(Data_folder, material_name)
    os.makedirs(Material_folder, exist_ok=True)

    for solvent_name in Scan_Data[material_name]:
        save_file_name = f'SCAN_{solvent_name}_{material_name}.txt'
        save_file_path = os.path.join(Material_folder, save_file_name)

        with open(save_file_path, 'w') as f:
            f.write("From fa scan:\n")
            f.write(f"{'Fixed'}   = {Scan_Data[material_name][solvent_name]['fscan-fa']['Fixed'][0]}\n")
            f.write(f"{'Opposit'} = {Scan_Data[material_name][solvent_name]['fscan-fa']['Fixed'][1]}\n\n")
            f.write(Scan_Data[material_name][solvent_name]['fscan-fa']['Data'].to_string(index=False))

            f.write("\n\n")

            f.write("From fb scan:\n")
            f.write(f"{'Fixed'}   = {Scan_Data[material_name][solvent_name]['fscan-fb']['Fixed'][0]}\n")
            f.write(f"{'Opposit'} = {Scan_Data[material_name][solvent_name]['fscan-fb']['Fixed'][1]}\n\n")
            f.write(Scan_Data[material_name][solvent_name]['fscan-fb']['Data'].to_string(index=False))

        print(f'Written data to file: {save_file_name} in folder: {Material_folder}')





