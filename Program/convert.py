import os
import json
import re
import pandas as pd


# S0_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/1-S0'
# S1_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/2-S1'
# Scan_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/3-scans'
#
# Solvent_in_filename = {'Acetone':'-acet',
#                        'DiMethylFormamide':'-dmfm',
#                        'DiMethylSulfoxide':'-dmso',
#                        'Methanol':'-meoh',
#                        'TetraHydroFuran':'-thfu',
#                        'Toluene':'-tolu'}
#
# Material_in_filename = {'BDP-nitrophenyl':'bdp-pnphen',
#                         'BDP-phenyl':'bdp-phenyl',
#                         'BDP-PP-phenyl':'ppp',
#                         'BDP-PP-phenyl_OMes_backward':'ppp-ome-out',
#                         'BDP-PP-phenyl_OMes_forward':'ppp-ome-in'}
#
# File_types = {'regular':['-optS0', '-tdS0'],
#               'scan':['-optS1', '-optR1']}

#bdp-phenyl-thfu-optS0.log

def Locate_Data_Files(folders, save_as_json=False):

    Data = {}

    for folder in folders:
        for file in os.listdir(folder):
            if file.endswith('.log'):
                file_path = os.path.join(folder, file)
                file_name = file.removesuffix('.log')

                parts = file_name.split('-')

                data_type = parts[-1]
                solvent_name = parts[-2]
                material_name = '-'.join(parts[:-2])

                if data_type not in list(Data.keys()):
                    Data[data_type] = {}

                if material_name not in list(Data[data_type].keys()):
                    Data[data_type][material_name] = {}

                Data[data_type][material_name][solvent_name] = file

    if save_as_json:
        with open("Data_files.json", "w") as file:
            json.dump(Data, file, indent=4)
        return Data
    else:
        return Data


#ppp-ome-in-meoh-S1-fscan-fa.log
#bdp-pnphen-dmfm-S1-fscan-fa-cont_1x32c.log
#bdp-phenyl-acet-S1-fscan-fa-cont.log
#ppp-ome-in-acet-S1-fscan-fb.log

def Locate_Scan_Files(folders, save_as_json=False):

    Data = {}

    for folder in folders:
        for file in os.listdir(folder):
            if file.endswith('.log'):
                file_path = os.path.join(folder, file)
                file_name = file.removesuffix('.log')

                parts_raw = file_name.split('-')

                parts = []
                for i in parts_raw:
                    if i =='fa' or i =='fb':
                        parts.append(i)
                        break
                    else:
                        parts.append(i)

                scan_direction = parts[-1]
                solvent_name = parts[-4]
                material_name = '-'.join(parts[:-4])

                if scan_direction not in list(Data.keys()):
                    Data[scan_direction] = {}

                if material_name not in list(Data[scan_direction].keys()):
                    Data[scan_direction][material_name] = {}

                Data[scan_direction][material_name][solvent_name] = file

    if save_as_json:
        with open("Scan_files.json", "w") as file:
            json.dump(Data, file, indent=4)
        return Data
    else:
        return Data

#
#
#
#             for material_name, material_abriviation in c.Material_to_fname.items():
#                 Material_folder = os.path.join(Data_folder, material_name)
#                 os.makedirs(Material_folder, exist_ok=True)
#
#                 if material_abriviation == 'ppp':
#                     carefull = ['ppp-ome-out', 'ppp-ome-in']
#                 else:
#                     carefull = ['99999999', '99999999']
#
#                 for material_solvent in Refined_Data.keys():
#
#                     if material_abriviation in material_solvent and not any(x in material_solvent for x in carefull):
#                         for solvent_name, solvent_abriviation in c.Solvent_to_fname.items():
#
#                             if solvent_abriviation in material_solvent:
#
#
#
#
#
# class Convert:
#     def __init__(self):
#         self.data = []





