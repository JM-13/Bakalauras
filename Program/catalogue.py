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


#ppp-ome-in-meoh-S1-fscan-fa.log
#bdp-pnphen-dmfm-S1-fscan-fa-cont_1x32c.log
#bdp-phenyl-acet-S1-fscan-fa-cont.log
#ppp-ome-in-acet-S1-fscan-fb.log

#['optS0', 'tdS0', 'optS1', 'optR1', 'fa', 'fb']

class Catalogue:
    def __init__(self, folders, filename_endings, save_as_json=False, save_location=''):

        self.folders = folders
        self.filename_endings = filename_endings
        self.save_as_json = save_as_json

        if save_location:
            self.save_location = save_location
        else:
            self.save_location = os.getenv("PWD", os.getcwd())


    def Clean_file_name(self, file, clean):

        filename_endings = self.filename_endings

        filename = file.removesuffix('.log')
        parts_raw = filename.split('-')

        if not clean:
            return parts_raw

        parts = []
        done = False
        for i in parts_raw:
            for fe in filename_endings:
                if i == fe:
                    done = True

            parts.append(i)
            if done:
                break

        return parts


    def Files(self, clean_name, solution_end, json_name):

        folders = self.folders
        save = self.save_as_json
        save_folder = self.save_location

        Data = {}

        for folder in folders:
            for file in os.listdir(folder):
                if file.endswith('.log'):

                    parts = self.Clean_file_name(file, clean_name)

                    file_type = parts[-1]
                    solute = '-'.join(parts[:solution_end])
                    solvent = parts[solution_end]

                    if file_type not in list(Data.keys()):
                        Data[file_type] = {}

                    if solute not in list(Data[file_type].keys()):
                        Data[file_type][solute] = {}

                    Data[file_type][solute][solvent] = file

        if save:
            file_location = os.path.join(save_folder, json_name)
            with open(file_location, "w") as file:
                json.dump(Data, file, indent=4)

        return Data


    def Data_Files(self, clean_name = False):

        solution_end = -2
        json_name = "Data_files.json"

        Data = self.Files(clean_name, solution_end, json_name)

        return Data


    def Scan_Files(self, clean_name = True):

        solution_end = -4
        json_name = "Scan_files.json"

        Data = self.Files(clean_name, solution_end, json_name)

        return Data







