import os
import json
import pandas as pd

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

                    file_path = os.path.join(folder, file)

                    parts = self.Clean_file_name(file, clean_name)

                    file_type = parts[-1]
                    solute = '-'.join(parts[:solution_end])
                    solvent = parts[solution_end]

                    if solute not in list(Data.keys()):
                        Data[solute] = {}

                    if solvent not in list(Data[solute].keys()):
                        Data[solute][solvent] = {}

                    Data[solute][solvent][file_type] = file_path

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







