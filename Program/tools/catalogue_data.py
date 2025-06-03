import os
import json

class Catalogue:
    def __init__(self, folders, filename_endings=[], save_as_json=False, save_location=''):

        self.folders = folders
        self.filename_endings = filename_endings
        self.save_as_json = save_as_json

        if save_location:
            self.save_location = save_location
        else:
            self.save_location = os.getenv("PWD", os.getcwd())


    def Clean_file_name(self, file, clean, suffix, sep):

        filename_endings = self.filename_endings

        filename = file.removesuffix(suffix)
        parts_raw = filename.split(sep)

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


    def Files(self, clean_name, solution_end, json_name, cat_order, suffix, sep, exact_locations=[]):

        #cat_order=[0 == file_type, 1 == solute, 2 == solvent] and the order will define the structure order
        #exact_locations=[file_type, solvent]

        folders = self.folders
        save = self.save_as_json
        save_folder = self.save_location

        Data = {}

        for folder in folders:
            for file in os.listdir(folder):
                if file.endswith(suffix):

                    file_path = os.path.join(folder, file)

                    parts = self.Clean_file_name(file, clean_name, suffix, sep)

                    if not exact_locations:
                        file_type = parts[-1]
                        solute = sep.join(parts[:solution_end])
                        solvent = parts[solution_end]

                    else:
                        file_type = parts[exact_locations[0]]
                        solute    = sep.join(parts[(exact_locations[1]+1):exact_locations[0]])
                        solvent   = parts[exact_locations[1]]


                    catagories = [file_type,
                                  solute,
                                  solvent]

                    cats = [catagories[cat_order[0]],
                            catagories[cat_order[1]],
                            catagories[cat_order[2]]]

                    if cats[0] not in Data:
                        Data[cats[0]] = {}

                    if cats[1] not in Data[cats[0]]:
                        Data[cats[0]][cats[1]] = {}

                    Data[cats[0]][cats[1]][cats[2]] = file_path

        if save:
            file_location = os.path.join(save_folder, json_name)
            with open(file_location, "w") as file:
                json.dump(Data, file, indent=4)

        return Data


    def Data_Files(self, cat_order=[1,2,0], clean_name = False, suffix='.log', sep='-'):

        solution_end = -2
        json_name = "Data_files.json"

        Data = self.Files(clean_name, solution_end, json_name, cat_order, suffix, sep)

        return Data


    def Scan_Files(self, cat_order=[1,2,0], clean_name = True, suffix='.log', sep='-'):

        solution_end = -4
        json_name = "Scan_files.json"

        Data = self.Files(clean_name, solution_end, json_name, cat_order, suffix, sep)

        return Data

    def Processed_Data_Files(self, cat_order, clean_name = True, suffix='.txt', sep='_'):

        solution_end = None
        json_name = "Processed_Data_files.json"

        Data = self.Files(clean_name, solution_end, json_name, cat_order, suffix, sep, exact_locations=[-1, 0])

        return Data







