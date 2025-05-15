import os
import pandas as pd
import io

class Retvieve:
    def __init__(self, foldername):
        self.material_name = os.path.basename(foldername)

        self.Data = {'optS0':{f'{self.material_name}':{}},
                     'tdS0': {f'{self.material_name}':{}},
                     'optS1':{f'{self.material_name}':{}},
                     'optR1':{f'{self.material_name}':{}},
                     'RAD':  {f'{self.material_name}':{}}
                    }

        self.Scan_Data = {f'{self.material_name}':{}}

        self.regular_file_names = []
        self.scan_file_names = []

        for file in os.listdir(foldername):
            full_file_path = os.path.join(foldername, file)
            if 'SCAN' not in file:
                self.regular_file_names.append(full_file_path)
            else:
                self.scan_file_names.append(full_file_path)


    def Regular_data(self):
        for file_path in self.regular_file_names:
            solvent_name = os.path.basename(file_path).split("_")[0]
            with open(file_path, "r") as file:
                for _ in range(4):
                    line = file.readline().strip()
                    key, value = line.split("=")
                    self.Data[key.strip()][self.material_name][solvent_name] = float(value.strip())

                file.readline()

                df = pd.read_csv(file, sep='\s+').convert_dtypes()

            self.Data['RAD'][self.material_name][solvent_name] = df

        return self.Data

    def Scan_data(self):
        for file_path in self.scan_file_names:
            solvent_name = os.path.basename(file_path).split("_")[1]
            self.Scan_Data[self.material_name][solvent_name] = {}

            with open(file_path, "r") as file:
                sections = file.read().strip().split("\n\n")

            df_fa = pd.read_csv(io.StringIO(sections[1]), sep='\s+').convert_dtypes()
            self.Scan_Data[self.material_name][solvent_name]['fa'] = df_fa

            df_fb = pd.read_csv(io.StringIO(sections[3]), sep='\s+').convert_dtypes()
            self.Scan_Data[self.material_name][solvent_name]['fb'] = df_fb

        return self.Scan_Data


