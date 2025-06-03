import os
import pandas as pd
import io

class Retvieve:
    def __init__(self, filepath):

        self.filepath = filepath
        self.Data = {'Energys':{},
                     'Coordinates':None}
        self.Scan_Data = {'fa':None,
                          'fb':None}


    def Regular_data(self):
        with open(self.filepath, "r") as file:
            for _ in range(4):
                line = file.readline().strip()
                key, value = line.split("=")
                self.Data['Energys'][key.strip()] = float(value.strip())

            file.readline()

            df = pd.read_csv(file, sep='\s+').convert_dtypes()

        self.Data['Coordinates'] = df

        return self.Data

    def Scan_data(self):

        with open(filepath, "r") as file:
            sections = file.read().strip().split("\n\n")

        df_fa = pd.read_csv(io.StringIO(sections[1]), sep='\s+').convert_dtypes()
        self.Scan_Data['fa'] = df_fa

        df_fb = pd.read_csv(io.StringIO(sections[3]), sep='\s+').convert_dtypes()
        self.Scan_Data['fb'] = df_fb

        return self.Scan_Data


