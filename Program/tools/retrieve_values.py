import os
import pandas as pd
import io

class Retvieve:
    def __init__(self, filepath=None, multiple=False):
        self.filepath = filepath
        self.multiple = multiple
        self.Data = {'Energys':{},
                     'Coordinates':[]}
        self.Scan_Data = {'fa':None,
                          'fb':None}

    def Regular_data(self, filepath=None):
        if not filepath:
            filepath = self.filepath

        with open(filepath, "r") as file:
            for _ in range(4):
                line = file.readline().strip()
                key, value = line.split("=")

                if self.multiple:
                    if not key.strip() in self.Data['Energys']:
                        self.Data['Energys'][key.strip()] = []

                    self.Data['Energys'][key.strip()] += [float(value.strip())]
                else:
                    self.Data['Energys'][key.strip()] = float(value.strip())

            file.readline()

            df = pd.read_csv(file, sep=r'\s+').convert_dtypes()

        if self.multiple:
            self.Data['Coordinates'] += [df]
        else:
            self.Data['Coordinates'] = df

        return self.Data


    def Scan_data(self, filepath):
        with open(filepath, "r") as file:
            sections = file.read().strip().split("\n\n")

        df_fa = pd.read_csv(io.StringIO(sections[1]), sep=r'\s+').convert_dtypes()
        self.Scan_Data['fa'] = df_fa

        df_fb = pd.read_csv(io.StringIO(sections[3]), sep=r'\s+').convert_dtypes()
        self.Scan_Data['fb'] = df_fb

        return self.Scan_Data

    def Return_multiple(self):
        return self.Data


