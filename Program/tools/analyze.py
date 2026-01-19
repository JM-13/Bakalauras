import os
import numpy as np
import pandas as pd
import math
import io
import time
import subprocess
import sys

import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, DotProduct, WhiteKernel
import warnings
from sklearn.exceptions import ConvergenceWarning

def degree_correction(df1, df2, scale=10000):
    # print(f'df1:\n{df1}')
    # print(f'\ndf2:\n{df2}')


    pi2 = 2 * math.pi

    a = df1.to_numpy(dtype=float)
    b = df2.to_numpy(dtype=float)

    mask1 = (a > 0) & (b < 0) & (np.abs(a - b) > 6)
    mask2 = (a < 0) & (b > 0) & (np.abs(b - a) > 6)

    a_corrected = np.where(mask1, np.round(a*scale - pi2*scale)/scale, a)
    b_corrected = np.where(mask2, np.round(b*scale - pi2*scale)/scale, b)

    # if (a != a_corrected).any():
    #     diff_num_a = np.sum(a != a_corrected)
    #     print(f'df1 was corrected {diff_num_a} times')
    # if (b != b_corrected).any():
    #     diff_num_b = np.sum(b != b_corrected)
    #     print(f'df2 was corrected {diff_num_b} times')

    df1_corrected = pd.DataFrame(a_corrected, columns=df1.columns)
    df2_corrected = pd.DataFrame(b_corrected, columns=df2.columns)

    return df1_corrected, df2_corrected

def root_mean_square(pd_series):
    rms = np.sqrt((pd_series**2).mean())

    return rms

def filter_by_atom_number(df, solute, bdp_central=""):
    if bdp_central == "":
        print(f"No bdp_central value given!\nbdp_central needs to be one of the following:\n{"A"} - All 21 atoms\n{"H"} - Exclude replaced Hydrogen atoms\n{"NoH"} - No Hydrogen atoms")
        sys.exit(1)
    elif bdp_central == "A":
        pass
    elif bdp_central == "H":
        pass
    elif bdp_central == "NoH":
        pass
    else:
        print(f"Incorect bdp_central value given: {bdp_central}\nbdp_central needs to be one of the following:\n{"A"} - All 21 atoms\n{"H"} - Exclude replaced Hydrogen atoms\n{"NoH"} - No Hydrogen atoms")
        sys.exit(1)


    atoms = {
                "A":{
                       'BDP-nitrophenyl':            [i+1 for i in range(21)],
                       'BDP-phenyl':                 [i+1 for i in range(21)],
                       'BDP-PP-phenyl':              [1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,28,29,30,41],
                       'BDP-PP-phenyl_OMes_backward':[1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,28,29,30,40],
                       'BDP-PP-phenyl_OMes_forward': [1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,28,29,30,40]
                    },

                "H":{
                       'BDP-nitrophenyl':            [1,2,3,4,5, 7, 8, 9,10,11,12,13,14,16,17,18,19,20,21],
                       'BDP-phenyl':                 [1,2,3,4,5, 7, 8, 9,10,11,12,13,14,16,17,18,19,20,21],
                       'BDP-PP-phenyl':              [1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,28,29],
                       'BDP-PP-phenyl_OMes_backward':[1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,28,29],
                       'BDP-PP-phenyl_OMes_forward': [1,2,3,4,5,11,12,13,14,15,16,17,18,19,20,21,22,28,29]
                    },

                "NoH":{
                       'BDP-nitrophenyl':            [1,2,3,5,7, 9,10,11,12,14,16,18,19,20,21],
                       'BDP-phenyl':                 [1,2,3,5,7, 9,10,11,12,14,16,18,19,20,21],
                       'BDP-PP-phenyl':              [1,2,3,4,5,11,12,13,14,15,16,17,18,19,20],
                       'BDP-PP-phenyl_OMes_backward':[1,2,3,4,5,11,12,13,14,15,16,17,18,19,20],
                       'BDP-PP-phenyl_OMes_forward': [1,2,3,4,5,11,12,13,14,15,16,17,18,19,20]
                      }
                }

    filter_df = pd.DataFrame()
    filter_df["nums"] = df["Definition"].str.findall(r'\d+').apply(lambda x: list(map(int, x)))

    allowed_set = set(atoms[bdp_central][solute])
    filter_df["valid"] = filter_df["nums"].apply(lambda group: all(x in allowed_set for x in group))

    # print(filter_df)
    # sys.exit(1)

    df_filtered = df[filter_df["valid"]]

    return df_filtered

def renumber_atoms_to_match(df1, df2, matching):
    filter_df = pd.DataFrame()
    filter_df["nums"] = df1["Definition"].str.findall(r'\d+').apply(lambda x: list(map(int, x)))
    filter_df["converted_nums"] = (filter_df["nums"].explode().map(matching).groupby(level=0).apply(list))

    df2["nums"] = df2["Definition"].str.findall(r'\d+').apply(lambda x: list(map(int, x)))
    df2_map = df2.set_index(df2["nums"].apply(lambda x: frozenset(x)))
    df2_sorted = df2_map.loc[filter_df["converted_nums"].apply(lambda x: frozenset(x))].reset_index(drop=True)
    df2_sorted.drop("nums", axis=1, inplace=True)

    return df2_sorted


class Retrieve:
    def __init__(self, multiple=False):

        self.multiple = multiple
        self.data = {'Energys':{},
                     'Coordinates':[]}

        self.data_scan = {'fa':[],
                          'fb':[]}

    def regular_data(self, filepath,
                     filter_cords=False, bdp_central="", solute="",
                     reorder_data=False, refrence_df=None, mapping={}):
        with open(filepath, "r") as file:
            for _ in range(4):

                line = file.readline().strip()
                key, value = line.split("=")
                if self.multiple:
                    if not key.strip() in self.data['Energys']:
                        self.data['Energys'][key.strip()] = []
                    self.data['Energys'][key.strip()] += [float(value.strip())]
                else:
                    self.data['Energys'][key.strip()] = float(value.strip())

            file.readline()

            df = pd.read_csv(file, sep=r'\s+').convert_dtypes()

        if filter_cords:
            df = filter_by_atom_number(df, solute, bdp_central)
        if reorder_data:
            df = renumber_atoms_to_match(refrence_df, df, mapping)


        if self.multiple:
            self.data['Coordinates'] += [df]
        else:
            self.data['Coordinates'] = df

        return df

    def return_regular_data(self):
        return self.data

    def scan_data(self, filepath):
        with open(filepath, "r") as file:
            sections = file.read().strip().split("\n\n")

        df_fa = pd.read_csv(io.StringIO(sections[1]), sep=r'\s+').convert_dtypes()
        df_fb = pd.read_csv(io.StringIO(sections[3]), sep=r'\s+').convert_dtypes()

        if self.multiple:
            self.data_scan['fa'] += [df_fa]
            self.data_scan['fb'] += [df_fb]


        else:
            self.data_scan['fa'] = df_fa
            self.data_scan['fb'] = df_fb

        return self.data_scan

    def return_scan_data(self):
        return self.data_scan

class Analyze:
    def __init__(self,
                 all_data_folder,
                 difference_function,
                 solutes,
                 solvents,
                 solutes_shorten,
                 solvents_shorten):

        self.All_Data_folder = all_data_folder
        self.difference_function = difference_function

        self.Solutes = solutes
        self.Solvents = solvents

        self.Short_Solutes = solutes_shorten
        self.Short_Solvents = solvents_shorten

        self.Energy_Energy_types = []
        self.Coordinate_Energy_types = []

        self.All_solute_energy_diffs = {}
        self.All_solute_coordinate_diffs = {}

        self.All_solute_energy_diffs_by_slv = {}
        self.All_solute_coordinate_diffs_by_slv = {}

        self.All_solvent_energy_diffs = {}
        self.All_solvent_coordinate_diffs = {}

        self.All_solvent_energy_diffs_by_slu = {}
        self.All_solvent_coordinate_diffs_by_slu = {}

    def _difference_calculator(self,
                               data,
                               short_names,
                               use_angstroms,
                               use_degrees,
                               save_energy_types=True):
        energy_difs = {}
        coordinate_difs = {}

        energy_energy_types = []
        coordinate_energy_types = []

        #ENERGIES
        for energy_type in data['Energys']:
            energy_difs[energy_type] = {}
            energy_energy_types.append(energy_type)

            for ssln, energy_1 in zip(short_names, data['Energys'][energy_type]):
                difference_column = []

                for energy_2 in data['Energys'][energy_type]:
                    difference = self.difference_function(energy_2, energy_1)
                    difference_column.append(difference)

                energy_difs[energy_type][ssln] = difference_column

            difference_df = pd.DataFrame(energy_difs[energy_type], index=short_names)
            energy_difs[energy_type] = difference_df

        #COORDINATES
        for ssln, coordinate_df in zip(short_names, data['Coordinates']):
            convert_a = True
            convert_d = True

            coord_df_R_1 = coordinate_df[coordinate_df['Name'].str.startswith('R')].reset_index(drop=True)
            coord_df_A_1 = coordinate_df[coordinate_df['Name'].str.startswith('A')].reset_index(drop=True)
            coord_df_D_1 = coordinate_df[coordinate_df['Name'].str.startswith('D')].reset_index(drop=True)

            coord_df_R_1.drop(coord_df_R_1.columns[[0, 1]], axis=1, inplace=True)
            coord_df_A_1.drop(coord_df_A_1.columns[[0, 1]], axis=1, inplace=True)
            coord_df_D_1.drop(coord_df_D_1.columns[[0, 1]], axis=1, inplace=True)

            refrence_df_D_1 = coord_df_D_1

            if not coordinate_energy_types:
                coordinate_energy_types = list(coord_df_R_1.columns)

            if not coordinate_difs:
                for energy_type in coordinate_energy_types:
                    coordinate_difs[energy_type] = {'MAX':{'R':{}, 'A':{}, 'D':{}},
                                                    'RMS':{'R':{}, 'A':{}, 'D':{}}}

            for energy_type in coordinate_difs:
                coordinate_difs[energy_type]['MAX']['R'][ssln] = []
                coordinate_difs[energy_type]['MAX']['A'][ssln] = []
                coordinate_difs[energy_type]['MAX']['D'][ssln] = []

                coordinate_difs[energy_type]['RMS']['R'][ssln] = []
                coordinate_difs[energy_type]['RMS']['A'][ssln] = []
                coordinate_difs[energy_type]['RMS']['D'][ssln] = []

            for coordinate_df2 in data['Coordinates']:
                coord_df_R_2 = coordinate_df2[coordinate_df2['Name'].str.startswith('R')].reset_index(drop=True)
                coord_df_A_2 = coordinate_df2[coordinate_df2['Name'].str.startswith('A')].reset_index(drop=True)
                coord_df_D_2 = coordinate_df2[coordinate_df2['Name'].str.startswith('D')].reset_index(drop=True)

                coord_df_R_2.drop(coord_df_R_2.columns[[0, 1]], axis=1, inplace=True)
                coord_df_A_2.drop(coord_df_A_2.columns[[0, 1]], axis=1, inplace=True)
                coord_df_D_2.drop(coord_df_D_2.columns[[0, 1]], axis=1, inplace=True)

                # print(f"{ssln} of df1:")
                coord_df_D_1, coord_df_D_2 = degree_correction(refrence_df_D_1, coord_df_D_2)

                if use_angstroms:
                    if convert_a:
                        coord_df_R_1 = coord_df_R_1 * 0.529177249
                        convert_a = False
                    coord_df_R_2 = coord_df_R_2 * 0.529177249

                if use_degrees:
                    if convert_d:
                        coord_df_A_1 = np.degrees(coord_df_A_1)
                        convert_d = False
                    coord_df_A_2 = np.degrees(coord_df_A_2)
                    coord_df_D_1 = np.degrees(coord_df_D_1)
                    coord_df_D_2 = np.degrees(coord_df_D_2)

                R_difference = self.difference_function(coord_df_R_2, coord_df_R_1)
                A_difference = self.difference_function(coord_df_A_2, coord_df_A_1)
                D_difference = self.difference_function(coord_df_D_2, coord_df_D_1)

                #MAX
                R_diff_MAX = R_difference.max().to_dict()
                A_diff_MAX = A_difference.max().to_dict()
                D_diff_MAX = D_difference.max().to_dict()

                #RMS
                R_diff_RMS = {col: root_mean_square(R_difference[col]) for col in R_difference}
                A_diff_RMS = {col: root_mean_square(A_difference[col]) for col in A_difference}
                D_diff_RMS = {col: root_mean_square(D_difference[col]) for col in D_difference}

                for energy_type in coordinate_difs:
                    coordinate_difs[energy_type]['MAX']['R'][ssln].append(R_diff_MAX[energy_type])
                    coordinate_difs[energy_type]['MAX']['A'][ssln].append(A_diff_MAX[energy_type])
                    coordinate_difs[energy_type]['MAX']['D'][ssln].append(D_diff_MAX[energy_type])

                    coordinate_difs[energy_type]['RMS']['R'][ssln].append(R_diff_RMS[energy_type])
                    coordinate_difs[energy_type]['RMS']['A'][ssln].append(A_diff_RMS[energy_type])
                    coordinate_difs[energy_type]['RMS']['D'][ssln].append(D_diff_RMS[energy_type])

        for energy_type in coordinate_difs:
            coordinate_difs[energy_type]['MAX']['R'] = pd.DataFrame(coordinate_difs[energy_type]['MAX']['R'], index=short_names)
            coordinate_difs[energy_type]['MAX']['A'] = pd.DataFrame(coordinate_difs[energy_type]['MAX']['A'], index=short_names)
            coordinate_difs[energy_type]['MAX']['D'] = pd.DataFrame(coordinate_difs[energy_type]['MAX']['D'], index=short_names)

            coordinate_difs[energy_type]['RMS']['R'] = pd.DataFrame(coordinate_difs[energy_type]['RMS']['R'], index=short_names)
            coordinate_difs[energy_type]['RMS']['A'] = pd.DataFrame(coordinate_difs[energy_type]['RMS']['A'], index=short_names)
            coordinate_difs[energy_type]['RMS']['D'] = pd.DataFrame(coordinate_difs[energy_type]['RMS']['D'], index=short_names)

        if save_energy_types:
            self.Energy_Energy_types = energy_energy_types
            self.Coordinate_Energy_types = coordinate_energy_types

        return energy_difs, coordinate_difs

    def solvent_differences(self, use_angstroms=False, use_degrees=False):
        solute_folders = {}
        for slu in self.Solutes:
            folder_path = os.path.join(self.All_Data_folder, slu)
            solute_folders[slu] = folder_path

        retrieved_data = {}
        for slu in solute_folders:
            solvent_data = Retrieve(multiple=True)
            for slv in self.Solvents:
                for file in os.listdir(solute_folders[slu]):
                    if file.endswith("_DATA.txt") and (slv in file):
                        solvent_file = os.path.join(solute_folders[slu], file)
                        solvent_data.regular_data(solvent_file)
                        break
            solvent_data = solvent_data.return_regular_data()
            retrieved_data[slu] = solvent_data

        for slu in retrieved_data:
            # print("")
            data = retrieved_data[slu]
            self.All_solute_energy_diffs[slu], self.All_solute_coordinate_diffs[slu] = self._difference_calculator(data, self.Short_Solvents, use_angstroms, use_degrees)

    def solute_differences(self, use_angstroms=False, use_degrees=False, bdp_central="", mapping={}):
        to_reorder = {'BDP-nitrophenyl':False,
                      'BDP-phenyl':False,
                      'BDP-PP-phenyl':True,
                      'BDP-PP-phenyl_OMes_backward':True,
                      'BDP-PP-phenyl_OMes_forward':True}

        solute_folders = {}
        for slu in self.Solutes:
            folder_path = os.path.join(self.All_Data_folder, slu)
            solute_folders[slu] = folder_path


        refrence_df = None

        retrieved_data = {}
        for slv in self.Solvents:
            solute_data = Retrieve(multiple=True)
            for slu in solute_folders:
                for file in os.listdir(solute_folders[slu]):
                    if file.endswith("_DATA.txt") and (slv in file):
                        solvent_file = os.path.join(solute_folders[slu], file)

                        if not to_reorder[slu]:
                            refrence_df = solute_data.regular_data(solvent_file,
                                                     filter_cords=True, bdp_central=bdp_central, solute=slu)
                        elif to_reorder[slu]:
                            solute_data.regular_data(solvent_file,
                                                     filter_cords=True, bdp_central=bdp_central, solute=slu,
                                                     reorder_data=True, refrence_df=refrence_df, mapping=mapping)

                        break
            solute_data = solute_data.return_regular_data()
            retrieved_data[slv] = solute_data

        for slv in retrieved_data:
            # print("")
            data = retrieved_data[slv]
            self.All_solvent_energy_diffs[slv], self.All_solvent_coordinate_diffs[slv] = self._difference_calculator(data, self.Short_Solutes, use_angstroms, use_degrees)

    def display_solvent_differences(self):
        pd.set_option('display.float_format', '{:.2e}'.format)
        for slu in self.All_solute_energy_diffs:
            for energy_type in self.All_solute_energy_diffs[slu]:
                print(f'\n\t\t\t\tSolute: {slu}\n')
                print(f'Energy: {energy_type}')
                print(f'{self.All_solute_energy_diffs[slu][energy_type]}\n\n')

                if energy_type in self.All_solute_coordinate_diffs[slu]:
                    for t in self.All_solute_coordinate_diffs[slu][energy_type]:
                        for coord in self.All_solute_coordinate_diffs[slu][energy_type][t]:
                            print(f'Coordinate: {coord} {t}')
                            print(f'{self.All_solute_coordinate_diffs[slu][energy_type][t][coord]}\n')
                        print('\n')
                else:
                    print('No coordinates were messured for this energy\n')

                input("Press Enter to continue...")
                os.system('clear')

    def display_solute_differences(self):
        pd.set_option('display.float_format', '{:.2e}'.format)
        for slv in self.All_solvent_energy_diffs:
            for energy_type in self.All_solvent_energy_diffs[slv]:
                print(f'\n\t\t\t\tSolvent: {slv}\n')
                print(f'Energy: {energy_type}')
                print(f'{self.All_solvent_energy_diffs[slv][energy_type]}\n\n')

                if energy_type in self.All_solvent_coordinate_diffs[slv]:
                    for t in self.All_solvent_coordinate_diffs[slv][energy_type]:
                        for coord in self.All_solvent_coordinate_diffs[slv][energy_type][t]:
                            print(f'Coordinate: {coord} {t}')
                            print(f'{self.All_solvent_coordinate_diffs[slv][energy_type][t][coord]}\n')
                        print('\n')
                else:
                    print('No coordinates were messured for this energy\n')

                input("Press Enter to continue...")
                os.system('clear')

    def _rearrange_differences(self, arrange_by="Solute"):

        if arrange_by == "Solute":
            sl_names = self.Solvents
            sl_names_short = self.Short_Solvents
            sl_sort_names = self.Solutes
            sl_sort_names_short = self.Short_Solutes
            all_energy_diffs = self.All_solute_energy_diffs
            all_coordinate_diffs = self.All_solute_coordinate_diffs

        elif arrange_by == "Solvent":
            sl_names = self.Solutes
            sl_names_short = self.Short_Solutes
            sl_sort_names = self.Solvents
            sl_sort_names_short = self.Short_Solvents
            all_energy_diffs = self.All_solvent_energy_diffs
            all_coordinate_diffs = self.All_solvent_coordinate_diffs


        All_energy_data = {sln: {} for sln in sl_names}
        for energy in self.Energy_Energy_types:
            for sln, ssln in zip(sl_names, sl_names_short):
                sl_columns = []
                for sls in sl_sort_names:
                    sl_columns.append(all_energy_diffs[sls][energy][ssln])
                sl_df = pd.concat(sl_columns, axis=1)
                sl_df.columns = sl_sort_names_short
                sl_df['Mean'] = sl_df.mean(axis=1)
                sl_df['Median'] = sl_df.drop(columns='Mean').median(axis=1)
                All_energy_data[sln][energy] = sl_df


        All_coordinate_data = {sln: {} for sln in sl_names}
        for sln in All_coordinate_data:
            All_coordinate_data[sln] = {e: {'MAX':{'R':None, 'A':None, 'D':None}, 'RMS':{'R':None, 'A':None, 'D':None}} for e in self.Coordinate_Energy_types}

        for energy in self.Coordinate_Energy_types:
            for sln, ssln in zip(sl_names, sl_names_short):
                for d_t in ['MAX', 'RMS']:
                    for c in ['R', 'A', 'D']:
                        sl_columns = []
                        for sls in sl_sort_names:
                            sl_columns.append(all_coordinate_diffs[sls][energy][d_t][c][ssln])
                        sl_df = pd.concat(sl_columns, axis=1)
                        sl_df.columns = sl_sort_names_short
                        sl_df['Mean'] = sl_df.mean(axis=1)
                        sl_df['Median'] = sl_df.drop(columns='Mean').median(axis=1)
                        All_coordinate_data[sln][energy][d_t][c] = sl_df


        if arrange_by == "Solute":
            self.All_solute_energy_diffs_by_slv = All_energy_data
            self.All_solute_coordinate_diffs_by_slv = All_coordinate_data

        elif arrange_by == "Solvent":
            self.All_solvent_energy_diffs_by_slu = All_energy_data
            self.All_solvent_coordinate_diffs_by_slu = All_coordinate_data

    def display_solvent_differences_by_solute(self):
        pd.set_option('display.float_format', '{:.2e}'.format)
        self._rearrange_differences(arrange_by="Solute")

        for slv in self.All_solute_energy_diffs_by_slv:
            for energy_type in self.All_solute_energy_diffs_by_slv[slv]:
                print(f'\n\t\t\t\tSolvent: {slv}\n')
                print(f'Energy: {energy_type}')
                print(f'{self.All_solute_energy_diffs_by_slv[slv][energy_type]}\n\n')

                if energy_type in self.All_solute_coordinate_diffs_by_slv[slv]:
                    for t in self.All_solute_coordinate_diffs_by_slv[slv][energy_type]:
                        for coord in self.All_solute_coordinate_diffs_by_slv[slv][energy_type][t]:
                            print(f'Coordinate: {coord} {t}')
                            print(f'{self.All_solute_coordinate_diffs_by_slv[slv][energy_type][t][coord]}\n')
                        print('\n')
                else:
                    print('No coordinates were messured for this energy\n')

                input("Press Enter to continue...")
                os.system('clear')

    def display_solute_differences_by_solvent(self):
        pd.set_option('display.float_format', '{:.2e}'.format)
        self._rearrange_differences(arrange_by="Solvent")

        for slu in self.All_solvent_energy_diffs_by_slu:
            for energy_type in self.All_solvent_energy_diffs_by_slu[slu]:
                print(f'\n\t\t\t\tSolute: {slu}\n')
                print(f'Energy: {energy_type}')
                print(f'{self.All_solvent_energy_diffs_by_slu[slu][energy_type]}\n\n')

                if energy_type in self.All_solvent_coordinate_diffs_by_slu[slu]:
                    for t in self.All_solvent_coordinate_diffs_by_slu[slu][energy_type]:
                        for coord in self.All_solvent_coordinate_diffs_by_slu[slu][energy_type][t]:
                            print(f'Coordinate: {coord} {t}')
                            print(f'{self.All_solvent_coordinate_diffs_by_slu[slu][energy_type][t][coord]}\n')
                        print('\n')
                else:
                    print('No coordinates were messured for this energy\n')

                input("Press Enter to continue...")
                os.system('clear')


    def generate_latex_results_document(self,
                                        file_name="",
                                        differences="",
                                        use_solvent_by_solute=False,
                                        use_solute_by_solvent=False,
                                        generate_pdf=True,
                                        clean_files=True):
        if differences == "":
            print(f"No differences value given\ndifferences value needs to be one of the following:\n{"Solute"} - differences between solutes\n{"Solvent"} - differences between solvents")
            sys.exit(1)
        elif differences == "Solvent":
            if not file_name:
                if use_solvent_by_solute:
                    file_name_d = "Differences_between_solvents_by_solute"
                else:
                    file_name_d = "Differences_between_solvents"
            pass
        elif differences == "Solute":
            if not file_name:
                if use_solute_by_solvent:
                    file_name_d = "Differences_between_solutes_by_solvent"
                else:
                    file_name_d = "Differences_between_solutes"
            pass
        else:
            print(f"Incorect differences value given: {differences}\ndifferences value needs to be one of the following:\n{"Solute"} - differences between solutes\n{"Solvent"} - differences between solvents")
            sys.exit(1)

        if use_solvent_by_solute:
            self._rearrange_differences(arrange_by="Solute")
        if use_solute_by_solvent:
            self._rearrange_differences(arrange_by="Solvent")

        if not file_name:
            file_name = file_name_d

        file_full_path = os.path.join(self.All_Data_folder, file_name+'.tex')

        latex_doc = r"""
\documentclass{article}
\usepackage{booktabs}
\usepackage{geometry}
\geometry{margin=0.5in}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{graphicx}
\usepackage{pdflscape}
\usepackage{adjustbox}

\begin{document}
"""
        for slv, slu in zip(self.Solvents, self.Solutes):
            if differences == "Solute":
                if use_solute_by_solvent:
                    page_caption = f"{slu} differences"
                    page_caption = page_caption.replace("_", r"\_")
                else:
                    page_caption = f"Differences in {slv}"
                    page_caption = page_caption.replace("_", r"\_")

            elif differences == "Solvent":
                if use_solvent_by_solute:
                    page_caption = f"{slv} differences"
                    page_caption = page_caption.replace("_", r"\_")
                else:
                    page_caption = f"{slu} solute"
                    page_caption = page_caption.replace("_", r"\_")

            latex_doc += rf"""
\begin{{landscape}}
\begin{{center}}
\LARGE \textbf{{{page_caption}}}
\end{{center}}
\vspace{{1em}}
\begin{{table}}[H]
\centering
\renewcommand{{\arraystretch}}{{1.3}}
\begin{{tabular}}{{ccc}}
"""
            # ---------- Energy tables ----------
            for i, energy in enumerate(self.Coordinate_Energy_types):
                if differences == "Solute":
                    if use_solute_by_solvent:
                        df_energy = self.All_solvent_energy_diffs_by_slu[slu][energy].copy()
                    else:
                        df_energy = self.All_solvent_energy_diffs[slv][energy].copy()
                elif differences == "Solvent":
                    if use_solvent_by_solute:
                        df_energy = self.All_solute_energy_diffs_by_slv[slv][energy].copy()
                    else:
                        df_energy = self.All_solute_energy_diffs[slu][energy].copy()

                df_energy = df_energy.map(lambda x: f"${x:.2e}$")
                table_body = df_energy.to_latex(escape=False)
                sub_caption = f"{energy}"

                latex_doc += rf"""
\begin{{minipage}}{{0.33\linewidth}}
\captionsetup{{skip=0.3em, position=top}}
\caption*{{{sub_caption}}}
\centering
\adjustbox{{max width=\linewidth}}{{{table_body}}}
\end{{minipage}}"""

                if i % 3 == 2:
                    latex_doc += r" \\" + "\n" + r"\vspace{0.8em}" + "\n"
                else:
                    latex_doc += " & "

            # ---------- RMS coordinate tables ----------
            for cord in ["R","A","D"]:
                for i, energy in enumerate(self.Coordinate_Energy_types):
                    if differences == "Solute":
                        if use_solute_by_solvent:
                            df_rms = self.All_solvent_coordinate_diffs_by_slu[slu][energy]['RMS'][cord].copy()
                        else:
                            df_rms = self.All_solvent_coordinate_diffs[slv][energy]['RMS'][cord].copy()
                    elif differences == "Solvent":
                        if use_solvent_by_solute:
                            df_rms = self.All_solute_coordinate_diffs_by_slv[slv][energy]['RMS'][cord].copy()
                        else:
                            df_rms = self.All_solute_coordinate_diffs[slu][energy]['RMS'][cord].copy()

                    df_rms = df_rms.map(lambda x: f"${x:.2e}$")
                    table_body = df_rms.to_latex(escape=False)
                    sub_caption = f"RMS {cord}"

                    latex_doc += rf"""
\begin{{minipage}}{{0.33\linewidth}}"""

                    if cord == 'R':
                        latex_doc += rf"""\vspace{{0.8em}}"""

                    latex_doc += rf"""
\captionsetup{{skip=0.3em, position=top}}
\caption*{{{sub_caption}}}
\centering
\adjustbox{{max width=\linewidth}}{{{table_body}}}
\end{{minipage}}"""

                    if i % 3 == 2:
                        latex_doc += r" \\" + "\n" + r"\vspace{0.8em}" + "\n"
                    else:
                        latex_doc += " & "

            latex_doc += r"""
\end{tabular}
\end{table}
\end{landscape}
\clearpage
"""
        latex_doc += r"\end{document}"

        with open(f"{file_full_path}", "w") as f:
            f.write(latex_doc)

        if generate_pdf:
            tex_file = file_full_path
            pdf_file = os.path.join(self.All_Data_folder, file_name + '.pdf')
            log_file = os.path.join(self.All_Data_folder, file_name + '.log')
            aux_file = os.path.join(self.All_Data_folder, file_name + '.aux')

            subprocess.run(
                            ["pdflatex", "-interaction=nonstopmode", tex_file],
                            cwd=self.All_Data_folder,
                            stdout=subprocess.DEVNULL)

            if clean_files:
                os.remove(tex_file)
                os.remove(log_file)
                os.remove(aux_file)
            print(f'PDF file created at location {pdf_file}')

        else:
            print(f'PDF file created at location {file_full_path}')

    def Scan_Data(self,
                  graph_config,
                  print_min_energy_data=True,
                  show_graphs=True,
                  save_combo_graphs=False,
                  save_single_graphs=False):

        warnings.filterwarnings("ignore", category=ConvergenceWarning)

        Scan_graph_folder = os.path.join(self.All_Data_folder, 'Scan_graphs')
        Scan_graph_single_graphs_folder = os.path.join(Scan_graph_folder, 'Separate_graphs')
        Scan_graph_single_graphs_folder_S0 = os.path.join(Scan_graph_single_graphs_folder, 'S0')
        Scan_graph_single_graphs_folder_S1 = os.path.join(Scan_graph_single_graphs_folder, 'S1')

        os.makedirs(Scan_graph_folder, exist_ok=True)
        os.makedirs(Scan_graph_single_graphs_folder, exist_ok=True)
        os.makedirs(Scan_graph_single_graphs_folder_S0, exist_ok=True)
        os.makedirs(Scan_graph_single_graphs_folder_S1, exist_ok=True)

        for num_slu, slu in enumerate(self.Solutes):
            folder_path = os.path.join(self.All_Data_folder, slu)

            scan_file_data = Retrieve(multiple=True)

            for slv in self.Solvents:
                for file in os.listdir(folder_path):
                    if file.endswith("_SCAN.txt") and (slv in file):
                        filepath = os.path.join(folder_path, file)

                        scan_file_data.scan_data(filepath)

            scan_file_data = scan_file_data.return_scan_data()

            min_S0 = 99999
            for num, slv in enumerate(self.Solvents):
                data = pd.concat([scan_file_data['fa'][num], scan_file_data['fb'][num]])
                if min_S0 > data['S0'].min():
                    min_S0 = data['S0'].min()

            fig_S0, ax_S0 =   plt.subplots()
            fig_S1, ax_S1 =   plt.subplots()

            y_value_ranges = {'S0':{'min': 99999,
                                    'max':-99999},
                            'S1':{'min': 99999,
                                    'max':-99999}}

            print(f'\nDoing {slu} scan data analasys')
            for num, slv in enumerate(self.Solvents):

                kwargs_line = {}
                kwargs_scatter = {}
                kwargs_location_line = {}
                kwargs_energy_line = {}

                for kwa in graph_config['Line']:
                    if graph_config['Line'][kwa]:
                        kwargs_line[kwa] = graph_config['Line'][kwa][num]


                for kwa in graph_config['Scatter']:
                    if graph_config['Scatter'][kwa]:
                        kwargs_scatter[kwa] = graph_config['Scatter'][kwa][num]


                for kwa in graph_config['Minima_location_dots']:
                    if graph_config['Minima_location_dots'][kwa]:
                        kwargs_location_line[kwa] = graph_config['Minima_location_dots'][kwa][num]


                for kwa in graph_config['Minima_energy_lines']:
                    if graph_config['Minima_energy_lines'][kwa]:
                        kwargs_energy_line[kwa] = graph_config['Minima_energy_lines'][kwa][num]

                fa_data = scan_file_data['fa'][num]
                fa_data['S0'] = fa_data['S0'] - min_S0
                fa_data['S1'] = fa_data['S1'] - min_S0

                fb_data = scan_file_data['fb'][num]
                fb_data['S0'] = fb_data['S0'] - min_S0
                fb_data['S1'] = fb_data['S1'] - min_S0

                data = pd.concat([fa_data, fb_data])

                degree_axis_values_average = np.array(data['D_average'])
                indices = np.argsort(degree_axis_values_average)

                degree_axis_values_fixed   = np.array(data['D_fixed'])[indices]
                degree_axis_values_average = np.array(data['D_average'])[indices]
                energy_axis_values_S0 = np.array(data['S0'])[indices]
                energy_axis_values_S1 = np.array(data['S1'])[indices]

                sorted_data = pd.DataFrame({
                    'D_fixed':degree_axis_values_fixed,
                    'D_average':degree_axis_values_average,
                    'S0':energy_axis_values_S0,
                    'S1':energy_axis_values_S1
                    })

                min_S0_row  = sorted_data.nsmallest(1, 'S0').to_dict()
                energy_S0  = round(list(min_S0_row['S0'].values())[0], 4)
                S0_location_fixed = round(list(min_S0_row['D_fixed'].values())[0], 4)
                S0_location_average = round(list(min_S0_row['D_average'].values())[0], 4)

                min_S1m_row = sorted_data.nsmallest(1, 'S1').to_dict()
                energy_S1m = round(list(min_S1m_row['S1'].values())[0], 4)
                S1m_location_fixed = round(list(min_S1m_row['D_fixed'].values())[0], 4)
                S1m_location_average = round(list(min_S1m_row['D_average'].values())[0], 4)

                max_S1_row  = sorted_data.nlargest(1, 'S1').to_dict()
                max_S1_index = list(max_S1_row['S1'].keys())[0]
                min_S1r_row = sorted_data.iloc[max_S1_index:].nsmallest(1, 'S1').to_dict()
                energy_S1r = round(list(min_S1r_row['S1'].values())[0], 4)
                S1r_location_fixed = round(list(min_S1r_row['D_fixed'].values())[0], 4)
                S1r_location_average = round(list(min_S1r_row['D_average'].values())[0], 4)

                if print_min_energy_data:
                    print(f'    Minimum energies in {slv}:')
                    print(f'        S0  : energy = {energy_S0}, fixed_cord = {S0_location_fixed}, average_cord = {S0_location_average}')
                    print(f'        S1m : energy = {energy_S1m}, fixed_cord = {S1m_location_fixed}, average_cord = {S1m_location_average}')
                    print(f'        S1r : energy = {energy_S1r}, fixed_cord = {S1r_location_fixed}, average_cord = {S1r_location_average}')

                #S0 plots /////////////////////////
                xdata_S0  = np.array(sorted_data['D_average'])
                x_plot_S0 = np.linspace(xdata_S0.min(), xdata_S0.max(), 500)
                ydata_S0  = np.array(sorted_data['S0'])

                RBF_bounds_S0 = RBF(length_scale_bounds=(graph_config['Fit_curvyness']['S0'], 50))
                kernel_S0 = RBF_bounds_S0 + WhiteKernel(noise_level_bounds=(1e-8, 1e5))
                gp_S0 = GaussianProcessRegressor(kernel=kernel_S0, normalize_y=True)
                gp_S0.fit(xdata_S0.reshape(-1,1), ydata_S0)

                y_mean_S0, y_std_S0 = gp_S0.predict(x_plot_S0.reshape(-1,1), return_std=True)

                if ydata_S0.min() < y_value_ranges['S0']['min']:
                    y_value_ranges['S0']['min'] = ydata_S0.min()

                if ydata_S0.max() > y_value_ranges['S0']['max']:
                    y_value_ranges['S0']['max'] = ydata_S0.max()

                ax_S0.scatter(xdata_S0, ydata_S0, **kwargs_scatter)
                ax_S0.plot(x_plot_S0, y_mean_S0, **kwargs_line)

                S0_location = list(min_S0_row['D_average'].values())[0]
                S0_energy   = list(min_S0_row['S0'].values())[0]
                ax_S0.scatter([S0_location], [S0_energy], **kwargs_location_line)
                # ax_S0.plot([graph_config['Xlimits'][0], S0_location], [S0_energy, S0_energy], **kwargs_energy_line)

                    #Single graph
                if save_single_graphs:
                    fig, ax = plt.subplots()
                    ax.scatter(xdata_S0, ydata_S0, **kwargs_scatter)
                    ax.plot(x_plot_S0, y_mean_S0, **kwargs_line)
                    ax.scatter([S0_location], [S0_energy], **kwargs_location_line)
                    ax.grid(True)
                    ax.legend(loc="upper left")
                    ax.set_title(f"{slu} S0 plane")
                    ax.set(xlabel = "Degrees",
                        ylabel = "Energy",
                        xlim=(graph_config['Xlimits'][0], graph_config['Xlimits'][1]),
                        xticks=np.arange(graph_config['Xlimits'][0], graph_config['Xlimits'][1]+1, 5)
                        )
                    fig_save_location = os.path.join(Scan_graph_single_graphs_folder_S0, f"{slu}_in_{slv}_S0_plane")
                    fig.savefig(fig_save_location, dpi=600, bbox_inches='tight')
                    plt.close(fig)
                    # print(f'    Saved {slu} in {slv} S0 plane graph to {fig_save_location}')

                #S1 plots /////////////////////////
                S1_limit = fb_data['S1'].max() #need smarter filter maybe?
                filterred_S1 = []
                filterred_Da = []
                for x,y in zip(sorted_data['D_average'], sorted_data['S1']):
                    if y <= S1_limit:
                        filterred_Da.append(x)
                        filterred_S1.append(y)
                filterred_S1 = np.array(filterred_S1)
                filterred_Da = np.array(filterred_Da)

                xdata_S1 = filterred_Da
                x_plot_S1 = np.linspace(xdata_S1.min(), xdata_S1.max(), 500)
                ydata_S1 = filterred_S1

                RBF_bounds_S1 = RBF(length_scale_bounds=(graph_config['Fit_curvyness']['S1'], 50))
                kernel_S1 = RBF_bounds_S1 + WhiteKernel(noise_level_bounds=(1e-8, 1e5))
                gp_S1 = GaussianProcessRegressor(kernel=kernel_S1, normalize_y=True)
                gp_S1.fit(xdata_S1.reshape(-1,1), ydata_S1)

                y_mean_S1, y_std_S1 = gp_S1.predict(x_plot_S1.reshape(-1,1), return_std=True)

                if sorted_data['S1'].min() < y_value_ranges['S1']['min']:
                    y_value_ranges['S1']['min'] = sorted_data['S1'].min()

                if sorted_data['S1'].max() > y_value_ranges['S1']['max']:
                    y_value_ranges['S1']['max'] = sorted_data['S1'].max()

                ax_S1.scatter(sorted_data['D_average'], sorted_data['S1'], **kwargs_scatter)
                ax_S1.plot(x_plot_S1, y_mean_S1, **kwargs_line)

                S1m_location = list(min_S1m_row['D_average'].values())[0]
                S1m_energy   = list(min_S1m_row['S1'].values())[0]
                ax_S1.scatter([S1m_location], [S1m_energy], **kwargs_location_line)
                # ax_S1.plot([graph_config['Xlimits'][0], S1m_location], [S1m_energy, S1m_energy], **kwargs_energy_line)

                S1r_location = list(min_S1r_row['D_average'].values())[0]
                S1r_energy   = list(min_S1r_row['S1'].values())[0]
                ax_S1.scatter([S1r_location], [S1r_energy], **kwargs_location_line)
                # ax_S1.plot([graph_config['Xlimits'][1], S1r_location], [S1r_energy, S1r_energy], **kwargs_energy_line)

                    #Single graph
                if save_single_graphs:
                    fig, ax = plt.subplots()
                    ax.scatter(sorted_data['D_average'], sorted_data['S1'], **kwargs_scatter)
                    ax.plot(x_plot_S1, y_mean_S1, **kwargs_line)
                    ax.scatter([S1m_location], [S1m_energy], **kwargs_location_line)
                    ax.scatter([S1r_location], [S1r_energy], **kwargs_location_line)
                    ax.grid(True)
                    ax.legend(loc="upper left")
                    ax.set_title(f"{slu} S1 plane")
                    ax.set(xlabel = "Degrees",
                        ylabel = "Energy",
                        xlim=(graph_config['Xlimits'][0], graph_config['Xlimits'][1]),
                        xticks=np.arange(graph_config['Xlimits'][0], graph_config['Xlimits'][1]+1, 5)
                        )
                    fig_save_location = os.path.join(Scan_graph_single_graphs_folder_S1, f"{slu}_in_{slv}_S1_plane")
                    fig.savefig(fig_save_location, dpi=600, bbox_inches='tight')
                    plt.close(fig)
                    # print(f'    Saved {slu} in {slv} S1 plane graph to {fig_save_location}')

            #S0 plot config
            y_start = graph_config['Ylimits']['S0']['start'][num_slu]
            y_stop  = graph_config['Ylimits']['S0']['stop'][num_slu]
            y_step  = graph_config['Ylimits']['S0']['step'][num_slu]
            yticks = np.arange(y_start, y_stop+y_step/10, y_step)

            ax_S0.grid(True)
            ax_S0.legend(loc="upper left")
            ax_S0.set_title(f"{slu} S0 plane")
            ax_S0.set(xlabel = "Degrees",
                    ylabel = "Energy",
                    xlim=(graph_config['Xlimits'][0], graph_config['Xlimits'][1]),
                    xticks=np.arange(graph_config['Xlimits'][0], graph_config['Xlimits'][1]+1, 5),
                    ylim=(y_start-y_step/2, y_stop),
                    yticks=yticks
                    )

            if save_combo_graphs:
                fig_S0_save_location = os.path.join(Scan_graph_folder, f"{slu}_S0_planes")
                fig_S0.savefig(fig_S0_save_location, dpi=600, bbox_inches='tight')
                print(f'Saved {slu} S0 planes graph to {fig_S0_save_location}')

            #S1 plot config
            y_start = graph_config['Ylimits']['S1']['start'][num_slu]
            y_stop  = graph_config['Ylimits']['S1']['stop'][num_slu]
            y_step  = graph_config['Ylimits']['S1']['step'][num_slu]
            yticks = np.arange(y_start, y_stop+y_step/10, y_step)

            ax_S1.grid(True)
            ax_S1.legend(loc="upper left")
            ax_S1.set_title(f"{slu} S1 plane")
            ax_S1.set(xlabel = "Degrees",
                    ylabel = "Energy",
                    xlim=(graph_config['Xlimits'][0], graph_config['Xlimits'][1]),
                    xticks=np.arange(graph_config['Xlimits'][0], graph_config['Xlimits'][1]+1, 5),
                    ylim=(y_start-y_step/2, y_stop),
                    yticks=yticks
                    )

            if save_combo_graphs:
                fig_S1_save_location = os.path.join(Scan_graph_folder, f"{slu}_S1_planes")
                fig_S1.savefig(fig_S1_save_location, dpi=600, bbox_inches='tight')
                print(f'Saved {slu} S1 planes graph to {fig_S1_save_location}')

            if show_graphs:
                plt.show()

            plt.close(fig_S0)
            plt.close(fig_S1)
