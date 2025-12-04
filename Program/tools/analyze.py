import os
import numpy as np
import pandas as pd
import math
import io
import time
import subprocess
import sys

def degree_correction(df1, df2, scale=10000):
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

    df_filtered = df[filter_df["valid"]]

    return df_filtered

class Retrieve:
    def __init__(self, multiple=False):

        self.multiple = multiple
        self.data = {'Energys':{},
                     'Coordinates':[]}

        self.scan_data = {'fa':None,
                          'fb':None}

    def regular_data(self, filepath, filter_cords=False, bdp_central="", solute=""):
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
        self.scan_data['fa'] = df_fa

        df_fb = pd.read_csv(io.StringIO(sections[3]), sep=r'\s+').convert_dtypes()
        self.scan_data['fb'] = df_fb

        return self.scan_data

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

    def solute_differences(self, use_angstroms=False, use_degrees=False, bdp_central=""):
        solute_folders = {}
        for slu in self.Solutes:
            folder_path = os.path.join(self.All_Data_folder, slu)
            solute_folders[slu] = folder_path

        retrieved_data = {}
        for slv in self.Solvents:
            solute_data = Retrieve(multiple=True)
            for slu in solute_folders:
                for file in os.listdir(solute_folders[slu]):
                    if file.endswith("_DATA.txt") and (slv in file):
                        solvent_file = os.path.join(solute_folders[slu], file)
                        solute_data.regular_data(solvent_file, filter_cords=True, bdp_central=bdp_central, solute=slu)
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

    def _rearrange_solvent_differences_by_solute(self):
        missing_variables = 0
        if not self.All_solute_energy_diffs:
            print("Missing All_solute_energy_diffs")
            missing_variables+=1

        if not self.All_solute_coordinate_diffs:
            print("Missing All_solute_coordinate_diffs")
            missing_variables+=1

        if not self.Energy_Energy_types:
            print("Missing Energy_Energy_types")
            missing_variables+=1

        if not self.Coordinate_Energy_types:
            print("Missing Coordinate_Energy_types")
            missing_variables+=1

        if missing_variables > 0:
            print(f"Missing {missing_variables} required variables for _rearrange_solvent_differences_by_solute to run\nPlease run everything in the correct order!")
            sys.exit(1)

        All_energy_data = {slv: {} for slv in self.Solvents}
        for energy in self.Energy_Energy_types:
            for slv, sslv in zip(self.Solvents, self.Short_Solvents):
                slv_columns = []
                for slu in self.Solutes:
                    slv_columns.append(self.All_solute_energy_diffs[slu][energy][sslv])
                slv_df = pd.concat(slv_columns, axis=1)
                slv_df.columns = self.Short_Solutes
                slv_df['Mean'] = slv_df.mean(axis=1)
                slv_df['Median'] = slv_df.drop(columns='Mean').median(axis=1)
                All_energy_data[slv][energy] = slv_df
        self.All_solute_energy_diffs_by_slv = All_energy_data

        All_coordinate_data = {slv: {} for slv in self.Solvents}
        for slv in All_coordinate_data:
            All_coordinate_data[slv] = {e: {'MAX':{'R':None, 'A':None, 'D':None}, 'RMS':{'R':None, 'A':None, 'D':None}} for e in self.Coordinate_Energy_types}

        for energy in self.Coordinate_Energy_types:
            for slv, sslv in zip(self.Solvents, self.Short_Solvents):
                for d_t in ['MAX', 'RMS']:
                    for c in ['R', 'A', 'D']:
                        slv_columns = []
                        for slu in self.Solutes:
                            slv_columns.append(self.All_solute_coordinate_diffs[slu][energy][d_t][c][sslv])
                        slv_df = pd.concat(slv_columns, axis=1)
                        slv_df.columns = self.Short_Solutes
                        slv_df['Mean'] = slv_df.mean(axis=1)
                        slv_df['Median'] = slv_df.drop(columns='Mean').median(axis=1)
                        All_coordinate_data[slv][energy][d_t][c] = slv_df
        self.All_solute_coordinate_diffs_by_slv = All_coordinate_data

    def display_solvent_differences_by_solute(self):
        pd.set_option('display.float_format', '{:.2e}'.format)
        self._rearrange_solvent_differences_by_solute()

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

    def generate_latex_results_document(self,
                                        file_name="",
                                        differences="",
                                        use_solvent_by_solute=True,
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
                file_name_d = "Differences_between_solutes"
            pass
        else:
            print(f"Incorect differences value given: {differences}\ndifferences value needs to be one of the following:\n{"Solute"} - differences between solutes\n{"Solvent"} - differences between solvents")
            sys.exit(1)

        if use_solvent_by_solute:
            self._rearrange_solvent_differences_by_solute()

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
                page_caption = f"Differences in {slv}"
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

