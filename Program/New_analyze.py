import os
import numpy as np
import pandas as pd
import math
import io
import time
import sys

def Degree_correction(df1, df2, scale=10000):
    pi2 = 2 * math.pi

    a = df1.to_numpy()
    b = df2.to_numpy()

    a = df1.to_numpy(dtype=float)
    b = df2.to_numpy(dtype=float)

    mask1 = (a > 0) & (b < 0) & (np.abs(a - b) > 6)
    mask2 = (a < 0) & (b > 0) & (np.abs(b - a) > 6)

    a_corrected = np.where(mask1, np.round(a*scale - pi2*scale)/scale, a)
    b_corrected = np.where(mask2, np.round(b*scale - pi2*scale)/scale, b)

    df1_corrected = pd.DataFrame(a_corrected, columns=df1.columns)
    df2_corrected = pd.DataFrame(b_corrected, columns=df2.columns)

    return df1_corrected, df2_corrected


def Root_Mean_Square(series):
    return np.sqrt((series**2).mean())

def Filter_cords_by_atom_number(df, solute, keep=""):
    if keep == "A":
        pass
    elif keep == "H":
        pass
    elif keep == "NoH":
        pass
    else:
        print(f"Incorect keep value given: {keep}\nkeep needs to be one of the following:\n{"A"} - All 21 atoms\n{"H"} - Exclude replaced Hydrogen atoms\n{"NoH"} - No Hydrogen atoms")
        sys.exit(1)


    Atoms = {
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

    allowed_set = set(Atoms[keep][solute])
    filter_df["valid"] = filter_df["nums"].apply(lambda group: all(x in allowed_set for x in group))

    df_filtered = df[filter_df["valid"]]

    return df_filtered

class Retvieve:
    def __init__(self, multiple=False):

        self.multiple = multiple
        self.Data = {'Energys':{},
                     'Coordinates':[]}

        self.Filtered_Coordinates = []

        self.Scan_Data = {'fa':None,
                          'fb':None}

    def Regular_data(self, filepath, filter_cords=False, keep="", solute=""):
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

            if filter_cords:
                df = Filter_cords_by_atom_number(df, solute, keep)

        if self.multiple:
            self.Data['Coordinates'] += [df]
        else:
            self.Data['Coordinates'] = df

        return df

    def Return_Data(self):
        return self.Data

    def Scan_data(self, filepath):
        with open(filepath, "r") as file:
            sections = file.read().strip().split("\n\n")

        df_fa = pd.read_csv(io.StringIO(sections[1]), sep=r'\s+').convert_dtypes()
        self.Scan_Data['fa'] = df_fa

        df_fb = pd.read_csv(io.StringIO(sections[3]), sep=r'\s+').convert_dtypes()
        self.Scan_Data['fb'] = df_fb

        return self.Scan_Data

class Difference_Between:
    def __init__(self,
                 All_Data_folder,
                 difference_function,
                 solutes, solvents,
                 solutes_shorten,
                 solvents_shorten):

        self.All_Data_folder = All_Data_folder
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

    def Difference_calculator(self, data, short_names_of_compared, save_energy_types=True):
        Energy_Energy_types = []
        Energy_difs = {}

        Coordinate_Energy_types = []
        Coordinate_difs = {}

        #ENERGIES
        for energy_type in data['Energys']:
            Energy_difs[energy_type] = {}
            Energy_Energy_types.append(energy_type)

            for ssl, energy_1 in zip(short_names_of_compared, data['Energys'][energy_type]):
                difference_column = []

                for energy_2 in data['Energys'][energy_type]:
                    difference = self.difference_function(energy_2, energy_1)
                    difference_column.append(difference)

                Energy_difs[energy_type][ssl] = difference_column

            difference_df = pd.DataFrame(Energy_difs[energy_type], index=short_names_of_compared)
            Energy_difs[energy_type] = difference_df

        #COORDINATES
        for ssl, coordinate_df in zip(short_names_of_compared, data['Coordinates']):

            coord_df_R_1 = coordinate_df[coordinate_df['Name'].str.startswith('R')].reset_index(drop=True)
            coord_df_A_1 = coordinate_df[coordinate_df['Name'].str.startswith('A')].reset_index(drop=True)
            coord_df_D_1 = coordinate_df[coordinate_df['Name'].str.startswith('D')].reset_index(drop=True)

            coord_df_R_1.drop(coord_df_R_1.columns[[0, 1]], axis=1, inplace=True)
            coord_df_A_1.drop(coord_df_A_1.columns[[0, 1]], axis=1, inplace=True)
            coord_df_D_1.drop(coord_df_D_1.columns[[0, 1]], axis=1, inplace=True)

            if not Coordinate_Energy_types:
                Coordinate_Energy_types = list(coord_df_R_1.columns)

            if not Coordinate_difs:
                for energy_type in Coordinate_Energy_types:
                    Coordinate_difs[energy_type] = {'MAX':{'R':{}, 'A':{}, 'D':{}},
                                                    'RMS':{'R':{}, 'A':{}, 'D':{}}}
            for energy_type in Coordinate_Energy_types:
                Coordinate_difs[energy_type]['MAX']['R'][ssl] = []
                Coordinate_difs[energy_type]['MAX']['A'][ssl] = []
                Coordinate_difs[energy_type]['MAX']['D'][ssl] = []

                Coordinate_difs[energy_type]['RMS']['R'][ssl] = []
                Coordinate_difs[energy_type]['RMS']['A'][ssl] = []
                Coordinate_difs[energy_type]['RMS']['D'][ssl] = []

            for coordinate_df2 in data['Coordinates']:
                coord_df_R_2 = coordinate_df2[coordinate_df2['Name'].str.startswith('R')].reset_index(drop=True)
                coord_df_A_2 = coordinate_df2[coordinate_df2['Name'].str.startswith('A')].reset_index(drop=True)
                coord_df_D_2 = coordinate_df2[coordinate_df2['Name'].str.startswith('D')].reset_index(drop=True)

                coord_df_R_2.drop(coord_df_R_2.columns[[0, 1]], axis=1, inplace=True)
                coord_df_A_2.drop(coord_df_A_2.columns[[0, 1]], axis=1, inplace=True)
                coord_df_D_2.drop(coord_df_D_2.columns[[0, 1]], axis=1, inplace=True)

                coord_df_D_1, coord_df_D_2 = Degree_correction(coord_df_D_1, coord_df_D_2)

                R_difference = coord_df_R_2.combine(coord_df_R_1, lambda s1, s2: s1.combine(s2, self.difference_function))
                A_difference = coord_df_A_2.combine(coord_df_A_1, lambda s1, s2: s1.combine(s2, self.difference_function))
                D_difference = coord_df_D_2.combine(coord_df_D_1, lambda s1, s2: s1.combine(s2, self.difference_function))

                #MAX
                R_diff_MAX = R_difference.max().to_dict()
                A_diff_MAX = A_difference.max().to_dict()
                D_diff_MAX = D_difference.max().to_dict()

                #RMS
                R_diff_RMS = {col: Root_Mean_Square(R_difference[col]) for col in R_difference}
                A_diff_RMS = {col: Root_Mean_Square(A_difference[col]) for col in A_difference}
                D_diff_RMS = {col: Root_Mean_Square(D_difference[col]) for col in D_difference}

                for energy_type in Coordinate_Energy_types:
                    Coordinate_difs[energy_type]['MAX']['R'][ssl].append(R_diff_MAX[energy_type])
                    Coordinate_difs[energy_type]['MAX']['A'][ssl].append(A_diff_MAX[energy_type])
                    Coordinate_difs[energy_type]['MAX']['D'][ssl].append(D_diff_MAX[energy_type])

                    Coordinate_difs[energy_type]['RMS']['R'][ssl].append(R_diff_RMS[energy_type])
                    Coordinate_difs[energy_type]['RMS']['A'][ssl].append(A_diff_RMS[energy_type])
                    Coordinate_difs[energy_type]['RMS']['D'][ssl].append(D_diff_RMS[energy_type])

        for energy_type in Coordinate_Energy_types:
            Coordinate_difs[energy_type]['MAX']['R'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['R'], index=short_names_of_compared)
            Coordinate_difs[energy_type]['MAX']['A'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['A'], index=short_names_of_compared)
            Coordinate_difs[energy_type]['MAX']['D'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['D'], index=short_names_of_compared)

            Coordinate_difs[energy_type]['RMS']['R'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['R'], index=short_names_of_compared)
            Coordinate_difs[energy_type]['RMS']['A'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['A'], index=short_names_of_compared)
            Coordinate_difs[energy_type]['RMS']['D'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['D'], index=short_names_of_compared)

        if save_energy_types:
            self.Energy_Energy_types = Energy_Energy_types
            self.Coordinate_Energy_types = Coordinate_Energy_types

        return Energy_difs, Coordinate_difs

    def Solvents_diff(self):
        solutes  = self.Solutes
        solvents = self.Solvents

        solute_folders = {}
        for slu in solutes:
            folder_path = os.path.join(self.All_Data_folder, slu)
            solute_folders[slu] = folder_path

        Retrieved_data = {}

        for slu in solute_folders:
            slv_data = Retvieve(multiple=True)
            for slv in solvents:
                for file in os.listdir(solute_folders[slu]):
                    if file.endswith("_DATA.txt") and (slv in file):
                        solvent_file = os.path.join(solute_folders[slu], file)
                        slv_data.Regular_data(solvent_file)
                        break
            slv_data = slv_data.Return_Data()
            Retrieved_data[slu] = slv_data

        for slu in Retrieved_data:
            data = Retrieved_data[slu]
            self.All_solute_energy_diffs[slu], self.All_solute_coordinate_diffs[slu] = self.Difference_calculator(data, self.Short_Solvents)

    def Solutes_diff(self, keep=""):
        solutes  = self.Solutes
        solvents = self.Solvents

        solute_folders = {}
        for slu in solutes:
            folder_path = os.path.join(self.All_Data_folder, slu)
            solute_folders[slu] = folder_path

        Retrieved_data = {}

        for slv in solvents:
            slu_data = Retvieve(multiple=True)
            for slu in solute_folders:
                for file in os.listdir(solute_folders[slu]):
                    if file.endswith("_DATA.txt") and (slv in file):
                        solvent_file = os.path.join(solute_folders[slu], file)
                        slu_data.Regular_data(solvent_file, filter_cords=True, keep=keep, solute=slu)
                        break
            slu_data = slu_data.Return_Data()
            Retrieved_data[slv] = slu_data

        for slv in Retrieved_data:
            data = Retrieved_data[slv]
            self.All_solvent_energy_diffs[slv], self.All_solvent_coordinate_diffs[slv] = self.Difference_calculator(data, self.Short_Solutes)

    def Display_in_neet_order(self):

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

    def Display_in_neet_order3(self):

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

    def Arrange_by_solute(self):

        All_energy_data = {slv: {} for slv in self.Solvents}
        for energy in self.Energy_Energy_types:
            for slv, sslv in zip(self.Solvents, self.Short_Solvents):
                slv_columns = []
                for slu in self.Solutes:
                    slv_column = self.All_solute_energy_diffs[slu][energy][sslv]
                    slv_columns.append(slv_column)
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
                            slv_column = self.All_solute_coordinate_diffs[slu][energy][d_t][c][sslv]
                            slv_columns.append(slv_column)
                        slv_df = pd.concat(slv_columns, axis=1)
                        slv_df.columns = self.Short_Solutes
                        slv_df['Mean'] = slv_df.mean(axis=1)
                        slv_df['Median'] = slv_df.drop(columns='Mean').median(axis=1)
                        All_coordinate_data[slv][energy][d_t][c] = slv_df

        self.All_solute_coordinate_diffs_by_slv = All_coordinate_data

    def Display_in_neet_order2(self):

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

    def Return_Energies(self):
        return self.All_solute_energy_diffs

    def Return_Coordinates(self):
        return self.All_solute_coordinate_diffs

    def Generate_Latex_document(self, file_name, differences=""):
        if differences == "Solute":
            pass
        elif differences == "Solvent":
            pass
        else:
            print(f"Incorect differences value given: {differences}\ndifferences needs to be one of the following:\n{"Solute"} - differences between solutes\n{"Solvent"} - differences between solvents")
            sys.exit(1)

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

        for slv in self.Solvents:
            page_caption = f"{slv} differences"

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

            # ---------- Row 1: Energy tables ----------
            for i, energy in enumerate(self.Coordinate_Energy_types):
                label = f"tab:{slv}_{energy}"
                if differences == "Solvent":
                    print(self.All_solute_energy_diffs_by_slv.keys())
                    df_energy = self.All_solute_energy_diffs_by_slv[slv][energy].copy()
                elif differences == "Solute":
                    df_energy = self.All_solvent_energy_diffs[slv][energy].copy()
                df_energy = df_energy.map(lambda x: f"${x:.2e}$" if isinstance(x, (float,int)) else x)
                table_body = df_energy.to_latex(escape=False)
                sub_caption = f"{energy}"

                latex_doc += rf"""
\begin{{minipage}}{{0.33\linewidth}}
\captionsetup{{skip=0.3em, position=top}}
\caption*{{{sub_caption}}}
\label{{{label}}}
\centering
\adjustbox{{max width=\linewidth}}{{{table_body}}}
\end{{minipage}}"""

                if i % 3 == 2:
                    latex_doc += r" \\" + "\n" + r"\vspace{0.8em}" + "\n"
                else:
                    latex_doc += " & "

            # ---------- Rows 2-4: RMS coordinate tables ----------
            for cord in ["R","A","D"]:
                for i, energy in enumerate(self.Coordinate_Energy_types):
                    label = f"tab:{slv}_RMS_{cord}_{energy}"
                    if differences == "Solvent":
                        df_rms = self.All_solute_coordinate_diffs_by_slv[slv][energy]['RMS'][cord].copy()
                    elif differences == "Solute":
                        df_rms = self.All_solvent_coordinate_diffs[slv][energy]['RMS'][cord].copy()
                    df_rms = df_rms.map(lambda x: f"${x:.2e}$" if isinstance(x,(float,int)) else x)
                    table_body = df_rms.to_latex(escape=False)
                    sub_caption = f"RMS {cord}"

                    latex_doc += rf"""
\begin{{minipage}}{{0.33\linewidth}}"""

                    if cord == 'R':
                        latex_doc += rf"""\vspace{{0.8em}}"""

                    latex_doc += rf"""
\captionsetup{{skip=0.3em, position=top}}
\caption*{{{sub_caption}}}
\label{{{label}}}
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

        with open(f"{file_name}", "w") as f:
            f.write(latex_doc)

        print("LATEX DONE!")









 # def Solvents_diff(self):
    #
    #     solutes  = self.Solutes
    #     solvents = self.Solvents
    #
    #     for slu in solutes:
    #         solute_folder = os.path.join(self.All_Data_folder, slu)
    #         solvent_data = Retvieve(multiple=True)
    #
    #         for slv in solvents:
    #             for file in os.listdir(solute_folder):
    #                 if file.endswith("_DATA.txt") and (slv in file):
    #                     solvent_file = os.path.join(solute_folder, file)
    #                     solvent_data.Regular_data(solvent_file)
    #                     break
    #
    #         solvent_data = solvent_data.Return_multiple()
    #
    #         #ENERGIES
    #         Energy_types = []
    #         Energy_difs = {}
    #         for energy_type in solvent_data['Energys']:
    #             Energy_difs[energy_type] = {}
    #             Energy_types.append(energy_type)
    #
    #             for slv, sslv, slv_energy_1 in zip(solvents,
    #                                                self.Short_Solvents,
    #                                                solvent_data['Energys'][energy_type]):
    #                 slv_difference_column = []
    #
    #                 for slv_energy_2 in solvent_data['Energys'][energy_type]:
    #                     slv_difference = self.difference_function(slv_energy_2, slv_energy_1)
    #                     slv_difference_column.append(slv_difference)
    #
    #                 Energy_difs[energy_type][sslv] = slv_difference_column
    #
    #             solvent_difference_df = pd.DataFrame(Energy_difs[energy_type], index=self.Short_Solvents)
    #             Energy_difs[energy_type] = solvent_difference_df
    #
    #         self.All_solute_energy_diffs[slu] = Energy_difs
    #         self.Energy_Energy_types = Energy_types
    #
    #         #COORDINATES
    #         Energy_types = []
    #         Coordinate_difs = {}
    #         for slv, sslv, slv_coordinate_df in zip(solvents,
    #                                                 self.Short_Solvents,
    #                                                 solvent_data['Coordinates']):
    #
    #             slv_coord_df_R_1 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('R')].reset_index(drop=True)
    #             slv_coord_df_A_1 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('A')].reset_index(drop=True)
    #             slv_coord_df_D_1 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('D')].reset_index(drop=True)
    #
    #             slv_coord_df_R_1.drop(slv_coord_df_R_1.columns[[0, 1]], axis=1, inplace=True)
    #             slv_coord_df_A_1.drop(slv_coord_df_A_1.columns[[0, 1]], axis=1, inplace=True)
    #             slv_coord_df_D_1.drop(slv_coord_df_D_1.columns[[0, 1]], axis=1, inplace=True)
    #
    #             if not Energy_types:
    #                 Energy_types = list(slv_coord_df_R_1.columns)
    #
    #             if not Coordinate_difs:
    #                 for energy_type in Energy_types:
    #                     Coordinate_difs[energy_type] = {'MAX':{'R':{}, 'A':{}, 'D':{}},
    #                                                     'RMS':{'R':{}, 'A':{}, 'D':{}}}
    #             for energy_type in Energy_types:
    #                 Coordinate_difs[energy_type]['MAX']['R'][sslv] = []
    #                 Coordinate_difs[energy_type]['MAX']['A'][sslv] = []
    #                 Coordinate_difs[energy_type]['MAX']['D'][sslv] = []
    #
    #                 Coordinate_difs[energy_type]['RMS']['R'][sslv] = []
    #                 Coordinate_difs[energy_type]['RMS']['A'][sslv] = []
    #                 Coordinate_difs[energy_type]['RMS']['D'][sslv] = []
    #
    #             for slv_coordinate_df2 in solvent_data['Coordinates']:
    #                 slv_coord_df_R_2 = slv_coordinate_df2[slv_coordinate_df2['Name'].str.startswith('R')].reset_index(drop=True)
    #                 slv_coord_df_A_2 = slv_coordinate_df2[slv_coordinate_df2['Name'].str.startswith('A')].reset_index(drop=True)
    #                 slv_coord_df_D_2 = slv_coordinate_df2[slv_coordinate_df2['Name'].str.startswith('D')].reset_index(drop=True)
    #
    #                 slv_coord_df_R_2.drop(slv_coord_df_R_2.columns[[0, 1]], axis=1, inplace=True)
    #                 slv_coord_df_A_2.drop(slv_coord_df_A_2.columns[[0, 1]], axis=1, inplace=True)
    #                 slv_coord_df_D_2.drop(slv_coord_df_D_2.columns[[0, 1]], axis=1, inplace=True)
    #
    #                 slv_coord_df_D_1, slv_coord_df_D_2 = Degree_correction(slv_coord_df_D_1, slv_coord_df_D_2)
    #
    #                 slv_R_difference = slv_coord_df_R_2.combine(slv_coord_df_R_1, lambda s1, s2: s1.combine(s2, self.difference_function))
    #                 slv_A_difference = slv_coord_df_A_2.combine(slv_coord_df_A_1, lambda s1, s2: s1.combine(s2, self.difference_function))
    #                 slv_D_difference = slv_coord_df_D_2.combine(slv_coord_df_D_1, lambda s1, s2: s1.combine(s2, self.difference_function))
    #
    #                 #MAX
    #                 slv_R_diff_MAX = slv_R_difference.max().to_dict()
    #                 slv_A_diff_MAX = slv_A_difference.max().to_dict()
    #                 slv_D_diff_MAX = slv_D_difference.max().to_dict()
    #
    #                 #RMS
    #                 slv_R_diff_RMS = {col: Root_Mean_Square(slv_R_difference[col]) for col in slv_R_difference}
    #                 slv_A_diff_RMS = {col: Root_Mean_Square(slv_A_difference[col]) for col in slv_A_difference}
    #                 slv_D_diff_RMS = {col: Root_Mean_Square(slv_D_difference[col]) for col in slv_D_difference}
    #
    #                 for energy_type in Energy_types:
    #                     Coordinate_difs[energy_type]['MAX']['R'][sslv].append(slv_R_diff_MAX[energy_type])
    #                     Coordinate_difs[energy_type]['MAX']['A'][sslv].append(slv_A_diff_MAX[energy_type])
    #                     Coordinate_difs[energy_type]['MAX']['D'][sslv].append(slv_D_diff_MAX[energy_type])
    #
    #                     Coordinate_difs[energy_type]['RMS']['R'][sslv].append(slv_R_diff_RMS[energy_type])
    #                     Coordinate_difs[energy_type]['RMS']['A'][sslv].append(slv_A_diff_RMS[energy_type])
    #                     Coordinate_difs[energy_type]['RMS']['D'][sslv].append(slv_D_diff_RMS[energy_type])
    #
    #         for energy_type in Energy_types:
    #             Coordinate_difs[energy_type]['MAX']['R'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['R'], index=self.Short_Solvents)
    #             Coordinate_difs[energy_type]['MAX']['A'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['A'], index=self.Short_Solvents)
    #             Coordinate_difs[energy_type]['MAX']['D'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['D'], index=self.Short_Solvents)
    #
    #             Coordinate_difs[energy_type]['RMS']['R'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['R'], index=self.Short_Solvents)
    #             Coordinate_difs[energy_type]['RMS']['A'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['A'], index=self.Short_Solvents)
    #             Coordinate_difs[energy_type]['RMS']['D'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['D'], index=self.Short_Solvents)
    #
    #         self.All_solute_coordinate_diffs[slu] = Coordinate_difs
    #         self.Coordinate_Energy_types = Energy_types
