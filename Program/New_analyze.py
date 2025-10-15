import os
import numpy as np
import pandas as pd
import math

def Degree_correction(df1, df2, scale=10000):
    pi2 = 2 * math.pi

    a = df1.to_numpy()
    b = df2.to_numpy()

    mask1 = (a > 0) & (b < 0) & (np.abs(a - b) > 6)
    mask2 = (a < 0) & (b > 0) & (np.abs(b - a) > 6)

    a_corrected = np.where(mask1, np.round(a*scale - pi2*scale)/scale, a)
    b_corrected = np.where(mask2, np.round(b*scale - pi2*scale)/scale, b)

    df1_corrected = pd.DataFrame(a_corrected, columns=df1.columns)
    df2_corrected = pd.DataFrame(b_corrected, columns=df2.columns)

    return df1, df2

def Root_Mean_Square(series):
    return np.sqrt((series**2).mean())


class Difference_Between:
    def __init__(self, All_Data_folder, difference_function):

        self.All_Data_folder = All_Data_folder
        self.difference_function = difference_function

    def Solvents(self, solutes, solvents):

        for slu in solutes:
            solute_folder = os.path.join(self.All_Data_folder, slu)
            solvent_data = Retvieve(multiple=True)

        for slv in solvents:
            for file in os.listdir(solute_folder):
                if file.endswith("_DATA.txt") and (slv in file):
                    solvent_file = os.path.join(solute_folder, file)
                    solvent_data.Regular_data(solvent_file)
                    break

        solvent_data = solvent_data.Return_multiple()

        #ENERGIES
        Energy_types = []
        Energy_difs = {}
        for energy_type in solvent_data['Energys']:
            Energy_difs[energy_type] = {}
            Energy_types.append[energy_type]

            for slv, slv_energy_value in zip(solvents, solvent_data['Energys'][energy_type]):
                slv_difference_column = []

                for energy in solvent_data['Energys'][energy_type]:
                    slv_difference = self.difference_function(energy, slv_energy_value)
                    slv_difference_column.append(slv_difference)

                Energy_difs[energy_type][slv] = slv_difference_column

            solvent_difference_df = pd.DataFrame(Energy_difs[energy_type]).reindex(solvents)
            Energy_difs[energy_type] = solvent_difference_df

        #COORDINATES
        Coordinate_difs = {}
        for energy_type in Energy_types:
            Coordinate_difs[energy_type] = {'MAX':{'R':{}, 'A':{}, 'D':{}},
                                            'RMS':{'R':{}, 'A':{}, 'D':{}}
                                            }
        for slv, slv_coordinate_df in zip(solvents, solvent_data['Coordinates']):

            for energy_type in Coordinate_difs:
                Coordinate_difs[energy_type]['MAX']['R'][slv] = []
                Coordinate_difs[energy_type]['MAX']['A'][slv] = []
                Coordinate_difs[energy_type]['MAX']['D'][slv] = []

                Coordinate_difs[energy_type]['RMS']['R'][slv] = []
                Coordinate_difs[energy_type]['RMS']['A'][slv] = []
                Coordinate_difs[energy_type]['RMS']['D'][slv] = []


            slv_coord_df_R_1 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('R')].reset_index(drop=True)
            slv_coord_df_A_1 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('A')].reset_index(drop=True)
            slv_coord_df_D_1 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('D')].reset_index(drop=True)

            slv_coord_df_R_1.drop(df.columns[[0, 1]], axis=1, inplace=True)
            slv_coord_df_A_1.drop(df.columns[[0, 1]], axis=1, inplace=True)
            slv_coord_df_D_1.drop(df.columns[[0, 1]], axis=1, inplace=True)

            for slv_coordinate_df in zip(solvents, solvent_data['Coordinates']):
                slv_coord_df_R_2 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('R')].reset_index(drop=True)
                slv_coord_df_A_2 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('A')].reset_index(drop=True)
                slv_coord_df_D_2 = slv_coordinate_df[slv_coordinate_df['Name'].str.startswith('D')].reset_index(drop=True)

                slv_coord_df_R_2.drop(df.columns[[0, 1]], axis=1, inplace=True)
                slv_coord_df_A_2.drop(df.columns[[0, 1]], axis=1, inplace=True)
                slv_coord_df_D_2.drop(df.columns[[0, 1]], axis=1, inplace=True)


                slv_coord_df_D_1, slv_coord_df_D_2 = Degree_correction(slv_coord_df_D_1, slv_coord_df_D_2)


                slv_R_difference = slv_coord_df_R_2.combine(slv_coord_df_R_1, difference_function)
                slv_A_difference = slv_coord_df_A_2.combine(slv_coord_df_A_1, difference_function)
                slv_D_difference = slv_coord_df_D_2.combine(slv_coord_df_D_1, difference_function)

                #MAX
                slv_R_diff_MAX = slv_R_difference.max().to_dict()
                slv_A_diff_MAX = slv_A_difference.max().to_dict()
                slv_D_diff_MAX = slv_D_difference.max().to_dict()

                #RMS
                slv_R_diff_RMS = {col: Root_Mean_Square(slv_R_difference[col]) for col in slv_R_difference}
                slv_A_diff_RMS = {col: Root_Mean_Square(slv_A_difference[col]) for col in slv_A_difference}
                slv_D_diff_RMS = {col: Root_Mean_Square(slv_D_difference[col]) for col in slv_D_difference}

                for energy_type in Coordinate_difs:
                    if energy_type in slv_coord_df_R_1:
                        Coordinate_difs[energy_type]['MAX']['R'][slv].append(slv_R_diff_MAX[energy_type])
                        Coordinate_difs[energy_type]['MAX']['A'][slv].append(slv_A_diff_MAX[energy_type])
                        Coordinate_difs[energy_type]['MAX']['D'][slv].append(slv_D_diff_MAX[energy_type])

                        Coordinate_difs[energy_type]['RMS']['R'][slv].append(slv_R_diff_RMS[energy_type])
                        Coordinate_difs[energy_type]['RMS']['A'][slv].append(slv_A_diff_RMS[energy_type])
                        Coordinate_difs[energy_type]['RMS']['D'][slv].append(slv_D_diff_RMS[energy_type])

        for energy_type in Coordinate_difs:
            Coordinate_difs[energy_type]['MAX']['R'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['R']).reindex(solvents)
            Coordinate_difs[energy_type]['MAX']['A'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['A']).reindex(solvents)
            Coordinate_difs[energy_type]['MAX']['D'] = pd.DataFrame(Coordinate_difs[energy_type]['MAX']['D']).reindex(solvents)

            Coordinate_difs[energy_type]['RMS']['R'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['R']).reindex(solvents)
            Coordinate_difs[energy_type]['RMS']['A'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['A']).reindex(solvents)
            Coordinate_difs[energy_type]['RMS']['D'] = pd.DataFrame(Coordinate_difs[energy_type]['RMS']['D']).reindex(solvents)

        return Energy_difs, Coordinate_difs

    def Solutes(self):
        return


