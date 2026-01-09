from analyze import Retrieve
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

Data_folder = '/home/deck/Desktop/Bak_Baig/Program/Data'

Solvent_to_shorten = {'Acetone':'Acetone',
                      'DiMethylFormamide':'DMF',
                      'DiMethylSulfoxide':'DMSO',
                      'Methanol':'Methanol',
                      'TetraHydroFuran':'THF'}

Solute_to_shorten = {'BDP-nitrophenyl':'BDP-NPh',
                     'BDP-phenyl':'BDP-Ph',
                     'BDP-PP-phenyl':'BDP-PPPh',
                     'BDP-PP-phenyl_OMes_backward':'OMe-out',
                     'BDP-PP-phenyl_OMes_forward':'OMe-in'}



Graph_settings = {'label': list(Solvent_to_shorten.values()),
                  'zorder':[2, 2, 2, 2, 2],
                  'alpha': [0.4, 0.4, 0.4, 0.4, 0.4],
                  'color': [],

                  'line_width':[],
                  'linestyle': [],

                  'facecolors':['blue', 'orange', 'green', 'red', 'purple'],
                  'edgecolors':[],
                  'dot_size':  [],
                  'marker':    []
                  }


def ScanAnalasys(All_data_folder, solvents, solutes, graph_config):

    for slu in solutes:
        folder_path = os.path.join(All_data_folder, slu)

        scan_file_data = Retrieve(multiple=True)

        for slv in solvents:
            for file in os.listdir(folder_path):
                if file.endswith("_SCAN.txt") and (slv in file):
                    filepath = os.path.join(folder_path, file)

                    scan_file_data.scan_data(filepath)

        scan_file_data = scan_file_data.return_scan_data()

        fig, ax = plt.subplots()
        ax.grid(True)
        ax.set_title(f"{slu} energy planes")

        for num, slv in enumerate(solvents):

            kwargs = {}
            for kwa in graph_config:
                if graph_config[kwa]:
                    kwargs[kwa] = graph_config[kwa][num]

            data = pd.concat([scan_file_data['fa'][num], scan_file_data['fb'][num]])

            degree_axis_values_average = np.array(data['D_average'])
            indices = np.argsort(degree_axis_values_average)

            degree_axis_values_fixed   = np.array(data['D_fixed'])[indices]
            degree_axis_values_average = np.array(data['D_average'])[indices]
            energy_axis_values_S0 = np.array(data['S0'])[indices]
            energy_axis_values_S1 = np.array(data['S1'])[indices]


            # ax.scatter(degree_axis_values_average, energy_axis_values_S0, **kwargs)
            ax.scatter(degree_axis_values_average, energy_axis_values_S1, **kwargs)

        ax.legend()
        plt.show()
        plt.close(fig)

# max_S1 = fb_data['S1'].max() #+ 0.0003
#
# filterred_fa_S1 = []
# filterred_fa_Da = []
# for x,y in zip(fa_data['D_average'], fa_data['S1']):
#     if y <= max_S1:
#         filterred_fa_Da.append(x)
#         filterred_fa_S1.append(y)


ScanAnalasys(Data_folder,
             list(Solvent_to_shorten.keys()),
             list(Solute_to_shorten.keys()),
             Graph_settings)






