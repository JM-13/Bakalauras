import os
import pandas as pd
import numpy as np
import math

from constants import *
from tools.catalogue_data import Catalogue
from tools.retrieve_values import Retvieve

Solute_order = ['BDP-nitrophenyl', 'BDP-phenyl', 'BDP-PP-phenyl', 'BDP-PP-phenyl_OMes_backward', 'BDP-PP-phenyl_OMes_forward']
Solvent_order = ['Acetone', 'DiMethylFormamide', 'DiMethylSulfoxide', 'Methanol', 'TetraHydroFuran']

energy_level_names = {'optS0':'S0',
                      'optS1':'S1',
                      'optR1':'R1',
                      'tdS0':'altS1'}


Solvent_order_short = []
for i in Solvent_order:
    Solvent_order_short.append(Solvent_to_shorten[i])

Solute_order_short = []
for i in Solute_order:
    Solute_order_short.append(Solute_to_shorten[i])

Folders = ['/home/deck/Desktop/Bak_Baig/Program/Data/BDP-nitrophenyl',
           '/home/deck/Desktop/Bak_Baig/Program/Data/BDP-phenyl',
           '/home/deck/Desktop/Bak_Baig/Program/Data/BDP-PP-phenyl',
           '/home/deck/Desktop/Bak_Baig/Program/Data/BDP-PP-phenyl_OMes_backward',
           '/home/deck/Desktop/Bak_Baig/Program/Data/BDP-PP-phenyl_OMes_forward']

Files_by_solute = Catalogue(Folders, save_as_json=False).Processed_Data_Files([0,1,2])
Files_by_solvent = Catalogue(Folders, save_as_json=False).Processed_Data_Files([0,2,1])

def find_RMS_value(num_list):

    number_of_vals = len(num_list)

    total_sq_sum = 0
    for v in num_list:
        total_sq_sum += v**2

    RMS = math.sqrt(total_sq_sum / number_of_vals)
    return RMS

def differencial_matrix(values, labels, decimal_places=None):

    diff_matrix = pd.DataFrame(index=labels, columns=labels, dtype=float)

    for i in range(len(values)):
        for j in range(len(values)):
            if i == j:
                diff = 0.0
            else:
                diff = abs(values[i] - values[j]) / max(abs(values[i]), abs(values[j])) * 100
            diff_matrix.iloc[i, j] = diff

    non_zero_values = diff_matrix[diff_matrix != 0].values.flatten()
    non_zero_values = non_zero_values[~np.isnan(non_zero_values)]  # remove NaNs just in case

    if not decimal_places:
        if len(non_zero_values) == 0:
            decimal_places = 0
        else:
            smallest = np.min(non_zero_values)
            decimal_places = max(0, -int(math.floor(math.log10(smallest))) + 1)

    def smart_format(x):
        rounded = round(x, decimal_places)
        return "0\\%" if rounded == 0 else f"{rounded:.{decimal_places}f}\\%"

    formatted_df = diff_matrix.apply(lambda col: col.map(smart_format))

    return formatted_df


def differencial_matrix_list(data, output, labels=None, decimal_places=None, is_D=False):

    #data is simple dict

    keys = list(data.keys())

    if not labels:
        labels = keys
    diff_matrix = pd.DataFrame(index=labels, columns=labels, dtype=float)

    # Compute max % difference between each pair of value lists
    for i in range(len(labels)):
        for j in range(len(labels)):
            if i == j:
                diff_matrix.iloc[i, j] = 0.0
            else:
                vals1 = np.array(data[keys[i]])
                vals2 = np.array(data[keys[j]])

                if is_D:
                    new_vals1 = vals1.copy()
                    new_vals2 = vals2.copy()

                    for idx in range(len(vals1)):
                        new_D = new_vals1[idx]
                        old_D = new_vals2[idx]
                        scale = 10

                        if new_D > 0 and old_D < 0:
                            if abs(new_D - old_D) > 6:
                                temp = round(new_D*scale - (2 * math.pi)*scale)
                                new_vals1[idx] = temp/scale
                        elif new_D < 0 and old_D > 0:
                            if abs(old_D - new_D) > 6:
                                temp = round(old_D*scale - (2 * math.pi)*scale)
                                new_vals2[idx] = temp/scale

                    vals1 = new_vals1
                    vals2 = new_vals2

                diffs = abs(vals1 - vals2) / np.maximum(abs(vals1), abs(vals2))

                if output == 'MAX':
                    diff = np.max(diffs * 100)

                elif output == 'RMS':
                    rel_diff_squared = diffs ** 2
                    diff = np.sqrt(np.mean(rel_diff_squared)) * 100

                diff_matrix.iloc[i, j] = diff

    non_zero_values = diff_matrix[diff_matrix != 0].values.flatten()
    non_zero_values = non_zero_values[~np.isnan(non_zero_values)]  # remove NaNs just in case

    if not decimal_places:
        if len(non_zero_values) == 0:
            decimal_places = 0
        else:
            smallest = np.min(non_zero_values)
            decimal_places = max(0, -int(math.floor(math.log10(smallest))) + 1)

    def smart_format(x):
        rounded = round(x, decimal_places)
        return "0\\%" if rounded == 0 else f"{rounded:.{decimal_places}f}\\%"

    formatted_df = diff_matrix.apply(lambda col: col.map(smart_format))

    return formatted_df









# === Build LaTeX document ===
latex_doc = r"""\documentclass{article}
\usepackage{booktabs}
\usepackage{geometry}
\geometry{margin=1in}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{graphicx}

\begin{document}
"""

for solute in Solute_order:
    short_solute = Solute_to_shorten[solute]


    Data = {'Energies':{'optS0':[],
                        'tdS0' :[],
                        'optS1':[],
                        'optR1':[]},

            'Coordinates':{'optS0':{'R':{},'A':{},'D':{}},
                           'optS1':{'R':{},'A':{},'D':{}},
                           'optR1':{'R':{},'A':{},'D':{}}}
           }

    for solvent in Solvent_order:
        dat = Retvieve(Files_by_solute['DATA'][solute][solvent]).Regular_data()
        for e in Data['Energies']:
            Data['Energies'][e].append(dat['Energys'][e])

        df = dat['Coordinates']
        df_R = df[df['Name'].str.startswith('R')]
        df_A = df[df['Name'].str.startswith('A')]
        df_D = df[df['Name'].str.startswith('D')]

        cords = {'R':df_R, 'A':df_A, 'D':df_D}

        for c in cords:
            Data['Coordinates']['optS0'][c][solvent] = list(cords[c]['optS0'])
            Data['Coordinates']['optS1'][c][solvent] = list(cords[c]['optS1'])
            Data['Coordinates']['optR1'][c][solvent] = list(cords[c]['optR1'])


    for e in Data['Energies']:

        e_name = energy_level_names[e]

        if e != 'tdS0':
            energy = np.array(Data['Energies'][e])
            diff_matrix = differencial_matrix(values=energy, labels=Solvent_order_short)

            # # Print the result
            # print(f'Energy {e} comparison between solvents:')
            # print(diff_matrix)
            # print('\n')

            table_body = diff_matrix.to_latex(index=True, escape=False, column_format="l" + "c"*len(diff_matrix.columns))
            caption = f"{short_solute} energy {e_name} difference between solvents"
            label = f"tab:{short_solute}_{e_name}"

            latex_doc += rf"""
  \begin{{table}}[H]
    \centering
    \small
    \caption{{{caption}}}
    {table_body}
  \end{{table}}
"""

            print(f'Comparing coordinate for energy value {e}:')

            global_caption = f"{short_solute} MAX and RMS coordinate differences between solvents for energy {e_name}"

            latex_data = {'R':{'MAX':{'table_body':None, 'caption':None, 'label':None}, 'RMS':{'table_body':None, 'caption':None, 'label':None}},
                          'A':{'MAX':{'table_body':None, 'caption':None, 'label':None}, 'RMS':{'table_body':None, 'caption':None, 'label':None}},
                          'D':{'MAX':{'table_body':None, 'caption':None, 'label':None}, 'RMS':{'table_body':None, 'caption':None, 'label':None}}}

            for c in Data['Coordinates'][e]:
                output = 'MAX'

                if c == 'D':
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solvent_order_short, output=output, is_D=True)
                else:
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solvent_order_short, output=output, is_D=False)

                # # Print the result
                # print(f'For energy {e} coordinate {c} comparison between solvents for {short_solute} taking the {output} difference value:')
                # print(diff_matrix)
                # print('\n')

                table_body_MAX = diff_matrix.to_latex(index=True, escape=False, column_format="l" + "c"*len(diff_matrix.columns))
                caption_MAX = f"MAX {c}"
                label_MAX = f"tab:{short_solute}_{e_name}_{c}_MAX"

                latex_data[c]['MAX']['table_body'] = table_body_MAX
                latex_data[c]['MAX']['caption'] = caption_MAX
                latex_data[c]['MAX']['label'] = label_MAX


                output = 'RMS'

                if c == 'D':
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solvent_order_short, output=output, is_D=True)
                else:
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solvent_order_short, output=output, is_D=False)

                # # Print the result
                # print(f'For energy {e} coordinate {c} comparison between solvents for {short_solute} taking the {output} difference value:')
                # print(diff_matrix)
                # print('\n')

                table_body_RMS = diff_matrix.to_latex(index=True, escape=False, column_format="l" + "c"*len(diff_matrix.columns))
                caption_RMS = f"RMS {c}"
                label_RMS = f"tab:{short_solute}_{e_name}_{c}_RMS"

                latex_data[c]['RMS']['table_body'] = table_body_RMS
                latex_data[c]['RMS']['caption'] = caption_RMS
                latex_data[c]['RMS']['label'] = label_RMS

            latex_doc += rf"""
  \begin{{table}}[H]
    \centering
    \caption{{{global_caption}}}

    % -------- Row 1 --------
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['R']['MAX']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['R']['MAX']['table_body']}
      }}
    \end{{subtable}}
    \hfill
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['R']['RMS']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['R']['RMS']['table_body']}
      }}
    \end{{subtable}}

    \vspace{{1em}}

    % -------- Row 2 --------
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['A']['MAX']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['A']['MAX']['table_body']}
      }}
    \end{{subtable}}
    \hfill
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['A']['RMS']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['A']['RMS']['table_body']}
      }}
    \end{{subtable}}

    \vspace{{1em}}

    % -------- Row 3 --------
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['D']['MAX']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['D']['MAX']['table_body']}
      }}
    \end{{subtable}}
    \hfill
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['D']['RMS']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['D']['RMS']['table_body']}
      }}
    \end{{subtable}}

  \end{{table}}
"""
            latex_doc += rf"""
  \clearpage
"""

latex_doc += r"\end{document}"

# print(latex_doc)
with open("All_solvent_differential_tables.tex", "w") as f:
    f.write(latex_doc)


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################



# === Build LaTeX document ===
latex_doc = r"""\documentclass{article}
\usepackage{booktabs}
\usepackage{geometry}
\geometry{margin=1in}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{graphicx}

\begin{document}
"""

for solvent in Solvent_order:
    short_solvent = Solvent_to_shorten[solvent]


    Data = {'Energies':{'optS0':[],
                        'tdS0' :[],
                        'optS1':[],
                        'optR1':[]},

            'Coordinates':{'optS0':{'R':{},'A':{},'D':{}},
                           'optS1':{'R':{},'A':{},'D':{}},
                           'optR1':{'R':{},'A':{},'D':{}}}
           }

    for solute in Solute_order:
        dat = Retvieve(Files_by_solvent['DATA'][solvent][solute]).Regular_data()
        for e in Data['Energies']:
            Data['Energies'][e].append(dat['Energys'][e])

        df = dat['Coordinates']
        df_R = df[df['Name'].str.startswith('R')]
        df_A = df[df['Name'].str.startswith('A')]
        df_D = df[df['Name'].str.startswith('D')]

        cords = {'R':df_R, 'A':df_A, 'D':df_D}

        for c in cords:
            Data['Coordinates']['optS0'][c][solute] = list(cords[c]['optS0'])
            Data['Coordinates']['optS1'][c][solute] = list(cords[c]['optS1'])
            Data['Coordinates']['optR1'][c][solute] = list(cords[c]['optR1'])


    for e in Data['Energies']:

        e_name = energy_level_names[e]

        if e != 'tdS0':
            energy = np.array(Data['Energies'][e])
            diff_matrix = differencial_matrix(values=energy, labels=Solute_order_short)

            # # Print the result
            # print(f'Energy {e} comparison between solvents:')
            # print(diff_matrix)
            # print('\n')

            table_body = diff_matrix.to_latex(index=True, escape=False, column_format="l" + "c"*len(diff_matrix.columns))
            caption = f"{short_solvent} energy {e_name} difference between solutes"
            label = f"tab:{short_solvent}_{e_name}"

            latex_doc += rf"""
  \begin{{table}}[H]
    \centering
    \small
    \caption{{{caption}}}
    {table_body}
  \end{{table}}
"""

            print(f'Comparing coordinate for energy value {e}:')

            global_caption = f"{short_solvent} MAX and RMS coordinate differences between solutes for energy {e_name}"

            latex_data = {'R':{'MAX':{'table_body':None, 'caption':None, 'label':None}, 'RMS':{'table_body':None, 'caption':None, 'label':None}},
                          'A':{'MAX':{'table_body':None, 'caption':None, 'label':None}, 'RMS':{'table_body':None, 'caption':None, 'label':None}},
                          'D':{'MAX':{'table_body':None, 'caption':None, 'label':None}, 'RMS':{'table_body':None, 'caption':None, 'label':None}}}

            for c in Data['Coordinates'][e]:
                output = 'MAX'

                if c == 'D':
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solute_order_short, output=output, is_D=True)
                else:
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solute_order_short, output=output, is_D=False)

                # # Print the result
                # print(f'For energy {e} coordinate {c} comparison between solvents for {short_solute} taking the {output} difference value:')
                # print(diff_matrix)
                # print('\n')

                table_body_MAX = diff_matrix.to_latex(index=True, escape=False, column_format="l" + "c"*len(diff_matrix.columns))
                caption_MAX = f"MAX {c} difference"
                label_MAX = f"tab:{short_solvent}_{e_name}_{c}_MAX"

                latex_data[c]['MAX']['table_body'] = table_body_MAX
                latex_data[c]['MAX']['caption'] = caption_MAX
                latex_data[c]['MAX']['label'] = label_MAX


                output = 'RMS'

                if c == 'D':
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solute_order_short, output=output, is_D=True)
                else:
                    diff_matrix = differencial_matrix_list(Data['Coordinates'][e][c], labels=Solute_order_short, output=output, is_D=False)

                # # Print the result
                # print(f'For energy {e} coordinate {c} comparison between solvents for {short_solute} taking the {output} difference value:')
                # print(diff_matrix)
                # print('\n')

                table_body_RMS = diff_matrix.to_latex(index=True, escape=False, column_format="l" + "c"*len(diff_matrix.columns))
                caption_RMS = f"RMS {c} difference"
                label_RMS = f"tab:{short_solvent}_{e_name}_{c}_RMS"

                latex_data[c]['RMS']['table_body'] = table_body_RMS
                latex_data[c]['RMS']['caption'] = caption_RMS
                latex_data[c]['RMS']['label'] = label_RMS

            latex_doc += rf"""
  \begin{{table}}[H]
    \centering
    \caption{{{global_caption}}}

    % -------- Row 1 --------
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['R']['MAX']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['R']['MAX']['table_body']}
      }}
    \end{{subtable}}
    \hfill
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['R']['RMS']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['R']['RMS']['table_body']}
      }}
    \end{{subtable}}

    \vspace{{1em}}

    % -------- Row 2 --------
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['A']['MAX']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['A']['MAX']['table_body']}
      }}
    \end{{subtable}}
    \hfill
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['A']['RMS']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['A']['RMS']['table_body']}
      }}
    \end{{subtable}}

    \vspace{{1em}}

    % -------- Row 3 --------
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['D']['MAX']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['D']['MAX']['table_body']}
      }}
    \end{{subtable}}
    \hfill
    \begin{{subtable}}[t]{{0.48\textwidth}}
      \centering
      \caption{{{latex_data['D']['RMS']['caption']}}}
      \resizebox{{\textwidth}}{{!}}{{
        {latex_data['D']['RMS']['table_body']}
      }}
    \end{{subtable}}

  \end{{table}}
"""
            latex_doc += rf"""
  \clearpage
"""

latex_doc += r"\end{document}"

# print(latex_doc)
with open("All_solute_differential_tables.tex", "w") as f:
    f.write(latex_doc)

