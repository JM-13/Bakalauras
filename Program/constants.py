from tools.calculation_functions import DifferenceFunction

# Folders of files
S0_folder   = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/1-S0'
S1_folder   = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/2-S1'
Scan_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/3-scans'

Save_optimized_coordinates = False #they are less accurate and double rounded
difference_function = DifferenceFunction.symmetric_percentage_difference


Filename_to_Solvent = {'acet':'Acetone',
                       'dmfm':'DiMethylFormamide',
                       'dmso':'DiMethylSulfoxide',
                       'meoh':'Methanol',
                       'thfu':'TetraHydroFuran',
                       'tolu':'Toluene'}

# #Measured in Polarity Index (P´)"
# Solvent_to_Polarity = {'Toluene':2.4,
#                        'TetraHydroFuran':4.0,
#                        'Acetone':5.1,
#                        'Methanol':5.1,
#                        'DiMethylFormamide':6.4,
#                        'DiMethylSulfoxide':7.2,
#                        }

Filename_to_Solute = {'bdp-pnphen' :'BDP-nitrophenyl',
                      'bdp-phenyl' :'BDP-phenyl',
                      'ppp'        :'BDP-PP-phenyl',
                      'ppp-ome-out':'BDP-PP-phenyl_OMes_backward',
                      'ppp-ome-in' :'BDP-PP-phenyl_OMes_forward'}

Solvent_to_shorten = {'Acetone':'Acetone',
                      'DiMethylFormamide':'DMF',
                      'DiMethylSulfoxide':'DMSO',
                      'Methanol':'Methanol',
                      'TetraHydroFuran':'THF'}

# Solute_to_shorten = {'BDP-nitrophenyl':'BDP-NPh',
#                      'BDP-phenyl':'BDP-Ph',
#                      'BDP-PP-phenyl':'BDP-PPPh',
#                      'BDP-PP-phenyl_OMes_backward':'BDP-PPPh-OMe-out',
#                      'BDP-PP-phenyl_OMes_forward':'BDP-PPPh-OMe-in'}

Solute_to_shorten = {'BDP-nitrophenyl':'BDP-NPh',
                     'BDP-phenyl':'BDP-Ph',
                     'BDP-PP-phenyl':'BDP-PPPh',
                     'BDP-PP-phenyl_OMes_backward':'OMe-out',
                     'BDP-PP-phenyl_OMes_forward':'OMe-in'}

#Atoms are numbered differently
#keys are:   BDP-NPh and BDP-Ph
#values are: BDP-PPPh, BDP-PPPh-OMe-out and BDP-PPPh-OMe-in
#4,6,8,13,15,17 are H. 6:41 and 15:30 replace H with C
Atom_number_conversion = { 1:19,
                           2:17,
                           3:16,
                           4:29,
                           5:15,
                           6:41,
                           7:14,
                           8:28,
                           9:13,
                          10:12,
                          11:1,
                          12:2,
                          13:21,
                          14:3,
                          15:30,
                          16:4,
                          17:22,
                          18:5,
                          19:18,
                          20:20,
                          21:11}

Scan_Graph_settings = {
    'Line':{   'label': [],
               'zorder':[3, 3, 3, 3, 3],
               'alpha': [0.8, 0.8, 0.8, 0.8, 0.8],

               'color': ['blue', 'orange', 'green', 'red', 'purple'],
               'linewidths':[],
               'linestyle': []
                        },

    'Scatter':{'label': list(Solvent_to_shorten.values()),
               'zorder':[4, 4, 4, 4, 4],
               'alpha': [0.8, 0.8, 0.8, 0.8, 0.8],

               'facecolors':['blue', 'orange', 'green', 'red', 'purple'],
               's':         [15,15,15,15,15],
               'edgecolors':['black', 'black', 'black', 'black', 'black'],
               'linewidths':[0.5,0.5,0.5,0.5,0.5],
               'marker':    []
                        },

    'Minima_location_dots':{'label': [],
                            'zorder':[5, 5, 5, 5, 5],
                            'alpha': [],

                            'facecolors':['white', 'white', 'white', 'white', 'white'],
                            's':         [15,15,15,15,15],
                            'edgecolors':['black', 'black', 'black', 'black', 'black'],
                            'linewidths':[0.5,0.5,0.5,0.5,0.5],
                            'marker':    []
                             },

    'Minima_energy_lines':{ 'label': [],
                            'zorder':[2, 2, 2, 2, 2],
                            'alpha': [],

                            'color': ['blue', 'orange', 'green', 'red', 'purple'],
                            'linewidths':[],
                            'linestyle': ['--','--','--','--','--']
                             },

    'Fit_curvyness':{'S0':14,
                     'S1':12},

    'Xlimits':[-60, 10],

    'Ylimits':{'S0':{'start':[0.000,0.000,0.000,0.000,0.000],
                     'stop' :[0.033,0.033,0.026,0.026,0.026],
                     'step' :[0.003,0.003,0.002,0.002,0.002]},

               'S1':{'start':[0.093,0.097,0.082,0.076,0.076],
                     'stop' :[0.102,0.106,0.094,0.089,0.089],
                     'step' :[0.001,0.001,0.001,0.001,0.001]}
              }
}
