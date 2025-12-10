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
