from tools.calculation_functions import DifferenceFunction

# Folders of files
S0_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/1-S0'
S1_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/2-S1'
Scan_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/3-scans'

Save_optimized_coordinates = False
Use_angstroms_and_degrees = False
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
