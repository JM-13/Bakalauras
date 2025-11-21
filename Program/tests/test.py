import pandas as pd
from extract import Extract
from catalogue import Catalogue


S0_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/1-S0'
S1_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/2-S1'
Scan_folder = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/3-scans'
Procesed_Data_Folder = '/home/deck/Desktop/Bak_Baig/Program/Data/'
D_endings = ['optS0', 'tdS0', 'optS1', 'optR1']
S_endings = ['fa', 'fb']

# Catalogue(folders = [S0_folder, S1_folder],
#           filename_endings = D_endings,
#           save_as_json = True,
#           save_location = Procesed_Data_Folder).Data_Files()
#
# Catalogue(folders = [Scan_folder],
#           filename_endings = S_endings,
#           save_as_json = True,
#           save_location = Procesed_Data_Folder).Scan_Files()


Extractor = Extract('/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/1-S0/bdp-phenyl-dmso-optS0.log')
Energys = Extractor.S0()
Coordinates = Extractor.RAD()
Optimized_Coordinates = Extractor.Optimized_RAD()

Coordinates_combined = {}

# print(len(Coordinates_combined))

print(Energys[-1])
print(Coordinates[-1])
print(Optimized_Coordinates[-1])
