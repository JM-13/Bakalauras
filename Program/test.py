import pandas as pd
from extract import Extract

file_path = '/home/deck/Desktop/Bak_Baig/(JV)/bdpolarity/1-S0/bdp-phenyl-dmfm-optS0.log'


Extractor = Extract(file_path)
Energys = Extractor.S0()
Coordinates = Extractor.RAD()
# Energys2 = Extractor.S1()
print(Coordinates)




