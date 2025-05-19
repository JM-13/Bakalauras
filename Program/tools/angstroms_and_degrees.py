import pandas as pd
import numpy as np

def Convert_measurement_type(row, rounding=True, rounding_digits=6):

    if 'R' in row.iloc[0]:
        converted_value = row.iloc[-1] * 0.529177249 #To angstroms
    else:
        converted_value = np.degrees(row.iloc[-1])   #To degrees

    if rounding:
        converted_value = round(converted_value, rounding_digits)

    return converted_value
