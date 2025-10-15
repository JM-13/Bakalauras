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

class Difference_function:
    def __init__(self):
        self.result = None

    def Symmetric_Percentage_Difference(a, b):
         difference = abs(a - b) / max(abs(a), abs(b)) * 100
         self.result = difference

         return difference
