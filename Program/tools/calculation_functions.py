import pandas as pd
import numpy as np

def convert_coordinate_measurement(row, rounding=True, rounding_digits=6):

    if 'R' in row.iloc[0]:
        converted_value = row.iloc[-1] * 0.529177249 #To angstroms
    else:
        converted_value = np.degrees(row.iloc[-1])   #To degrees

    if rounding:
        converted_value = round(converted_value, rounding_digits)

    return converted_value

class DifferenceFunction:
    @staticmethod
    def symmetric_percentage_difference(a, b):
        a_abs = np.abs(a)
        b_abs = np.abs(b)

        denom = np.maximum(a_abs, b_abs)
        denom = np.where(denom == 0, 1, denom)

        return np.abs(a - b) / denom * 100
