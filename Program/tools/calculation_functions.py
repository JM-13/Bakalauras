import pandas as pd
import numpy as np

class DifferenceFunction:
    @staticmethod
    def symmetric_percentage_difference(a, b):
        a_abs = np.abs(a)
        b_abs = np.abs(b)

        denom = np.maximum(a_abs, b_abs)
        denom = np.where(denom == 0, 1, denom)

        return np.abs(a - b) / denom * 100
