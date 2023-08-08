from sys import argv
from itertools import permutations
from itertools import combinations
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def Getinfo(filename, alt):
    Costfile = pd.read_excel(filename, sheet_name = "costs")
    altcosts = Costfile.loc[int(alt)].values
    Impactfile = pd.read_excel(filename, sheet_name = "impactECI")
    altimpact = Impactfile.loc[int(alt)].values
    Impact2file = pd.read_excel(filename, sheet_name = "impactPP")
    altimpact2 = Impact2file.loc[int(alt)].values
    print(altcosts[0], altcosts[1], altcosts[2], altcosts[3], altcosts[4], altcosts[5], altcosts[6], altcosts[7], altcosts[8], np.sum(altcosts) - altcosts[8], np.sum(altcosts), altimpact[0], altimpact[1], altimpact[2], altimpact[3], altimpact[4], np.sum(altimpact), altimpact2[0], altimpact2[1], np.sum(altimpact2))

if __name__ == "__main__":
    Getinfo(argv[1], argv[2])