import sys
from sys import argv
import getopt
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None  # default='warn'

def picture(datafile):
    database = pd.read_excel(datafile)
    database.drop(database.index[database.Length < 1], inplace = True)
    H000200 = []
    H200400 = []
    H400600 = []
    H6001000 = []

    labelz = ["h = 0-200mm", "h = 200-400 mm", "h = 400-600 mm", "h = 600-1000 mm"]

    for i in range(len(database)):
        h = int(database.h.iloc[i])
        l = database.Length.iloc[i]
        if h <= 200:
            H000200.append(l)
        elif h > 200 and h <= 400:
            H200400.append(l)
        elif h > 400 and h <= 600:
            H400600.append(l)
        else:
            H6001000.append(l)    

    arr = [H000200, H200400, H400600, H6001000]

    plt.figure(figsize = (12,5))
    plt.hist(arr, label = labelz, bins = np.arange(0, 16,1))
    plt.grid()
    plt.xticks(np.arange(0, 16,1))
    plt.legend()
    plt.title("Availability of structural elements of a certain structural height and length")
    plt.ylabel(f"Occurance in dataset (#)")
    plt.xlabel(f"Element length (bins) --> e.g. 1 to 2 meters")
    plt.savefig(f"Initialviewofelements", bbox_inches = 'tight')

if __name__ == "__main__":
    # myfunc(argv)
    picture(argv[1])