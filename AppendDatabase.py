import sys
from sys import argv
import getopt
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

def Append(FileName, fy, validated, Trader, identifier, Toxic, Cday, Ctest, Crecalc, Crecalcp, Cclean, Year, SSP):
    data = pd.read_excel(FileName, skiprows = 10)
    data.drop(data.index[data.Length < 1], inplace = True)
    NewElems = pd.read_excel("NewElements.xlsx")
    database = pd.read_excel("Database.xlsx")
    df = pd.DataFrame(index = None, columns = ["Identifier", "Cross-section", "Length", "Fy", "G", "A", "h", "Iy", "Iz", "W", "bucmaj", "bucmin", "TestCosts", "StoreCosts", "ToxicFabCosts"])
    df["Cross-section"] = data.Type
    df["Length"] = data.Length
    df["Identifier"] = identifier
    
    for i in range(len(df)):
        df.G.iloc[i] = NewElems.Weight.iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.A.iloc[i] = NewElems.Area.iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.Iy.iloc[i] = NewElems["Second moment of area y"].iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.Iz.iloc[i] = NewElems["Second moment of area z"].iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.W.iloc[i] = NewElems["Plastic section modulus"].iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.bucmaj.iloc[i] = NewElems["Buckling about major axis y-y"].iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.bucmin.iloc[i] = NewElems["Buckling about minor axis z-z"].iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
        df.h.iloc[i] = NewElems.Height.iloc[NewElems.index[df["Cross-section"].iloc[i] == NewElems.Profile]].to_list()[0]
    Wtot = np.sum(df.Length * df.G)
    
    ##Code for determining test costs scenario 1a
    testunits = len(df["Cross-section"].unique())
    Elemcount = len(df)
    Tc = (testunits/20) * float(Cday) + testunits * float(Ctest) 
    Tcpe = Tc / Elemcount
    
    if validated == "Yes":
        YS = int(fy)
        df.TestCosts = 0
    else:
        if Year == "1955-1972":
            YS = 235
        elif Year == "1990-1997":
            YS = 235
        elif Year == "1997-2005":
            YS = 235
        elif Year == "2005-present":
            YS = 235
        else:
            YS = 200
        df.TestCosts = Tcpe
    
    if Trader == "Yes":
        df.StoreCosts = (0.075 + 0.35 * (float(SSP) + 0.20)) * df.G * df.Length
    else:
        df.StoreCosts = 0
            
            
    df.Fy = YS
    
    if Toxic == "Yes":
        df.ToxicFabCosts = float(Cclean) * 4 * df.G * df.Length
    else:
        df.ToxicFabCosts = 0

    database = database.append(df, ignore_index = True)
    database.to_excel("Database.xlsx", index = False)
    print(np.round(Wtot / 1000, 1))

if __name__ == "__main__":
    # myfunc(argv)
    Append(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12], argv[13])