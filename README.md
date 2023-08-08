# MsC_Thesis_DonorSteel
Codes and application developed during my Master's Thesis in Civil Engineering.


ZIP-FILE
-------------------------------------
The zip-file contains the tool SteelIT. To be able to use the tool you will need to install python and rhino.
For the python installation it is important to include the following packages:
- numpy
- pandas
- itertools

The tool will ask for the path to your python executable. Normally this path is: C:\Users\USER\AppData\Local\Microsoft\WindowsApps\python.exe

Rhino is used with grasshopper. Grasshopper is included in Rhino 7. In addition the Karamba3D plug-in is required if the user wants to visualise the solution that are generated with the tool.

The zip itself contains a 'SteelIT' folder and a shortcut to the executable. To simply run the tool, run the executable. To access the project files or solution files head to: SteelIT\SteelIT\SteelIT\bin\Debug\net6.0-windows. Here the following files can be accessed.
Excel files from the generated solutions
Donor steel database
Graphs generated by the tool 

OTHER FILES
--------------------------------------
The other files are the source codes from the too. A short description of all filenames is provided.

- AppenDatabase.py: the python script used to create the database of donor elements within the tool.
- DesignCreator.gh: the grasshopper script that allows the user to create an office building geometry.
- GetCostImpactDetails.py: python script to show the cost and impact breakdown in the tool.
- InitialView.py: Creates the visualisation of both the database and used elements in a reclaimed steel design.
- TheCodeV2.py: This code performs the structural analysis, assignment optimisation and cost/impact assessment of the designs.
