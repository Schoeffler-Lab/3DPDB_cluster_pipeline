## This script was written by Alyce Fields, Larissa Cortes Morales, and Allyn Schoeffler
# Loyola Univeristy New Orleans

#edit this section for pdb directory and desired color for histogram
user_path_to_pdbs1 = "/Users/larissacortes/Desktop/SchoefflerLab/Scripts/1qtm/Mesophiles/Pro/"
user_color = 'blue'

#Importing the PDB package from biopandas
from re import A
from biopandas.pdb import PandasPdb as ppdb
#import glob for user parameters
import glob
#need this to do the math 
import numpy as np 
import pandas as pd
import sys

#formatting user input parameters
path_to_pdbs1 = user_path_to_pdbs1 + "/*.pdb"
pdb_files1 = glob.glob(path_to_pdbs1)

print ('calculating...')

def distanceArrayMaker (pdbFile, cutoff):
    #count rows of DNA dataframe created above to be used in loop below 
    DNAcount = len(DNA_df.index)

    #Importing the PDB package from biopanda
    from biopandas.pdb import PandasPdb as ppdb

    #creating empty arrays to append to later 
    distanceArray = []
    filteredAllDistance = []
    coordArray = []

    #for loop to read coordinates from hetatm dataframe created above. Goes through the loop based on the number of rows counted and appends
    #the coordinates (coord) to empty array coordArray created above
    for n in range(DNAcount):
        x = DNA_df.at[n, "x_coord"]
        y = DNA_df.at[n, "y_coord"]
        z = DNA_df.at[n, "z_coord"]
        coord = [x,y,z]
        coordArray.append(coord)
    #reading pdb file that was created in seperate code named: 
    p_1 = ppdb().read_pdb(pdbFile)
    #Biopandas function that does distance math from DNA coords to all the coordinates from the p_1 pdb that was created with other code  
    for individualCoord in coordArray: 
        distances = p_1.distance(xyz=individualCoord, records=('ATOM',))
    #Appending distances to empty array alldistance so that it can be manipulated 
        distanceArray.append(distances)
    #for loop that appends distances under the cutoff to a new array called filteredAllDistance   
    for distanceSet in distanceArray:
        for element in distanceSet: 
            if (element < cutoff):
                filteredAllDistance.append(element)
    return filteredAllDistance

### Parameters below must be adjusted to match the desired ligand

#Getting the biopandas pdb to fetch the PDB file from the Protein Data Bank 
ppdb = ppdb().fetch_pdb('1qtm')
#creating DNA dataframe 
DA_df = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name']=='DA']
DG_df = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name']=='DG']
DC_df = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name']=='DC']
DT_df = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name']=='DT']
frames = [DA_df, DG_df, DC_df, DT_df]
DNA_df = pd.concat(frames)

set_of_distances1 = []
set_of_colors1 = []
#for loop to avoid adding filenames
for fileName in pdb_files1:
    Distances1 = distanceArrayMaker (fileName, 5)
    set_of_distances1.append (Distances1)
    set_of_colors1.append(user_color)
    
Array_of_Distances1 = np.concatenate((set_of_distances1), axis=0)

#needed for histogram 
import matplotlib.pyplot as plt 


#defining function for histogram
plt.hist (Array_of_Distances1, color=set_of_colors1)
plt.show ()
