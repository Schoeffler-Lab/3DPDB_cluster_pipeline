## This script was written by Alyce Fields, Larissa Cortes Morales, and Allyn Schoeffler
# Loyola Univeristy New Orleans

from posixpath import split
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
import glob

#Below, specify where output files will be stored
io_path = "/enter/the/absolute/path/to/your/directory/here

#Below, specify where the threaded structures are located
user_path_to_pdbs = "/enter/the/absolute/path/to/your/threaded/structures/here"

path_to_pdbs = user_path_to_pdbs + "/*.pdb"

###End of user-specified parameterss ###

#Creates a list of PDBs in the thread directory
pdb_files = glob.glob(path_to_pdbs)
#Class built into biopython to select certain atoms of certain residues 
#Here there are all 21 of the most common amino acids named with their 3 letter code followed by Select (example:ThrSelect for threonine)
#These modules also filter to the most likely interacting atom; for example, the amino group of Lysine

class ArgSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "ARG": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "NH1" or atom.get_name()=="NH2": 
            return 1

class GlySelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "GLY": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class ProSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "PRO": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class LysSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "LYS": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "NZ": 
            return 1
class LysArg (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "LYS" or residue.get_resname()== "ARG":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "NH1" or atom.get_name()=="NH2" or atom.get_name()=="NZ":
            return 1 
        
class AsnSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "ASN":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OD1" or atom.get_name()=="ND2": 
            return 1 

class HisSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "HIS":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "ND1" or atom.get_name()=="NE2":
            return 1 

class AspSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "ASP":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OD1" or atom.get_name()=="OD2":
            return 1 

class GluSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "GLU":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OE1" or atom.get_name()=="OE2":
            return 1 

class SerSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "SER":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OG":
            return 1 

class ThrSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "THR":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OG1":
            return 1
class GlnSelect (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "GLN":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OE1" or atom.get_name()=="NE2":
            return 1 

class AlaSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "ALA": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class ValSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "VAL": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class IleSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "ILE": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class LeuSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "LEU": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class MetSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "MET": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class PheSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "PHD": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class TyrSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "TYR": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class TrpSelect (Select):
    def accept_residue(self, residue): 
        if residue.get_resname ()== "TRP": 
            return 1 
    def accept_atom(self, atom): 
        if atom.get_name()== "CA": 
            return 1

class GluAsp (Select): 
    def accept_residue(self,residue):
        if residue.get_resname()== "GLU" or residue.get_resname()== "ASP":
            return 1
    def accept_atom(self, atom):
        if atom.get_name()== "OE1" or atom.get_name()=="OE2" or atom.get_name()== "OD1" or atom.get_name()=="OD2":
            return 1 

#The loop below creates PDBs filteres for particular amino acid types and atoms. These are used in the next stage of the analysis

for pdb in pdb_files:
    pdb_path_list = pdb.split ('/')
    pdb_name_list = pdb_path_list [-1] .split ('.')
    pdb_name = pdb_name_list [0]
    parser= PDBParser (PERMISSIVE=1)
    structure = parser.get_structure (pdb_name, open (pdb))
    io = PDBIO()
    io.set_structure(structure)
    io.save( pdb + "Arg.pdb", ArgSelect())
    io.save( pdb + "Gly.pdb", GlySelect())
    io.save( pdb + "Pro.pdb", ProSelect())
    io.save( pdb + "Lys.pdb", LysSelect())
    io.save( pdb + "LysArg.pdb", LysArg())
    io.save( pdb + "Asn.pdb", AsnSelect())
    io.save( pdb + "His.pdb", HisSelect())
    io.save( pdb + "Asp.pdb", AspSelect())
    io.save( pdb + "Glu.pdb", GluSelect())
    io.save( pdb + "Ser.pdb", SerSelect())
    io.save( pdb + "Thr.pdb", ThrSelect())
    io.save( pdb + "Gln.pdb", GlnSelect())
    io.save( pdb + "Ala.pdb", AlaSelect())
    io.save( pdb + "Val.pdb", ValSelect())
    io.save( pdb + "Ile.pdb", IleSelect())
    io.save( pdb + "Leu.pdb", LeuSelect())
    io.save( pdb + "Met.pdb", MetSelect())
    io.save( pdb + "Phe.pdb", PheSelect())
    io.save( pdb + "Tyr.pdb", TyrSelect())
    io.save( pdb + "Trp.pdb", TrpSelect())
    io.save( pdb + "GluAsp.pdb", GluAsp())
