#!/usr/bin/env python
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdqueries
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import molzip
from rdkit.Chem import RWMol
from rdkit.Chem import rdDepictor
from rdkit import Geometry
rdDepictor.SetPreferCoordGen(True)
from rdkit.Chem import rdmolops
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import itertools
import rdkit
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import pandas as pd
import random
import argparse 

def main(molecules, cores):

	# Import molecules
	mols = [mol for mol in Chem.SmilesMolSupplier(molecules ,delimiter='\t', titleLine=True)]


	# Import core scaffolds
	scaffolds = [mol for mol in Chem.SmilesMolSupplier(cores ,delimiter=' ', titleLine=False)]

	# R_group decomposition
	groups,unmatched = rdRGroupDecomposition.RGroupDecompose(scaffolds, mols, asSmiles=False, asRows=False)
	#len(unmatched)

	keis = len(groups.keys()) # First key is a core the rest are R groups

	# Get R groups but keep duplicates
	# Remove R groups that have more than 1 dummy atom
	r = {}
	for i in range(1, keis):
		query_atom = rdqueries.AtomNumEqualsQueryAtom(0)
		r["r{0}".format(i)] = [x for x in groups['R'f'{i}'] if len(x.GetAtomsMatchingQuery(query_atom))==1]

	#Create a list with all R groups
	rgs_list = list(r.values()) #list of a list
	rgroups_list = [item for sublist in rgs_list for item in sublist] #flat list

	rgroups_smis = []
	for mol in rgroups_list:
		smi = Chem.MolToSmiles(mol)
		rgroups_smis.append(smi)


	# Make a dataframe with column R1 using flat list and rename all dummy atom to position 1
	rgs_df = pd.DataFrame(rgroups_smis, columns = ["R1"])
	rgs_df['R1'] = rgs_df['R1'].str.replace(":[1-9]+", ":1", regex = True)
	#Remove R grups that contains Silicium
	rgs_df = rgs_df[rgs_df['R1'].str.contains("Si")==False]
	#Remove R grups that begin with double bonds
	rgs_df = rgs_df[rgs_df["R1"].str.contains("=\[\*")==False]


	# Save df to check the groups
	rgs_df.to_csv('./r_groups.smi', sep = '\t', index = False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = 'Collect R-groups from molecules', epilog = "Good Luck !")
	parser.add_argument('-m', '--molecules', help = "The first input file is the r group collection (smi)", required=True)
	parser.add_argument('-c', '--cores', help = "The second input file is the list of cores (smi)", required=True)
	
	args = parser.parse_args()
	sys.exit(main(args.molecules, args.cores))



