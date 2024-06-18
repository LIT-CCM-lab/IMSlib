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
import regex as re
import pandas as pd
import random
import numpy
import itertools
import string
import argparse 

import timeit
start = timeit.default_timer()


def counter(smiles):
    """This function counts number of atoms in smiles"""
    letter_list = string.ascii_uppercase + string.ascii_lowercase
    letter_count = 0
    for i in smiles:
        if i in letter_list:
            letter_count += 1
    return letter_count


def convert_to_mols(rgs):
    """
    This function converts smiles to mols
    """
    mols = []
    for smi in rgs:
        mol = Chem.MolFromSmiles(smi)
        mols.append(mol)
    return mols


def combine_groups(mols, n_samples):
    """This function gives specified number of combinations of R groups"""
    return [[random.choice(mols[key]) for key in mols.keys()] for i in range(n_samples)]



def is_r_type(smi, r_type):
    """This function tells the type of R groups --> aromaric, bulky, hydrophilic, etc."""
    if smi in r_type:
        return True
    else:
        return False


def count_functional_groups(order, group_type):
    """This function counts the number of specified interactions in one combination"""
    i = 0
    for mol in order:
        if is_r_type(mol, group_type)==True:
            i+=1
    return i


def neutralize_atoms(mol):
    """This function neutralizes charged molecules by adding and/or removing hydrogens"""
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def find_substitution_points(smile):
    """This function finds the number of substitution points in core"""
    numbers = re.findall(r'\[\*:(\d+)\]', smile)
    numbers = [int(num) for num in numbers]
    
    if numbers:
        highest_number = max(numbers)
    else:
        highest_number = None
    return highest_number

def find_max_subs_points_all(smiles):
    """This function finds the highest number of substitution points across the list of cores"""
    highest_numbers = []
    
    for smile in smiles:
        highest_number = find_substitution_points(smile)
        if highest_number is not None:
            highest_numbers.append(highest_number)
            
    if highest_numbers:
        highest_overall = max(highest_numbers)
    else:
        highest_overall = None
    
    return highest_overall


def main(rgs, cores, num_combinations):


    # # Import R-groups
    df = pd.read_csv(rgs, sep = " ")


    canonical = []
    for smi in df['R1']:
        can = Chem.CanonSmiles(smi, useChiral=0)
        canonical.append(can)
    df['R1'] = canonical


    # # Import core scaffolds
    cores = [mol for mol in Chem.SmilesMolSupplier(cores, titleLine=False)]


    #Get number of substitution points on the biggest core 
    core_smiles = []
    for core in cores:
        smile = Chem.MolToSmiles(core)
        core_smiles.append(smile)
        
    subs_points = find_max_subs_points_all(core_smiles)


    #Multiply R groups 
    for i in range(1, int(subs_points)+1): #number of possible supstitutions 
        df['R'f'{i}'] = df["R1"]
        df['R'f'{i}'] = df['R'f'{i}'].str.replace(':[1-int(subs_points)]', ':'f'{i}', regex = True)


    rgs = df.to_dict('list')

    # Define types of R groups
    hydrogens = []
    aromatic = []
    halogens = []
    bulky = []
    heterocyclic = []
    hydrophilic = [] # HBAs and HBDs
    alkyl = [] #alkyl, alkenyl, alkynyl, cicloalkyl
    sulphur = []
    small = []
    for key in rgs.keys():
        #print(rgs[key])
        ar = [smi for smi in rgs[key] if "c" in smi 
              and counter(smi) < 12
              and "s" not in smi
              and "n" not in smi
              and "o" not in smi]
        aromatic.append(ar) # Nested list - must flatten it 
        aromatic_flat = list(numpy.concatenate(aromatic).flat) 
        
        pattern = Chem.MolFromSmarts('*@[n,N,o,O,s,S]')
        for smiles in rgs[key]:
            m = Chem.MolFromSmiles(smiles)
            if m.HasSubstructMatch(pattern):
                het = Chem.MolToSmiles(m)
                if counter(het) < 12:
                    heterocyclic.append(het)
        heterocyclic_flat = heterocyclic

        h = [smi for smi in rgs[key] if "[H]" in smi]
        hydrogens.append(h)
        hydrogens_flat = list(numpy.concatenate(hydrogens).flat)

        cl = [smi for smi in rgs[key] if "Cl[" in smi
              or "]Cl" in smi
              or "F[" in smi
              or "]F" in smi
              or "]Br" in smi
              or "Br[" in smi
              or "I[" in smi
              or "]I" in smi]
        halogens.append(cl)
        halogens_flat = list(numpy.concatenate(halogens).flat)

        big = [smi for smi in rgs[key] if counter(smi) >= 12]
        bulky.append(big)
        bulky_flat = list(numpy.concatenate(bulky).flat)

        o = [smi for smi in rgs[key] if "O" in smi]
        n = [smi for smi in rgs[key] if "N" in smi]
        hydrophilic.extend([o,n])
        hydrophilic_flat = list(numpy.concatenate(hydrophilic).flat)
        
        s = [smi for smi in rgs[key] if "S" in smi or "s" in smi]
        sulphur.append(s)
        sulphur_flat = list(numpy.concatenate(sulphur).flat)

        r = [smi for smis in rgs.values() for smi in smis if not re.search('[a-zA-BD-Z]', smi)]
        alkyl.append(r)
        alkyl_flat = list(numpy.concatenate(alkyl).flat)
        
        sm = [smi for smi in rgs[key] if counter(smi) <= 4]
        small.append(sm)
        small_flat = list(numpy.concatenate(small).flat)


    # Set the number of combinations
    combo = combine_groups(rgs, num_combinations)


    combo_new = []
    if subs_points <= 4:
        for combination in combo:
            if (count_functional_groups(combination, aromatic_flat) <= 2
                and count_functional_groups(combination, hydrogens_flat) <= 1
                and count_functional_groups(combination, bulky_flat) <= 1
                and count_functional_groups(combination, halogens_flat) <= 2
                and count_functional_groups(combination, heterocyclic_flat) <= 2
                and count_functional_groups(combination, hydrophilic_flat) <= 2
                and count_functional_groups(combination, sulphur_flat) <= 2
                and count_functional_groups(combination, alkyl_flat) <= 3):
                combo_new.append(combination)
                
    if subs_points > 4:
        for combination in combo:
            if (count_functional_groups(combination, aromatic_flat) <= 2
                and count_functional_groups(combination, hydrogens_flat) >= 2
                and count_functional_groups(combination, hydrogens_flat) <= 4
                and count_functional_groups(combination, bulky_flat) <= 2
                and count_functional_groups(combination, halogens_flat) <= 2
                and count_functional_groups(combination, heterocyclic_flat) <= 2
                and count_functional_groups(combination, hydrophilic_flat) <= 3
                and count_functional_groups(combination, sulphur_flat) <= 2
                and count_functional_groups(combination, alkyl_flat) <= 3):
                combo_new.append(combination)


    # Remove duplicates
    combo_new.sort()
    final_combo = list(l for l, _ in itertools.groupby(combo_new))

    # Convert smiles to molecules
    r_combos = []
    for i in range(len(final_combo)):
        r_combos.append(convert_to_mols(final_combo[i]))

    combo_tpl = [tuple(x) for x in r_combos]

    # # Make new molecules

    products = []
    for core in cores:
        for tpl in combo_tpl:
            tm = Chem.RWMol(core)
            for r in tpl:
                tm.InsertMol(r)
            prod = Chem.molzip(tm)
            if prod is not None:
                # Convert to smiles 
                neutral = neutralize_atoms(prod)
                mols = Chem.MolToSmiles(neutral)
                for fragment in mols.split("."):
                    if ":" not in fragment:
                        canonical_fragment = Chem.CanonSmiles(fragment, useChiral=0)
                        products.append(canonical_fragment)

    # Remove duplicates
    smis_set = list(set(products))
    print(len(smis_set))


    # # Save as .smi
    with open("./products.smi", "a") as f:
         for line in smis_set:
                f.write(f"{line}\n")


    stop = timeit.default_timer()
    print('Time: ', stop - start) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Generate molecules from cores and R-groups', epilog = "Good Luck !")
    parser.add_argument('-r', '--rgroups', help = "The first input file is the r group collection (smi)", required=True)
    parser.add_argument('-c', '--cores', help = "The second input file is the list of cores (smi)", required=True)
    parser.add_argument('-n', '--num_combinations', type=int, help = "The number of R-groups combinations", required=True)
    args = parser.parse_args()

    sys.exit(main(args.rgroups, args.cores, args.num_combinations))
