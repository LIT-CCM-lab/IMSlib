---
title: Scaffold Decorator
author: Teodora Djikic-Stojsic
date: 2024-06-14
---

# Scaffold Decorator

Scaffold decorator consists of two python scripts: one for the collection of R-groups and other one for generation of new molecules. 

## Requirements
* Python 3
* RDKit version 2022.03.5 (or higher)
* Regex
* Pandas
* NumPy

## Folder Structure
scaffold_decorator

├── collect_r_groups

    ├── collect_r_groups.py
    ├── example_molecules.smi
    └── r_groups_output.smi
   
└── generate_new_molecules

    ├── make_new_molecules.py
    ├── r_groups_unique.smi
    └── example_scaffolds.smi
    
└── other_files
          
      ├── bioinfo_db.smi
      ├── chembl_clinical.smi

## Step 1: R-Group Collection
* Run the script ‘collect_r_groups.py’, providing molecules and scaffolds as an input to generate the ‘r_groups.smi’ file
```sh
python collect_r_groups.py -m example_molecules.smi -c example_core_scaffolds.smi
```
* Input Files:
    - example_molecules.smi: Contains existing molecules in SMILES format.
    - example_core_scaffolds.smi: Contains core scaffolds in SMILES format.

* Output:
    - r_groups.smi: Contains the R-groups identified from the decomposition

## Step 2: Generation of New Molecules
* Run the script ‘make_new_molecules.py’, providing the scaffolds with the substitution points, a collection of R groups, and the desired number of R group combinations to generate a file ‘products.smi’ that contains novel molecules
```sh
python make_new_molecules.py -r r_groups_unique.smi -c example_scaffolds.smi -n 1000
```
* Input Files:
    - r_groups.smi: Contains the R-groups for decoration.
    - example_scaffolds.smi: Contains scaffolds to be decorated.
    - Make sure that scaffolds 2-4 substitution points and > 4 substitution points are separated
    - Number of substitution points and number of R-group combinations should be provided.
* Output:
     - new_molecules.smi: Contains the newly generated molecules
## Note
* Use R-groups with multiplicities to ensure the good distribution of different R-group types
* Make sure to filter the molecules after generating them. The filtering rules are provided in Supplementary Materials 5.

## Summary
* Ensure that your environment has Python 3, Regex, Pandas, NumPy, and RDKit 2022.03.5 (or higher).
* Follow the directory structure and place input files accordingly.
* Run the provided scripts to collect R-groups and generate new molecules.

