#!/usr/bin/env python3
"""
pdbfixer_cleaning_local.py

Given a folder containing a single PDB file, this script will:
 1. Change into the directory ./<folder_name>/.
 2. Find the single PDB file within it.
 3. Rename that PDB file to <folder_name>.pdb.
 4. Strip out all heterogens except water.
 5. Build any missing residues/atoms.
 6. Add all hydrogens at pH 7.0.
 7. Write out <folder_name>_cleaned.pdb in that folder.

Usage:
    # Assuming you have a folder named 'my_protein' with one .pdb file inside
    python3 pdbfixer_cleaning_local.py my_protein

Requirements:
    conda install -c conda-forge pdbfixer openmm
"""

import os
import sys
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def preprocess(input_pdb, output_pdb, target_pH=7.0):
    """Applies PDBFixer cleaning operations."""
    print(f"[Preprocess] Loading {input_pdb}")
    fixer = PDBFixer(filename=input_pdb)
    print("[Preprocess] Stripping heterogens (keeping waters)…")
    fixer.removeHeterogens(keepWater=True)
    print("[Preprocess] Finding/building missing residues & atoms…")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    print(f"[Preprocess] Adding hydrogens at pH {target_pH}…")
    fixer.addMissingHydrogens(pH=target_pH)
    print(f"[Preprocess] Writing cleaned PDB → {output_pdb}")
    with open(output_pdb, 'w') as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out)

def main():
    parser = argparse.ArgumentParser(description="Find, rename, & preprocess a single PDB file in a given folder.")
    parser.add_argument('folder_name', help="Name of the directory containing a single PDB file.")
    args = parser.parse_args()

    folder_name = args.folder_name.lower()

    # 1) Check if folder exists and change into it
    if not os.path.isdir(folder_name):
        print(f"❌ Error: Directory not found at './{folder_name}/'")
        sys.exit(1)
    os.chdir(folder_name)
    print(f"[Setup] Changed directory to ./{folder_name}/")

    # 2) Find the single PDB file
    pdb_files = [f for f in os.listdir('.') if f.lower().endswith('.pdb')]

    if len(pdb_files) == 0:
        print(f"❌ Error: No PDB files (.pdb) found in '{folder_name}'.")
        sys.exit(1)
    elif len(pdb_files) > 1:
        print(f"❌ Error: More than one PDB file found in '{folder_name}'. Please ensure there is only one.")
        print("Found:", pdb_files)
        sys.exit(1)
    
    original_pdb_filename = pdb_files[0]
    target_pdb_filename = f"{folder_name}.pdb"
    
    # 3) Rename the PDB file to match the folder name
    if original_pdb_filename != target_pdb_filename:
        print(f"[Setup] Renaming '{original_pdb_filename}' → '{target_pdb_filename}'")
        os.rename(original_pdb_filename, target_pdb_filename)
    else:
        print(f"[Setup] PDB file is already correctly named: '{target_pdb_filename}'")

    cleaned_pdb = f"{folder_name}_cleaned.pdb"

    # 4-7) Preprocess the renamed PDB file
    preprocess(target_pdb_filename, cleaned_pdb)

    # Change back to the original directory to be clean
    os.chdir('..')
    
    print("\n✅ All done!")
    print(f"→ Check directory ./{folder_name}/ for {target_pdb_filename} and {cleaned_pdb}")

if __name__ == "__main__":
    main()