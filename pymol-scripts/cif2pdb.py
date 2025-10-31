#!/usr/bin/env python3
"""
CIF to PDB Converter

Converts CIF (Crystallographic Information File) format files to PDB format using PyMOL.
Automatically detects crystal structures and applies symmetry expansion when available.
For predicted structures (e.g., AlphaFold), saves the original structure without symmetry operations.

Author: PyMOL Script Repository
License: Open Source
"""

import sys
import glob
import os
import pymol2

def main():
    if len(sys.argv) != 3 or sys.argv[1] in ("-h","--help"):
        print("""
CIF to PDB Converter

Usage:
    python cif2pdb.py <input_pattern> <output_dir>

Arguments:
    input_pattern    Glob pattern for CIF files (e.g., "*.cif" or "/path/*.cif")
    output_dir       Directory to save converted PDB files

Features:
    - Batch conversion of multiple CIF files
    - Automatic symmetry expansion for crystal structures
    - Handles predicted structures (AlphaFold) without symmetry warnings
    - Creates output directory if it doesn't exist

Examples:
    python cif2pdb.py "*.cif" ./pdb_output/
    python cif2pdb.py "~/data/*.cif" ~/converted/
""")
        return

    input_pattern = sys.argv[1]
    input_pattern = os.path.expanduser(input_pattern)
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    cif_files = glob.glob(input_pattern)

    print(f"Found {len(cif_files)} CIF files")

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        for cif_file in cif_files:
            name = os.path.splitext(os.path.basename(cif_file))[0]
            output_file = os.path.join(output_dir, name+".pdb")
            print(f"Converting {cif_file} -> {output_file}")
            cmd.load(cif_file, name)
            
            # Check if symmetry info exists
            if cmd.get_symmetry(name):
                cmd.symexp("sym", name, name, 10)
                cmd.save(output_file, "sym*")
            else:
                cmd.save(output_file, name)
            
            cmd.reinitialize()

if __name__ == "__main__":
    main()
