# !/usr/bin/env python
# File name: itp_to_top.py
# Author: Matt Robinson, Jorgensen Lab @ Yale
# Email: matthew.robinson@yale.edu
# Date created: 01/31/2017
# Python Version: 3.6

"""
Description:
Converts .itp files from LigParGen into .top Amber files using ParmEd

Usage:
python itp_to_top.py -i initial_file.itp -f final_file.top 

Requirements:
Preferably Anaconda python with following modules:
    ParmEd
    
"""
#from rdkit.Chem import rdMolAlign
import parmed as pmd
import argparse
import os

def main():

    parser = argparse.ArgumentParser(
        prog='make_single_topology.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        Converts .itp files from LigParGen into .top Amber files using ParmEd
         
        @author: Matt Robinson, Jorgensen Lab @ Yale
        @email:  matthew.robinson@yale.edu

        Usage: 
        python itp_to_top.py -i initial_file.itp 

        REQUIREMENTS:
        Preferably Anaconda python with following modules
            ParmEd
        """
    )
    parser.add_argument(
        "-i", "--itp", help="initial .itp filename",required=True)

    args = parser.parse_args()

    make_top_file(args.itp)

    run_parmed(args.itp)


def make_top_file(itp_filename):
    '''
    Converts .itp file to .top file
    '''
    with open(itp_filename) as itp_file:
        itp_data = itp_file.readlines()

    ligand_resname = itp_filename[:-4]
    top_filename = itp_filename[:-4]+'.top'

    # find relevant indicies to parse only important part of z-matrices
    line_index = 0
    for line in itp_data:
        if '[ atomtypes ]' in line:
            atomtypes_index = line_index
            break
        line_index = line_index + 1

    # include the necessary [ defaults ] info
    defaults_list = \
    [
        '[ defaults ]\n',
        '; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n',
        '1               3               yes             0.5     0.5\n',
        '\n'
    ]

    itp_data_with_atomtypes = itp_data[:atomtypes_index] + defaults_list + itp_data[atomtypes_index:]

    # NOTE: the name in 'moleculetype' in .itp file must match name in ' molecules' section
    ending_list = \
    [
        '\n',
        '[ system ]\n',
        '; Name\n',
        'OPLSAA UNK\n',
        '\n',
        '[ molecules ]\n',
        '; Compound        #mols\n',
        'UNK' + '     1\n'
        #ligand_resname + '     1\n'
    ]

    top_data = itp_data_with_atomtypes + ending_list

    with open(top_filename,'w') as top_file:
        for line in top_data:
            top_file.write(line)

def run_parmed(itp_filename):

    top_filename = itp_filename[:-4]+'.top'
    # pdb_filename = itp_filename[:-4]+'.pdb'
    gro_filename = itp_filename[:-4]+'.gro'
    parm_filename = itp_filename[:-4]+'.parm7'
    rst_filename = itp_filename[:-4]+'.rst7'

    # gmx_top = pmd.load_file(top_filename,xyz=pdb_filename)
    gmx_top = pmd.load_file(top_filename,xyz=gro_filename)
    gmx_top.save(parm_filename, format='amber',overwrite=True)
    gmx_top.save(rst_filename, format='rst7',overwrite=True)


if __name__ == "__main__":
    main()
