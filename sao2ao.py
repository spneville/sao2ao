#!/usr/bin/env python
"""
A simple program to convert a Turbomole mos file
from the SAO basis to the AO basis
"""
import sys
import os.path
from os import path
import errno
import os
import aobas
import readbas
import readinfsao
import writemos


def get_arg(kwd, default):
    """Gets a keyword argument from sys.argv if present, otherwise sets
    the value to a given default."""
    args = sys.argv
    if kwd in args:
        return args[args.index(kwd)+1]
    else:
        return default


def read_args():
    """
    Reads the names of the mo and sao files if they have been
    provided by the user. Else, returns the default names
    """
    
    # Default filenames
    mofile='mos'
    saofile='infsao'

    # User-provided filenames
    mofile  = get_arg('-mos', mofile)
    saofile = get_arg('-sao', saofile)
    
    return mofile, saofile


def check_files_exist(mofile, saofile):
    """
    Checks that the mos and infsao files exist
    """
    filelist=[mofile, saofile, 'basis', 'coord']
    for filename in filelist:
        if not path.exists(filename):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),filename)


if __name__ == '__main__':
    # Read the command line arguments
    mofile, saofile = read_args()

    # Check that all required files exist
    check_files_exist(mofile, saofile)

    # Number of atoms
    natm = readbas.get_natm()

    # Create the atom classes
    atoms=[]
    for i in range(natm):
        atoms.append(aobas.Atom())

    # Get the atom labels
    readbas.get_atom_labels(atoms, natm)
    
    # Read the AO basis information
    readbas.read_basis_file(atoms, natm)

    # Read the SAO contraction coefficients
    saos=[]
    readinfsao.read_infsao_file(saofile, atoms, natm, saos)
    
    # Write the MOs in terms of the AOs to file
    writemos.write_mos(mofile, atoms, natm, saos)
