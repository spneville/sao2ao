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

def check_files_exist():
    """
    Checks that the mos and infsao files exist
    """
    filelist=['mos', 'infsao', 'basis', 'coord']
    for filename in filelist:
        if not path.exists(filename):
            raise FileNotFoundError(errno.ENOENT,
                                    os.strerror(errno.ENOENT),filename)


if __name__ == '__main__':
    # Check that all required files exist
    check_files_exist()

    # Number of atoms
    natm = readbas.get_natm()

    # Create the atom classes
    atoms=[]
    for i in range(natm):
        atoms.append(aobas.Atom())

    # Get the atom labels
    readbas.get_atom_labels(atoms,natm)
    
    # Read the AO basis information
    readbas.read_basis_file(atoms,natm)

    # Read the SAO contraction coefficients
    saos=[]
    readinfsao.read_infsao_file(atoms,natm,saos)
    
