"""
Parsing of the basis file
"""
import sys
import os
import numpy as np
import aobas

def get_natm():
    """
    Parses the coord file to determine the number of atoms
    """

    # Read in the coord file
    with open('coord', 'r') as coord:
        coord_lines = coord.readlines()

    # Determine the no. atoms
    natm=0
    for line in coord_lines[1:]:
        if '$' not in line:
            natm+=1
        
    return natm


def get_atom_labels(atoms,natm):
    """
    Parses the coord file to determine the atom labels
    """

    # Read in the coord file
    with open('coord', 'r') as coord:
        coord_lines = coord.readlines()

    # Determine the atom labels
    for i in range(0,natm):
        line=coord_lines[i+1].split()
        lbl=line[len(line)-1]
        atoms[i].set_label(lbl)
    
def read_basis_file(atoms,natm):
    """
    Parses the basis file, filling in the atom objects as we go
    """

    # Read in the basis file
    with open('basis', 'r') as basis:
        basis_lines = basis.readlines()

    # Get to the start of the $basis section
    i = 0
    while basis_lines[i].split()[0] != '$basis':
        i += 1
    i += 1
    
    # Parse the basis information
    sec_start = False
    while True:

        line = basis_lines[i].split()

        if line[0] == '$end':
            # stepping on $end line breaks us out of parser
            break
        elif len(line) == 0:
            # ignore blank lines
            i += 1
        elif line[0][0] == '#':
            # ignore comment lines
            i += 1
        elif line[0] == '*':
            # if first string is star, either
            # beginning or ending basis section
            sec_start = not sec_start
            i += 1
        elif sec_start:
            # if starting section, first line
            # is atom and basis set line
            a_sym = line[0]
            b_set = line[1]
            a_list=[]
            for indx in range(natm):
                if a_sym == atoms[indx].label:
                    a_list.append(indx)
            i += 1
        else:
            # if we get this far, we're parsing the basis set!
            nprim, ang_sym = line
            ang_mom        = aobas.shell_label.index(ang_sym)
            nfunc          = aobas.shell_nfunc[ang_mom]
            for j in range(nfunc):
                lbl=aobas.ao_order[ang_mom][j]
                for atom in a_list:
                    atoms[atom].add_ao_label(lbl)
            for j in range(int(nprim)):
                i += 1
                exp, coef = basis_lines[i].split()

            i += 1

    # Set the total number of AOs and the AO index offsets
    nao=atoms[0].naos
    for n in range(1,natm):
        atoms[n].set_offset(nao)
        nao+=atoms[n].naos
    aobas.nao_tot=nao
