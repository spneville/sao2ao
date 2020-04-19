"""
Structures and functions to store the AO basis.
Note that we will only support up to f-functions for the moment.
"""
import os
import numpy as np

# Shell labels
shell_label  = ['s', 'p', 'd', 'f']

# Number of AOs per shell for each angular momentum quantum number
shell_nfunc  = [1, 3, 5, 7]

# Turbomole ordering of the AOs within each shell
ao_order     = [['s',  '',    '',    '',    '',    '',    ''   ],
                ['px', 'py',  'pz',  '',    '',    '',    ''   ],
                ['d0', 'd1a', 'd1b', 'd2a', 'd2b', '',    ''   ],
                ['f0', 'f1a', 'f1b', 'f2a', 'f2b', 'f3a', 'f3b']]

# Total number of AOs
nao_tot      = 0

# Irreps
nirrep       = 0
irrep_labels = []

class Atom:
    """
    Class to hold all the information about the AOs associated
    with a single atom
    """

    def __init__(self):
        self.offset = 0
        self.nbas   = 0
        self.label  = ''
        self.aolbls = []
        self.naos   = 0
        
    def set_offset(self, i):
        """
        Set the offset for the AO indices for the atom
        """
        self.offset = i

    def set_label(self, lbl):
        """
        Set the atom label 
        """
        self.label = lbl

    def add_ao_label(self, lbl):
        """
        Add the label of an AO for the atom
        """
        # Number of times the AO label has already
        # occured for the atoms
        k=1+sum(lbl in s for s in self.aolbls)
        # Add the AO label to the list
        self.aolbls.append(str(k)+' '+lbl)
        # Update the number of AOs
        self.naos+=1

        
class SAO:
    """
    Class to hold the information about the SAOs for a single irrep
    """

    def __init__(self):
        self.nsao   = 0
        self.label  = ''
        self.ao2sao = []
                
    def set_label(self, lbl):
        """
        Set the irrep label
        """
        self.label = lbl

    def set_nsao(self, n):
        """
        Set the number of SAOs
        """
        self.nsao = n

    def add_column(self, column):
        """
        Appends a column to the AO-to-SAO transformation matrix
        """
        self.ao2sao.append(column)
