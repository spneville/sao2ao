"""
Writes an Turbomole-formatted MOs file in terms of the AO basis
"""
import sys
import os
import numpy as np
import aobas

def write_mos(atoms,natm,saos):
    """
    Simultaneous reading of the original mos file and writing of
    the new mos file
    """
    
    print("HERE")
