"""
Parsing of the infsao file
"""
import sys
import os
import numpy as np
import copy
import aobas

def read_infsao_file(saofile, atoms, natm, saos):
    """
    Parses the infsao file
    """
    
    #
    # Read in the infsao file
    #
    with open(saofile, 'r') as infsao:
        infsao_lines = infsao.readlines()

    #
    # (1) Initialise the SAO objects
    #
    n=-1
    for line in infsao_lines:
        if 'representation' in line:
            # Update the number of irreps
            aobas.nirrep+=1
            n+=1
            # Keywords on the line
            keywords=line.split()
            # Append the irrep label to the list
            aobas.irrep_labels.append(keywords[1])
            # Create an SAO class for the irrep
            saos.append(aobas.SAO())            
            # Set the irrep label and number of SAOs for the irrep
            saos[n].set_label(keywords[1])
            saos[n].set_nsao(int(keywords[5]))

    #
    # (2) Read in the AO-to-SAO transformations
    #
    # Initialise the irrep and line counters
    irrep = -1
    i     = -1

    # Loop over lines in infsao
    while True:

        # Increment the line counter
        i += 1
        
        # Are we at the next irrep?
        if 'SAO-index ' in infsao_lines[i]:
            # Increment te irrep counter
            irrep += 1

            # Initialise the SAO counter
            nsao = -1

            # Loop over SAOs
            while True:

                # Increment the line counter
                i += 1

                # Split the line into keywords
                keywords=infsao_lines[i].split()
                
                # Are we at the next SAO?
                if len(keywords) == 7:
                    # Increment the SAO counter
                    nsao += 1

                    # Initialise the current column in the
                    # AO-to-SAO transformation matrix
                    column = np.zeros(aobas.nao_tot)
                    
                    # Read in the SAO coefficients
                    atnum = int(keywords[2])-1
                    aolbl = str(keywords[4])+' '+str(keywords[5])
                    aonum = atoms[atnum].offset+atoms[atnum].aolbls.index(aolbl)
                    column[aonum]=float(keywords[6])
                    while True:
                        i += 1
                        keywords1=infsao_lines[i].split()
                        if len(keywords1) == 5:
                            atnum = int(keywords1[0])-1
                            aolbl = str(keywords1[2])+' '+str(keywords1[3])
                            aonum = atoms[atnum].offset+atoms[atnum].aolbls.index(aolbl)
                            column[aonum]=float(keywords1[4])
                        else:
                            # Break out of the SAO coeff. loop
                            i -= 1
                            break

                    # Add the current current column in the
                    # AO-to-SAO transformation matrix
                    saos[irrep].add_column(column)
                                        
                # Break out of the SAO loop if we have passed the
                # last SAO for the current irrep
                if nsao+1 == saos[irrep].nsao:
                    break
                
        # Break out of the main loop of we have parsed the entry for
        # the last irrep
        if irrep+1 == aobas.nirrep:
            break
            
