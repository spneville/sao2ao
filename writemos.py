"""
Writes an Turbomole-formatted MOs file in terms of the AO basis
"""
import sys
import os
import numpy as np
import math
import aobas

def write_mos(mofile, atoms, natm, saos):
    """
    Simultaneous reading of the original mos file and writing of
    the new mos file
    """

    #
    # Construct the AO-to-SAO transformation matrices
    # for all irreps
    #
    # Maximum number of SAOs across all irreps
    ns=[]
    for i in range(aobas.nirrep):
        ns.append(saos[i].nsao)
    maxnsao = max(ns)
    # Construct the transformation matrices
    trafo = np.zeros((aobas.nao_tot,maxnsao,aobas.nirrep))
    for i in range(aobas.nirrep):
        for j in range(saos[i].nsao):
            for k in range(aobas.nao_tot):
                trafo[k][j][i] = saos[i].ao2sao[j][k]

    #
    # Open the output file
    #
    outfile = open('mos.ao', 'w')

    #
    # Read in the mos file
    #
    with open(mofile, 'r') as mosfile:
        mos = mosfile.readlines()

    #
    # Simulataneously read the mos file, transform the MOs to the
    # AO basis and write the mos.ao file
    #
    wid = 20
    il = -1
    while True:

        # Increment the line index
        il += 1

        # Split the line into keywords
        keywords=mos[il].split()

        # Are we at the start of another MO entry?
        if 'eigenvalue' in mos[il]:
            # Irrep index
            irrep=aobas.irrep_labels.index(keywords[1])

            # Number of SAOs
            nsaos=saos[irrep].nsao

            # Write the MO section header to the mos.ao file
            # with nsao set to the total number of AOs
            ilbl=mos[il].rfind('=')
            string = mos[il][:ilbl+1]
            outfile.write(string+str(aobas.nao_tot)+'\n')

            # Convert the MO from the SAO basis to the
            # AO basis
            nlines = math.ceil(nsaos/4)
            raw_orbs = []
            vec = np.zeros(nsaos)
            for k in range(nlines):
                il += 1
                mo_row = mos[il].lower().replace('d', 'e')
                n_cf = int(len(mo_row) // wid)
                for j in range(n_cf):
                    raw_orbs.append(float(mo_row[j*wid:(j+1)*wid]))
            for j in range(len(raw_orbs)):
                vec[j] = raw_orbs[j]
            trafo1  = np.zeros((aobas.nao_tot,nsaos))
            new_vec = np.zeros(aobas.nao_tot)
            for j in range(nsaos):
                for k in range(aobas.nao_tot):
                    trafo1[k][j] = trafo[k][j][irrep]
            new_vec = np.matmul(trafo1, vec)
            
            # Write the MO in terms of the AO basis to mos.ao
            nlines = math.ceil(aobas.nao_tot/4)
            nleft = aobas.nao_tot
            count = 0
            for k in range(nlines):

                # No. coefficients on this line
                n_mo = min(4,nleft)

                # Indices of the first and last coefficients
                # on this line
                iindx = count+1
                findx = count+n_mo

                # Write the coefficient string
                string=''
                for indx in range(iindx-1,findx):
                    # Hack to get around the need to use Fortran 20d.14 formatting,
                    # which doesn't seem to be possible using Python formatting...
                    if (new_vec[indx] < 0):
                        string1 = '{:{width}.{prec}e}'.format(-new_vec[indx], width=20, prec=14)
                        pieces = string1.split('.')
                        string1 = '-.'+pieces[0]+pieces[1]
                        pieces = string1.split('e')
                        string1 = pieces[0][0:len(pieces[0])-1]+'D'+'{0:+03d}'.format(int(pieces[1])+1)
                    else:
                        string1 = '{:{width}.{prec}e}'.format(new_vec[indx], width=20, prec=14)
                        pieces = string1.split('e')
                        string1 = pieces[0]+'D'+'{0:+03d}'.format(int(pieces[1]))
                    string = string + string1
                
                # Write the coefficient string to mos.ao
                outfile.write(string+'\n')
                
                # Update the counters
                nleft -= 4
                count += n_mo
                
        else:
            # Write the line to the mos.ao file
            outfile.write(mos[il])
            
        # Break out of the loop if the end of the file is reached
        if keywords[0] == '$end':
            break
        
    #
    # Close the output file
    #
    outfile.close()
