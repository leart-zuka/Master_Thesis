# -*- coding: utf-8 -*-
"""
Created on Fri Sep 05 13:45:19 2014

@author: Manuel

Function to automatically generate the state space and all operators for an atom-cavity system

Usage:

from CavityOperators import operators
a,S,astate = operators(photondim,atomdim)

where photondim and atomdim are lists of integers with each integer specifiying an n-dimensional atom/mode
e.g.

phondim = [3,3]
atomdim = [5]

specifies a five-state atom and two modes up to two photons
"""



from qutip import *
import numpy

def field_operators(photondim,atomdim):
    modes = numpy.size(photondim)
    ops = []
    for i in range(0,numpy.size(atomdim)):
        ops.append(qeye(atomdim[i]))
    
    idatoms = tensor(ops)
    
    a = []
    for i in range(0,modes):
        ops = []
        ops.append(idatoms)
        for j in range(0,modes):
            if j==i:
                ops.append(destroy(photondim[j]))
            else:
                ops.append(qeye(photondim[j]))
        a.append(tensor(ops))

    ops = []

    
    return a
    
def atom_operators(photondim, atomdim):
    ops = []    
    for i in range(0,numpy.size(photondim)):
        ops.append(qeye(photondim[i]))
    
    idfields = tensor(ops)
    
    astates = []
    for i in range(0,numpy.size(atomdim)):
        tempstate = []
        for j in range(0,atomdim[i]):
            tempstate.append(basis(atomdim[i],j))
        astates.append(tempstate)
    
    S = []   
    for k in range(0,numpy.size(atomdim)):
        tempS = [[None] * atomdim[k] for i in range(atomdim[k])]
        for i in range(0,atomdim[k]):
            for j in range(0,atomdim[k]):
                ops = []
                for l in range(0,numpy.size(atomdim)):
                    if l==k:
                        ops.append(astates[k][i] * astates[k][j].dag())
                    else:
                        ops.append(qeye(atomdim[l]))
                ops.append(idfields)
                tempS[i][j] = tensor(ops)
        S.append(tempS)
    return astates,S
    
def operators(photondim, atomdim):
    """The function builds the operators for a cavity-atom system with (potentially) multiple atoms and multiple modes
    photondim must be a list of integers, which describe at which occupation state to truncate the Hilbert space for each mode.
    atomdim must be a list of integers, describing how many levels each atom has.
    Example: operators([4,4],[3,3]) creates a system with two modes, which can be occupied up to n=3 and two three-level atoms.
    a returns a list of destruction field operators (e.g. a[0] destroys a photon in mode 0)
    S returns a list of two-dimensional lists of operators for each atom. For example S[0][1][2] couples level 2 to level 1 for the 0th atom.
    atomstates is a list of lists, describing the basis states for each atom. atomstates[0][1] is the state 1 of atom 0.
    """
    a = field_operators(photondim,atomdim)
    atomstates, S = atom_operators(photondim,atomdim)
    return a,S,atomstates
    
        
        
    