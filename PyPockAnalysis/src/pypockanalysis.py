import MDAnalysis as md
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import os, sys

def extract_reference(tpr,trr):
    doc=''' This function just extracts the first frame of a trr trajectory in 
    order to align the rest of the trajectory to it'''
    myMD=md.Universe(tpr,trr)
    atomGroup=myMD.select_atoms('all')
    #extract the reference frame
    nameOut='ref.gro'
    md.Writer(nameOut).write(atomGroup)
    return nameOut

def align(tpr,trr):
    nameOut=extract_reference(tpr,trr)
    myMD=md.Universe(tpr,trr)
    ref=md.Universe(tpr,nameOut)
    rms_fit_trj(myMD,ref,select='name CA')
        

#########################################
# MAIN FUNCITON                         #
#########################################

tpr=sys.argv[1]
trr=sys.argv[2]
myMD=align(tpr,trr)

    