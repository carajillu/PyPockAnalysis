'''
MODULE  COORDINATES.
1. Input files:
   - Topology (tpr): GROMACS tpr file.
   - Trajectory (trr): GROMACS trr or xtc file.
   - Reference (nameRef): pdb or gro file.

2. Functions:
   2.1. align(tpr,trr,ref): 
   aligns trr to ref using backbone. Produces an aligned trajectory file and 
   returns its name.
                         
   2.2 split(tpr,nameAli,resToIgnore,ext): 
   splits a trajectory (usually aligned, but NOT mandatory) in individual frames
   with the desired extension. It ignores HOH, WAT, SOL, NA nad CL by default
   
   '''

import os
import MDAnalysis as md
from MDAnalysis.analysis.align import *

def align(tpr,trr,nameRef):
    '''
    Aligns trr to ref using backbone. Produces an aligned
    trajectory file and returns its name.
    '''
    myMD=md.Universe(tpr,trr)
    ref=md.Universe(tpr,nameRef)
    nameAli=trr.split('.')[0]+'_ali.trr'
    rms_fit_trj(myMD,ref,select='protein and name CA',filename=nameAli)
    return nameAli

def split(tpr,nameAli,resToIgnore,ext):
    
    myMD=md.Universe(tpr,nameAli)
    selection='not resname '+' and not resname '.join(resToIgnore)
    atomGroup=myMD.select_atoms(selection)
    dirName='split'+ext.upper()
    mkdir='mkdir '+dirName
    os.system(mkdir)
    listName=dirName+'.txt'
    pdbList=[]
    fileout=open(listName,'w')
    i=0
    for ts in myMD.trajectory:
        PDBOut=dirName+'/'+nameAli.split('_ali')[0]+'_'+'0'*(5-len(str(i)))+str(i)+'.'+ext
        fileout.write(PDBOut+'\n')
        pdbList.append(PDBOut)
        md.Writer(PDBOut).write(atomGroup)
        print 'wrote file '+PDBOut
        i=i+1
    fileout.close()
    return pdbList