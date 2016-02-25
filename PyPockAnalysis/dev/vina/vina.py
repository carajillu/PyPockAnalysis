'''
MODULE VINA
Runs docking calculations using vina in all the MD snapshots. It can use a
predefined binding site specified in the input file (vina_ensemble) or calculate
it from the coordinates of a pocket calculated by fpocket (vina_fpocket).
1. Input files:
   - pdbList (txt): if split is not selected, the names of the aligned pdb
                    files will be read from this file.
2. Functions:
   2.1. gen_pdbqt: uses AutoDock Tools to generate pdbqt files from all the 
                  MD snapshots. If the system includes cofactors, it has to be
                  tricked to include its residue name. Otherwise it will be
                  ignored.
   2.2. fpocket_box: This funtion makes a vina suitable input file from a 
                     binding site computed by fpocket. NOT recommended for 
                     normal docking, as fpocket does not always find the 
                     correct binding site.
   2.3. vina_ensemble: This funtion docks the ligands in all the MD snapshots but
                      in a binding site defined a priori. 
   2.4. get_vina_scores: extracts the scores of the best pose and stores them
                      in a dictionary. The dictionary has a key for each ligand
                      and a list with all the scores of this ligand.
   2.5. get_vina_rmsd: the same as get_vina_scores but with the ligand RMSD.
   2.6. vina_plot: extracts the info from get_vina_scores and get_vina_rmsd and 
                   outputs it to a file that we can plot.
    
'''
import os, math, numpy
import MDAnalysis as md
from MDAnalysis.analysis.align import *

def gen_pdbqt(pdbList):
    '''
    I am NOT doing this with MDAnalysis because Vina doesn't take the pdbqt
    files produced by it. To make it work with cofactors, you have to trick
    the file MoleculePreparation.py (in AutoDockTools dir) to accept the
    resname of the cofactor.
    I think that it understands the cofactor topology but I am not sure.
    Use this at your own risk.
    '''
    ligprep='/usr/local/MGLTools-1.5.6/bin/pythonsh '+\
            '/usr/local/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/'+\
            'Utilities24/prepare_receptor4.py'
    os.system('mkdir splitPDBQT')
    pdbqtList=[]
    fileout=open('pdbqtList.txt','w')
    for snapIn in pdbList:
        snapName=snapIn.split('/')[1].split('.')[0]
        snapOut='splitPDBQT/'+snapName+'.pdbqt'
        command=ligprep+' -r '+snapIn+' -o '+snapOut
        print command
        os.system(command)
        pdbqtList.append(snapOut)
        fileout.write(snapOut+'\n')
    fileout.close()
    return pdbqtList

def fpocket_box():
    '''
    This funtion makes a vina suitable input file from a binding site 
    computed by fpocket. NOT recommended for normal docking, as fpocket
    does not always find the correct binding site
    '''
    pocket='dry_apo_out/pockets/pocket0_vert.pqr'
    pocketIn=open(pocket,'r')
    x=[]
    y=[]
    z=[]
    for line in pocketIn:
        if line.startswith('ATOM'):
            line=line.split()
            x.append(float(line[5]))
            y.append(float(line[6]))
            z.append(float(line[7]))
        else:
            continue
    center=[numpy.mean(x),numpy.mean(y),numpy.mean(z)]
    size=[abs(max(x)-min(x)),abs(max(y)-min(y)),abs(max(z)-min(z))]
    pocketIn.close()
    return center, size

def vina_ensemble(pdbqtList,ligands,center,size):
    for ligand in ligands:
       lig=ligand.split('.')[0]
       command='mkdir splitPDBQT/'+lig
       os.system(command)

       # Generate vina input file
       vinaPar='ligand = '+ligand+'\n'+\
               'center_x = '+str(center[0])+'\n'+\
               'center_y = '+str(center[1])+'\n'+\
               'center_z = '+str(center[2])+'\n'+\
               'size_x = '+str(size[0])+'\n'+\
               'size_y = '+str(size[1])+'\n'+\
               'size_z = '+str(size[2])+'\n'+\
               'num_modes = 10'
       vinaConfName='splitPDBQT/'+lig+'/vina.txt'
       vinaConf=open(vinaConfName,'w')
       vinaConf.write(vinaPar)
       vinaConf.close()
    
       # Run Docking     
       for frame in pdbqtList:
           print "docking "+ligand+" in frame "+frame
          # dirName=lig+'/'+frame.split('.')[0]+'_docked'
           dirName=frame.split('/')[0]+'/'+lig+'/'+frame.split('/')[1].split('.')[0]+'_docked'
           mkdir='mkdir '+dirName
           print mkdir
           os.system(mkdir)
           docking='vina --config '+vinaConfName+\
                   ' --receptor '+frame+' '+\
                   '--out '+dirName+'/out.pdbqt '+\
                   '--log '+dirName+'/log.txt '
           print docking
           os.system(docking)
    return 0

def get_vina_scores(ligands,pdbqtList):
    scores={}
    for ligand in ligands:
        lig=ligand.split('.')[0]
        scores[lig]=[]
        for snapshot in pdbqtList:
            nameIn='splitPDBQT/'+ligand.split('.')[0]+'/'+\
                    snapshot.split('/')[1].split('.')[0]+'_docked/log.txt'
            #print nameIn
            filein=open(nameIn,'r')
            energy=float(filein.readlines()[25][12:19])
            print "getting docking score for snapshot: "+nameIn
            scores[lig].append(energy)    
            filein.close()
    return scores

def get_vina_rmsd(ligands,pdbqtList):
    rmsd={}
    for ligand in ligands:
        lig=ligand.split('.')[0]
        refStruct=lig+'_refStruct.pdb' # might or might not be a copy of ref
        reftpr=lig+'_refTop.tpr' # might ot might not be a copy of tpr
        rmsd[lig]=[]
        for snapshot in pdbqtList:
            nameIn='splitPDBQT/'+ligand.split('.')[0]+'/'+\
                    snapshot.split('/')[1].split('.')[0]+'_docked/out.pdbqt'
            nameMod1='splitPDBQT/'+ligand.split('.')[0]+'/'+\
                    snapshot.split('/')[1].split('.')[0]+'_docked/model1.pdbqt'
            filein=open(nameIn,'r')
            fileout=open(nameMod1,'w')
            for line in filein:
                if 'MODEL 2' in line:
                    break
                fileout.write(line)
            fileout.close()
            filein.close()
            dockPose=md.Universe(nameMod1)
            dockAtoms=dockPose.select_atoms('not name H*') #
            mobile=dockAtoms.coordinates()
            expPose=md.Universe(reftpr,refStruct)
            expAtoms=expPose.select_atoms('(resname '+lig+') and (not name H*)')
            static=expAtoms.coordinates()
            print "calculating RMSD for snapshot: "+nameIn
            rmsd[lig].append(rms.rmsd(static,mobile))
    return rmsd

def vina_plot(scores,rmsd,ligands):
    for ligand in ligands:
        lig=ligand.split('.')[0]
        nameOut=lig+'_dockAnalysis.out'
        fileout=open(nameOut,'w')
        fileout.write('Score RMSD\n')
        for i in range(0,len(rmsd[lig])):
            line=str(scores[lig][i])+' '+str(rmsd[lig][i])+'\n'
            fileout.write(line)            
    return 0