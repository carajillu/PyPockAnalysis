import MDAnalysis as md
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import MDAnalysis.analysis.hbonds as hb
import os, sys, math, numpy

#########################################


#########################################
# SECTION DOCKING: vina                       #
#########################################
def genPDBQT(pdbList):
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
        

def fpocket_docking_site(nameRef):
    doc='''
        This funtion makes a vina suitable input file from a binding site 
        computed by fpocket. NOT recommended for normal docking, as fpocket
        does not always find the correct binding site
        '''
    pocket=nameRef.split('.')[0]+'_out/pockets/pocket0_vert.pqr'
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

def vina_plot(scores,rmsd,lignads):
    for ligand in ligands:
        lig=ligand.split('.')[0]
        nameOut=lig+'_dockAnalysis.out'
        fileout=open(nameOut,'w')
        fileout.write('Snapshot Score RMSD\n')
        for i in range(0,len(rmsd[lig])):
            line=str(scores[lig][i])+' '+str(rmsd[lig][i])+'\n'
            fileout.write(line)            
    return 0  
            
    

#########################################
# SECTION HBONDS                        #
#########################################
def hbond_calc(resOfInt,tpr,trr):
    for residue in resOfInt:
        myMD=md.Universe(tpr,trr)
        resid='resid '+residue+' and not backbone'
        print resid
        h=hb.HydrogenBondAnalysis(myMD,resid, 'all', distance=3.0, angle=120.0)
        h.run()
        Hbonds=h.timeseries
        nameOut='hbonds_'+residue+'.txt'
        fileout=open(nameOut,'w')
        for frame in Hbonds:
            line=[frame[0][2].split(':')[0]]
            for bond in frame:
                line.append(bond[3].split(':')[0])
            fileout.write(' '.join(line)+'\n')
    return 0
         
#########################################
# SECTION SOLVATION                     #
#########################################



#########################################
# MAIN FUNCITON                         #
#########################################

tpr=sys.argv[1]
trr=sys.argv[2]
nameRef=sys.argv[3]
resToIgnore=['SOL','WAT','HOH','NA','CL'] # maybe you don't need to ignore the ligand...
ligands=[]
lig=''
while lig!='STOP':
   lig=raw_input('Enter the name of the ligand PDBQT file or STOP to finish: ')
   if lig!='STOP':
      ligands.append(lig)
      resToIgnore.append(lig.split('.')[0])

if '-split' in sys.argv:
   from coordinates import coordinates
   nameAli=coordinates.align(tpr,trr,nameRef)
   pdbList=coordinates.split(tpr,nameAli,resToIgnore,'pdb')
else:
   pdbList=[]
   filein=open('splitPDB.txt','r')
   for line in filein:
      pdbList.append(line.strip('\n'))

if '-hbond' in sys.argv:
    resOfInt=[]
    entry=''
    while entry!='STOP':
        entry=raw_input('Enter a residue number (they begin at 0) or STOP to finish: ')
        if entry!='STOP':
            resOfInt.append(entry)  
    hbond_calc(resOfInt,tpr,trr)

if '-fpocket' in sys.argv:
    from pockets import pockets
    pockets.runFPocket(pdbList)
    pockets.fpocket_analysis(pdbList,trr)

if '-mdpocket' in sys.argv:
    from pockets import pockets
    pockets.runMDPocket(pdbList,nameRef,resToIgnore)

if '-vina' in sys.argv:
    print 'Defining binding site. Units are angstrom:'
    center=[]
    center.append(raw_input("Enter the x component of the center: "))
    center.append(raw_input("Enter the y component of the center: "))
    center.append(raw_input("Enter the z component of the center: "))
    
    size=[]
    size.append(raw_input("give me the size in the x axis: "))
    size.append(raw_input("give me the size in the y axis: "))
    size.append(raw_input("give me the size in the z axis: "))
    
    pdbqtList=genPDBQT(pdbList)
    vina_ensemble(pdbqtList,ligands,center,size)
    scores=get_vina_scores(ligands,pdbqtList)
    rmsd=get_vina_rmsd(ligands,pdbqtList)
    
if '-fpocket_docking' in sys.argv:
    if '-fpocket' not in sys.argv:
       command='fpocket -f '+nameRef
       os.system(command)
    pdbqtList=genPDBQT(pdbList)
    center,size=fpocket_docking_site(nameRef)
    vina_ensemble(pdbqtList,ligands,center,size)

if '-vina_analysis' in sys.argv:
    pdbqtList=[]
    filein=open('pdbqtList.txt','r')
    for line in filein:
        pdbqtList.append(line)
    scores=get_vina_scores(ligands,pdbqtList)
    rmsd=get_vina_rmsd(ligands,pdbqtList)
    vina_plot(scores,rmsd,ligands)