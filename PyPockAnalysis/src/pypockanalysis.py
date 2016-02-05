import MDAnalysis as md
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import MDAnalysis.analysis.hbonds as hb
import os, sys, math, numpy

###########################
# SECTION COORDINATES     #
###########################

def extract_reference(tpr,trr):
    doc=''' This function just extracts the first frame of a trr trajectory'''
    myMD=md.Universe(tpr,trr)
    atomGroup=myMD.select_atoms('all')
    #extract the reference frame
    nameOut='ref.gro'
    md.Writer(nameOut).write(atomGroup)
    return nameOut

def align(tpr,trr):
    doc='''this function aligns the whole trajectory to the fist frame'''
    nameRef=extract_reference(tpr,trr)
    myMD=md.Universe(tpr,trr)
    ref=md.Universe(tpr,nameRef)
    nameAli=trr.split('.')[0]+'_ali.trr'
    rms_fit_trj(myMD,ref,select='protein and name CA',filename=nameAli)
    return nameAli

def splitPDB(tpr,nameAli,resToIgnore,ext):
    doc='''this function extracts all the frames of the trajectory in the 
    chosen format'''
    myMD=md.Universe(tpr,nameAli)
    selection='not resname '+' and not resname '.join(resToIgnore)
    atomGroup=myMD.select_atoms(selection)
    dirName='split'+ext.upper()
    mkdir='mkdir '+dirName
    os.system(mkdir)
    listName=dirName+'.txt'
    frameList=[]
    fileout=open(listName,'w')
    i=0
    for ts in myMD.trajectory:
        PDBOut=dirName+'/'+trr.split('.')[0]+'_'+'0'*(5-len(str(i)))+str(i)+'.'+ext
        fileout.write(PDBOut+'\n')
        frameList.append(PDBOut)
        md.Writer(PDBOut).write(atomGroup)
        print 'wrote file '+PDBOut
        i=i+1
    fileout.close()
    return frameList

#########################################
# SECTION POCKET HUNTING                #
#########################################

def runFPocket(pdbList):
    for snapshot in pdbList:
        print 'Running Fpocket on file :'+snapshot
        command='fpocket -f '+snapshot
        os.system(command)
    return 0

def runMDPocket(pdbList):
    if '-fpocket' not in sys.argv:
       command='fpocket -f '+pdbList[0]
       os.system(command)
    alphaSph=pdbList[0].split('.')[0]+'_out/pockets/pocket0_vert.pqr'
    print "Running mdpocket starting at file: ",pdbList[0]
    command='mdpocket -L splitPDB.txt -f '+alphaSph
    print command
    os.system(command)
    return 0

def Fpocket_analysis(pdbList,trr):
    # Initializing pocket properties file
    propName=trr.split('.')[0]+'_druggability.txt'
    propOut=open(propName,'w')
    propOut.write('Pocket_Score Drug_Score Num_Vornoi Mean_alpha_radius Mean_alpha_SA Mean_Bf Hydrophob_score Polar_score Vol_score Real_Vol Char_score Loc_hydroph_dens_score Num_apol_spheres Prop_apol_sph\n')
    # Initializing pocket coordinates file
    pocketName=trr.split('.')[0]+'_pocket.pqr'
    pocketOut=open(pocketName,'w')
    i=0
    for snapshot in pdbList:
       props=[]
       SnapPocketName=snapshot.split('.')[0]+'_out/pockets/pocket0_vert.pqr'
       snapPocketIn=open(SnapPocketName,'r')
       print 'Opened file '+SnapPocketName+' for analysis'
       pocketOut.write('MODEL '+str(i)+'\n')
       i=i+1
       for line in snapPocketIn:
           if 'HEADER' in line and ' - ' in line:
               props.append((line.split(':')[1]).strip('\n').strip(' '))
           elif 'ATOM' in line:
               pocketOut.write(line)
       pocketOut.write('TER\nENDMDL\n\n')
       props=(' '.join(props))+'\n'
       propOut.write(props)
       snapPocketIn.close()
    propOut.close()
    pocketOut.close()
    return 0


#########################################
# SECTION DOCKING: vina                       #
#########################################
def genPDBQT(pdbList):
    '''
    I am NOT doing this with MDAnalysis because Vina doesn't take the pdbqt
    files produced by it. To make it work with cofactors, you have to trick
    the file MoleculePreparation.py (in AutoDockTools dir) to accept the
    resname of the cofactor
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
    return pdbqtList
        

def vina_ensemble(pdbqtList,pockFrame,ligand):
    # define box center and size
    pocket=pockFrame.split('.')[0]+'_out/pockets/pocket0_vert.pqr'
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
    
    # Generate vina input file
    vinaPar='ligand = '+ligand+'\n'+\
            'center_x = '+str(center[0])+'\n'+\
            'center_y = '+str(center[1])+'\n'+\
            'center_z = '+str(center[2])+'\n'+\
            'size_x = '+str(size[0])+'\n'+\
            'size_y = '+str(size[1])+'\n'+\
            'size_z = '+str(size[2])+'\n'+\
            'num_modes = 10'
    vinaConf=open('vina.txt','w')
    vinaConf.write(vinaPar)
    vinaConf.close()
    
    # Run Docking     
    for frame in pdbqtList:
        print "docking "+ligand+" in frame "+frame
        dirName=frame.split('.')[0]+'_docked'
        mkdir='mkdir '+dirName
        os.system(mkdir)
        docking='vina --config vina.txt '+\
                '--receptor '+frame+' '+\
                '--out '+dirName+'/out.pdbqt '+\
                '--log '+dirName+'/log.txt '
        print docking
        os.system(docking)
    return 0

#########################################
# SECTION HBONDS                        #
#########################################
def HbondCalc(resOfInt,tpr,trr):
    for residue in resOfInt:
        myMD=md.Universe(tpr,trr)
        resid='resid '+residue+' and not backbone'
        print resid
        h=hb.HydrogenBondAnalysis(myMD,resid, 'all', distance=3.0, angle=120.0)
        h.run()
        Hbonds=h.timeseries
        nameOut='hbonds_'+residue+'.txt'
        fileout=open(nameOut,'w')
        hbList=[]
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
ligand=sys.argv[3]
resToIgnore=['SOL','NA','CL',ligand.split('.')[0]] # maybe you don't need to ignore the ligand...

if '-split' in sys.argv:
   nameAli=align(tpr,trr)
   pdbList=splitPDB(tpr,nameAli,resToIgnore,'pdb')
else:
   frameList=[]
   filein=open('splitPDB.txt','r')
   for line in filein:
      frameList.append(line.strip('\n'))

if '-hbond' in sys.argv:
    resOfInt=[]
    entry=''
    while entry!='STOP':
        entry=raw_input('Enter a residue number (they begin at 0) or STOP to finish: ')
        if entry!='STOP':
            resOfInt.append(entry)  
        HbondCalc(resOfInt,tpr,trr)

if '-fpocket' in sys.argv:
    runFPocket(pdbList)
    Fpocket_analysis(pdbList,trr)

if '-mdpocket' in sys.argv:
    runMDPocket(pdbList)
    
if '-docking' in sys.argv:
    if '-fpocket' not in sys.argv:
       command='fpocket -f '+pdbList[0]
       os.system(command)
    pdbqtList=genPDBQT(pdbList)
    vina_ensemble(pdbqtList,pdbList[0],ligand)
    










#HbondCalc(resOfInt,tpr,trr)

#nameAli=align(tpr,trr)
#pdbList=splitPDB(tpr,nameAli,resToIgnore,'pdb')
#runFPocket(pdbList)
#Fpocket_analysis(pdbList,trr)
#runMDPocket(pdbList)

#pdbqtList=splitPDB(tpr,nameAli,resToIgnore,'pdbqt')
#pdbqtList=genPDBQT(pdbList)
'''
for snapshot in mol2List:
    pdbqtName=snapshot.split('.')[0]+'.pdbqt'
    conversion='/usr/local/MGLTools-1.5.6/bin/pythonsh '+\
               '/usr/local/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py '+\
               '-r '+snapshot+' -o '+pdbqtName
    pdbqtList.append(pdbqtName)'''
    
#vina_ensemble(pdbqtList,pdbList[0],ligand)
