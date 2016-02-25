'''
MODULE POCKETS
1. Input Files:
   - Reference (pdb): pdb of the reference structure
   - pdbList (txt):   if split is not selected, the names of the aligned pdb
                      files will be read from this file.

2. Functions:
   2.1. run_fpocket:  runs fpocket to independently calculate pockets for
                      each file in pdbList (or pdbList.txt)
   2.2. ref_pocket: strips the reference structure of anything that is in
                    resToIgnore and runs fpocket on the result.
   2.3. run_mdpocket: runs mdpocket on all the snapshots in pdbList using as
                      a reference a pocket calculated from the X-Ray structure.
   2.4. fpocket_analysis: extracts the information from the top scored pocket of
                          each snapshot and puts all the properties in a file.
'''
import os

def run_fpocket(pdbList):
    for snapshot in pdbList:
        print 'Running Fpocket on file :'+snapshot
        command='fpocket -f '+snapshot
        os.system(command)
    return 0

def ref_pocket(nameRef,resToIgnore):
    #delete waters and ligands from nameRef
    filein=open(nameRef,'r')
    fileout=open('dry_apo.pdb','w')
    for line in filein:
        if 'ATOM' in line or 'HETATM' in line:
            if line.split()[3] in resToIgnore:
                continue
        fileout.write(line)
    filein.close()
    fileout.close()
    command='fpocket -f dry_apo.pdb'
    os.system(command)
    return 0

def run_mdpocket(pdbList):
    alphaSph='dry_apo_out/pockets/pocket0_vert.pqr'
    print "Running mdpocket starting at file: ",pdbList[0]
    command='mdpocket -L splitPDB.txt -f '+alphaSph
    print command
    os.system(command)
    return 0

def fpocket_analysis(pdbList,trr):
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