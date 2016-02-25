import sys

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
    from atom_contacts import atom_contacts
    resOfInt=[]
    entry=''
    while entry!='STOP':
        entry=raw_input('Enter a residue number (they begin at 0) or STOP to finish: ')
        if entry!='STOP':
            resOfInt.append(entry)  
    atom_contacts.hbond_calc(resOfInt,tpr,trr)

if '-fpocket' in sys.argv:
    from pockets import pockets
    pockets.run_fpocket(pdbList)
    pockets.fpocket_analysis(pdbList,trr)

if '-mdpocket' in sys.argv:
    from pockets import pockets
    pockets.ref_pocket(nameRef,resToIgnore)
    pockets.run_mdpocket(pdbList)

if '-vina' in sys.argv:
    from vina import vina
    print 'Defining binding site. Units are angstrom:'
    center=[]
    center.append(raw_input("Enter the x component of the center: "))
    center.append(raw_input("Enter the y component of the center: "))
    center.append(raw_input("Enter the z component of the center: "))
    
    size=[]
    size.append(raw_input("give me the size in the x axis: "))
    size.append(raw_input("give me the size in the y axis: "))
    size.append(raw_input("give me the size in the z axis: "))
    
    pdbqtList=vina.gen_pdbqt(pdbList)
    vina.vina_ensemble(pdbqtList,ligands,center,size)
    scores=vina.get_vina_scores(ligands,pdbqtList)
    rmsd=vina.get_vina_rmsd(ligands,pdbqtList)
    vina.vina_plot(scores,rmsd,ligands)
    
if '-fpocket_docking' in sys.argv:
    from vina import vina
    if '-mdpocket' not in sys.argv:
       from pockets import pockets
       pockets.ref_pocket(nameRef,resToIgnore)
    pdbqtList=vina.gen_pdbqt(pdbList)
    center,size=vina.fpocket_box()
    vina.vina_ensemble(pdbqtList,ligands,center,size)
    scores=vina.get_vina_scores(ligands,pdbqtList)
    rmsd=vina.get_vina_rmsd(ligands,pdbqtList)
    vina.vina_plot(scores,rmsd,ligands)

if '-vina_analysis' in sys.argv:
    from vina import vina
    pdbqtList=[]
    filein=open('pdbqtList.txt','r')
    for line in filein:
        pdbqtList.append(line)
    scores=vina.get_vina_scores(ligands,pdbqtList)
    rmsd=vina.get_vina_rmsd(ligands,pdbqtList)
    vina.vina_plot(scores,rmsd,ligands)