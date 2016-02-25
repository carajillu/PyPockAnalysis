import os
import MDAnalysis as md
import MDAnalysis.analysis.hbonds as hb

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
