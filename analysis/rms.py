import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF
import os

# This script is used to calculate the RMSD and RMSF for YiiP

def rms_analysis(u, ref, folder):
    HELIX1 = "resid 11-37"          
    HELIX2 = "resid 38-72"
    HELIX3 = "resid 73-112"          
    HELIX4 = "resid 113-146"          
    HELIX5 = "resid 147-176"          
    HELIX6 = "resid 177-211"          
    HELIX7 = "resid 212-228" 
    HELIX8 = "resid 256-278"          
    SHEET1 = "resid 229-242"                                              
    SHEET2 = "resid 243-255"                 
    SHEET3 = "resid 279-288"                 

    TM = "(({}) or ({}) or ({}) or ({}) or ({}) or ({}))".format(HELIX1, HELIX2, HELIX3, HELIX4, HELIX5, HELIX6)
    CT = "(({}) or ({}) or ({}) or ({}) or ({}))".format(HELIX7, HELIX8, SHEET1, SHEET2, SHEET3)


    #RMSD

    R = RMSD(u, ref,
             select="name CA and resid 11-288",             # superimpose on whole backbone of the whole protein
             groupselections=["name CA and resid 11-288"])
    R.run(verbose=True)
    np.save('{}rmsd_sys_whole'.format(folder), R.rmsd)

    R_TM = RMSD(u, ref,
               select="name CA and {}".format(TM),             # superimpose on whole backbone of the whole protein
               groupselections=["name CA and resid 11-288"])                                   
    R_TM.run(verbose=True)
    np.save('{}rmsd_sys_TM'.format(folder), R_TM.rmsd)

    R_CT = RMSD(u, ref,
               select="name CA and {} ".format(CT),             # superimpose on whole backbone of the whole protein
               groupselections=["name CA and resid 11-288"])
    R_CT.run(verbose=True)
    np.save('{}rmsd_sys_CT'.format(folder), R_CT.rmsd)

    #RMSF

    #whole
    prealigner = align.AlignTraj(u, u, select="protein",
                                 in_memory=True).run()
    protein = u.select_atoms("protein")

    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(
                ref_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select="protein and name CA", in_memory=True).run()
    calphas = protein.select_atoms("name CA")

    rmsfer = RMSF(calphas, verbose=True).run()
    np.save('{}rmsf_whole'.format(folder), rmsfer.rmsf)

    #CT
    prealigner = align.AlignTraj(u, u, select="protein and {}".format(CT),
                                 in_memory=True).run()
    protein = u.select_atoms("protein")

    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(
                ref_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select="protein and name CA and {}".format(CT), in_memory=True).run()
    calphas = protein.select_atoms("name CA")

    rmsfer = RMSF(calphas, verbose=True).run()
    np.save('{}rmsf_CT'.format(folder), rmsfer.rmsf)

    #TM
    prealigner = align.AlignTraj(u, u, select="protein and {}".format(TM),
                                 in_memory=True).run()
    protein = u.select_atoms("protein")

    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(
                ref_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select="protein and name CA and {}".format(TM), in_memory=True).run()
    calphas = protein.select_atoms("name CA")

    rmsfer = RMSF(calphas, verbose=True).run()
    np.save('{}rmsf_TM'.format(folder), rmsfer.rmsf)

for i in range(3):
    u = mda.Universe('md.gro', 'md{}.xtc'.format(i))
    ref = mda.Universe('5vrf.pdb')
    folder = 'md{}/'.format(i)
    rms_analysis(u, ref, folder)
