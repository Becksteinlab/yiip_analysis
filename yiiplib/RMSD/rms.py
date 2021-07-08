import MDAnalysis as mda
import numpy as np
import xarray as xr
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

def cal_rmsds(u, ref, filename, whole_sel):
    R = RMSD(u, ref,
             select=whole_sel,
             groupselections=[whole_sel])
    R.run(verbose=True)

    R_TM = RMSD(u, ref,
                select="{} and {}".format(whole_sel, u.tmd_selection_string()),
                groupselections=[whole_sel]) 
    R_TM.run(verbose=True)

    R_CT = RMSD(u, ref,
                select="{} and {}".format(whole_sel, u.ctd_selection_string()),
                groupselections=[whole_sel])
    R_CT.run(verbose=True)
    RMSDs=[R.rmsd[:,2:], R_CT.rmsd[:,2:], R_TM.rmsd[:,2:]]
    da = xr.DataArray(RMSDs, dims=['alignment', 'time', 'sel'])
    da = da.assign_coords(alignment=['whole', 'ctd', 'tmd'],
                          time=R.rmsd[:,1],
                          sel=['fit', 'whole'])
    da.to_netcdf(filename)

def cal_rmsfs(u, ref, filename):
    RMSFs = []
    #whole
    prealigner = align.AlignTraj(u, ref, select="name CA",
                                 in_memory=True).run()
    protein = u.select_atoms("protein")

    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(
                ref_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select="name CA", in_memory=True).run()
    calphas = protein.select_atoms("name CA")

    rmsfer = RMSF(calphas, verbose=True).run()
    
    n_res = int(len(rmsfer.rmsf)/2)
    resids = calphas.resids[:n_res]

    RMSFs.append([rmsfer.rmsf[:n_res], rmsfer.rmsf[n_res:]])

    #CTD
    prealigner = align.AlignTraj(u, ref, select='name CA and ' + u.ctd_selection_string(),
                                 in_memory=True).run()
    protein = u.select_atoms("protein")

    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(
                ref_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select='name CA and ' + u.ctd_selection_string(), in_memory=True).run()
    calphas = protein.select_atoms("name CA")

    rmsfer = RMSF(calphas, verbose=True).run()

    RMSFs.append([rmsfer.rmsf[:n_res], rmsfer.rmsf[n_res:]])

    #TMD
    prealigner = align.AlignTraj(u, u, select='name CA and ' + u.tmd_selection_string(),
                                 in_memory=True).run()
    protein = u.select_atoms("protein")

    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    reference = mda.Merge(protein).load_new(
                ref_coordinates[:, None, :], order="afc")
    aligner = align.AlignTraj(u, reference, select='name CA and ' + u.tmd_selection_string(), in_memory=True).run()
    calphas = protein.select_atoms("name CA")

    rmsfer = RMSF(calphas, verbose=True).run()

    RMSFs.append([rmsfer.rmsf[:n_res], rmsfer.rmsf[n_res:]])

    da = xr.DataArray(RMSFs, dims=['alignment', 'protomer', 'resid'])
    da = da.assign_coords(alignment=['whole', 'ctd', 'tmd'],
                          protomer=['A', 'B'],
                          resid=resids)
    da.to_netcdf(filename)
