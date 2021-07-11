import MDAnalysis as mda
import numpy as np
import xarray as xr
from MDAnalysis.analysis.rms import RMSD
from pmda import rdf

def cal_site_rmsd(u, ref, filename):
    siteA = 'resid 47 51 159 155'
    siteB = 'resid 70 73 77'
    siteC = 'resid 285 263 287'
    siteD = 'resid 234 235 250 287'
    
    AA = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid A".format(siteA))
    AA.run(verbose=True)
    
    AB = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid B".format(siteA))
    AB.run(verbose=True)    
    
    BA = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid A".format(siteB))
    BA.run(verbose=True)
    
    BB = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid B".format(siteB))
    BB.run(verbose=True)
    
    CA = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid A".format(siteC))
    CA.run(verbose=True)
    
    CB = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid B".format(siteC))
    CB.run(verbose=True)    
    
    DA = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid A".format(siteD))
    DA.run(verbose=True)
    
    DB = RMSD(u, ref,
             select="{} and not backbone and not prop mass < 10.0 and segid B".format(siteD))
    DB.run(verbose=True)
    
    data = [[AA.rmsd[:,2], AB.rmsd[:,2]], [BA.rmsd[:,2], BB.rmsd[:,2]],
            [CA.rmsd[:,2], CB.rmsd[:,2]], [DA.rmsd[:,2], DB.rmsd[:,2]]]
    da = xr.DataArray(data, dims=['site', 'protomer', 'time'])
    da = da.assign_coords(site=['A', 'B', 'C1', 'C2'],
                          protomer=['A', 'B'],
                          time=AA.rmsd[:,1])
    da.to_netcdf(filename)

def cal_rdf_A(zn_ids, u, n, filename):
    siteA = u.select_atoms('(name OD1 OD2 and resid 47 51 159) or (name NE2 and resid 155)')
    AA_atom = siteA.select_atoms('segid A')
    AB_atom = siteA.select_atoms('segid B')
    sitename = ["{}{}{}".format(mda.lib.util.convert_aa_code(atom.resname), atom.resid, atom.name) for atom in AA_atom]
    AA = [u.select_atoms('bynum {}'.format(zn_ids[0]+1)), AA_atom]
    AB = [u.select_atoms('bynum {}'.format(zn_ids[1]+1)), AB_atom]
    RDF = rdf.InterRDF_s(u, [AA, AB], nbins=300, range=(0, 3.0))
    RDF.run(n_blocks=n, n_jobs=n)

    da = xr.DataArray([RDF.rdf[0][0], RDF.rdf[1][0]], dims=['ZN', 'site', 'distance'])
    da = da.assign_coords(ZN=zn_ids,
                          site=sitename,
                          distance=RDF.bins)
    da.to_netcdf(filename)
