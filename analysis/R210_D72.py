import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance
from tqdm import tqdm

# This script is used to calculate the D72-R210 distances for YiiP

def link_loop_AB(u):
    R210A = u.select_atoms('resid 210 and name CB and segid A')
    R210B = u.select_atoms('resid 210 and name CB and segid B')
    D72A = u.select_atoms('resid 72 and name CA and segid A')
    D72B = u.select_atoms('resid 72 and name CA and segid B')


    distances = []
    for ts in tqdm(u.trajectory):
        pairsAB, distAB = capped_distance(R210A.positions,
                                          D72B.positions,
                                          30,
                                          box=u.dimensions)
        pairsBA, distBA = capped_distance(R210B.positions,
                                          D72A.positions,
                                          30,
                                          box=u.dimensions)
        if len(distAB) == 0:
            distAB=[30]
        if len(distBA) == 0:
            distBA=[30]

        distances.append([np.min(distAB), np.min(distBA)])
    return np.array(distances)


for i in range(3):
    u = mda.Universe('md.gro', 'md{}.xtc'.format(i), in_memory=True)
    ref = mda.Universe('5vrf.pdb')
    res=link_loop_AB(u)
    np.save('md{}_210_72_AB'.format(i), res)