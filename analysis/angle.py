import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from scipy.spatial.transform import Rotation
from glob import glob

# This script is used to calculate the TMD-CTD angle for YiiP

ref = mda.Universe('5vrf.pdb')
ref0 = ref.select_atoms('name CA and resid 212-288').positions - ref.select_atoms('name CA and resid 212-288').center_of_mass()

def matrix2angle(matrix):
    # Function to convert rotation matrix to rotation angle.
    angles=[]
    vectors=[]
    for i in matrix:
        phi = np.arccos((np.trace(i)-1)/2)
        vector = [(i[2,1]-i[1,2])/2/np.sin(phi), (i[0,2]-i[2,0])/2/np.sin(phi), (i[1,0]-i[0,1])/2/np.sin(phi),]
        angles.append(np.rad2deg(phi))
        vectors.append(vector)
    return angles

if __name__ == "__main__":
    for i in range(3):
        u = mda.Universe('md.gro', '/md{}.xtc'.format(i), in_memory=True)
        aligner = align.AlignTraj(u, ref, select="(name CA and ((resid 79-109) or (resid 178-208))) and (segid BP1 or segid AP1 or segid A or segid B)", in_memory=True).run()
        CT = u.select_atoms('(name CA and resid 212-288) and (segid BP1 or segid AP1 or segid A or segid B)')
        matrix = []
        angles = []
        for ts in u.trajectory:
            mobile = CT.positions - CT.center_of_mass()
            R, rmsd = align.rotation_matrix(mobile, ref0)
            matrix.append(R)

        matrix = np.array(matrix)
        angles = matrix2angle(matrix)

        np.save('matrix'.format(i), matrix)
        np.save('angles'.format(i), angles)
