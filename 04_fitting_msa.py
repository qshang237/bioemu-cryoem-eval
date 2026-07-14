"""
04_fitting_msa.py
=================
Real-space model-map evaluation for MSA subsampling ensemble.
Reads {protein}/{protein}_unrelaxed_rank_*.pdb from msa_results/{protein}/
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import glob
import mdtraj as md
from shared_fitting import run_analysis
import numpy as np

MSA_DIR = 'msa_results'
OUT_DIR = 'results/msa'


def load_msa(protein):
    pdb_dir = f'{MSA_DIR}/{protein}'
    if not os.path.isdir(pdb_dir):
        return None
    pattern = f'{pdb_dir}/{protein}_unrelaxed_rank_*.pdb'
    pdbs = sorted(glob.glob(pattern))
    if not pdbs:
        return None
    print(f'  Found {len(pdbs)} PDB files for {protein}')

    # use first PDB topology as reference
    ref = md.load(pdbs[0])
    n_atoms = ref.n_atoms
    
    xyz_list = [ref.xyz[0]]
    for p in pdbs[1:]:
        t = md.load(p)
        if t.n_atoms == n_atoms:
            xyz_list.append(t.xyz[0])
    
    xyz = np.stack(xyz_list, axis=0)  # (n_frames, n_atoms, 3)
    traj = md.Trajectory(xyz, ref.topology)
    print(f'  Kept {traj.n_frames} frames')
    return traj


if __name__ == '__main__':
    run_analysis(load_msa, model_name='msa', out_dir=OUT_DIR)
