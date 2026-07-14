"""
04_fitting_alphaflow.py
========================
Real-space model-map evaluation for AlphaFlow ensemble.
Reads {protein}.pdb from alphaflow_results/ -- a single multi-MODEL PDB
file per protein (each MODEL block is one sampled conformation, up to 500
per protein). mdtraj loads multi-MODEL PDBs natively as a multi-frame
trajectory, so no special-casing is needed vs the BioEmu/MSA loaders.
"""
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import time
import mdtraj as md
from shared_fitting import run_analysis

ALPHAFLOW_DIR = 'alphaflow_results'
OUT_DIR       = 'results/alphaflow'


def load_alphaflow(protein):
    pdb_path = f'{ALPHAFLOW_DIR}/{protein}.pdb'
    if not os.path.exists(pdb_path):
        return None
    print(f'  Loading {pdb_path} (multi-MODEL PDB, mdtraj parses this in '
         f'pure Python -- large files can take 30s-2min, this is normal, '
         f'not stuck)...', flush=True)
    t0 = time.time()
    traj = md.load(pdb_path)
    print(f'  Loaded {traj.n_frames} frames from {pdb_path} '
         f'({time.time()-t0:.0f}s)')
    return traj


if __name__ == '__main__':
    run_analysis(load_alphaflow, model_name='alphaflow', out_dir=OUT_DIR)
