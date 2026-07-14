"""
04_fitting_bioemu.py
====================
Real-space model-map evaluation for BioEmu ensemble.
Reads topology.pdb + samples.xtc from bioemu_results/{protein}/
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import mdtraj as md
from shared_fitting import run_analysis

BIOEMU_DIR = 'bioemu_results'
OUT_DIR    = 'results/bioemu'


def load_bioemu(protein):
    xtc_path = f'{BIOEMU_DIR}/{protein}/samples.xtc'
    top_path = f'{BIOEMU_DIR}/{protein}/topology.pdb'
    if not os.path.exists(xtc_path):
        return None
    return md.load(xtc_path, top=top_path)


if __name__ == '__main__':
    run_analysis(load_bioemu, model_name='bioemu', out_dir=OUT_DIR)