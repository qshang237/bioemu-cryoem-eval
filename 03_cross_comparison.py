"""
00_cross_comparison.py
=======================
Standalone experimental cross-comparison matrix: for each protein, how well
does each deposited experimental structure (PDB_i) fit into each
experimental map (map_j)?

This answers "how structurally different are the experimental states
themselves" -- the baseline against which the model-fitting results
(04_fitting_bioemu.py / 04_fitting_msa.py / ...) should be read. If two
experimental states are nearly identical (high off-diagonal CCmask), no
model can be expected to tell them apart, and that's fine/expected. If they
are very different (low off-diagonal CCmask), a good model should show
different density ratios for the two states.

Uses real CCmask: Phenix-style masked Pearson correlation between a
B-factor-weighted Gaussian-blurred model density (rho_calc, via gemmi's
DensityCalculatorE) and the experimental map (rho_obs), restricted to a
mask around the model atoms (gemmi SolventMasker-style sphere mask).
This is bounded and self-fit naturally comes out close to (not forced to)
1.0 -- unlike the old compute_cross_matrix() in shared_fitting.py, which
reused the no-B-factor interpolation-ratio method meant for BioEmu
conformations and hard-coded the diagonal to 1.0. That was a bug.

Independent of any predicted ensemble -- only needs the experimental PDBs +
maps + dataset table, so it only needs to be run once (not once per model).
"""
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import mdtraj as md
import gemmi
import matplotlib.pyplot as plt

from shared_fitting import DATASET, download_map, download_pdb, superpose_to_target

OUT_DIR       = 'results/cross_comparison'
MASK_RADIUS   = 2.5  # Angstrom, standard atom-centred mask radius for CCmask
SELF_FIT_MIN  = 0.1  # a map's own self-fit CCmask must clear this to be used
                     # as the normalising denominator (guards against noisy/
                     # near-zero self-fits blowing up the ratio)
os.makedirs(OUT_DIR, exist_ok=True)


def ccmask(struct_i_path, target_pdb_path, target_map_grid, d_min,
          radius=MASK_RADIUS):
    """
    Real CCmask: Pearson correlation, restricted to a mask around the model
    atoms, between a Gaussian-blurred model density (generated from the
    model's own B-factors/occupancies via gemmi.DensityCalculatorE) and the
    experimental map. struct_i is superposed onto target_pdb's frame
    (common Calpha residues) before blurring, so it's evaluated in the
    target map's coordinate system.
    Returns (cc, n_mask_points).
    """
    st = gemmi.read_structure(struct_i_path)
    st.setup_entities()
    traj_i = md.load(struct_i_path)
    superpose_to_target(traj_i, target_pdb_path)

    gemmi_names  = [cra.atom.name for cra in st[0].all()]
    mdtraj_names = [a.name for a in traj_i.topology.atoms]
    if gemmi_names != mdtraj_names:
        raise ValueError(
            f'atom order mismatch between gemmi ({len(gemmi_names)}) and '
            f'mdtraj ({len(mdtraj_names)}) parse of {struct_i_path}')

    xyz_A = traj_i.xyz[0] * 10.0  # nm -> Angstrom
    for cra, pos in zip(st[0].all(), xyz_A):
        cra.atom.pos = gemmi.Position(float(pos[0]), float(pos[1]), float(pos[2]))

    dc = gemmi.DensityCalculatorE()
    dc.d_min = float(d_min)
    dc.rate  = 1.5
    st.cell          = target_map_grid.unit_cell
    st.spacegroup_hm = 'P 1'
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    rho_calc = dc.grid

    mask = gemmi.Int8Grid(target_map_grid.nu, target_map_grid.nv, target_map_grid.nw)
    mask.set_unit_cell(target_map_grid.unit_cell)
    mask.spacegroup = target_map_grid.spacegroup
    for cra in st[0].all():
        mask.set_points_around(cra.atom.pos, radius=radius, value=1)
    idx = np.argwhere(np.array(mask, copy=False) > 0)

    obs_arr  = np.array(target_map_grid, copy=False)
    obs_vals = obs_arr[idx[:, 0], idx[:, 1], idx[:, 2]]
    calc_vals = np.empty(len(idx), dtype=np.float64)
    for k in range(len(idx)):
        pos = target_map_grid.point_to_position(
            target_map_grid.get_point(int(idx[k, 0]), int(idx[k, 1]), int(idx[k, 2])))
        calc_vals[k] = rho_calc.interpolate_value(pos)

    if len(idx) < 10 or np.std(obs_vals) == 0 or np.std(calc_vals) == 0:
        return float('nan'), len(idx)
    return float(np.corrcoef(obs_vals, calc_vals)[0, 1]), len(idx)


def run():
    df = pd.read_excel(DATASET)
    all_rows = []

    for protein in df['protein_label'].unique():
        prot_rows   = df[df['protein_label'] == protein]
        emdb_ids    = prot_rows['emdb_id'].tolist()
        pdb_ids     = prot_rows['pdb_id'].tolist()
        resolutions = prot_rows['resolution'].tolist()
        n = len(emdb_ids)
        print(f'\n{"="*60}\n{protein}  ({n} experimental states)\n{"="*60}')

        pdb_paths, map_grids = {}, {}
        ok = True
        for eid, pid in zip(emdb_ids, pdb_ids):
            pdb_path = download_pdb(pid)
            map_path = download_map(eid)
            if pdb_path is None or map_path is None:
                print(f'  Missing PDB/map for {eid}/{pid}, skipping {protein}')
                ok = False
                break
            pdb_paths[eid] = pdb_path
            ccp4 = gemmi.read_ccp4_map(map_path)
            ccp4.setup(0.0)
            map_grids[eid] = ccp4.grid
        if not ok:
            continue

        cross = np.full((n, n), np.nan)
        for i, eid_i in enumerate(emdb_ids):
            for j, eid_j in enumerate(emdb_ids):
                try:
                    cc, npts = ccmask(pdb_paths[eid_i], pdb_paths[eid_j],
                                      map_grids[eid_j], d_min=resolutions[j])
                except Exception as e:
                    print(f'  {eid_i} -> {eid_j}: FAILED ({e})')
                    cc, npts = float('nan'), 0
                cross[i, j] = cc

        # Normalise each column j by that map's own self-fit (cross[j,j]),
        # i.e. the same raw/cc_ref convention used everywhere else in this
        # project (compute_score_matrix etc). This maps every diagonal entry
        # to exactly 1.0 -- not hard-coded, it falls out of self_fit/self_fit
        # -- and makes cross-protein comparisons meaningful: raw CCmask
        # ceilings differ per protein/map (resolution, particle count, model
        # quality), so "0.27" for one protein and "0.60" for another aren't
        # directly comparable, but "27% of achievable fit" and "60% of
        # achievable fit" are.
        self_fit = np.diag(cross).copy()
        bad = self_fit < SELF_FIT_MIN
        if bad.any():
            print(f'  WARNING: self-fit below {SELF_FIT_MIN} for '
                 f'{[emdb_ids[k] for k in np.where(bad)[0]]}; '
                 f'their columns will not be normalised (kept as NaN)')
            self_fit[bad] = np.nan
        cross_norm = cross / self_fit[np.newaxis, :]

        for i, eid_i in enumerate(emdb_ids):
            for j, eid_j in enumerate(emdb_ids):
                tag = ' (self-fit)' if i == j else ''
                print(f'  PDB[{pdb_ids[i]}] -> map[{eid_j}]: '
                     f'CCmask={cross[i,j]:.4f}  normalised={cross_norm[i,j]:.4f}{tag}')
                all_rows.append(dict(protein=protein, pdb_i=pdb_ids[i],
                                     emdb_i=eid_i, emdb_j=eid_j,
                                     ccmask=cross[i, j],
                                     ccmask_normalized=cross_norm[i, j]))

        short_ids = [e.replace('EMD-', '') for e in emdb_ids]
        fig, ax = plt.subplots(figsize=(max(4, n * 1.2), max(3.5, n * 1.0)))
        im = ax.imshow(cross_norm, cmap='RdYlGn', vmin=0, vmax=1.0, aspect='auto')
        ax.set_xticks(range(n)); ax.set_xticklabels(short_ids, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(n)); ax.set_yticklabels(short_ids, fontsize=9)
        ax.set_xlabel('Experimental map'); ax.set_ylabel('Deposited PDB (superposed)')
        ax.set_title(f'{protein}: Cross-comparison\n'
                     f'(CCmask normalised by each map\'s own self-fit)')
        plt.colorbar(im, ax=ax, label='Normalised CCmask (1.0 = self-fit)', shrink=0.85)
        for i in range(n):
            for j in range(n):
                val = cross_norm[i, j]
                txt = f'{val:.2f}' if not np.isnan(val) else 'n/a'
                ax.text(j, i, txt, ha='center', va='center', fontsize=8)
        fig.tight_layout()
        fig.savefig(f'{OUT_DIR}/{protein}_cross.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Figure -> {OUT_DIR}/{protein}_cross.png')

    pd.DataFrame(all_rows).to_excel(f'{OUT_DIR}/cross_comparison_all.xlsx', index=False)
    print(f'\nResults saved to {OUT_DIR}/')


if __name__ == '__main__':
    run()
