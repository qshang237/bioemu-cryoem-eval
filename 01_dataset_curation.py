
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
 
import os
import pickle
import urllib.request
import pandas as pd
import mdtraj as md
from emdb.client import EMDB
 
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs('data', exist_ok=True)
os.makedirs('pdb_structures', exist_ok=True)
 
DF_CACHE = 'data/emdb_parsed_df.csv'

if os.path.exists(DF_CACHE):
    print(f'Loading parsed dataframe from cache: {DF_CACHE}')
    df = pd.read_csv(DF_CACHE, dtype=str).fillna('')
    print(f'Loaded {len(df)} entries from cache')
else:
    print('Querying EMDB...')
    client = EMDB()
    query = (
        '* AND sample_type:"protein"'
        ' AND structure_determination_method:"singleparticle"'
        ' AND resolution:[0 TO 4}'
        ' AND assembly_molecular_weight:{0 TO 80000]'
    )
    results = client.search(query)
    print(f'EMDB query returned: {len(results)} entries')

    filter_sequence = []
    for i, entry in enumerate(results):
        if i % 100 == 0:
            print(f'  Processing {i}/{len(results)}...')
        try:
            macros = entry.sample['macromolecule_list']['macromolecule']
            seq = macros[0]['sequence']['string']
            if 50 < len(seq) < 500:
                filter_sequence.append(entry)
        except:
            continue

    print(f'Sequence length 50-500: {len(filter_sequence)} entries')

    rows = []
    for entry in filter_sequence:
        try:
            macro = entry.sample['macromolecule_list']['macromolecule'][0]
            seq = macro['sequence']['string'].replace(' ', '')

            uniprot_id = ''
            try:
                uniprot_id = macro['sequence']['external_references'][0]['valueOf_']
            except:
                pass

            pdb_id = ''
            try:
                pdb_id = entry.related_pdb_ids[0]['pdb_id']
            except:
                pass

            entry_name = ''
            try:
                entry_name = entry.title
            except:
                try:
                    entry_name = macro['name']['valueOf_']
                except:
                    pass

            rows.append({
                'emdb_id': entry.id,
                'pdb_id': pdb_id,
                'name': entry_name,
                'seq_len': len(seq),
                'resolution': entry.resolution,
                'uniprot_id': uniprot_id,
                'sequence': seq,
            })
        except Exception as e:
            continue

    df = pd.DataFrame(rows)
    df.to_csv(DF_CACHE, index=False)
    print(f'Saved parsed dataframe to cache: {DF_CACHE}')
 
print(f'Parsed {len(df)} entries')
print(f'  with PDB: {(df["pdb_id"] != "").sum()}')
print(f'  with UniProt: {(df["uniprot_id"] != "").sum()}')
 
df_has_pdb = df[df['pdb_id'] != ''].copy()
print(f'Keeping entries with PDB: {len(df_has_pdb)}')
 
print('Downloading PDB/CIF files (skipping existing)...')
for _, row in df_has_pdb.iterrows():
    pdb_id = row['pdb_id']
    
    if os.path.exists(f'pdb_structures/{pdb_id}.pdb') or \
       os.path.exists(f'pdb_structures/{pdb_id}.cif'):
        continue
    
    downloaded = False
    for fmt in ['pdb', 'cif']:
        url = f'https://files.rcsb.org/download/{pdb_id}.{fmt}'
        filename = f'pdb_structures/{pdb_id}.{fmt}'
        try:
            urllib.request.urlretrieve(url, filename)
            print(f'  Downloaded {pdb_id}.{fmt}')
            downloaded = True
            break
        except:
            continue
    
    if not downloaded:
        print(f'  FAILED {pdb_id}: both pdb and cif unavailable')
 
single_chain_ids = []
failed_download = []

for _, row in df_has_pdb.iterrows():
    pdb_id = row['pdb_id']
    
    filename = None
    for fmt in ['pdb', 'cif']:
        candidate = f'pdb_structures/{pdb_id}.{fmt}'
        if os.path.exists(candidate):
            filename = candidate
            break
    
    if filename is None:
        failed_download.append(pdb_id)
        continue
    
    try:
        ref = md.load(filename)
        if ref.n_chains == 1:
            single_chain_ids.append(pdb_id)
    except Exception as e:
        print(f'  Could not load {pdb_id}: {e}')

print(f'Single chain (verified from PDB/CIF): {len(single_chain_ids)}')
print(f'Download failed (excluded): {len(failed_download)}')
if failed_download:
    print(f'Failed IDs: {failed_download}')
 
df_single = df_has_pdb[df_has_pdb['pdb_id'].isin(single_chain_ids)].copy()
print(f'Single chain (verified from PDB): {len(df_single)}')
 
df_single = df_single.copy()
df_single['has_his_tag'] = df_single['sequence'].str.contains('HHHHHH', na=False)
df_nohtag = df_single[~df_single['has_his_tag']].copy()
print(f'No His-tag: {len(df_nohtag)}')
 
uniprot_counts = df_nohtag[df_nohtag['uniprot_id'] != ''].groupby('uniprot_id').size()
multi_state_uniprots = uniprot_counts[uniprot_counts >= 2].index.tolist()
print(f'UniProt IDs with >= 2 entries: {len(multi_state_uniprots)}')
 
df_multi = df_nohtag[df_nohtag['uniprot_id'].isin(multi_state_uniprots)].copy()
 
also_no_uniprot = df_nohtag[df_nohtag['uniprot_id'] == '']
name_groups = also_no_uniprot.groupby('name').size()
multi_by_name = name_groups[name_groups >= 2].index.tolist()
df_multi_nouni = also_no_uniprot[also_no_uniprot['name'].isin(multi_by_name)]
 
df_candidates = pd.concat([df_multi, df_multi_nouni], ignore_index=True)
print(f'Total candidate entries (same protein >= 2 entries): {len(df_candidates)}')
 
df_candidates = df_candidates.sort_values(['uniprot_id', 'name', 'resolution'])
df_candidates.to_excel('data/dataset_multistate_candidates_v3.xlsx', index=False)
print(f'Saved to data/dataset_multistate_candidates_v3.xlsx')
 
print('\n--- Summary by UniProt ID ---')
for uid, group in df_candidates[df_candidates['uniprot_id'] != ''].groupby('uniprot_id'):
    print(f'\nUniProt: {uid} ({len(group)} entries)')
    for _, row in group.iterrows():
       print(f'  {row["emdb_id"]} | {row["pdb_id"]:8s} | {float(row["resolution"]):.2f}A | {row["name"][:80]}')
print('\n--- Entries without UniProt, grouped by name ---')
for name, group in df_candidates[df_candidates['uniprot_id'] == ''].groupby('name'):
    print(f'\n{name[:60]} ({len(group)} entries)')
    for _, row in group.iterrows():
        print(f'  {row["emdb_id"]} | {row["pdb_id"]:8s} | {float(row["resolution"]):.2f}A')

