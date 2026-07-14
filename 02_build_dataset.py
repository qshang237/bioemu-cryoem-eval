import pandas as pd

df = pd.read_csv('data/emdb_parsed_df.csv', dtype=str).fillna('')

selected_emdb_ids = [
    'EMD-67847', 'EMD-64144', 'EMD-64147',
    'EMD-66194', 'EMD-66660', 'EMD-62584',
    'EMD-44530', 'EMD-44529',
    'EMD-60057', 'EMD-61839',
    'EMD-28651', 'EMD-28653',
    'EMD-64490', 'EMD-61397',
]

df_selected = df[df['emdb_id'].isin(selected_emdb_ids)].copy()

gltph_seq = df_selected[df_selected['emdb_id'] == 'EMD-44530']['sequence'].values[0]

manual_entry = pd.DataFrame([{
    'emdb_id': 'EMD-21986',
    'pdb_id': '6X12',
    'name': 'Inward-facing Apo-open state of the glutamate transporter homologue GltPh',
    'seq_len': len(gltph_seq),
    'resolution': '3.50',
    'uniprot_id': 'O59010',
    'sequence': gltph_seq,
}])

df_dataset = pd.concat([df_selected, manual_entry], ignore_index=True)

protein_labels = {
    'EMD-67847': 'SLC37A4', 'EMD-64144': 'SLC37A4', 'EMD-64147': 'SLC37A4',
    'EMD-66194': 'SLC37A4', 'EMD-66660': 'SLC37A4', 'EMD-62584': 'SLC37A4',
    'EMD-44530': 'GltPh', 'EMD-44529': 'GltPh', 'EMD-21986': 'GltPh',
    'EMD-60057': 'GPR4', 'EMD-61839': 'GPR4',
    'EMD-28651': 'SPNS2', 'EMD-28653': 'SPNS2',
    'EMD-64490': 'AuxinTransporter', 'EMD-61397': 'AuxinTransporter',
}
df_dataset['protein_label'] = df_dataset['emdb_id'].map(protein_labels)

df_dataset = df_dataset.sort_values(['protein_label', 'resolution'])
df_dataset = df_dataset[['emdb_id', 'pdb_id', 'protein_label', 'name',
                          'uniprot_id', 'resolution', 'seq_len', 'sequence']]

df_dataset.to_excel('data/dataset_multistate_v1.xlsx', index=False)
print(f'Saved {len(df_dataset)} entries to data/dataset_multistate_v1.xlsx')
print(df_dataset[['emdb_id', 'pdb_id', 'protein_label', 'resolution', 'seq_len']].to_string())

seq_file = 'data/sequences_v1.txt'
unique_seqs = df_dataset.drop_duplicates('protein_label')[['protein_label', 'sequence']]
with open(seq_file, 'w') as f:
    for _, row in unique_seqs.iterrows():
        f.write(f"{row['protein_label']} {row['sequence']}\n")
print(f'Saved {len(unique_seqs)} sequences to {seq_file}')

# Save colabfold format CSV for MSA subsampling
colabfold_csv = 'data/colabfold_input.csv'
unique_seqs_cf = df_dataset.drop_duplicates('protein_label')[['protein_label', 'sequence']]
unique_seqs_cf.columns = ['id', 'sequence']
unique_seqs_cf.to_csv(colabfold_csv, index=False)
print(f'Saved colabfold input to {colabfold_csv}')