import pandas as pd
#pbase and pset


# === Configuration ===
NUM_MODULES = 7
GRN_TEMPLATE = '/Users/irmaskaljo/codedev/masters_thesis/Data/Enrichment/Yeast_0.01/Raw/select_pbase1_r5_fpkm_6400_net_m_{}_0.01'
REACTOME_TEMPLATE = '/Users/irmaskaljo/codedev/masters_thesis/Data/Enrichment/Yeast_0.01/results/Reactome/select_pset1_r5_fpkm_6400_net_m_{}_0.01_al_1'
KEGG_TEMPLATE = '/Users/irmaskaljo/codedev/masters_thesis/Data/Enrichment/Yeast_0.01/results/KEGG/select_pset1_r5_fpkm_6400_net_m_{}_0.01_al_2'
OUTPUT_GRN = '/Users/irmaskaljo/codedev/masters_thesis/Data_preprocessing/Yeast_Enrichment/scr_size_10/processed_grn_yeast_0.01.csv'
OUTPUT_ENRICH = '/Users/irmaskaljo/codedev/masters_thesis/Data_preprocessing/Yeast_Enrichment/scr_size_10/processed_enrichment_0.01.csv'


# === Step 1: Preprocess GRNs ===
grn_list = []

for module in range(1, NUM_MODULES + 1):
    df = pd.read_csv(GRN_TEMPLATE.format(module), sep='\t', index_col=0)
    df = df.reset_index().rename(columns={'index': 'Reg'})
    df[['Regulator', 'TargetGene']] = df['Reg'].str.split('_', expand=True)
    
    df_long = df[['Reg', 'Regulator', 'TargetGene', f'weight_{module}']].copy()
    df_long.rename(columns={f'weight_{module}': 'Weight'}, inplace=True)
    df_long['Module'] = module
    grn_list.append(df_long)

grn_df = pd.concat(grn_list, ignore_index=True)
grn_df.to_csv(OUTPUT_GRN, index=False)

# === Step 2: Combine Reactome & KEGG Enrichment ===
enrich_list = []

def parse_enrichment(path, module, source):
    df = pd.read_csv(path, sep='\t')
    records = []
    for _, row in df.iterrows():
        genes = row['Gene_Set'].split(',')
        for gene in genes:
            records.append({
                'Gene': gene,
                'Path_ID': row['Path_ID'],
                'Path_Description': row['Path_Description'],
                'Module': module,
                'Source': source
            })
    return pd.DataFrame(records)

for module in range(1, NUM_MODULES + 1):
    reactome_df = parse_enrichment(REACTOME_TEMPLATE.format(module), module, 'Reactome')
    kegg_df = parse_enrichment(KEGG_TEMPLATE.format(module), module, 'KEGG')
    enrich_list.extend([reactome_df, kegg_df])

enrich_df = pd.concat(enrich_list, ignore_index=True)
enrich_df.to_csv(OUTPUT_ENRICH, index=False)

print("Preprocessing complete. Files saved:")
print(f"- {OUTPUT_GRN}")
print(f"- {OUTPUT_ENRICH}")
