#%%
import pandas as pd
import pybedtools
import os
#%%
motif_folder = '/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/homer_known_motif_hg38'
#%%
entries = os.listdir(motif_folder)
# Filter out directories, only keep files
motif_files = [entry for entry in entries if (os.path.isfile(os.path.join(motif_folder, entry))) and (entry != 'homer.KnownMotifs.hg38.191020.bed')]
# %%
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
subfamily = 'THE1C'
coord_file = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/coord_internal_id/{subfamily}.txt'
coord_df=pd.read_csv(coord_file, sep='\t', header= None)
coord_bed=pybedtools.BedTool.from_dataframe(coord_df)
#%%
motif_file = f'{motif_folder}/{motif_files[0]}'
motif_df = pd.read_csv(motif_file, sep='\t', header=None)
motif_bed = pybedtools.BedTool.from_dataframe(motif_df)
# %%
intersect_df=coord_bed.intersect(motif_bed, c=True, s=True).to_dataframe()
# %%
intersect_df[intersect_df['thickStart']>0].shape[0]
# %% loop
import sys
sys.path.append('/home/pc575/rds/rds-mi339-kzfps/users/pakkanan/phd_project_development/dev/packaging_dir/ma_mapper/test/te_age')
import config_hg38 as config
subfamily = 'THE1C'
coord_file = f'/rds/project/rds-XrHDlpCeVDg/users/pakkanan/phd_project_development/data/_mapper_output/hg38_repeatmasker/coord_internal_id/{subfamily}.txt'
coord_df=pd.read_csv(coord_file, sep='\t', header= None)
coord_bed=pybedtools.BedTool.from_dataframe(coord_df)
overlap_count = []
for motif_file in motif_files:
    print(motif_file)
    motif_filepath = f'{motif_folder}/{motif_file}'
    motif_df = pd.read_csv(motif_filepath, sep='\t', header=None)
    motif_bed = pybedtools.BedTool.from_dataframe(motif_df)
    intersect_df=coord_bed.intersect(motif_bed, c=True, s=True).to_dataframe()
    overlap_count.append(intersect_df[intersect_df['thickStart']>0].shape[0])
#%%
motif_count_df=pd.DataFrame({'motif':motif_files,'overlap_count':overlap_count})
motif_count_df['motif'] = motif_count_df['motif'].str.replace('.bed', '', regex=False)
# %%
motif_count_df.to_csv(f'{motif_folder}/{subfamily}_count.txt',sep='\t')
# %%
motif_count_df.sort_values('overlap_count', ascending=False)
# %%
motif_count_df[motif_count_df['overlap_count']>2000]
# %%
import re
def extract_motif_name(motif):
    # Match everything before the first parenthesis or colon
    match = re.match(r'^[^\(]+', motif)
    if match:
        # Remove everything after the first colon
        if ':' in name:
            name = name.split(':')[0]
        return match.group(0).strip()
    return motif

# Apply the function to the 'motif' column
motif_count_df['motif_name'] = motif_count_df['motif'].apply(extract_motif_name)
# %%
for idx, row in motif_count_df.iterrows():
    print(row['motif_name'])
# %%
import pyjaspar 
import requests

# Initialize JASPAR database
jaspar=pyjaspar.jaspardb(release='JASPAR2024')
def motif_to_tf(motif_name):
    # Search for the motif in JASPAR
    results = jaspar.search(motif_name)
    tf_names = [tf.name for tf in results]
    return tf_names
def get_target_genes(tf_name):
    # Use Enrichr to get target genes for a TF
    url = "https://maayanlab.cloud/Enrichr/addList"
    response = requests.post(url, files={"list": (None, tf_name)})
    if response.status_code == 200:
        list_id = response.json()["userListId"]
        enrich_url = f"https://maayanlab.cloud/Enrichr/enrich?dataset={list_id}&backgroundType=ChEA"
        enrich_response = requests.get(enrich_url)
        if enrich_response.status_code == 200:
            data = enrich_response.json().get("ChEA_2021", [])
            genes = []
            for item in data:
                genes.extend(item.get("genes", []))
            return list(set(genes))  # Remove duplicates
    return []
#%%
#%%
from pyjaspar import jaspardb
jdb_obj = jaspardb(release='JASPAR2024')
motifs = jdb_obj.fetch_motifs_by_name('Pitx1:Ebox')
#%%
#Map motifs to transcription factors
motif_count_df['transcription_factors'] = motif_count_df['motif_name'].apply(motif_to_tf)
#%%
# Flatten the list of TFs and get unique TFs
all_tfs = set(tf for sublist in motif_count_df['transcription_factors'] for tf in sublist)

# Get target genes for each TF
tf_to_genes = {tf: get_target_genes(tf) for tf in all_tfs}

# Map TFs to target genes in the DataFrame
motif_count_df['target_genes'] = motif_count_df['transcription_factors'].apply(lambda tf_list: [tf_to_genes.get(tf, []) for tf in tf_list])

# Print the DataFrame with target genes
print(df)