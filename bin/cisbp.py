###############################################################################
# Description:
#   This script reads a JSON file containing PLT gene-specific cisBP URLs,
#   fetches TF tables from those URLs, and filters them for relevant organisms
#   as specified in the JSON. It then applies organism-specific ID conversions
#   (e.g., for Oryza, Solanum, Glycine, and Zea) using the provided mapping files
#   (RAP-MSU, MaizeMine). The final table of PLT-related TFs is saved as a TSV.
###############################################################################

import pandas as pd
import numpy as np
import argparse
import json

def parse_table(url, org_list, genid):
    """
    Fetches HTML tables from a given URL, extracts the relevant table 
    (4th table on the page), and filters out unwanted columns and rows 
    while labeling each row with the given genid (PLT ID). 
    Keeps only entries matching the target organisms (org_list).
    """
    # Read all tables from the URL
    tables_on_page = pd.read_html(url)
    # Extract the 4th table 
    table = tables_on_page[4]
    # Assign genid to the first column
    table[table.columns[0]] = np.repeat(genid, len(table.index))
    # Drop unneeded columns
    table.drop(table.columns[[3,5,6]], axis=1, inplace=True)
    # Drop the first two rows (usually header/info rows in the table)
    table.drop([0,1], axis=0, inplace=True)
    # Rename columns
    table.columns = ['Ref', 'Organism', 'GeneID', 'SR_Score']
    # Keep only rows where 'Organism' is in org_list
    table = table.loc[[i in org_list for i in table['Organism']]]
    return table

# Command-line argument parsing
parser = argparse.ArgumentParser(
    description= 'This script filters a transfac database based on an orthologs table and modifies gene IDs.'
)
parser.add_argument('-j', '--cisbpjson', required=True, type=str, help='JSON file with cisBP urls.')
parser.add_argument('-r', '--rapmsu', required=True, type=str, help='RAP-MSU file (ID-conversion table).')
parser.add_argument('-z', '--zeaids', required=True, type=str, help='gene_to_alias_ids file from MaizeMine.')
parser.add_argument('-o', '--output', required=True, type=str, help='PLT related TFs output table file.')

args = parser.parse_args()

# Load JSON with cisBP URLs and a list of organisms to keep
with open(args.cisbpjson, 'r') as f:
    org_list, plt_urls = json.load(f).values()

# Initialize an empty DataFrame for storing merged TF data
tfs_plt_related = pd.DataFrame()

# Loop over each PLT entry and fetch/parse the table from the cisBP URL
for plt, plt_url in plt_urls.items():
    plt_table = parse_table(url=plt_url, genid=plt, org_list=org_list)
    tfs_plt_related = pd.concat([tfs_plt_related, plt_table], axis=0)

# Read the RAP-MSU ID mapping table for Oryza
rap_msu = pd.read_csv(
    args.rapmsu, 
    sep='\t', 
    header=None, 
    names=['orig', 'mapped'], 
    comment='#', 
    keep_default_na=False,
    na_values=['NaN']
)

# Remove version suffix after the first dot from 'mapped' column
mapped = []
for mapped_id in rap_msu.mapped:
    mapped_id = mapped_id.split('.', maxsplit=1)[0]
    mapped.append(mapped_id)
rap_msu.mapped = mapped

# Convert original IDs to uppercase for consistent matching
rap_msu.orig = rap_msu.orig.str.upper()

# Replace 'Oryza sativa' label with 'Oryza sativa Japonica Group'
tfs_plt_related.loc[tfs_plt_related.Organism == 'Oryza sativa', 'Organism'] = 'Oryza sativa Japonica Group'

# For each Oryza gene, map the GeneID using the RAP-MSU table
for gen in tfs_plt_related.loc[tfs_plt_related['Organism'] == 'Oryza sativa Japonica Group', 'GeneID']:
    tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = rap_msu.loc[rap_msu.orig == gen, 'mapped'].values[0]

# For Solanum lycopersicum, remove the version suffix from the GeneID
for gen in tfs_plt_related.GeneID[tfs_plt_related.Organism == 'Solanum lycopersicum']:
    tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = gen.split('.', maxsplit=1)[0]

# For Glycine max, convert GLYMA prefix to Glyma.
for gen in tfs_plt_related.GeneID[tfs_plt_related.Organism == 'Glycine max']:
    tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = gen.replace('GLYMA', 'Glyma.')

# Read the maize ID mapping file for Zea mays
zea_map = pd.read_csv(
    args.zeaids, 
    sep='\t', 
    header=None, 
    comment='#',
    names=['newID', 'newGenome', 'olderID', 'olderGenome']
)

# For Zea mays, map older IDs to the new ones. If no match, drop the row.
for gen in tfs_plt_related.GeneID[tfs_plt_related.Organism == 'Zea mays']:
    if gen in zea_map.olderID.values:
        tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = zea_map.newID[zea_map.olderID == gen].values[0]
    else:
        tfs_plt_related.drop(tfs_plt_related.index[tfs_plt_related.GeneID == gen], axis=0, inplace=True)

# Write the final merged table to a TSV
tfs_plt_related.to_csv(args.output, sep='\t', index=False)
