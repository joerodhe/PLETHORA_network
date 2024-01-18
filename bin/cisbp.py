import pandas as pd
import numpy as np
import argparse
import json

def parse_table(url, org_list, genid):
    tables_on_page = pd.read_html(url)
    table = tables_on_page[4]
    table[table.columns[0]] = np.repeat(genid, len(table.index))
    table.drop(table.columns[[3,5,6]], axis=1, inplace = True)
    table.drop([0,1], axis=0, inplace = True)
    table.columns = ['Ref', 'Organism', 'GeneID', 'SR_Score']
    table = table.loc[[i in org_list for i in table['Organism']]]
    return table

parser=argparse.ArgumentParser(description= 'This script filters a transfac database based on a orthologs table and plot the logos.')

parser.add_argument('-j', '--cisbpjson', required=True, type=str, help='JSON file with cisBP urls.')
parser.add_argument('-r', '--rapmsu', required=True, type=str, help='RAP-MSU file (ID-conversion table).')
parser.add_argument('-z', '--zeaids', required=True, type=str, help='gene_to_alias_ids file from MaizeMine.')
parser.add_argument('-o', '--output', required=True, type=str, help='PLT related TFs output table file.')

args=parser.parse_args()

with open(args.cisbpjson, 'r') as f:
    org_list, plt_urls = json.load(f).values()

tfs_plt_related = pd.DataFrame()
for plt, plt_url in plt_urls.items():
    plt_table = parse_table(url=plt_url, genid=plt, org_list=org_list)
    tfs_plt_related = pd.concat([tfs_plt_related, plt_table], axis=0)

rap_msu = pd.read_csv(
    args.rapmsu, 
    sep='\t', 
    header=None, 
    names=['orig', 'mapped'], 
    comment='#', 
    keep_default_na=False,
    na_values=['NaN'])

mapped = []
for mapped_id in rap_msu.mapped:
    mapped_id = mapped_id.split('.', maxsplit=1)[0]
    mapped.append(mapped_id)
rap_msu.mapped = mapped
rap_msu.orig = rap_msu.orig.str.upper()

tfs_plt_related.loc[tfs_plt_related.Organism == 'Oryza sativa', 'Organism'] = 'Oryza sativa Japonica Group'

for gen in tfs_plt_related.loc[tfs_plt_related['Organism'] == 'Oryza sativa Japonica Group', 'GeneID']:
    tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = rap_msu.loc[rap_msu.orig == gen, 'mapped'].values[0]

for gen in tfs_plt_related.GeneID[tfs_plt_related.Organism == 'Solanum lycopersicum']:
    tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = gen.split('.', maxsplit=1)[0]

for gen in tfs_plt_related.GeneID[tfs_plt_related.Organism == 'Glycine max']:
    tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = gen.replace('GLYMA', 'Glyma.')

zea_map = pd.read_csv(args.zeaids, sep='\t', header=None, comment='#',
                      names=['newID', 'newGenome', 'olderID', 'olderGenome'])

for gen in tfs_plt_related.GeneID[tfs_plt_related.Organism == 'Zea mays']:
    if gen in zea_map.olderID.values:
        tfs_plt_related.loc[tfs_plt_related.GeneID == gen, 'GeneID'] = zea_map.newID[zea_map.olderID == gen].values[0]
    else:
        tfs_plt_related.drop(tfs_plt_related.index[tfs_plt_related.GeneID == gen], axis=0, inplace=True)

tfs_plt_related.to_csv(args.output, sep='\t', index=False)
