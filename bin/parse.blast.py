import pandas as pd
import argparse
from os.path import join

parser=argparse.ArgumentParser(description= 'This script joins (merge) blastp tables.')

parser.add_argument('-o', '--organisms', required=True, type=str, help='Organism list (comma separated)')
parser.add_argument('-f', '--files', required=True, type=str, help='BLASTp/orthologs table files (comma separated)')
parser.add_argument('-r', '--reference', required=True, type=str, help='Reference organism.')
parser.add_argument('-u', '--output', required=True, type=str, help='Output file.')

args=parser.parse_args()

files = args.files.split(",")
organisms = args.organisms.split(",")

files = dict(map(lambda i,j : (i,j) , organisms, files))

del files[args.reference]

tair_ids = []
for file in files.values():
    orts_df = pd.read_csv(file, sep="\t", usecols=["subject_id"])
    tair_ids.append(orts_df.subject_id.to_list())

tair_genes = list(set.intersection(*map(set,tair_ids)))

orthologs = pd.DataFrame({args.reference: tair_genes})
orthologs.set_index(args.reference, drop=False, inplace=True)

for org,file in files.items():
    tmp = pd.read_csv(file, sep='\t', usecols=['query_id','subject_id'])
    orts_tmp = tmp.drop_duplicates('subject_id').set_index('subject_id', drop=False)
    orts_tmp = orts_tmp.loc[tair_genes, 'query_id']
    orts_tmp = pd.DataFrame({org: orts_tmp})
    orthologs = pd.concat([orthologs, orts_tmp], axis=1)

orthologs.set_axis(['Node'+str(i).zfill(5) for i in range(len(orthologs.index))], inplace=True)

orthologs.to_csv(args.output, sep='\t')

#for organism in x.columns:
#    if organism != args.reference:
#        with open(join("output", organism, "PLTs", "degs.plt.compendium.txt"), "w") as f:
#            f.write("\n".join(x[organism]))