###############################################################################
# Description:
#   This script merges multiple BLASTp (or ortholog) tables so that all 
#   organisms share a common reference organism. The intersection of subject IDs
#   from each BLASTp table (excluding the reference) is taken to find genes 
#   present across all organisms, and those genes are mapped to their ortholog 
#   IDs in each other organism. The final output is a table where each row 
#   corresponds to a reference gene and each column corresponds to an organism's 
#   ortholog ID.
###############################################################################

import pandas as pd
import argparse
from os.path import join

###############################################################################
# Set up argument parsing for command-line usage
###############################################################################
parser = argparse.ArgumentParser(description='This script joins (merges) multiple BLASTp tables on a common set of genes.')
parser.add_argument('-o', '--organisms', required=True, type=str, 
                    help='Comma-separated list of organisms.')
parser.add_argument('-f', '--files', required=True, type=str, 
                    help='Comma-separated list of BLASTp/orthologs table file paths.')
parser.add_argument('-r', '--reference', required=True, type=str, 
                    help='Reference organism name.')
parser.add_argument('-u', '--output', required=True, type=str, 
                    help='Path to the output file.')

args = parser.parse_args()

###############################################################################
# Parse the list of organisms and matching file paths
###############################################################################
files = args.files.split(",")         # e.g. ["path1.tsv","path2.tsv","path3.tsv"]
organisms = args.organisms.split(",") # e.g. ["Org1","Org2","Org3"]

# Create a dictionary mapping organism -> file
# e.g. { "Org1": "path1.tsv", "Org2": "path2.tsv", "Org3": "path3.tsv" }
files = dict(map(lambda i,j : (i,j), organisms, files))

# Remove the reference organism from the dictionary so we only have other organisms
del files[args.reference]

###############################################################################
# Gather subject (reference) gene IDs from each file, then find the intersection
###############################################################################
tair_ids = []
for file in files.values():
    # Read only the "subject_id" column, which corresponds to reference genes
    orts_df = pd.read_csv(file, sep="\t", usecols=["subject_id"])
    tair_ids.append(orts_df.subject_id.to_list())

# Intersect the set of reference genes across all non-reference organisms
tair_genes = list(set.intersection(*map(set, tair_ids)))

###############################################################################
# Initialize a DataFrame using the reference genes as rows
###############################################################################
orthologs = pd.DataFrame({args.reference: tair_genes})
orthologs.set_index(args.reference, drop=False, inplace=True)

###############################################################################
# For each other organism, look up ortholog IDs that map to those reference genes
###############################################################################
for org, file in files.items():
    # Each file has a "query_id" (organism gene) and "subject_id" (reference gene)
    tmp = pd.read_csv(file, sep='\t', usecols=['query_id', 'subject_id'])
    # Remove duplicates on subject_id, then index by subject_id
    orts_tmp = tmp.drop_duplicates('subject_id').set_index('subject_id', drop=False)
    
    # Filter only those subject (reference) IDs in the intersection
    orts_tmp = orts_tmp.loc[tair_genes, 'query_id']
    
    # Convert this Series to a DataFrame with a single column named after the org
    orts_tmp = pd.DataFrame({org: orts_tmp})
    
    # Merge into the main orthologs table
    orthologs = pd.concat([orthologs, orts_tmp], axis=1)

###############################################################################
# Rename the rows to "Node000xx" style for convenience (optional step)
###############################################################################
orthologs.set_axis(['Node'+str(i).zfill(5) for i in range(len(orthologs.index))], inplace=True)

###############################################################################
# Write the final merged table to the specified output file
###############################################################################
orthologs.to_csv(args.output, sep='\t')

# For example usage, once you have "orthologs", you could loop over columns and 
# extract ortholog sets for further processing, but that code is commented out 
# here as it was in the original script:

# for organism in x.columns:
#     if organism != args.reference:
#         with open(join("output", organism, "PLTs", "degs.plt.compendium.txt"), "w") as f:
#             f.write("\n".join(x[organism]))
