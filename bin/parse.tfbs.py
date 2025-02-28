###############################################################################
# Description:
#   This script processes a TRANSFAC-format motif file and filters motifs 
#   based on a table of TFs associated with PLTs. It also adjusts gene IDs 
#   according to organism-specific mappings (e.g., removing versions, mapping 
#   older IDs to newer ones). After filtering the motifs, the script writes the 
#   filtered set back to a TRANSFAC file. If a set of reference motifs is 
#   needed (for the reference organism), it merges them into the TF table. 
#   Finally, each filtered motif is converted into a WebLogo (SVG) file.
###############################################################################

import argparse
from Bio import motifs
import subprocess
import pandas as pd
from os.path import join
import numpy as np

# Disable chained assignment warnings in pandas
pd.options.mode.chained_assignment = None

###############################################################################
# Function: run_weblogo
# ---------------------
# Description:
#   Runs the WebLogo command-line tool to produce a motif logo (in SVG format) 
#   from a given motif in TRANSFAC format.
#
# Parameters:
#   input_data : A string representation of the motif in TRANSFAC format
#   id         : A label for the motif that will appear in the logo title
#
# Returns:
#   A string containing the SVG data for the generated WebLogo
###############################################################################
def run_weblogo(input_data, id):
    command = f"weblogo -F svg -c classic -t {id} --title-fontsize 10"
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    output, _ = process.communicate(input_data.encode())
    return output.decode().strip()

###############################################################################
# Parse command-line arguments
###############################################################################
parser = argparse.ArgumentParser(
    description='This script filters a transfac database based on an orthologs table and plots logos.')
parser.add_argument('-t', '--tfs', required=True, type=str, 
                    help='File path to the table of TFs related to PLTs.')
parser.add_argument('-m', '--transfac', required=True, type=str, 
                    help='File path to the TRANSFAC input file (motifs).')
parser.add_argument('-p', '--organism', required=True, type=str, 
                    help='Organism name (with underscores if needed).')
parser.add_argument('-r', '--ref', required=True, type=str, 
                    help='Reference organism name.')
parser.add_argument('-i', '--tfsRef', required=True, type=str, 
                    help='Reference organism cisBP file.')
parser.add_argument('-z', '--zeaids', required=True, type=str, 
                    help='MaizeMine gene_to_alias_ids file for Zea mays.')
parser.add_argument('-u', '--rapmsu', required=True, type=str, 
                    help='RAP-MSU file (ID-conversion table) for Oryza.')
parser.add_argument('-o', '--output', required=True, type=str, 
                    help='Path to output the filtered TRANSFAC file.')
parser.add_argument('-d', '--logos', required=True, type=str, 
                    help='Output directory to save motif logos (SVG).')

args = parser.parse_args()

###############################################################################
# Step 1: Read the input TRANSFAC file with Bio.motifs
###############################################################################
with open(args.transfac, 'r') as transfacFile:
    motifs_transfac = motifs.parse(transfacFile, 'TRANSFAC', strict=False)

###############################################################################
# Step 2: Read the table of TFs associated with PLTs (related_tfs)
###############################################################################
related_tfs = pd.read_csv(args.tfs, sep='\t')

###############################################################################
# Step 3: Organism-specific gene ID transformations
###############################################################################
# If the organism is Solanum lycopersicum, strip version suffix from each ID
if args.organism == 'Solanum_lycopersicum':
    for i, motif_obj in enumerate(motifs_transfac):
        ids = motif_obj['ID'].split(',')
        ids = [id.split('.', maxsplit=1)[0] for id in ids]
        motifs_transfac[i]['ID'] = ','.join(ids)

# If the organism is Glycine max, convert 'GLYMA' prefixes to 'Glyma.'
if args.organism == 'Glycine_max':
    for i, motif_obj in enumerate(motifs_transfac):
        ids = motif_obj['ID'].split(',')
        ids = [id.replace('GLYMA', 'Glyma.') for id in ids]
        motifs_transfac[i]['ID'] = ','.join(ids)

# If the organism is Zea mays, map older IDs to the new ones using zeaids
if args.organism == 'Zea_mays':
    zea_map = pd.read_csv(args.zeaids, sep='\t', header=None, comment='#',
                          names=['newID', 'newGenome', 'olderID', 'olderGenome'])
    for i, motif_obj in enumerate(motifs_transfac):
        ids = motif_obj['ID'].split(',')
        genes = []
        for gen in ids:
            if gen in zea_map.olderID.values:
                genes.append(zea_map.newID[zea_map.olderID == gen].values[0])
        motifs_transfac[i]['ID'] = ','.join(genes)

# If the organism is Oryza sativa Japonica Group, map IDs using RAP-MSU file
if args.organism == 'Oryza_sativa_Japonica_Group':
    rap_msu = pd.read_csv(args.rapmsu, sep='\t', header=None, names=['orig', 'mapped'], comment='#')
    # Truncate version suffix from mapped IDs
    mapped = []
    for mapped_id in rap_msu.mapped:
        mapped_id = mapped_id.split('.', maxsplit=1)[0]
        mapped.append(mapped_id)
    rap_msu.mapped = mapped
    # Convert original IDs to uppercase for matching
    rap_msu.orig = rap_msu.orig.str.upper()
    
    for i, motif_obj in enumerate(motifs_transfac):
        ids = motif_obj['ID'].split(',')
        genes = []
        for gen in ids:
            if gen in rap_msu.orig.values:
                genes.append(rap_msu.loc[rap_msu.orig == gen, 'mapped'].values[0])
        motifs_transfac[i]['ID'] = ','.join(genes)

###############################################################################
# Step 4: Filter motifs to keep only those that have at least one ID matching 
#         the relevant TFs
###############################################################################
motifs_filtered = []

# If the current organism is NOT the reference
if args.organism != args.ref:
    for motif_obj in motifs_transfac:
        ids = motif_obj['ID'].split(',')
        # Keep if any ID is in the 'GeneID' column of related_tfs
        if any([True for id in ids if id in related_tfs.GeneID.values]):
            motifs_filtered.append(motif_obj)

# If the current organism IS the reference
else:
    for motif_obj in motifs_transfac:
        ids = motif_obj['ID'].split(',')
        # Keep if any ID is in the 'Ref' column of related_tfs
        if any([True for id in ids if id in related_tfs.Ref.values]):
            motifs_filtered.append(motif_obj)

###############################################################################
# Step 5: Write filtered motifs back to a TRANSFAC file
###############################################################################
with open(args.output, 'w') as transfacFile:
    transfacFile.write(motifs.write(motifs_filtered, 'transfac'))

###############################################################################
# Step 6: If 'refMotif' isn't in the related_tfs table yet, annotate reference 
#         motifs by reading the reference's cisBP file.
###############################################################################
if 'refMotif' not in related_tfs:
    with open(args.tfsRef, 'r') as transfacFile:
        ref_motifs = motifs.parse(transfacFile, 'TRANSFAC', strict=False)
    
    motifs_list = []
    # For each reference gene in related_tfs, find matching motif AC from ref_motifs
    for tf in related_tfs.Ref.values:
        motif_ids = [motif_obj['AC'] for motif_obj in ref_motifs if tf in motif_obj['ID'].split(',')]
        if len(motif_ids) > 0:
            motif_ids = ','.join(motif_ids)
            motifs_list.append(motif_ids)
        else:
            motifs_list.append(None)
    
    related_tfs['refMotif'] = motifs_list

###############################################################################
# Step 7: If this is not the reference organism, annotate 'relatedMotif' 
#         for each gene if not already present
###############################################################################
if args.organism != args.ref:
    if 'relatedMotif' not in related_tfs:
        related_tfs['relatedMotif'] = np.repeat(None, len(related_tfs.index))

    # Build a dictionary from motif AC -> list of gene IDs
    ids_dict = {}
    for motif_obj in motifs_filtered:
        ids = motif_obj['ID'].split(',')
        ids_dict[motif_obj['AC']] = ids

    # For each gene in related_tfs, find all motifs that match
    for i, gen in enumerate(related_tfs.GeneID.values):
        motifs_related = [ac for ac, idlist in ids_dict.items() if gen in idlist]
        if len(motifs_related) > 0:
            related_tfs.relatedMotif[i] = ','.join(motifs_related)

###############################################################################
# Step 8: Write the updated related_tfs table back to the same file
###############################################################################
related_tfs.to_csv(args.tfs, sep='\t', index=False)

###############################################################################
# Step 9: For each filtered motif, run WebLogo to generate an SVG logo
###############################################################################
for motif_obj in motifs_filtered:
    # Convert motif to TRANSFAC format and pass into run_weblogo
    motif_logo = run_weblogo(motif_obj.format('transfac'), motif_obj['AC'])
    # Save the resulting SVG file in the specified logos folder
    with open(join(args.logos, (motif_obj['AC']+'.svg')), 'w') as f:
        f.write(motif_logo)
