import argparse
from Bio import motifs
import subprocess
import pandas as pd
from os.path import join
import numpy as np

pd.options.mode.chained_assignment = None

def run_weblogo(input_data, id):
    command = f"weblogo -F svg -c classic -t {id} --title-fontsize 10"
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    output, _ = process.communicate(input_data.encode())
    return output.decode().strip()

parser=argparse.ArgumentParser(description= 'This script filters a transfac database based on a orthologs table and plot the logos.')

parser.add_argument('-t', '--tfs', required=True, type=str, help='TFs related to PLTs.')
parser.add_argument('-m', '--transfac', required=True, type=str, help='Transfac input file.')
parser.add_argument('-p', '--organism', required=True, type=str, help='Organism name (separated by underscore).')
parser.add_argument('-r', '--ref', required=True, type=str, help='Reference organism.')
parser.add_argument('-i', '--tfsRef', required=True, type=str, help='Reference organism cisBP file.')
parser.add_argument('-z', '--zeaids', required=True, type=str, help='gene_to_alias_ids file from MaizeMine.')
parser.add_argument('-u', '--rapmsu', required=True, type=str, help='RAP-MSU file (ID-conversion table).')
parser.add_argument('-o', '--output', required=True, type=str, help='Transfac output file.')
parser.add_argument('-d', '--logos', required=True, type=str, help='Output dir to save logos.')

args=parser.parse_args()

with open(args.transfac, 'r') as transfacFile:
    motifs_transfac = motifs.parse(transfacFile, 'TRANSFAC', strict=False)

related_tfs = pd.read_csv(args.tfs, sep='\t')

if args.organism == 'Solanum lycopersicum':
    for i,motif in enumerate(motifs_transfac):
        ids = motif['ID'].split(',')
        ids = [id.split('.', maxsplit=1)[0] for id in ids]
        motifs_transfac[i]['ID'] = ','.join(ids)

if args.organism == 'Glycine max':
    for i,motif in enumerate(motifs_transfac):
        ids = motif['ID'].split(',')
        ids = [id.replace('GLYMA', 'Glyma.') for id in ids]
        motifs_transfac[i]['ID'] = ','.join(ids)

if args.organism == 'Zea mays':
    zea_map = pd.read_csv(args.zeaids, sep='\t', header=None, comment='#',
                      names=['newID', 'newGenome', 'olderID', 'olderGenome'])
    for i,motif in enumerate(motifs_transfac):
        ids = motif['ID'].split(',')
        genes = []
        for gen in ids:
            if gen in zea_map.olderID.values:
                genes.append(zea_map.newID[zea_map.olderID == gen].values[0])
        motifs_transfac[i]['ID'] = ','.join(genes)

if args.organism == 'Oryza sativa Japonica Group':
    rap_msu = pd.read_csv(args.rapmsu, sep='\t', header=None, names=['orig', 'mapped'], comment='#')
    mapped = []
    for mapped_id in rap_msu.mapped:
        mapped_id = mapped_id.split('.', maxsplit=1)[0]
        mapped.append(mapped_id)
    rap_msu.mapped = mapped
    rap_msu.orig = rap_msu.orig.str.upper()
    for i,motif in enumerate(motifs_transfac):
        ids = motif['ID'].split(',')
        genes = []
        for gen in ids:
            if gen in rap_msu.orig.values:
                genes.append(rap_msu.loc[rap_msu.orig == gen, 'mapped'].values[0])
        motifs_transfac[i]['ID'] = ','.join(genes)

motifs_filtered = []
if args.organism != args.ref:
    for motif in motifs_transfac:
        ids = motif['ID'].split(',')
        if any([True for id in ids if id in related_tfs.GeneID.values]):
            motifs_filtered.append(motif)
if args.organism == args.ref:
    for motif in motifs_transfac:
        ids = motif['ID'].split(',')
        if any([True for id in ids if id in related_tfs.Ref.values]):
            motifs_filtered.append(motif)

with open(args.output, 'w') as transfacFile:
    transfacFile.write(motifs.write(motifs_filtered, 'transfac'))

if 'refMotif' not in related_tfs:
    with open(args.tfsRef, 'r') as transfacFile:
        ref_motifs = motifs.parse(transfacFile, 'TRANSFAC', strict=False)
    motifs_list = []
    for tf in related_tfs.Ref.values:
        motif_ids = [motif['AC'] for motif in ref_motifs if tf in motif['ID'].split(',')]
        if len(motif_ids) > 0:
            motif_ids = ','.join(motif_ids)
            motifs_list.append(motif_ids)
        else:
            motifs_list.append(None)
    related_tfs['refMotif'] = motifs_list

if args.organism != args.ref:
    if 'relatedMotif' not in related_tfs:
        related_tfs['relatedMotif'] = np.repeat(None, len(related_tfs.index))

    ids_dict = {}
    for motif in motifs_filtered:
        ids = motif['ID'].split(',')
        ids_dict[motif['AC']] = ids

    for i,gen in enumerate(related_tfs.GeneID.values):
        motifs_related = [ac for ac,idlist in ids_dict.items() if gen in idlist]
        if len(motifs_related) > 0:
            related_tfs.relatedMotif[i] = ','.join(motifs_related)

related_tfs.to_csv(args.tfs, sep='\t', index=False)

for motif in motifs_filtered:
    motif_logo = run_weblogo(motif.format('transfac'), motif['AC'])
    with open(join(args.logos, (motif['AC']+'.svg')), 'w') as f:
        f.write(motif_logo)
