#!/usr/bin/env python

from __future__ import print_function

import gffutils

import sys
import re

if len(sys.argv) < 3:
    print("Usage: gfftools.py pasa_gff transdecoder_predict_bed protein_db_fasta")
    sys.exit(0)

prot=sys.argv[3]

genes = dict()
with open(prot) as inh:
    for line in inh:
        if line.startswith(">"):
            try:
                gene = line.split("GN=")[1].split(" ")[0]
            except:
                gene = line.split("|")[2].split(" ")[0]
            pid = line.split("|")[1]
            genes[pid] = gene

bed=sys.argv[2]

mapper = dict()
with open(bed) as inh:
    for line in inh:
        cols = line.strip().split("\t")
        if len(cols)<3:
            continue
        name = cols[0]
        attrb = cols[3]

        score = attrb.split(",")
        if len(score) > 2:
            mapper[attrb.split(";")[0][3:].split(".p")[0]]=score[2].split("|")[0]

def transform(f):
    f.featuretype = "exon"
    _id = f.attributes
    # f.attributes['exon_id'] = _id['ID'][0]
    f.attributes['gene_id'] = ["G" + _id['Target'][0].split("_")[1]]
    f.attributes['transcript_id'] = ["T" + _id['Target'][0].split("_")[1]]
    f.attributes.pop("ID", None)
    return f

db = gffutils.create_db(sys.argv[1], ':memory:',
                        disable_infer_transcripts=False,
                        disable_infer_genes=False,
                        id_spec={'gene': 'gene_id', 'transcript': 'transcript_id' },
                        gtf_transcript_key='transcript_id', gtf_gene_key='gene_id', transform=transform)

gene_cache = ""
updated = list()
for f in db.all_features(order_by=('seqid', 'attributes', 'start')):
    node = f.attributes['Target'][0].split(" ")[0]
    f.attributes['Target'] = node
    if gene_cache != f.attributes['gene_id'][0]:
        gene_cache = f.attributes['gene_id'][0]
        exon_idx = 1
    if node in mapper:
        f.attributes['protein_id'] = mapper[node]
    f.attributes['gene_name'] = node
    if mapper[node] in genes:
        f.attributes['gene_name'] = genes[mapper[node]]
    f.attributes['exon_id'] = "E" + f.attributes['gene_id'][0] + "." + str(exon_idx)
    f.attributes['exon_number'] = [str(exon_idx)]
    exon_idx += 1
    updated.append(f)

dbu = gffutils.create_db(updated, ':memory:',
                        disable_infer_transcripts=False,
                        disable_infer_genes=False,
                        id_spec={'gene': 'gene_id', 'transcript': 'transcript_id' },
                        gtf_transcript_key='transcript_id', gtf_gene_key='gene_id')

with open("tmp.gff", 'w') as outh:
    for f in dbu.all_features(order_by=('seqid', 'attributes', 'start')):
        print(f, file=outh)   

fixed = list()
with open("tmp.gff") as inh:
    for line in inh:
        cols = line.split("\t")
        cols[8] = cols[8].replace("=", " ")
        cols[8] = cols[8].replace(";", "; ").strip()
        fixed.append("\t".join(cols))

dbu = gffutils.create_db("\n".join(fixed), ':memory:',
                        from_string = True,
                        disable_infer_transcripts=False,
                        disable_infer_genes=False,
                        id_spec={'gene': 'gene_id', 'transcript': 'transcript_id' },
                        gtf_transcript_key='transcript_id', gtf_gene_key='gene_id')

for f in dbu.all_features(order_by=('seqid', 'attributes', 'start')):
    f.dialect['quoted GFF2 values'] = True
    print(f)