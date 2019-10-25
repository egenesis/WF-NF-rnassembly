#!/usr/bin/env python
# gffread -M --cluster-only -T  -F --keep-genes --keep-exon-attrs
from __future__ import print_function
from collections import Counter

import gffutils

import sys
import re

if len(sys.argv) < 2:
    print("Usage: gfftools.py full.gtf")
    sys.exit(0)


def fix(number):
    if number > 100:
        return round(number/100)
    elif number > 10:
        return round(number/10)
    return number

def clean_tx(tx, cleaned):
    thesame = 0
    if not cleaned:
        thesame = 0
    else:
        for tx_kept in cleaned:
            # compare exons coverage
            es = [ "%s-%s" % (fix(e.start), fix(e.end))for e in tx["exons"] ]
            es_kept = [ "%s-%s" % (fix(e.start), fix(e.end)) for e in tx_kept["exons"] ]
            # print([len(es), len(es_kept)])
            if len(es) == len(es_kept):
                for e_pair in es:
                    # check start and end in second digit
                    # print([e_pair, es_kept])
                    if e_pair in es_kept:
                        thesame = 1
                    else:
                        thesame = 0
                        break
            else:
               thesame = 0
            # detected being the same in this iteration,then stop   
            # print("the same is %s" % thesame)
            if thesame > 0:
                # print("skip")
                return cleaned
    if thesame == 0:
        # print("the same is %s, then add %s" % (thesame, tx["source"]) )
        cleaned.append(tx)
    return cleaned

def solve_structure(fs, chr, start, end, rloc):
    if not fs:
        return None
    attrb = "; ".join(["locus \"RLOC_%s\"" % rloc]) + ";"
    line = "\t".join([chr, "guessed", "gene", str(start), ".", str(end), ".", attrb])
    gene = gffutils.feature.feature_from_line(line)

    fs.append(gene)
    db = gffutils.create_db(fs, ':memory:',
                        disable_infer_transcripts=True,
                        disable_infer_genes=True,
                        id_spec={'gene': 'locus', 'transcript': 'transcript_id' },
                        gtf_transcript_key='transcript_id', gtf_gene_key='locus')
    txs = dict()
    possible_genes = Counter()
    ensembl_id = None
    cleaned = list()
    total = 0
    for f in db.iter_by_parent_childs("transcript", order_by=("transcript_id", "start")):
        ne = len(f[1:])
        total += 1
        size = sum([abs(e.end-e.start) for e in f[1:]])
        if 'gene_name' in f[0].attributes:
            possible_genes[f[0].attributes['gene_name'][0]] += 1
        if 'ref_gene_id' in f[0].attributes:
            ensembl_id = [f[0].attributes['ref_gene_id'][0]]
        ctx = {'source': f[0], 'size': size, 'ne': ne, 'exons': f[2:], 'id': f[0].id}

        if ctx["source"].source == "ensembl":
            cleaned.append(ctx)
        elif ctx["source"].source == "StringTie" and ctx["source"].attributes["transcript_id"][0].find("ENS") > -1:
            cleaned.append(ctx)
        else:
            txs[f[0].attributes['transcript_id'][0]] = ctx
    
    gene_name = possible_genes.most_common(1)[0][0]
    txs_sorted = sorted(txs, key=lambda x: (txs[x]['size'], txs[x]['ne']), reverse=True)
    for tx in txs_sorted:
        # print(tx)
        ctx = txs[tx]
        cleaned=clean_tx(ctx, cleaned)
        # print("how many in cleaned %s out of %s" % (len(cleaned), total))

    # print(gene)
    for tx in cleaned:
        # if ensembl_id:
        #     tx["source"].attributes['gene_id'] = ensembl_id   
        # tx["source"].attributes['gene_name'] = gene_name
        # print(tx["source"])
        # print(tx["id"])
        for e in db.children(tx["id"], order_by="start"):
            # print(e)
        # for e in tx["exons"]:
            e.attributes['gene_name'] = gene_name
            if ensembl_id:
                e.attributes['gene_id'] = ensembl_id
            e.attributes['locus'] = "RLOC_%s" % rloc
            print(e)

current  = ""
rloc = ""
chr = ""
features = []
start = 0
end = 0
with open(sys.argv[1]) as inh:
    for line in inh:
        if line.find("RLOC_") > -1:
            # TODO when life is better and re.compile(r'RLOC_\d+') works this will be nicer
            rloc = line.split("RLOC_")[1][:-3]
        f = gffutils.feature.feature_from_line(line)
        if current != rloc:
            # resolve gene
            solve_structure(features, chr, start, end, rloc)
            # when is done:
            current = rloc
            chr = f.seqid
            start = f.start
            end = f.end
            # add gene to database?
            features = [f]
        else:
            if f.start < start:
                start = f.start
            if f.end > end:
                end = f.end
            features.append(f)
    solve_structure(features, chr, start, end, rloc)
                

sys.exit(0)
