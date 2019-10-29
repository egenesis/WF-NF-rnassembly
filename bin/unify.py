#!/usr/bin/env python
# gffread -M --cluster-only -T  -F --keep-genes --keep-exon-attrs
from __future__ import print_function
from collections import Counter
import pandas as pd 

import gffutils

import sys
import re

import logging
import argparse
import progressbar

parser = argparse.ArgumentParser(description='Process GTF from different sources.')
parser.add_argument('-d', action='store_true', help='turn on debugging messages')
parser.add_argument('--gtf', help='GTF file', required=True)
args = parser.parse_args()

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("unify")
logger.setLevel(logging.INFO)

if args.d:
    logger.setLevel(logging.DEBUG)
    logger.debug("Debug is on.")

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
            logger.debug("### The exon numbers for each tx is: %s" % [len(es), len(es_kept)])
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
            logger.debug("### the same is %s" % thesame)
            if thesame > 0:
                logger.debug("### skip this tx")
                return cleaned
    if thesame == 0:
        logger.debug("### the same is %s, then add %s" % (thesame, tx["source"]) )
        cleaned.append(tx)
    return cleaned

def calc_size(exons):
    return sum([abs(e.start - e.end) for e in exons])

def solve_structure(fs, chr, start, end, strand, rloc):
    if not fs:
        return None
    attrb = "; ".join(["locus \"RLOC_%s\"" % rloc,
                       "gene_id \"RLOC_%s\"" % rloc,
                       "gene_name \"RLOC_%s\"" % rloc,
                       "gene_biotype \"predicted\""]) + ";"
    line = "\t".join([chr, "egenesis", "gene", str(start), str(end), ".", strand ,".", attrb])
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
    gene_name = None
    in_ensembl = False
    cleaned = list()
    total = 0
    stats = dict()
    for f in db.iter_by_parent_childs("transcript", order_by=("transcript_id", "start")):
        ne = len(f[1:])
        total += 1
        size = sum([abs(e.end-e.start) for e in f[1:]])
        if 'gene_name' in f[0].attributes:
            possible_genes[f[0].attributes['gene_name'][0]] += 1
        if 'ref_gene_id' in f[0].attributes:
            ensembl_id = [f[0].attributes['ref_gene_id'][0]]
        
        stats[f[0].id] = {'locus': rloc, 'tx': f[0].id, 'source': f[0].source, 'nexons': len(f[1:]), 'size': calc_size(f[1:]), 'added': 0, 'gene_name': 'na'}

        ctx = {'source': f[0], 'size': size, 'ne': ne, 'exons': f[2:], 'id': f[0].id}
        
        if ctx["source"].source == "ensembl":
            logger.debug("## This tx %s comes from ensembl, added" % f[0].id)
            cleaned.append(ctx)
            in_ensembl = True
        elif ctx["source"].source == "StringTie" and ctx["source"].attributes["transcript_id"][0].find("ENS") > -1:
            logger.debug("## This tx %s comes from stringTie and ensembl, added" % f[0].id)
            cleaned.append(ctx)
            # stats[f[0].id]['added'] = 1
            in_ensembl = True
        else:
            logger.debug("## This tx %s is novel, to be analyzed" % f[0].id)
            txs[f[0].attributes['transcript_id'][0]] = ctx
    
    if in_ensembl:
        logger.debug("## This gene is in ensembl, so skipping")
        return [tx for tx in list(stats.values())]

    if possible_genes:
        gene_name = possible_genes.most_common(1)[0][0]

    txs_sorted = sorted(txs, key=lambda x: (txs[x]['size'], txs[x]['ne']), reverse=True)
    for tx in txs_sorted:
        logger.debug("## tx being analyzed: %s" % tx)
        ctx = txs[tx]
        cleaned = clean_tx(ctx, cleaned)
        logger.debug("## how many in cleaned %s out of %s" % (len(cleaned), total))

    if gene_name:
        gene.attributes['gene_name'] =  gene_name
    if ensembl_id:
        gene.attributes['gene_id'] = ensembl_id
    print(gene)
    for tx in cleaned:           
        if tx["id"] in stats:
            stats[tx["id"]]['added'] = 1
            if gene_name:
                stats[tx["id"]]['gene_name'] = gene_name
        for e in db.children(tx["id"], order_by="start"):
            if gene_name:
                e.attributes['gene_name'] = gene_name
            if ensembl_id:
                e.attributes['gene_id'] = ensembl_id
            e.attributes['gene_biotype'] = "predicted"
            e.attributes['locus'] = "RLOC_%s" % rloc
            print(e)
    return [tx for tx in list(stats.values())]

read = 0
current  = ""
rloc = ""
chr = ""
features = []
start = 0
end = 0
stats = list()
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
with open(args.gtf) as inh:
    import datetime
    today = datetime.date.today()
    print("#!eGenesis: %s" % today)
    for line in inh:
        read += 1
        bar.update(read)
        if line.startswith("#"):
            continue
        if line.find("RLOC_") > -1:
            # TODO when life is better and re.compile(r'RLOC_\d+') works this will be nicer
            rloc = line.split("RLOC_")[1][:-3]
        f = gffutils.feature.feature_from_line(line)
        # logger.debug(f)
        if current != rloc:
            # resolve gene
            if current: 
                stats.extend(solve_structure(features, chr, start, end, strand, rloc))
            # when is done:
            current = rloc
            chr = f.seqid
            start = f.start
            end = f.end
            strand = f.strand
            # add gene to database?
            features = [f]
        else:
            if f.start < start:
                start = f.start
            if f.end > end:
                end = f.end
            features.append(f)
    stats.extend(solve_structure(features, chr, start, end, strand, rloc))

df = pd.DataFrame(stats) 
              
df.to_csv("unify.stats.txt")