#!/usr/bin/python
import sys
import re

fmt=sys.argv[1]

keep = set()
pattern = re.compile('score=[+-][0-9]+\.[0-9]+,(\w+)|')
with open(fmt) as inh:
    for line in inh:
        cols = line.strip().split("\t")
        if len(cols)<3:
            continue
        name = cols[0]
        attrb = cols[3]

        score = attrb.split(",")

        if len(score) > 2:
            print(name)
