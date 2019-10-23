#!/usr/bin/env python

import sys

fmt=sys.argv[1]

keep = set()
with open(fmt) as inh:
    for line in inh:
        cols = line.strip().split("\t")
        name = cols[0]
        size = int(name.split("_")[3])
        hit = int(cols[3])
        if size * 1.0 > hit * 0.9:
            keep.add(name)

for name in keep:
    print(name)