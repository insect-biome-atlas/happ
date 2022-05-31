#!/usr/bin/env python

def sum_counts(f, ids=[]):
    counts = {}
    with open(f, 'r') as fhin:
        for i, line in enumerate(fhin):
            if i == 0:
                continue
            asv = line.rsplit()[0]
            if len(ids)>0:
                if asv not in ids:
                    continue
            summed_count = sum([int(x) for x in line.rsplit()[1:]])
            counts[asv] = summed_count
    return counts