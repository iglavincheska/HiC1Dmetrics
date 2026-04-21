#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys


def resolve_chr_length(genometable, chrom):
    """Resolve chromosome length for both '4' and 'chr4' style genome tables."""
    candidates = [chrom]
    if str(chrom).startswith("chr"):
        candidates.append(str(chrom).replace("chr", "", 1))
    else:
        candidates.append("chr" + str(chrom))

    for key in candidates:
        if key in genometable.index:
            value = genometable.loc[key, 1]
            # If duplicated chromosome labels exist, pick the first row.
            if hasattr(value, "iloc"):
                value = value.iloc[0]
            return int(value)

    raise KeyError(
        "Chromosome '{}' not found in genome table. Tried: {}".format(
            chrom, ", ".join(candidates)
        )
    )

def parse_argv():
    usage = 'Usage: \n    python {} <inputfile> <outputfile> <genometable> <chr> <resolution> [--help]'.format(__file__)
    arguments = sys.argv
    if len(arguments) == 1:
        print(usage)
        exit()
    # ファイル自身を指す最初の引数を除去
    arguments.pop(0)
    # - で始まるoption
    options = [option for option in arguments if option.startswith('-')]

    if '-h' in options or '--help' in options:
        print(usage)

    return arguments

#def getfilename(i, j):
#    return dir + "/" + str(res) + "/chr" + str(i) + "-chr" + str(j) + "/" + ntype + ".matrix.gz"

if __name__ == '__main__':
    arguments = parse_argv()
    inputfile = arguments[0]
    outputfile = arguments[1]
    gtfile = arguments[2]
    chr = arguments[3]
    resolution = int(arguments[4])

    genometable = pd.read_csv(gtfile, delimiter='\t', index_col=[0], header=None)
    chrlen = resolve_chr_length(genometable, chr)
    binlen = int(chrlen / resolution) + 1

    arr = np.zeros((binlen, binlen))
    #    d = pd.read_csv(inputfile, delimiter='\t', header=None, index_col=[0,1])
    d = pd.read_csv(inputfile, delimiter='\t', header=None)
    d = d.set_index([0,1])
    d = d.iloc[:,0]
    d = d.unstack(fill_value=0)
    index = np.arange(binlen) * resolution
    d = d.reindex(index, columns=index, fill_value=0)
    d.index.name = None
    d.columns.name = None

    triu = np.triu(d)
    array = triu + triu.T - np.diag(np.diag(triu))
    df = pd.DataFrame(array, index=d.index, columns=d.columns)

    df.to_csv(outputfile, sep='\t', compression='gzip')
