#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys


def _norm_chr(value):
    text = str(value).strip()
    return text.lower().replace("chr", "", 1) if text.lower().startswith("chr") else text.lower()


def resolve_chr_length(genometable, chrom):
    """Resolve chromosome length with robust matching for varied chromosome naming."""
    chrom_text = str(chrom).strip()
    candidates = [chrom_text]
    if chrom_text.lower().startswith("chr"):
        candidates.append(chrom_text[3:])
    else:
        candidates.append("chr" + chrom_text)

    # 1) Exact string match first.
    for key in candidates:
        if key in genometable.index:
            value = genometable.loc[key, 1]
            if hasattr(value, "iloc"):
                value = value.iloc[0]
            return int(float(value))

    # 2) Case-insensitive + optional 'chr' prefix normalization.
    norm_target = _norm_chr(chrom_text)
    for idx_value in genometable.index:
        if _norm_chr(idx_value) == norm_target:
            value = genometable.loc[idx_value, 1]
            if hasattr(value, "iloc"):
                value = value.iloc[0]
            return int(float(value))

    raise KeyError(
        "Chromosome '{}' not found in genome table. Tried: {}".format(
            chrom_text, ", ".join(candidates)
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

    genometable = pd.read_csv(gtfile, delimiter='\t', index_col=[0], header=None, dtype={0: str})
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
