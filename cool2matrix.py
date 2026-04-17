from optparse import OptionParser
import cooler
import pandas as pd

desc="Transfer cool to hic."

parser = OptionParser(description=desc)

parser.add_option("-i", "--input", action="store", type="string",
                  dest="input", help="Input cool or mcool file.", metavar="<file>")

parser.add_option("-r", "--resolution", action="store", type="int",
                  dest="res", help="Resolution for the output hic file.", metavar="<int>")

parser.add_option("-o", "--output", action="store", type="string",
                  dest="output", help="Output matrix file.", metavar="<file>")

(opt, args) = parser.parse_args()
file =  opt.input
resolution = opt.res
out = opt.output

print("Input file: %s" % file)
print("Resolution: %s" % resolution)

#==========================

file_type = file.split('.')[len(file.split('.'))-1]
if file_type == 'cool':
    c = cooler.Cooler(file)
if file_type == 'mcool':
    c = cooler.Cooler(file+'::resolutions/'+str(resolution))

chrom = c.chroms()[0:len(c.chroms())]
chrom = chrom[~chrom.name.str.contains('_|M')]
chrom = chrom.sort_values('name')

import os

if os.path.exists(out):
    os.remove(out)

sel = c.matrix(balance=False, as_pixels=True, join=True, sparse=False)
chroms = list(chrom.name)

for a_idx, chrA in enumerate(chroms):
    for chrB in chroms[a_idx:]:
        # chrA==chrB gives cis; chrA!=chrB gives trans
        hic_df = sel.fetch(chrA, chrB)

        if hic_df.empty:
            continue

        # keep only needed columns
        hic_df = hic_df.loc[:, ["chrom1", "start1", "chrom2", "start2", "count"]].copy()

        # add dummy juicer fields
        hic_df["str1"] = 0
        hic_df["frag1"] = 0
        hic_df["str2"] = 0
        hic_df["frag2"] = 1

        # order for juicer short-with-score style
        hic_df = hic_df[["str1", "chrom1", "start1", "frag1",
                         "str2", "chrom2", "start2", "frag2", "count"]]

        hic_df.to_csv(out, sep="\t", mode="a", index=False, header=False)
    
#===========================
