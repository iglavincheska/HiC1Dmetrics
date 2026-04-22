import pandas as pd
import numpy as np
import re

ipo_if_file = "/Users/iglavinchesk/WP2/HiC1Dmetrics/IPO323_WT_IFregions_merged.bed"
ir_if_file = "/Users/iglavinchesk/WP2/HiC1Dmetrics/IR0126B_WT_IFregions_merged.bed"
ipo_gene_file = "/Users/iglavinchesk/WP1/FINAL/data/z.tritici.IP0323.reannot.gff3"
ir_gene_file = "/Users/iglavinchesk/WP2/data/00_genomes/IR0126B.genes.gff3"
ortho_file = "/Users/iglavinchesk/WP1/FINAL/data/orthoproteins.tsv"
out_prefix = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts"

CHR_MIN = 1
CHR_MAX = 13


def chrom_num_from_name(chrom):
    s = str(chrom).strip().lower()
    m = re.search(r"(?:^|[._-])chr(?:omosome)?[._-]?(\d+)$", s)
    if m:
        return int(m.group(1))
    m = re.search(r"^(?:chr(?:omosome)?)?[._-]?(\d+)$", s)
    if m:
        return int(m.group(1))
    return np.nan


def normalize_chr(chrom):
    n = chrom_num_from_name(chrom)
    if pd.isna(n):
        return str(chrom).strip()
    return str(int(n))


def load_ifregions(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["chromosome", "start", "end"])
    df["chromosome"] = df["chromosome"].astype(str).map(normalize_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["chr_num"] = df["chromosome"].map(chrom_num_from_name)
    df = df[df["chr_num"].between(CHR_MIN, CHR_MAX)].drop(columns=["chr_num"])
    return df.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)


def parse_attrs(s):
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def load_genes_ipo_reannot(path):
    cols = ["chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, comment="#")
    g = df[df["type"] == "gene"].copy()
    attrs = g["attributes"].map(parse_attrs)
    g["gene_id"] = attrs.map(lambda x: x.get("ID") or x.get("Name") or x.get("locus_tag") or "")
    g["chromosome"] = g["chromosome"].astype(str).map(normalize_chr)
    g["start"] = pd.to_numeric(g["start"], errors="coerce")
    g["end"] = pd.to_numeric(g["end"], errors="coerce")
    g = g.dropna(subset=["start", "end"])
    g["chr_num"] = g["chromosome"].map(chrom_num_from_name)
    g = g[g["chr_num"].between(CHR_MIN, CHR_MAX)].drop(columns=["chr_num"])
    g = g[g["gene_id"] != ""].drop_duplicates(["gene_id", "chromosome", "start", "end"])
    return g[["gene_id", "chromosome", "start", "end"]].sort_values(["chromosome", "start", "end"])


def load_genes_ir(path):
    cols = ["chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, comment="#")
    g = df[df["type"] == "mRNA"].copy()
    attrs = g["attributes"].map(parse_attrs)
    g["gene_id"] = attrs.map(lambda x: x.get("geneID") or x.get("Parent") or x.get("ID") or "")
    g["chromosome"] = g["chromosome"].astype(str).map(normalize_chr)
    g["start"] = pd.to_numeric(g["start"], errors="coerce")
    g["end"] = pd.to_numeric(g["end"], errors="coerce")
    g = g.dropna(subset=["start", "end"])
    g["chr_num"] = g["chromosome"].map(chrom_num_from_name)
    g = g[g["chr_num"].between(CHR_MIN, CHR_MAX)].drop(columns=["chr_num"])
    g = g[g["gene_id"] != ""].drop_duplicates(["gene_id", "chromosome", "start", "end"])
    return g[["gene_id", "chromosome", "start", "end"]].sort_values(["chromosome", "start", "end"])


def overlap_genes_in_if(genes, ifr):
    out = []
    if_by_chr = {c: df[["start", "end"]].to_numpy() for c, df in ifr.groupby("chromosome")}
    for row in genes.itertuples(index=False):
        regs = if_by_chr.get(row.chromosome)
        if regs is None:
            continue
        gs, ge = row.start, row.end
        hit = False
        for rs, re_ in regs:
            if gs < re_ and ge > rs:
                hit = True
                break
        if hit:
            out.append(row)
    if not out:
        return pd.DataFrame(columns=genes.columns)
    return pd.DataFrame(out, columns=genes.columns)


def clean_id(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    if s.startswith("ZtIPO323_"):
        s = s.split(".", 1)[0]
    return s


ipo_if = load_ifregions(ipo_if_file)
ir_if = load_ifregions(ir_if_file)
ipo_genes = load_genes_ipo_reannot(ipo_gene_file)
ir_genes = load_genes_ir(ir_gene_file)

ipo_if_genes = overlap_genes_in_if(ipo_genes, ipo_if)
ir_if_genes = overlap_genes_in_if(ir_genes, ir_if)

orth = pd.read_csv(ortho_file, sep="\t", dtype=str)
orth = orth[["IPO323", "IR01_26b"]].copy()
orth["IPO323"] = orth["IPO323"].map(clean_id)
orth["IR01_26b"] = orth["IR01_26b"].map(clean_id)
orth = orth.dropna(subset=["IPO323", "IR01_26b"]).drop_duplicates()

ipo2ir = orth.groupby("IPO323")["IR01_26b"].apply(lambda s: sorted(set(s))).to_dict()
ir2ipo = orth.groupby("IR01_26b")["IPO323"].apply(lambda s: sorted(set(s))).to_dict()

ipo_if_set = set(ipo_if_genes["gene_id"])
ir_if_set = set(ir_if_genes["gene_id"])

ipo_with_ortho = {g for g in ipo_if_set if g in ipo2ir}
ir_with_ortho = {g for g in ir_if_set if g in ir2ipo}

ipo_with_if_ortho = []
for g in sorted(ipo_if_set):
    partners = ipo2ir.get(g, [])
    in_if = [p for p in partners if p in ir_if_set]
    if in_if:
        ipo_with_if_ortho.append((g, ";".join(in_if)))

pair_rows = []
for g, plist in ipo_with_if_ortho:
    for p in plist.split(";"):
        pair_rows.append((g, p))
pairs = pd.DataFrame(pair_rows, columns=["IPO323_gene", "IR01_26b_gene"]).drop_duplicates()

ipo_report = pd.DataFrame({"IPO323_gene": sorted(ipo_if_set)})
ipo_report["has_any_ortholog"] = ipo_report["IPO323_gene"].isin(ipo_with_ortho)
ipo_report["ortholog_in_IR_IF"] = ipo_report["IPO323_gene"].isin(set(pairs["IPO323_gene"])) if len(pairs) else False
ipo_report["IR01_26b_orthologs"] = ipo_report["IPO323_gene"].map(lambda x: ";".join(ipo2ir.get(x, [])))

ir_report = pd.DataFrame({"IR01_26b_gene": sorted(ir_if_set)})
ir_report["has_any_ortholog"] = ir_report["IR01_26b_gene"].isin(ir_with_ortho)
ir_report["ortholog_in_IPO_IF"] = ir_report["IR01_26b_gene"].isin(set(pairs["IR01_26b_gene"])) if len(pairs) else False
ir_report["IPO323_orthologs"] = ir_report["IR01_26b_gene"].map(lambda x: ";".join(ir2ipo.get(x, [])))

ipo_out = f"{out_prefix}/IPO323_IF_genes_ortholog_check.tsv"
ir_out = f"{out_prefix}/IR0126B_IF_genes_ortholog_check.tsv"
pair_out = f"{out_prefix}/IPO323_IR0126B_IF_ortholog_pairs.tsv"
ipo_report.to_csv(ipo_out, sep="\t", index=False)
ir_report.to_csv(ir_out, sep="\t", index=False)
pairs.to_csv(pair_out, sep="\t", index=False)

print("IPO323 IF genes total:", len(ipo_if_set))
print("IPO323 IF genes with any IR ortholog:", len(ipo_with_ortho))
print("IPO323 IF genes with ortholog also in IR IF:", pairs["IPO323_gene"].nunique() if len(pairs) else 0)
print("IR0126B IF genes total:", len(ir_if_set))
print("IR0126B IF genes with any IPO ortholog:", len(ir_with_ortho))
print("IR0126B IF genes with ortholog also in IPO IF:", pairs["IR01_26b_gene"].nunique() if len(pairs) else 0)
print("Ortholog pairs where both genes are in IF:", len(pairs))
print("out1", ipo_out)
print("out2", ir_out)
print("out3", pair_out)
