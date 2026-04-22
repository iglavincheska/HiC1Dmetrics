#%%
import pandas as pd
import numpy as np
import re

# =========================
# INPUTS
# =========================
ipo_if_file = "/Users/iglavinchesk/WP2/HiC1Dmetrics/IPO323_WT_IFregions_merged.bed"
ir_if_file = "/Users/iglavinchesk/WP2/HiC1Dmetrics/IR0126B_WT_IFregions_merged.bed"

ipo_te_file = "/Users/iglavinchesk/WP2/data/00_genomes/IPO323.TE_2023.gff"
ir_te_file = "/Users/iglavinchesk/WP2/data/00_genomes/IR0126B.TE_2023.gff"

ipo_gene_file = "/Users/iglavinchesk/WP2/data/00_genomes/IPO323.genes.gff3"
ir_gene_file = "/Users/iglavinchesk/WP2/data/00_genomes/IR0126B.genes.gff3"

ipo_gt_file = "/Users/iglavinchesk/packages/3D-genome-tools/IPO323_genome.txt"
ir_gt_file = "/Users/iglavinchesk/WP2/data/00_genomes/IR0126b.genome"

CHR_MIN = 1
CHR_MAX = 13

# =========================
# HELPERS
# =========================
def chrom_num_from_name(chrom):
    s = str(chrom).strip()
    s_low = s.lower()

    # Prefer a terminal chr token to avoid mixing digits from strain IDs
    # e.g. IR01_26b.chr_1 -> 1 (not 1261).
    m = re.search(r"(?:^|[._-])chr(?:omosome)?[._-]?(\d+)$", s_low)
    if m:
        return int(m.group(1))

    # Fallback for simple names like 1 / chr1 / chromosome1.
    m = re.search(r"^(?:chr(?:omosome)?)?[._-]?(\d+)$", s_low)
    if m:
        return int(m.group(1))

    return np.nan

def normalize_chr(chrom):
    # Keep a standard representation for matching between files
    n = chrom_num_from_name(chrom)
    if np.isnan(n):
        return str(chrom).strip()
    return str(int(n))  # use plain numeric chromosome labels

def load_ifregions(path):
    # merged bed: chr, start, end
    df = pd.read_csv(path, sep="\t", header=None, names=["chromosome", "start", "end"])
    df["chromosome"] = df["chromosome"].astype(str).map(normalize_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["chr_num"] = df["chromosome"].map(chrom_num_from_name)
    df = df[df["chr_num"].between(CHR_MIN, CHR_MAX)]
    df = df.drop(columns=["chr_num"])
    df["length"] = df["end"] - df["start"]
    return df.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)

def load_tes(path):
    te_cols = ["chromosome", "source", "type", "start", "end", "score", "strand", "frame", "attributes"]
    df = pd.read_csv(path, sep="\t", header=None, names=te_cols, comment="#")
    df["chromosome"] = df["chromosome"].astype(str).map(normalize_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["chr_num"] = df["chromosome"].map(chrom_num_from_name)
    df = df[df["chr_num"].between(CHR_MIN, CHR_MAX)]
    df = df.drop(columns=["chr_num"])
    df["length"] = df["end"] - df["start"]
    return df.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)

def _parse_gff_attributes(attr_text):
    d = {}
    for item in str(attr_text).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def _infer_gene_id(attrs, feature_type):
    if feature_type == "gene":
        return attrs.get("ID") or attrs.get("Name")
    return attrs.get("geneID") or attrs.get("Parent") or attrs.get("ID")

def load_protein_coding_genes(path):
    cols = ["chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, comment="#")
    df["type"] = df["type"].astype(str)

    # Prefer gene features if present; otherwise use mRNA as gene proxies.
    if (df["type"] == "gene").any():
        features = df[df["type"] == "gene"].copy()
    else:
        features = df[df["type"] == "mRNA"].copy()

    parsed_attrs = features["attributes"].map(_parse_gff_attributes)

    # Keep explicit protein-coding entries when biotype is provided;
    # otherwise keep all selected gene-like features.
    biotype = parsed_attrs.map(lambda x: (x.get("biotype") or x.get("gene_biotype") or "").lower())
    has_biotype = biotype != ""
    features = features[(~has_biotype) | (biotype == "protein_coding")].copy()
    parsed_attrs = features["attributes"].map(_parse_gff_attributes)

    features["gene_id"] = [
        _infer_gene_id(attrs, ftype)
        for attrs, ftype in zip(parsed_attrs, features["type"])
    ]
    features["gene_id"] = features["gene_id"].fillna("")
    features = features[features["gene_id"] != ""].copy()

    features["chromosome"] = features["chromosome"].astype(str).map(normalize_chr)
    features["start"] = pd.to_numeric(features["start"], errors="coerce")
    features["end"] = pd.to_numeric(features["end"], errors="coerce")
    features = features.dropna(subset=["start", "end"])
    features["chr_num"] = features["chromosome"].map(chrom_num_from_name)
    features = features[features["chr_num"].between(CHR_MIN, CHR_MAX)]
    features = features.drop(columns=["chr_num"])
    features["length"] = features["end"] - features["start"]

    # De-duplicate multi-isoform proxy entries by gene_id and coordinates.
    features = features.drop_duplicates(subset=["chromosome", "start", "end", "gene_id"])
    return features.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)

def genome_size_from_gt(path):
    gt = pd.read_csv(path, sep="\t", header=None, names=["chromosome", "length"])
    gt["chromosome"] = gt["chromosome"].astype(str).map(normalize_chr)
    gt["chr_num"] = gt["chromosome"].map(chrom_num_from_name)
    gt = gt[gt["chr_num"].between(CHR_MIN, CHR_MAX)]
    return pd.to_numeric(gt["length"], errors="coerce").fillna(0).sum()

def summarize_ifregions(df, genome_size, name):
    print("\n" + "="*65)
    print(f"IFregions summary for {name}")
    print("="*65)

    n_regions = len(df)
    total_bp = df["length"].sum()
    mean_bp = df["length"].mean() if n_regions else 0
    median_bp = df["length"].median() if n_regions else 0
    frac = (total_bp / genome_size * 100) if genome_size > 0 else np.nan

    print(f"Number of IF regions: {n_regions:,}")
    print(f"Total IF-region length: {int(total_bp):,} bp")
    print(f"Mean IF-region size: {mean_bp:,.0f} bp")
    print(f"Median IF-region size: {median_bp:,.0f} bp")
    print(f"Genome covered by IF regions: {frac:.2f}%")

def overlap_stats(if_df, te_df, name):
    print("\n" + "-"*65)
    print(f"TE overlap with IFregions: {name}")
    print("-"*65)

    total_tes = len(te_df)
    total_te_bp = te_df["length"].sum()

    tes_overlapping = 0
    overlap_bp_total = 0

    # per chromosome to reduce comparisons
    for chrom in sorted(set(te_df["chromosome"]).intersection(set(if_df["chromosome"]))):
        tes_chr = te_df[te_df["chromosome"] == chrom]
        if_chr = if_df[if_df["chromosome"] == chrom]

        # simple interval loop; fine for moderate-size sets
        for _, te in tes_chr.iterrows():
            te_s, te_e = te["start"], te["end"]
            te_overlaps_any = False
            te_overlap_bp = 0

            for _, reg in if_chr.iterrows():
                rs, re = reg["start"], reg["end"]
                if te_s < re and te_e > rs:
                    te_overlaps_any = True
                    ov = min(te_e, re) - max(te_s, rs)
                    if ov > 0:
                        te_overlap_bp += ov

            if te_overlaps_any:
                tes_overlapping += 1
            overlap_bp_total += te_overlap_bp

    pct_tes = (tes_overlapping / total_tes * 100) if total_tes else np.nan
    pct_bp = (overlap_bp_total / total_te_bp * 100) if total_te_bp else np.nan

    print(f"Total TEs: {total_tes:,}")
    print(f"TEs overlapping IF regions: {tes_overlapping:,}")
    print(f"% of TEs overlapping IF regions: {pct_tes:.2f}%")
    print(f"Total TE bp: {int(total_te_bp):,}")
    print(f"TE bp overlapping IF regions: {int(overlap_bp_total):,}")
    print(f"% of TE bp overlapping IF regions: {pct_bp:.2f}%")

    return {
        "total_tes": total_tes,
        "overlap_tes": tes_overlapping,
        "pct_tes": pct_tes,
        "total_te_bp": total_te_bp,
        "overlap_te_bp": overlap_bp_total,
        "pct_te_bp": pct_bp
    }

def gene_overlap_stats(if_df, gene_df, name):
    print("\n" + "-"*65)
    print(f"Protein-coding gene overlap with IFregions: {name}")
    print("-"*65)

    total_genes = len(gene_df)
    total_gene_bp = gene_df["length"].sum()

    genes_overlapping = 0
    overlap_bp_total = 0

    for chrom in sorted(set(gene_df["chromosome"]).intersection(set(if_df["chromosome"]))):
        genes_chr = gene_df[gene_df["chromosome"] == chrom]
        if_chr = if_df[if_df["chromosome"] == chrom]

        for _, gene in genes_chr.iterrows():
            g_s, g_e = gene["start"], gene["end"]
            gene_overlaps_any = False
            gene_overlap_bp = 0

            for _, reg in if_chr.iterrows():
                rs, re = reg["start"], reg["end"]
                if g_s < re and g_e > rs:
                    gene_overlaps_any = True
                    ov = min(g_e, re) - max(g_s, rs)
                    if ov > 0:
                        gene_overlap_bp += ov

            if gene_overlaps_any:
                genes_overlapping += 1
            overlap_bp_total += gene_overlap_bp

    pct_genes = (genes_overlapping / total_genes * 100) if total_genes else np.nan
    pct_bp = (overlap_bp_total / total_gene_bp * 100) if total_gene_bp else np.nan

    print(f"Total protein-coding genes: {total_genes:,}")
    print(f"Genes overlapping IF regions: {genes_overlapping:,}")
    print(f"% genes overlapping IF regions: {pct_genes:.2f}%")
    print(f"Total gene bp: {int(total_gene_bp):,}")
    print(f"Gene bp overlapping IF regions: {int(overlap_bp_total):,}")
    print(f"% gene bp overlapping IF regions: {pct_bp:.2f}%")

    return {
        "total_genes": total_genes,
        "overlap_genes": genes_overlapping,
        "pct_genes": pct_genes,
        "total_gene_bp": total_gene_bp,
        "overlap_gene_bp": overlap_bp_total,
        "pct_gene_bp": pct_bp
    }

# =========================
# LOAD
# =========================
ipo_if = load_ifregions(ipo_if_file)
ir_if = load_ifregions(ir_if_file)

ipo_tes = load_tes(ipo_te_file)
ir_tes = load_tes(ir_te_file)

ipo_genes = load_protein_coding_genes(ipo_gene_file)
ir_genes = load_protein_coding_genes(ir_gene_file)

ipo_genome = genome_size_from_gt(ipo_gt_file)
ir_genome = genome_size_from_gt(ir_gt_file)

# =========================
# SUMMARIES
# =========================
summarize_ifregions(ipo_if, ipo_genome, "IPO323_WT")
summarize_ifregions(ir_if, ir_genome, "IR0126B_WT")

ipo_te_stats = overlap_stats(ipo_if, ipo_tes, "IPO323_WT")
ir_te_stats = overlap_stats(ir_if, ir_tes, "IR0126B_WT")

ipo_gene_stats = gene_overlap_stats(ipo_if, ipo_genes, "IPO323_WT")
ir_gene_stats = gene_overlap_stats(ir_if, ir_genes, "IR0126B_WT")

# =========================
# COMPARISON
# =========================
print("\n" + "="*65)
print("Comparison: IPO323_WT vs IR0126B_WT")
print("="*65)

ipo_if_bp = ipo_if["length"].sum()
ir_if_bp = ir_if["length"].sum()

print(f"IF-region bp IPO323_WT: {int(ipo_if_bp):,}")
print(f"IF-region bp IR0126B_WT: {int(ir_if_bp):,}")
print(f"Difference (abs): {int(abs(ipo_if_bp - ir_if_bp)):,} bp")

print(f"\n% TE bp overlapping IF regions")
print(f"IPO323_WT: {ipo_te_stats['pct_te_bp']:.2f}%")
print(f"IR0126B_WT: {ir_te_stats['pct_te_bp']:.2f}%")

print(f"\n% protein-coding genes overlapping IF regions")
print(f"IPO323_WT: {ipo_gene_stats['pct_genes']:.2f}%")
print(f"IR0126B_WT: {ir_gene_stats['pct_genes']:.2f}%")