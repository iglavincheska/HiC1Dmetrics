#!/usr/bin/env python3
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Inputs
IPO_IF_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/IPO323_WT_IFregions_merged.bed"
IR_IF_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/IR0126B_WT_IFregions_merged.bed"
IPO_GENE_FILE = "/Users/iglavinchesk/WP1/FINAL/data/z.tritici.IP0323.reannot.gff3"
IR_GENE_FILE = "/Users/iglavinchesk/WP2/data/00_genomes/IR0126B.genes.gff3"
ORTHO_FILE = "/Users/iglavinchesk/WP1/FINAL/data/orthoproteins.tsv"
IPO_EXPR_FILE = "/Users/iglavinchesk/WP2/6_RNAseq/data/RPKM/IPO323_control_wt.bedgraph"
IR_EXPR_FILE = "/Users/iglavinchesk/WP2/6_RNAseq/data/RPKM/IR26b_control_wt.bed"

OUT_DIR = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/if_distance_ortholog_expression"
PLOT_DIR = os.path.join(OUT_DIR, "plots")

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


def parse_attrs(s):
    d = {}
    for item in str(s).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def load_ifregions(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["chromosome", "start", "end"])
    df["chromosome"] = df["chromosome"].astype(str).map(normalize_chr)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["chr_num"] = df["chromosome"].map(chrom_num_from_name)
    df = df[df["chr_num"].between(CHR_MIN, CHR_MAX)].drop(columns=["chr_num"])
    return df.sort_values(["chromosome", "start", "end"]).reset_index(drop=True)


def load_ipo_genes(path):
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
    return g[["gene_id", "chromosome", "start", "end"]].reset_index(drop=True)


def load_ir_genes(path):
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
    return g[["gene_id", "chromosome", "start", "end"]].reset_index(drop=True)


def load_expression(path, gene_col=5, expr_col=4):
    # chrom, start, end, value1, value2, gene_id, feature
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    out = pd.DataFrame()
    out["gene_id"] = df.iloc[:, gene_col].astype(str)
    out["expr"] = pd.to_numeric(df.iloc[:, expr_col], errors="coerce")
    out = out.dropna(subset=["gene_id", "expr"])
    out = out.groupby("gene_id", as_index=False)["expr"].median()
    return out


def clean_ipo_id(x):
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    if s.startswith("ZtIPO323_"):
        return s.split(".", 1)[0]
    return s


def clean_ir_id(x):
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return None
    return s


def interval_distance(gs, ge, rs, re_):
    if gs < re_ and ge > rs:
        return 0
    if ge <= rs:
        return rs - ge
    return gs - re_


def annotate_distance_to_if(genes, if_regions):
    regs_by_chr = {c: df[["start", "end"]].to_numpy() for c, df in if_regions.groupby("chromosome")}
    dists = []
    for row in genes.itertuples(index=False):
        regs = regs_by_chr.get(row.chromosome)
        if regs is None or len(regs) == 0:
            dists.append(np.nan)
            continue
        gs, ge = row.start, row.end
        best = None
        for rs, re_ in regs:
            d = interval_distance(gs, ge, rs, re_)
            if best is None or d < best:
                best = d
                if best == 0:
                    break
        dists.append(best)
    out = genes.copy()
    out["distance_to_if_bp"] = dists
    return out


def distance_bin(d):
    if pd.isna(d):
        return "NA"
    d = float(d)
    if d == 0:
        return "overlap"
    if d <= 5000:
        return "0-5kb"
    if d <= 20000:
        return "5-20kb"
    if d <= 100000:
        return "20-100kb"
    return ">100kb"


def plot_expression_by_distance(pair_df, out_png):
    order = ["overlap", "0-5kb", "5-20kb", "20-100kb", ">100kb"]

    ipo_data = []
    ir_data = []
    for b in order:
        ipo_vals = np.log1p(pd.to_numeric(pair_df.loc[pair_df["ipo_distance_bin"] == b, "IPO_expr"], errors="coerce").dropna())
        ir_vals = np.log1p(pd.to_numeric(pair_df.loc[pair_df["ir_distance_bin"] == b, "IR_expr"], errors="coerce").dropna())
        ipo_data.append(ipo_vals.values)
        ir_data.append(ir_vals.values)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    axes[0].boxplot(ipo_data, tick_labels=order, showfliers=False)
    axes[0].set_title("IPO323 ortholog genes")
    axes[0].set_xlabel("Distance to nearest IPO323 IF region")
    axes[0].set_ylabel("log1p(Expression)")
    axes[0].grid(axis="y", alpha=0.2)

    axes[1].boxplot(ir_data, tick_labels=order, showfliers=False)
    axes[1].set_title("IR0126B ortholog genes")
    axes[1].set_xlabel("Distance to nearest IR0126B IF region")
    axes[1].grid(axis="y", alpha=0.2)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_pair_difference_by_distance(pair_df, out_png):
    order = ["overlap", "0-5kb", "5-20kb", "20-100kb", ">100kb"]
    tmp = pair_df.copy()
    tmp["log2_fc_ipo_over_ir"] = np.log2((tmp["IPO_expr"] + 1.0) / (tmp["IR_expr"] + 1.0))

    data_ipo = [pd.to_numeric(tmp.loc[tmp["ipo_distance_bin"] == b, "log2_fc_ipo_over_ir"], errors="coerce").dropna().values for b in order]
    data_ir = [pd.to_numeric(tmp.loc[tmp["ir_distance_bin"] == b, "log2_fc_ipo_over_ir"], errors="coerce").dropna().values for b in order]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    axes[0].boxplot(data_ipo, tick_labels=order, showfliers=False)
    axes[0].axhline(0, linestyle="--", color="gray", linewidth=1)
    axes[0].set_title("log2(IPO/IR) by IPO323 IF distance")
    axes[0].set_xlabel("IPO323 distance bin")
    axes[0].set_ylabel("log2((IPO_expr+1)/(IR_expr+1))")
    axes[0].grid(axis="y", alpha=0.2)

    axes[1].boxplot(data_ir, tick_labels=order, showfliers=False)
    axes[1].axhline(0, linestyle="--", color="gray", linewidth=1)
    axes[1].set_title("log2(IPO/IR) by IR0126B IF distance")
    axes[1].set_xlabel("IR0126B distance bin")
    axes[1].grid(axis="y", alpha=0.2)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)

    ipo_if = load_ifregions(IPO_IF_FILE)
    ir_if = load_ifregions(IR_IF_FILE)

    ipo_genes = load_ipo_genes(IPO_GENE_FILE)
    ir_genes = load_ir_genes(IR_GENE_FILE)

    ipo_expr = load_expression(IPO_EXPR_FILE).rename(columns={"gene_id": "IPO323", "expr": "IPO_expr"})
    ir_expr = load_expression(IR_EXPR_FILE).rename(columns={"gene_id": "IR01_26b", "expr": "IR_expr"})

    ipo_dist = annotate_distance_to_if(ipo_genes, ipo_if).rename(columns={"gene_id": "IPO323", "distance_to_if_bp": "ipo_distance_bp"})
    ir_dist = annotate_distance_to_if(ir_genes, ir_if).rename(columns={"gene_id": "IR01_26b", "distance_to_if_bp": "ir_distance_bp"})

    ipo_dist["IPO323"] = ipo_dist["IPO323"].map(clean_ipo_id)
    ir_dist["IR01_26b"] = ir_dist["IR01_26b"].map(clean_ir_id)

    orth = pd.read_csv(ORTHO_FILE, sep="\t", dtype=str)[["IPO323", "IR01_26b"]].copy()
    orth["IPO323"] = orth["IPO323"].map(clean_ipo_id)
    orth["IR01_26b"] = orth["IR01_26b"].map(clean_ir_id)
    orth = orth.dropna(subset=["IPO323", "IR01_26b"]).drop_duplicates()

    pair = orth.merge(ipo_dist[["IPO323", "ipo_distance_bp"]], on="IPO323", how="left")
    pair = pair.merge(ir_dist[["IR01_26b", "ir_distance_bp"]], on="IR01_26b", how="left")
    pair = pair.merge(ipo_expr, on="IPO323", how="left")
    pair = pair.merge(ir_expr, on="IR01_26b", how="left")

    pair["ipo_distance_bin"] = pair["ipo_distance_bp"].map(distance_bin)
    pair["ir_distance_bin"] = pair["ir_distance_bp"].map(distance_bin)
    pair["log2_fc_ipo_over_ir"] = np.log2((pair["IPO_expr"] + 1.0) / (pair["IR_expr"] + 1.0))

    pair_out = os.path.join(OUT_DIR, "ortholog_pairs_expression_if_distance.tsv")
    pair.to_csv(pair_out, sep="\t", index=False)

    # Long-format summary by strain and distance bin
    ipo_summary = pair.groupby("ipo_distance_bin", dropna=False).agg(
        n_pairs=("IPO323", "size"),
        median_expr=("IPO_expr", "median"),
        mean_expr=("IPO_expr", "mean"),
        median_log2_fc=("log2_fc_ipo_over_ir", "median"),
    ).reset_index().rename(columns={"ipo_distance_bin": "distance_bin"})
    ipo_summary["strain"] = "IPO323"

    ir_summary = pair.groupby("ir_distance_bin", dropna=False).agg(
        n_pairs=("IR01_26b", "size"),
        median_expr=("IR_expr", "median"),
        mean_expr=("IR_expr", "mean"),
        median_log2_fc=("log2_fc_ipo_over_ir", "median"),
    ).reset_index().rename(columns={"ir_distance_bin": "distance_bin"})
    ir_summary["strain"] = "IR0126B"

    summary = pd.concat([ipo_summary, ir_summary], ignore_index=True)
    summary_out = os.path.join(OUT_DIR, "expression_by_if_distance_summary.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)

    # Matched-bin summary where ortholog genes are in same distance class in both genomes
    matched = pair[pair["ipo_distance_bin"] == pair["ir_distance_bin"]].copy()
    matched_summary = matched.groupby("ipo_distance_bin", dropna=False).agg(
        n_pairs=("IPO323", "size"),
        ipo_median_expr=("IPO_expr", "median"),
        ir_median_expr=("IR_expr", "median"),
        median_log2_fc=("log2_fc_ipo_over_ir", "median"),
    ).reset_index().rename(columns={"ipo_distance_bin": "distance_bin"})
    matched_out = os.path.join(OUT_DIR, "matched_distance_bin_ortholog_expression_summary.tsv")
    matched_summary.to_csv(matched_out, sep="\t", index=False)

    # Plots
    p1 = os.path.join(PLOT_DIR, "expression_by_if_distance_bins.png")
    p2 = os.path.join(PLOT_DIR, "ortholog_log2fc_by_if_distance_bins.png")
    plot_expression_by_distance(pair, p1)
    plot_pair_difference_by_distance(pair, p2)

    print("Saved:")
    print(pair_out)
    print(summary_out)
    print(matched_out)
    print(p1)
    print(p2)


if __name__ == "__main__":
    main()
