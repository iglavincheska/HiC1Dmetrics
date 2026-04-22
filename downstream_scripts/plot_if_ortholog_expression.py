import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

IPO_ONLY_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/ortholog_pairs_IF_only_IPO323_with_expression.tsv"
IR_ONLY_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/ortholog_pairs_IF_only_IR0126B_with_expression.tsv"
SUMMARY_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/ortholog_IF_one_strain_expression_summary.tsv"
OUT_DIR = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/plots"


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def clean_series(s):
    s = pd.to_numeric(s, errors="coerce")
    return s.replace([np.inf, -np.inf], np.nan).dropna()


def plot_boxplots(ipo_df, ir_df, out_path):
    ipo_if = clean_series(ipo_df["IPO_expr_value2"])
    ipo_partner = clean_series(ipo_df["IR_expr_value2"])
    ir_if = clean_series(ir_df["IR_expr_value2"])
    ir_partner = clean_series(ir_df["IPO_expr_value2"])

    data = [
        np.log1p(ipo_if),
        np.log1p(ipo_partner),
        np.log1p(ir_if),
        np.log1p(ir_partner),
    ]
    labels = [
        "IPO323 IF-only genes\n(IPO323 expr)",
        "IPO323 IF-only orthologs\n(IR0126B expr)",
        "IR0126B IF-only genes\n(IR0126B expr)",
        "IR0126B IF-only orthologs\n(IPO323 expr)",
    ]

    fig, ax = plt.subplots(figsize=(11, 6))
    box = ax.boxplot(data, patch_artist=True, tick_labels=labels, showfliers=False)
    colors = ["#1b9e77", "#66a61e", "#d95f02", "#7570b3"]
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)

    ax.set_ylabel("log1p(Expression value2)")
    ax.set_title("Expression of orthologs when IF occurs in one strain only")
    ax.grid(axis="y", alpha=0.2)
    plt.xticks(rotation=12, ha="right")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_paired_scatter(ipo_df, ir_df, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=False, sharey=False)

    # Panel 1: IPO IF-only pairs
    x1 = np.log1p(pd.to_numeric(ipo_df["IR_expr_value2"], errors="coerce"))
    y1 = np.log1p(pd.to_numeric(ipo_df["IPO_expr_value2"], errors="coerce"))
    m1 = x1.notna() & y1.notna()
    x1 = x1[m1]
    y1 = y1[m1]
    axes[0].scatter(x1, y1, s=18, alpha=0.55, color="#1b9e77")
    mn1 = min(x1.min(), y1.min()) if len(x1) else 0
    mx1 = max(x1.max(), y1.max()) if len(x1) else 1
    axes[0].plot([mn1, mx1], [mn1, mx1], "--", color="gray", linewidth=1)
    axes[0].set_title("IPO323 IF-only ortholog pairs")
    axes[0].set_xlabel("IR0126B ortholog expression (log1p)")
    axes[0].set_ylabel("IPO323 IF-gene expression (log1p)")
    axes[0].grid(alpha=0.2)

    # Panel 2: IR IF-only pairs
    x2 = np.log1p(pd.to_numeric(ir_df["IPO_expr_value2"], errors="coerce"))
    y2 = np.log1p(pd.to_numeric(ir_df["IR_expr_value2"], errors="coerce"))
    m2 = x2.notna() & y2.notna()
    x2 = x2[m2]
    y2 = y2[m2]
    axes[1].scatter(x2, y2, s=18, alpha=0.55, color="#d95f02")
    mn2 = min(x2.min(), y2.min()) if len(x2) else 0
    mx2 = max(x2.max(), y2.max()) if len(x2) else 1
    axes[1].plot([mn2, mx2], [mn2, mx2], "--", color="gray", linewidth=1)
    axes[1].set_title("IR0126B IF-only ortholog pairs")
    axes[1].set_xlabel("IPO323 ortholog expression (log1p)")
    axes[1].set_ylabel("IR0126B IF-gene expression (log1p)")
    axes[1].grid(alpha=0.2)

    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_summary_bars(summary_df, out_path):
    labels = summary_df["group"].tolist()
    ipo_m = pd.to_numeric(summary_df["IPO_expr_value2_median"], errors="coerce").values
    ir_m = pd.to_numeric(summary_df["IR_expr_value2_median"], errors="coerce").values

    x = np.arange(len(labels))
    w = 0.35

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - w / 2, ipo_m, w, label="IPO323 median expr", color="#4c78a8", alpha=0.85)
    ax.bar(x + w / 2, ir_m, w, label="IR0126B median expr", color="#f58518", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Median expression value2")
    ax.set_title("Median expression by IF-only ortholog group")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    ensure_dir(OUT_DIR)
    ipo_df = pd.read_csv(IPO_ONLY_FILE, sep="\t")
    ir_df = pd.read_csv(IR_ONLY_FILE, sep="\t")
    summary_df = pd.read_csv(SUMMARY_FILE, sep="\t")

    plot_boxplots(ipo_df, ir_df, os.path.join(OUT_DIR, "if_only_ortholog_expression_boxplots.png"))
    plot_paired_scatter(ipo_df, ir_df, os.path.join(OUT_DIR, "if_only_ortholog_expression_paired_scatter.png"))
    plot_summary_bars(summary_df, os.path.join(OUT_DIR, "if_only_ortholog_expression_median_bars.png"))

    print("Saved plots:")
    print(os.path.join(OUT_DIR, "if_only_ortholog_expression_boxplots.png"))
    print(os.path.join(OUT_DIR, "if_only_ortholog_expression_paired_scatter.png"))
    print(os.path.join(OUT_DIR, "if_only_ortholog_expression_median_bars.png"))


if __name__ == "__main__":
    main()
