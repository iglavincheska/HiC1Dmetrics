import pandas as pd

ORTHO_FILE = "/Users/iglavinchesk/WP1/FINAL/data/orthoproteins.tsv"
IPO_IF_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/IPO323_IF_genes_ortholog_check.tsv"
IR_IF_FILE = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts/IR0126B_IF_genes_ortholog_check.tsv"
IPO_EXPR_FILE = "/Users/iglavinchesk/WP2/6_RNAseq/data/RPKM/IPO323_control_wt.bedgraph"
IR_EXPR_FILE = "/Users/iglavinchesk/WP2/6_RNAseq/data/RPKM/IR26b_control_wt.bed"
OUT_DIR = "/Users/iglavinchesk/WP2/HiC1Dmetrics/downstream_scripts"


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


def load_expression(path):
    # Expected columns: chrom, start, end, value1, value2, gene_id, feature
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    if df.shape[1] < 6:
        raise ValueError(f"Unexpected expression format in {path}: {df.shape[1]} columns")

    out = pd.DataFrame()
    out["gene_id"] = df.iloc[:, 5].astype(str)
    out["expr_value1"] = pd.to_numeric(df.iloc[:, 3], errors="coerce")
    out["expr_value2"] = pd.to_numeric(df.iloc[:, 4], errors="coerce")
    return out.dropna(subset=["gene_id"]).drop_duplicates(subset=["gene_id"]) 


orth = pd.read_csv(ORTHO_FILE, sep="\t", dtype=str)[["IPO323", "IR01_26b"]].copy()
orth["IPO323"] = orth["IPO323"].map(clean_ipo_id)
orth["IR01_26b"] = orth["IR01_26b"].map(clean_ir_id)
orth = orth.dropna(subset=["IPO323", "IR01_26b"]).drop_duplicates()

ipo_if = pd.read_csv(IPO_IF_FILE, sep="\t")
ir_if = pd.read_csv(IR_IF_FILE, sep="\t")

ipo_if_set = set(ipo_if["IPO323_gene"].astype(str))
ir_if_set = set(ir_if["IR01_26b_gene"].astype(str))

pairs = orth.copy()
pairs["ipo_in_if"] = pairs["IPO323"].isin(ipo_if_set)
pairs["ir_in_if"] = pairs["IR01_26b"].isin(ir_if_set)
pairs["if_xor"] = pairs["ipo_in_if"] ^ pairs["ir_in_if"]

ipo_only_if = pairs[(pairs["ipo_in_if"]) & (~pairs["ir_in_if"])].copy()
ir_only_if = pairs[(pairs["ir_in_if"]) & (~pairs["ipo_in_if"])].copy()

ipo_expr = load_expression(IPO_EXPR_FILE).rename(
    columns={"gene_id": "IPO323", "expr_value1": "IPO_expr_value1", "expr_value2": "IPO_expr_value2"}
)
ir_expr = load_expression(IR_EXPR_FILE).rename(
    columns={"gene_id": "IR01_26b", "expr_value1": "IR_expr_value1", "expr_value2": "IR_expr_value2"}
)

merged = pairs.merge(ipo_expr, on="IPO323", how="left").merge(ir_expr, on="IR01_26b", how="left")
merged_ipo_only = merged[(merged["ipo_in_if"]) & (~merged["ir_in_if"])].copy()
merged_ir_only = merged[(merged["ir_in_if"]) & (~merged["ipo_in_if"])].copy()

merged_ipo_only.to_csv(f"{OUT_DIR}/ortholog_pairs_IF_only_IPO323_with_expression.tsv", sep="\t", index=False)
merged_ir_only.to_csv(f"{OUT_DIR}/ortholog_pairs_IF_only_IR0126B_with_expression.tsv", sep="\t", index=False)

summary = pd.DataFrame([
    {
        "group": "IPO323_in_IF_only",
        "ortholog_pairs": len(merged_ipo_only),
        "unique_IPO323_genes": merged_ipo_only["IPO323"].nunique(),
        "unique_IR0126B_genes": merged_ipo_only["IR01_26b"].nunique(),
        "IPO_expr_value2_median": merged_ipo_only["IPO_expr_value2"].median(),
        "IR_expr_value2_median": merged_ipo_only["IR_expr_value2"].median(),
    },
    {
        "group": "IR0126B_in_IF_only",
        "ortholog_pairs": len(merged_ir_only),
        "unique_IPO323_genes": merged_ir_only["IPO323"].nunique(),
        "unique_IR0126B_genes": merged_ir_only["IR01_26b"].nunique(),
        "IPO_expr_value2_median": merged_ir_only["IPO_expr_value2"].median(),
        "IR_expr_value2_median": merged_ir_only["IR_expr_value2"].median(),
    },
])
summary.to_csv(f"{OUT_DIR}/ortholog_IF_one_strain_expression_summary.tsv", sep="\t", index=False)

print("Total ortholog pairs:", len(pairs))
print("Pairs with IF in IPO323 only:", len(merged_ipo_only))
print("Pairs with IF in IR0126B only:", len(merged_ir_only))
print("\nSummary")
print(summary.to_string(index=False))
print("\nOutputs:")
print(f"{OUT_DIR}/ortholog_pairs_IF_only_IPO323_with_expression.tsv")
print(f"{OUT_DIR}/ortholog_pairs_IF_only_IR0126B_with_expression.tsv")
print(f"{OUT_DIR}/ortholog_IF_one_strain_expression_summary.tsv")
