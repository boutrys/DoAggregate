"""Summarise the DoAggregate run on the public 1000 Genomes exome example.

Reads the DoAggregate output folder produced on data/big_example
(26 CEU individuals as "cases" vs 25 GBR individuals as "controls") and the
small test dataset (data/small_example), and writes the tables, figure
and numbers used in the Results section of the bioRxiv manuscript.

Usage (from the repository root):
    python3 article/scripts/01_summarise_public_data_results.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BIG = ROOT / "data" / "big_example"
RES = BIG / "results"
SMALL = ROOT / "data" / "small_example"
OUT = ROOT / "article"
LEVELS = {"gene": "Gene", "over_represented_pathway": "Over-represented pathway",
          "over_represented_GO": "Over-represented GO term"}
Q = 0.05

numbers = {}

# Cohort
cases = [l.strip() for l in open(BIG / "cases.txt") if l.strip()]
controls = [l.strip() for l in open(BIG / "controls.txt") if l.strip()]
numbers["n_cases"], numbers["n_controls"] = len(cases), len(controls)

# Global information per level of analysis (file ID 35)
glob = pd.read_excel(RES / "general" / "global_info_analysis.xlsx").set_index("name_of_experiment")
numbers["global_info"] = glob.T.to_dict()

# Over-representation summary (file ID 34)
ora = pd.read_excel(RES / "general" / "Summary_results_over_representation.xlsx")
numbers["ora_terms_total"] = len(ora)
numbers["ora_terms_by_database"] = ora["Database"].value_counts().to_dict()

rows, sig_tables = [], []
for folder, label in LEVELS.items():
    gi = pd.read_excel(RES / folder / "general_info.xlsx")
    sig = gi[gi["Excalibur_qvalue"] <= Q].copy()
    case_driven = sig[sig["Total_nbr_variant_in_patients"] >= sig["Total_nbr_variant_in_controls"]]
    best_only = gi[(gi["best_qvalue"] <= Q)]
    rows.append({
        "Level of analysis": label,
        "Regions defined": int(glob.loc[folder, "Nbr_region"]),
        "Regions tested": int(glob.loc[folder, "Nbr_region_analyzed"]),
        "Regions not tested (<2 unique variants)": int(glob.loc[folder, "Nbr_region_not_analyzed"]),
        "Median aggregation tests applied per region": float(gi["Nbr_test_applied"].median()),
        "Excalibur q <= 0.05": len(sig),
        "...of which more variants in cases": len(case_driven),
        "Best single test q <= 0.05": len(best_only),
        "Lowest Excalibur q": float(gi["Excalibur_qvalue"].min()),
    })
    numbers[f"{folder}_best_test_counts"] = best_only["best_Stat_test"].value_counts().to_dict()
    numbers[f"{folder}_tests_applied_range"] = [int(gi["Nbr_test_applied"].min()), int(gi["Nbr_test_applied"].max())]
    numbers[f"{folder}_median_excalibur_q"] = float(gi["Excalibur_qvalue"].median())
    sig.insert(0, "Level of analysis", label)
    sig_tables.append(sig)

table2 = pd.DataFrame(rows)
table3 = pd.concat(sig_tables, ignore_index=True)
keep = ["Level of analysis", "name_genetic_region", "Excalibur_qvalue", "best_qvalue", "best_Stat_test",
        "Nbr_test_applied", "Nbr_patient_with_mutation", "Nbr_control_with_mutation",
        "Total_nbr_variant_in_patients", "Total_nbr_variant_in_controls", "Nbr_unique_variant"]
table3 = table3[keep]

# Small example (GVCF input with coverage-based cleaning)
clean = pd.read_excel(SMALL / "impact_of_cleaning.xlsx").set_index("Info")
numbers["small_example_cleaning"] = clean.to_dict()
sglob = pd.read_excel(SMALL / "results" / "general" / "global_info_analysis.xlsx").set_index("name_of_experiment")
numbers["small_example_global_info"] = sglob.T.to_dict()

with pd.ExcelWriter(OUT / "tables" / "Table1_summary_by_level.xlsx") as w:
    table2.to_excel(w, index=False)
table2.to_csv(OUT / "tables" / "Table1_summary_by_level.csv", index=False)
with pd.ExcelWriter(OUT / "tables" / "Table2_significant_regions.xlsx") as w:
    table3.to_excel(w, index=False)
table3.to_csv(OUT / "tables" / "Table2_significant_regions.csv", index=False)
with open(OUT / "tables" / "results_numbers.json", "w") as f:
    json.dump(numbers, f, indent=2, default=str)

# Figure 2
gene = pd.read_excel(RES / "gene" / "general_info.xlsx").sort_values("Excalibur_qvalue")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

ax = axes[0]
labels = ["Gene", "Pathway", "GO term"]
defined = table2["Regions defined"].values
tested = table2["Regions tested"].values
sigs = table2["Excalibur q <= 0.05"].values
x = np.arange(3)
ax.bar(x - 0.27, defined, 0.27, label="Defined", color="#9db4c0")
ax.bar(x, tested, 0.27, label="Tested", color="#5c6b73")
ax.bar(x + 0.27, sigs, 0.27, label="Excalibur q ≤ 0.05", color="#c2563a")
ax.set_yscale("symlog", linthresh=1)
ax.set_xticks(x, labels)
ax.set_ylabel("Number of genetic regions")
for xi, v in zip(x + 0.27, sigs):
    ax.text(xi, max(v, 0) + 0.3, str(v), ha="center", va="bottom", fontsize=8)
ax.legend(frameon=False, fontsize=8)
ax.set_title("a  Regions per level of analysis", loc="left", fontsize=10)

ax = axes[1]
y = -np.log10(gene["Excalibur_qvalue"].clip(lower=1e-300))
ax.plot(np.arange(1, len(y) + 1), y, ".", ms=3, color="#5c6b73")
ax.axhline(-np.log10(Q), ls="--", lw=0.8, color="#c2563a")
top = gene.head(5)
for i, (_, r) in enumerate(top.iterrows()):
    ax.annotate(r["name_genetic_region"], (i + 1, -np.log10(r["Excalibur_qvalue"])),
                xytext=(40, -12 * i), textcoords="offset points", fontsize=7,
                arrowprops=dict(arrowstyle="-", lw=0.4))
ax.set_xscale("log")
ax.set_xlabel("Gene rank")
ax.set_ylabel("-log10(Excalibur q-value)")
ax.set_title("b  Gene-level ranking", loc="left", fontsize=10)

ax = axes[2]
sigg = gene[gene["Excalibur_qvalue"] <= Q]
ind = np.arange(len(sigg))
ax.barh(ind - 0.2, sigg["Total_nbr_variant_in_patients"], 0.4, label="Cases (CEU)", color="#c2563a")
ax.barh(ind + 0.2, sigg["Total_nbr_variant_in_controls"], 0.4, label="Controls (GBR)", color="#5c6b73")
ax.set_yticks(ind, sigg["name_genetic_region"])
ax.invert_yaxis()
ax.set_xlabel("Total rare variant calls")
ax.legend(frameon=False, fontsize=8)
ax.set_title("c  Significant genes: variant burden", loc="left", fontsize=10)

for a in axes:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "figures" / "Figure2_public_data.png", dpi=300)
fig.savefig(OUT / "figures" / "Figure2_public_data.pdf")
print(table2.to_string())
print(table3.to_string())
print(json.dumps({k: v for k, v in numbers.items() if "global" not in k and "cleaning" not in k}, indent=1, default=str))
print(clean.to_string())
