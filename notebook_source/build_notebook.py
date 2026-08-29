#!/usr/bin/env python3
"""
build_notebook.py  -  Generate RHT_Adult_Analysis.ipynb.

Single optical section analysis, WT vs Opn4 KO. Reports the outliers-removed
statistics (green sections with low SNR or tissue damage excluded), and shows
figures for the full dataset followed by the outliers-removed dataset.

Reads:
    data/stats_single_section_outliers_removed.mat   (the reported statistics)
    data/figures/single_section_all_data/            (full-dataset figures)
    data/figures/single_section_outliers_removed/    (outliers-removed figures)

Run:
    python build_notebook.py
"""

import os
import sys

import numpy as np

try:
    import h5py
except ImportError:
    sys.exit("h5py is required (pip install -r requirements-dev.txt)")

try:
    import nbformat as nbf
except ImportError:
    sys.exit("nbformat is required (pip install -r requirements-dev.txt)")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
DATA = os.path.join(REPO, "data")
FIGROOT = os.path.join(DATA, "figures")

STATS_FILE = "stats_single_section_outliers_removed.mat"

FULL_RUN = "single_section_all_data"          # figures only
CURATED_RUN = "single_section_outliers_removed"  # figures + stats
CURATED_SECTIONS = ["01_Linear_Stretch", "Aligned_Profiles", "05_Heatmaps"]

CHANNELS = ["Red", "Green"]
HALVES = ["negative", "positive"]

# Section order + headings for the full-dataset figure set.
SECTION_ORDER = [
    ("01_Linear_Stretch",       "Image montages"),
    ("02_Histograms",           "Aligned intensity-profile montages"),
    ("Before_After_Comparison", "Profiles before vs after alignment"),
    ("QC_Plots",                "Alignment QC"),
    ("Aligned_Profiles",        "Aligned WT vs KO profiles (normalized)"),
    ("Statistical_Analysis",    "Per-animal statistics (normalized)"),
    ("Outliers",                "Largest deviations from group mean"),
    ("05_Heatmaps",             "Plasma difference heatmaps with average images"),
]
CURATED_HEADINGS = {
    "01_Linear_Stretch": "Image montages",
    "Aligned_Profiles":  "Aligned WT vs KO profiles (normalized)",
    "05_Heatmaps":       "Plasma difference heatmaps with average images",
}


def _num(h, path):
    return float(np.array(h[path]).squeeze())


def read_stats(mat_path):
    """Per-channel normalized stats + Ns for one run."""
    out = {"n": {}, "stats": {}}
    with h5py.File(mat_path, "r") as h:
        for ch in CHANNELS:
            b = f"data/statsData/Linear/{ch}/normalized"
            out["n"][ch] = (int(_num(h, b + "/num_WT")), int(_num(h, b + "/num_KO")))
            out["stats"][ch] = {}
            for half in HALVES:
                out["stats"][ch][half] = {
                    "p": _num(h, f"{b}/p_{half}"),
                    "power": _num(h, f"{b}/power_{half}"),
                    "d": _num(h, f"{b}/cohens_d_{half}"),
                }
    return out


def fmt_sig(p, alpha=0.05):
    return "**significant**" if p < alpha else "n.s."


def stats_table_md(data, channel):
    n = data["n"][channel]
    rows = [
        f"#### {channel} channel (n = {n[0]} WT vs {n[1]} KO)",
        "",
        "| SCN half | p | power | Cohen's d | result |",
        "|---|---|---|---|---|",
    ]
    for half in HALVES:
        s = data["stats"][channel][half]
        rows.append(f"| {half} | {s['p']:.4f} | {s['power']:.2f} | {s['d']:+.2f} | {fmt_sig(s['p'])} |")
    rows.append("")
    return "\n".join(rows)


def find_pngs(run_dir, section_substr):
    base = os.path.join(FIGROOT, run_dir)
    if not os.path.isdir(base):
        return []
    hits = []
    for dirpath, _dirs, files in os.walk(base):
        for fn in files:
            if fn.lower().endswith(".png") and section_substr.lower() in os.path.join(dirpath, fn).lower():
                hits.append(os.path.relpath(os.path.join(dirpath, fn), REPO))
    return sorted(hits)


def gallery_cell(paths):
    listing = "[\n    " + ",\n    ".join(repr(p) for p in paths) + ",\n]"
    return nbf.v4.new_code_cell(
        "from IPython.display import Image, display, Markdown\n"
        f"_figs = {listing}\n"
        "for _p in _figs:\n"
        "    display(Markdown('*' + _p.split('/')[-1] + '*'))\n"
        "    display(Image(filename=_p))\n"
    )


def build():
    data_path = os.path.join(DATA, STATS_FILE)
    if not os.path.exists(data_path):
        sys.exit(f"Missing stats file: {data_path}")
    data = read_stats(data_path)

    cells = []
    def md(t):
        cells.append(nbf.v4.new_markdown_cell(t))

    md(
        "# Adult RHT analysis - RHT axon distribution in the SCN (WT vs *Opn4* KO)\n\n"
        "1D spatial-profile analysis of retinohypothalamic tract (RHT) axon signal in "
        "the suprachiasmatic nucleus (SCN), wild-type vs melanopsin (*Opn4*) knockout. "
        "Single optical sections. Dual-eye CTB tracing; the **red** and **green** "
        "channels report the two eyes' projections within the same section."
    )

    md(
        "## Result\n\n"
        "No spatial-pattern difference in RHT innervation between WT and KO, in either "
        "channel. The spatial-pattern (normalized) profiles are the analysis metric."
    )

    md(
        "## Statistics\n\n"
        "Single optical section, outliers removed (see below). Spatial-pattern "
        "(normalized) profiles; WT vs KO by Mann-Whitney U with Cohen's d and post-hoc "
        "power."
    )
    for ch in CHANNELS:
        md(stats_table_md(data, ch))
    md(
        "Red and green both show no spatial-pattern difference between WT and KO "
        "(non-significant, small effect sizes)."
    )

    md(
        "## Quality control and exclusions\n\n"
        "Green sections were assessed by eye and excluded for low signal-to-noise or "
        "tissue damage. The exclusions follow the data structure and standard practice "
        "in the field.\n\n"
        "**Excluded (green):** WT 927, 928, 929, 930, 931, 932; KO 937. "
        "**Excluded (red):** KO 937.\n\n"
        "The green (CTB-488) channel has low signal-to-noise; many green sections were "
        "rejected on quality. The montages below show the section quality directly."
    )

    # Full dataset figures
    md("## Figures: full dataset")
    for section, heading in SECTION_ORDER:
        paths = find_pngs(FULL_RUN, section)
        if paths:
            md(f"**{heading}**")
            cells.append(gallery_cell(paths))

    # Outliers-removed figures (curated)
    md("## Figures: outliers removed")
    for section in CURATED_SECTIONS:
        paths = find_pngs(CURATED_RUN, section)
        if paths:
            md(f"**{CURATED_HEADINGS[section]}**")
            cells.append(gallery_cell(paths))

    md(
        "## Notes\n\n"
        "- A non-significant test indicates no difference detected, not proof of "
        "equivalence.\n"
        "- Red gives the same no-difference result with and without exclusions."
    )

    md("---\n*Generated by `notebook_source/build_notebook.py` from the committed "
       "`data/` stats and figures. Code MIT; data CC BY 4.0. See `LICENSE.md`.*")

    nb = nbf.v4.new_notebook()
    nb["cells"] = cells
    return nb


def main():
    nb = build()
    path = os.path.join(REPO, "RHT_Adult_Analysis.ipynb")
    nbf.write(nb, path)
    print("wrote", os.path.relpath(path))


if __name__ == "__main__":
    main()
