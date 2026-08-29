#!/usr/bin/env python3
"""
build_notebook.py  -  Generate RHT_Adult_Analysis.ipynb.

Single optical section analysis, WT vs Opn4 KO. Reports the outliers-removed
statistics (green sections with low SNR or tissue damage excluded), and shows
figures for the full dataset followed by the outliers-removed dataset.

The figures are EMBEDDED directly into the notebook as base64 markdown
attachments, so the notebook renders fully on GitHub, nbviewer, or a fresh
JupyterLab session WITHOUT the user running any cell and without depending on
the working directory. The embedded copies are downscaled/recompressed (the
full-resolution PNGs remain in data/figures/); this keeps the notebook well
under GitHub's 25 MB render limit.

Reads:
    data/stats_single_section_outliers_removed.mat   (the reported statistics)
    data/figures/single_section_all_data/            (full-dataset figures)
    data/figures/single_section_outliers_removed/    (outliers-removed figures)

Run:
    python build_notebook.py
"""

import base64
import io
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

try:
    from PIL import Image
except ImportError:
    sys.exit("pillow is required (pip install -r requirements-dev.txt)")

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

# Embedded-image settings. Full-res originals are left untouched in
# data/figures/; only the copies baked into the notebook are shrunk.
MAX_EMBED_WIDTH = 1600   # px; wider figures are downscaled to this width
PNG_OPTIMIZE = True

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


def _encode_png(rel_path):
    """Load a PNG, downscale to MAX_EMBED_WIDTH, return base64 PNG bytes."""
    abs_path = os.path.join(REPO, rel_path)
    with Image.open(abs_path) as im:
        im.load()
        if im.width > MAX_EMBED_WIDTH:
            new_h = round(im.height * MAX_EMBED_WIDTH / im.width)
            im = im.resize((MAX_EMBED_WIDTH, new_h), Image.LANCZOS)
        # Palette/greyscale-with-alpha may not save cleanly; normalise those.
        if im.mode in ("P", "LA"):
            im = im.convert("RGBA")
        buf = io.BytesIO()
        im.save(buf, format="PNG", optimize=PNG_OPTIMIZE)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def embed_gallery_cell(paths):
    """A single markdown cell that displays every figure in `paths`, with the
    image bytes embedded as cell attachments (no execution, no file access)."""
    attachments = {}
    body = []
    for rel in paths:
        name = os.path.basename(rel)
        attachments[name] = {"image/png": _encode_png(rel)}
        body.append(f"*{name}*\n\n![{name}](attachment:{name})")
    cell = nbf.v4.new_markdown_cell("\n\n".join(body))
    cell["attachments"] = attachments
    return cell


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
            cells.append(embed_gallery_cell(paths))

    # Outliers-removed figures (curated)
    md("## Figures: outliers removed")
    for section in CURATED_SECTIONS:
        paths = find_pngs(CURATED_RUN, section)
        if paths:
            md(f"**{CURATED_HEADINGS[section]}**")
            cells.append(embed_gallery_cell(paths))

    md(
        "## Notes\n\n"
        "- A non-significant test indicates no difference detected, not proof of "
        "equivalence.\n"
        "- Red gives the same no-difference result with and without exclusions.\n"
        "- Figures shown inline are downscaled for display; the full-resolution PNGs "
        "are in `data/figures/`."
    )

    md("---\n*Generated by `notebook_source/build_notebook.py` from the committed "
       "`data/` stats and figures. Figures are embedded (base64) so the notebook "
       "renders without execution. Code MIT; data CC BY 4.0. See `LICENSE.md`.*")

    nb = nbf.v4.new_notebook()
    nb["cells"] = cells
    return nb


def main():
    nb = build()
    path = os.path.join(REPO, "RHT_Adult_Analysis.ipynb")
    nbf.write(nb, path)
    size_mb = os.path.getsize(path) / 1e6
    print(f"wrote {os.path.relpath(path)} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
