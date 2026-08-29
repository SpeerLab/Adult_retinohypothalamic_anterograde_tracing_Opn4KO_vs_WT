# Adult RHT analysis

### RHT axon distribution in the SCN: wild-type vs *Opn4* knockout

Comparative analysis of retinohypothalamic tract (RHT) innervation of the
suprachiasmatic nucleus (SCN) in adult wild-type and melanopsin (*Opn4*) knockout
mice. Each animal received a dual-eye injection of spectrally distinct CTB tracers,
so the **red** and **green** channels report the projections of the two eyes
independently within the same section. Each aligned section is collapsed to a 1D
mediolateral intensity profile; WT and KO are compared per SCN half.

**Design:** WT vs KO, adult mice, one coronal SCN section per animal, two
fluorescence channels. Analyzed as single optical sections. Included N (QC-excluded
primary analysis):

| Channel | WT | KO |
|---|---|---|
| Red   | 11 | 7 |
| Green | 5  | 7 |

## Result

No spatial-pattern difference in RHT innervation between WT and KO, in either channel
(red and green both non-significant, small effect sizes).

Analysis is on single optical sections; the reported statistics are the dataset with
outliers removed (green sections excluded for low SNR or tissue damage). The
spatial-pattern (normalized) profiles are the analysis metric.

## Start here

Open `RHT_Adult_Analysis.ipynb`. It reports the primary single-section QC-excluded
statistics and shows the full single-section image set plus a curated single-section
set.

To build the notebook from the committed stats files:

```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt
python notebook_source/build_notebook.py    # writes both notebooks
```

The notebooks read only from `data/` and perform no image processing.

### Figures

The notebook embeds the PNG figures produced by the MATLAB pipeline, committed under
`data/figures/`. The full single-section dataset is shown in full; the QC-excluded run
is shown as a curated set (montages, aligned profiles, heatmaps). Re-running
`build_notebook.py` regenerates the notebook from the committed stats and figures.

## Coordinate convention

After alignment, **x = 0 is the SCN midline**; negative x is the left half of the
image, positive x the right half.

Left and right are *image coordinates*. Ipsilateral and contralateral are defined
relative to the *injected eye*. Because the two eyes were injected with different
tracers, each eye's ipsilateral SCN lobe lies on its own side, so the mapping from
image half to laterality is channel-specific:

| Channel | Ipsilateral half | Contralateral half |
|---------|------------------|--------------------|
| Red     | x < 0            | x > 0              |
| Green   | x > 0            | x < 0              |

This mapping is declared once - `config.ipsiHalf` in the MATLAB pipeline - and every
figure, table and statistic derives its laterality labels from it.

## Quality control

Green sections were assessed by eye and excluded for low signal-to-noise or tissue
damage. The exclusions follow the data structure and standard practice in the field.

- **Green excluded:** WT 927, 928, 929, 930, 931, 932; KO 937
- **Red excluded:** KO 937

The green (CTB-488) channel has low signal-to-noise, and many green sections were
rejected on quality. The notebook displays the section montages so the section
quality can be seen directly.

## Contents

```
├── RHT_Adult_Analysis.ipynb           Results notebook
├── data/                              Committed stats + figures read by the notebook
│   ├── stats_single_section_all_data.mat
│   ├── stats_single_section_outliers_removed.mat
│   └── figures/                       PNG figures embedded in the notebook
├── matlab/                            The MATLAB pipeline that produced the stats
│   ├── adult_RHT_analysis_single_section_outliers_removed.m
│   └── adult_RHT_analysis_single_section_ALL_DATA.m
├── notebook_source/build_notebook.py  Regenerates the notebook from data/
├── requirements.txt                   To run the notebook
└── requirements-dev.txt               To regenerate the notebook
```

## The pipeline

Each MATLAB script runs the same three steps and differs only in whether QC
exclusions are applied (`config.exclude`):

1. **Normalization** - linear histogram stretch (1st-99th percentile) per genotype
   and channel.
2. **1D profiling** - each image collapsed to a mediolateral intensity profile
   (80 px bin for alignment, 20 px bin for analysis).
3. **Alignment and statistics** - the midline valley is detected and each section
   aligned to a common origin; each animal contributes one mean value per SCN half;
   WT vs KO compared by Mann-Whitney U with Cohen's d and post-hoc power.

The `_ALL_DATA` and `_outliers_removed` scripts are identical except for
`config.exclude` (and the output folder), so the two runs are directly comparable.
Input/output paths in the scripts use a generic `C:\Users\user\RHT_project\...` root;
set them to your own data location before running.

## Citation

See `CITATION.cff`. Code is MIT licensed; the contents of `data/` are CC BY 4.0.
See `LICENSE.md`.
