# Dynamic imbalances in cell-type-specific striatal ensembles reflect learned coupling between trajectory representations and locomotor dynamics
doi: https://doi.org/10.1101/2024.10.29.620847

MATLAB code that generates the figure panels and the statistics / source-data
workbook for the paper's main Figures 1-6 and Extended Data Figures 2–9

## Requirements

- MATLAB (tested on R2023b)
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox
- Parallel Computing Toolbox

## Usage

1. Open `runAllFigures.m`.
2. In the **Inputs** section, set `dataRoot` to the processed-data folder and
   `cfg.outDir` to where outputs should go. All other data paths are built
   from `dataRoot`; edit them individually if your layout differs.
3. Optionally trim `figsToRun` to generate a subset.
4. Run the script.

The script adds the repository folders to the MATLAB path, checks that every
input exists, starts a thread-based parallel pool and runs each figure in turn.

## Outputs

- `<outDir>/<Figure>/` — the panels of each figure
- `<outDir>/PaperStats_report.xlsx` — one tab per figure with the plotted
  values, statistics and sample sizes for each panel

## Repository layout

| Folder / file | Contents |
|---|---|
| `runAllFigures.m` | Entry point; all inputs are set here |
| `figures/` | One function per figure (`Figure1`–`Figure5`, `Sup2`–`Sup9`), each taking the `cfg` struct |
| `figTools/` | Figure-specific data loading and selection (two-track, infinite-track, novel/familiar datasets; field width) |
| `reportTools/` | Writing statistics and source data to the workbook |
| `commonFunction/` | Shared analysis and plotting functions (table building, binning, triggered averages, mixed-model comparisons, plotting) |

## Inputs

| `cfg` field | Used by | Content |
|---|---|---|
| `track1w`, `gui` | Fig 1–2, Sup 2–5 | Linear track (1 world) and spontaneous locomotion sessions |
| `speedTuning1w`, `speedTuningGui` | Fig 1, Sup 2 | Speed-tuning correlations |
| `roiTrack`, `roiGui` | Sup 3 | Summary ROI folders (cells per cohort) |
| `prctActive1w`, `prctActiveGui` | Sup 4 | Percent-active-cell data |
| `tuned2w`, `sc2w`, `switchInfo` | Fig 3, Sup 6–7 | Two-track tuning subtypes, single-cell data, session switch sheet |
| `trackNovFam` | Fig 4, Sup 8 | All linear-track sessions for novel vs. familiar |
| `infPop`, `infTuned`, `infSc` | Fig 5, Sup 9 | Infinite track: population, tuning subtypes, single-cell data |
