# Ray-etal-2026

Code and data for Ray et al. (2026), *Observed stratospheric mean age decrease consistent with circulation acceleration*, Nature Geoscience, https://doi.org/10.1038/s41561-026-02011-3.

All analysis code is written in IDL (tested with IDL 9.1). The merged aircraft and balloon data are stored as IDL save files; the data behind every figure are written as NetCDF.

## Repository layout

```
code/
  run_all.pro       runs the analysis steps in order, then draws every figure
  steps/            the analysis, one program per step (see below)
  figures/          one plotting routine per figure (fig1.pro ... supp_fig10.pro) and make_figures.pro
  lib/              helper routines (data and output paths, NetCDF read/write, smoothing, gridding)
  original/         the original analysis programs, kept for reference (see below)
data/               data files, in the folder structure the code expects
  MANIFEST.md       every data file, with size, checksum, source and which code uses it
docs/figures/       the published main-text and Extended Data figures (PDF)
output/             written by the code (not in the repository)
  intermediate/     NetCDF files passed between steps
  figure_data/      one NetCDF file per figure
  figures/          PDF and PNG of every figure
```

## Running the code

1. Put `code/` on the IDL path:
   ```idl
   !path = expand_path('+/path/to/Ray-etal-2026/code') + ':' + !path
   ```
2. Gather the data. Small files are in `data/`. Larger files are on the NOAA CSL data site (https://csl.noaa.gov/groups/csl8/modeldata/) and third-party files come from their original sources; `data/MANIFEST.md` lists where each file comes from. Place them under `data/` at the paths given in the manifest, or point the code at another directory:
   ```idl
   setenv, 'MEAN_AGE_DATA_DIR=/path/to/Data/'
   setenv, 'MEAN_AGE_OUTPUT_DIR=/path/to/output/'   ; optional, defaults to output/ in this repository
   ```
3. Run everything (about 10 minutes):
   ```idl
   run_all                 ; all steps and all figures
   run_all, from=4         ; rerun from step 4 onward, using the outputs of steps 1-3
   make_figures            ; only redraw the figures from output/figure_data/*.nc
   ```

## Analysis steps

Each step is a single IDL procedure in `code/steps/`. It reads the original data files and the NetCDF files written by earlier steps, and writes NetCDF.

| Step | Program | What it does | Figure data written |
|---|---|---|---|
| 1 | `step1_equivalent_latitude` | N2O equivalent latitude for every in situ measurement; gridded 1990s and 2010s-2020s fields | |
| 2 | `step2_seasonal_grids` | seasonal 1990s equivalent latitude-altitude grids | |
| 3 | `step3_time_series_and_trends` | latitude-adjusted time series in altitude and N2O/CH4 bins, 1990s mean age vs tracer relationships, trend profiles | Supp. Figs 3, 5 |
| 4 | `step4_balloon_profiles` | balloon profile quality control, binning and 1970s-2000s flask averages | Extended Data Figs 2, 3 |
| 5 | `step5_observed_trends` | combined aircraft and balloon time series and trends, 1990s and 2020s relationships | Extended Data Figs 4, 5; Supp. Figs 8, 9, 10 |
| 6 | `step6_models_and_main_figures` | ACE-FTS, CCMI-2022 and Tropical Leaky Pipe comparisons | Figs 1-4; Extended Data Figs 1, 6; Supp. Fig. 4 |
| 7 | `step7_normalization_and_periods` | normalization and trend-period comparisons | Supp. Figs 1, 2, 6, 7 |

Each figure's data is in `output/figure_data/<name>.nc` (variables carry units and descriptions; trends are stored per year and plotted per decade) and is drawn by `code/figures/<name>.pro`:

| Figure | Data file | Plot routine |
|---|---|---|
| Figs 1-4 | `fig1.nc` ... `fig4.nc` | `fig1.pro` ... `fig4.pro` |
| Extended Data Fig. 1 | `ed_fig1a.nc`, `ed_fig1b.nc` | `ed_fig1.pro` |
| Extended Data Figs 2, 3 | `ed_fig2_3.nc` | `ed_fig2.pro`, `ed_fig3.pro` |
| Extended Data Figs 4-6 | `ed_fig4.nc` ... `ed_fig6.nc` | `ed_fig4.pro` ... `ed_fig6.pro` |
| Supplementary Figs 1-10 | `supp_fig1.nc` ... `supp_fig10.nc` | `supp_fig1.pro` ... `supp_fig10.pro` |

The published figures were finished in other software (several panels in matplotlib), so the IDL versions match their content but not their exact styling.

## Original programs

`code/original/` holds the four programs used for the paper (`check_elats.pro`, `age_time_series.pro`, `make_engel_merge_file.pro`, `n2o_elats.pro`) with only their data paths fixed, plus the code that writes their figure data and `run_original.pro`, which runs them in order. The programs in `code/steps/` were produced from them by keeping only the code that contributes to the paper's figures and splitting it at natural hand-off points. Running both on the same data gives identical figure data. The original programs also need some files that the steps do not (see the "Used by" column in `data/MANIFEST.md`).

The balloon profile file `Balloon/Engel_profiles_v2025.sav` is an input to step 4. It was built by `code/original/age_time_series.pro` (keyword `/BUILD_PROFILES`) from the cryo-flask, OMS, WAS and AirCore profile data.

## Data sources

The merged aircraft and balloon files were built from the original mission data:

- NASA ER-2 missions (ASHOE-MAESA, STRAT, POLARIS, SOLVE) and OMS balloon flights: NASA ESPO archive, https://espoarchive.nasa.gov/archive/browse
- DCOTSS: https://www-air.larc.nasa.gov/missions/dcotss/index.html
- SABRE: https://csl.noaa.gov/projects/sabre/data.html
- NOAA GML AirCore: https://gml.noaa.gov/ccgg/arc/?id=144
- Surface trace gases: NOAA GML Greenhouse Gas Marine Boundary Layer Reference (https://gml.noaa.gov/ccgg/mbl/) and HATS combined data (https://gml.noaa.gov/hats/)

Third-party data used in the paper but not redistributed here:

- CCMI-2022 model output: CEDA / BADC archive, https://data.ceda.ac.uk/badc/ccmi/data/post-cmip6/ccmi-2022
- ACE-FTS satellite data: https://databace.scisat.ca/ (registration required)
- Cryo-flask balloon data (Engel et al.), WACCM, CLaMS and other collaborator-provided files: see `data/MANIFEST.md`

Figure data files that contain values derived from these third-party data are produced locally by the code and are not part of this repository.

## License

The code is released under the MIT License (see `LICENSE`). Data files are provided for reproducibility of the paper; please cite the paper and the original data sources when using them.

## Citation

Ray et al. (2026). Observed stratospheric mean age decrease consistent with circulation acceleration. *Nature Geoscience*. https://doi.org/10.1038/s41561-026-02011-3
