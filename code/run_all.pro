;+
; RUN_ALL [, FROM=n] [, TO=n] [, /FIGURES_ONLY]
;
; Produces the data for every figure in Ray et al. (2026) and draws the figures.
; Each step reads the original data files (mean_age_data_dir()) and/or the NetCDF
; files written by earlier steps (ma_output_dir()/intermediate/), and writes the
; data behind its figures to ma_output_dir()/figure_data/.
;
;   1  step1_equivalent_latitude        N2O equivalent latitudes and gridded fields
;   2  step2_seasonal_grids             seasonal 1990s mean age and N2O grids
;   3  step3_time_series_and_trends     gridded time series, relationships, trend profiles  (Supp. Figs 3, 5)
;   4  step4_balloon_profiles           balloon profile quality control, flask averages     (ED Figs 2, 3)
;   5  step5_observed_trends            combined time series and trends                     (ED Figs 4, 5; Supp. Figs 8-10)
;   6  step6_models_and_main_figures    ACE-FTS, CCMI-2022 and TLP comparisons              (Figs 1-4; ED Figs 1, 6; Supp. Fig. 4)
;   7  step7_normalization_and_periods  normalization and period comparisons                (Supp. Figs 1, 2, 6, 7)
;   8  make_figures                     all figures from figure_data/*.nc
;
; Keywords:
;   FROM, TO      run only steps FROM..TO (default 1..8); later steps need the outputs of earlier ones
;   FIGURES_ONLY  same as FROM=8
;-
pro run_all, FROM=from, TO=to, FIGURES_ONLY=figures_only

  compile_opt idl2

  steps = ['step1_equivalent_latitude','step2_seasonal_grids','step3_time_series_and_trends','step4_balloon_profiles', $
    'step5_observed_trends','step6_models_and_main_figures','step7_normalization_and_periods','make_figures']
  if n_elements(from) eq 0 then from = 1
  if n_elements(to) eq 0 then to = n_elements(steps)
  if keyword_set(figures_only) then from = n_elements(steps)

  t0 = systime(1)
  print, 'Data directory:   ', mean_age_data_dir()
  print, 'Output directory: ', ma_output_dir()
  for i = from-1, to-1 do begin
    t1 = systime(1)
    print, '' & print, '=== ', steps[i], '  (', systime(), ')'
    call_procedure, steps[i]
    print, steps[i], ' took ', string((systime(1)-t1)/60., format='(f6.1)'), ' minutes'
  endfor
  print, '' & print, 'run_all finished in ', string((systime(1)-t0)/60., format='(f6.1)'), ' minutes'

end
