;+
; RUN_ORIGINAL [, /CACHED] [, /BUILD_PROFILES] [, /FIGURES_ONLY]
;
; Runs the full analysis for Ray et al. (2026) in order, then draws the figures.
;
;   1. check_elats, stage='elats'        N2O equivalent latitudes and gridded fields
;   2. make_balloon_seas_grid            seasonal 1990s mean age and N2O grids
;   3. check_elats, stage='time_series'  gridded time series, relationships and trends
;                                        (+ data for Supp. Figs 3 and 5)
;   4. age_time_series                   main-text and Extended Data figure data
;                                        (+ data for Supp. Figs 4, 8, 9 and 10)
;   5. make_engel_merge_file             data for Supp. Figs 1, 2, 6 and 7
;   6. make_figures                      all figures from output/figure_data/*.nc
;
; Intermediate files are written to the data directory (mean_age_data_dir());
; NetCDF figure data and figures go to the output directory (ma_output_dir()).
;
; Keywords:
;   CACHED          skip steps 1-3 and use the intermediate files already in the data directory
;   BUILD_PROFILES  also rebuild the balloon profile file Engel_profiles_v2025.sav in step 4
;                   (needs the third-party cryo-flask data)
;   FIGURES_ONLY    only run step 6
;-
pro run_all_step, n, name
  print, '' & print, '=== Step ', strtrim(n,2), ': ', name, '  (', systime(), ')'
end

pro run_original, CACHED=cached, BUILD_PROFILES=build_profiles, FIGURES_ONLY=figures_only

  compile_opt idl2

  t0 = systime(1)
  print, 'Data directory:   ', mean_age_data_dir()
  print, 'Output directory: ', ma_output_dir()

  if ~keyword_set(figures_only) then begin
    if ~keyword_set(cached) then begin
      run_all_step, 1, 'check_elats (equivalent latitudes)' & check_elats, stage='elats'
      run_all_step, 2, 'make_balloon_seas_grid'             & make_balloon_seas_grid
      run_all_step, 3, 'check_elats (time series, trends)'  & check_elats, stage='time_series'
    endif
    run_all_step, 4, 'age_time_series' & age_time_series, build_profiles=build_profiles
    run_all_step, 5, 'make_engel_merge_file' & make_engel_merge_file
  endif
  run_all_step, 6, 'make_figures' & make_figures

  print, 'run_original finished in ', string((systime(1)-t0)/60., format='(f6.1)'), ' minutes'

end

