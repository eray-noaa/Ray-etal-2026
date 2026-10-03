; Included at the end of age_time_series.pro with @export_figure_data_supp_ats.
;
; Writes the data plotted in Supplementary Figures 4, 8, 9 and 10 to NetCDF in
; <output>/figure_data/. (In the plotting code above these figures carry their
; pre-revision numbers 2, 4, 5 and 6.)

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()

; ---------- Supplementary Figure 4: 1990s gridded mean age and N2O, 20 km latitude distributions
z20__ = 28
v__ = orderedhash()
v__['equiv_lat'] = ma_var(lat_grid_a, ['lat'], 'degrees_north', 'N2O equivalent latitude bin centres')
v__['alt'] = ma_var(alt_grid, ['alt'], 'km', 'altitude bin centres')
v__['tropopause_alt_on_grid'] = ma_var(tp_alt_grid, ['lat'], 'km', 'climatological tropopause altitude on the equivalent latitude grid')
v__['age_1990s'] = ma_var(age_combo_grid_early, ['stat','lat','alt'], 'years', '1990s gridded mean age', $
  statistics='index 0 mean; index 2 number of seasons sampled (full-size symbols where >= 3, small where 1-2)')
v__['n2o_norm_1990s'] = ma_var(n2o_norm_grid_early, ['stat','lat','alt'], '', '1990s gridded normalized N2O', $
  statistics='index 0 mean; index 2 number of seasons sampled')
v__['merra2_lat'] = ma_var(lat_tp, ['lat_tp'], 'degrees_north')
v__['merra2_tropopause_alt'] = ma_var(tpause_alt, ['lat_tp'], 'km', 'MERRA-2 climatological tropopause altitude (magenta line)')
v__['alt_20km'] = ma_var(alt_grid[z20__], '', 'km', 'altitude of panels c and d')
v__['ccmi_lat'] = ma_var(lat_c, ['lat_ccmi'], 'degrees_north')
v__['ccmi_age_20km_1990s'] = ma_var(reform(age_c_1990s_zobs[*,z20__,*]), ['lat_ccmi','member'], 'years', 'CCMI-2022 refD1 1990s mean age at 20 km, each model member')
v__['ccmi_n2o_norm_20km_1990s'] = ma_var(reform(n2o_norm_c_1990s_zobs[*,z20__,*]), ['lat_ccmi','member'], '', $
  'CCMI-2022 refD1 1990s normalized N2O at 20 km, each model member (the figure omits members 0-2 and 15)')
ma_nc_write, fig_dir__+'supp_fig4.nc', v__, global=hash('title', 'Supplementary Figure 4: 1990s gridded mean age and N2O')

; ---------- Supplementary Figure 8: mean age vs normalized N2O, this study and earlier relationships
; (plotting code for this figure is in make_engel_merge_file.pro; all inputs are available here)
n2o_ppb_engel__ = 313. * n2o_norm_grid
v__ = orderedhash()
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['age_1990s'] = ma_var(bins_a0_n2o_norm_age_nh_coarse, ['norm'], 'years', 'NH mean age vs normalized N2O, 1990s in situ (this study)')
v__['age_2020s'] = ma_var(bins_a2_n2o_norm_age_nh_trend, ['norm'], 'years', 'NH mean age vs normalized N2O, 2020s in situ (this study)')
v__['n2o_norm_grid_fine'] = ma_var(n2o_norm_grid, ['norm_fine'], '', 'normalized N2O grid of the published relationships')
v__['age_andrews2001'] = ma_var(n2o_age_a, ['norm_fine'], 'years', 'Andrews et al. (2001) ER-2 relationship, converted to normalized N2O')
v__['age_engel2002'] = ma_var(6.03 - 0.0136*n2o_ppb_engel__ + 8.5892e-5*n2o_ppb_engel__^2 - 3.38e-7*n2o_ppb_engel__^3, ['norm_fine'], 'years', $
  'Engel et al. (2002) polynomial mean age vs N2O, evaluated at N2O = 313 ppb x normalized N2O')
ma_nc_write, fig_dir__+'supp_fig8.nc', v__, global=hash('title', 'Supplementary Figure 8: mean age vs normalized N2O')

; ---------- Supplementary Figure 9: surface CO2 and the 6 Aug 1994 ASHOE flight
i_ft__ = 29
v__ = orderedhash()
v__['years_surface'] = ma_var(years_sm+2./12., ['time'], 'year', 'decimal year, shifted by the 2-month delay used in the figure')
v__['co2_mlo_smo'] = ma_var(co2_sm, ['time'], 'ppmv', 'Mauna Loa and Samoa average CO2')
v__['co2_mlo_smo_smooth'] = ma_var(co2_sm_smooth, ['time'], 'ppmv', 'smoothed Mauna Loa and Samoa average CO2')
v__['steady_growth_time'] = ma_var(float(co2_steady_growth_time), ['endpoint'], 'year', 'constant growth rate CO2 used in a previous study')
v__['steady_growth_co2'] = ma_var(co2_steady_growth, ['endpoint'], 'ppmv')
v__['flight_date'] = ma_var(flight_track[i_ft__].date, '', '', 'ER-2 flight date (YYYYMMDD)')
v__['flight_date_float'] = ma_var(flight_track[i_ft__].date_float, '', 'year')
v__['flight_co2_adj'] = ma_var(*flight_track[i_ft__].co2_adj, ['obs'], 'ppmv', 'flight CO2 adjusted for CH4 oxidation')
v__['flight_co2_archive_adj'] = ma_var(*flight_track[i_ft__].co2 - (1.725-*flight_track[i_ft__].ch4/1e3), ['obs'], 'ppmv', $
  'flight CO2 with the archived CH4 adjustment (CO2 - (1.725 - CH4/1000))')
v__['flight_co2_age_archive'] = ma_var(*flight_track[i_ft__].co2_age, ['obs'], 'years', 'archived mean age')
v__['flight_co2_lag_age'] = ma_var(*flight_track[i_ft__].co2_lag_age, ['obs'], 'years', 'simple lag-technique mean age (this study)')
ma_nc_write, fig_dir__+'supp_fig9.nc', v__, global=hash('title', 'Supplementary Figure 9: CO2 time series and mean age example')

; ---------- Supplementary Figure 10: archived minus updated CO2 mean age vs N2O, by 1990s mission
mi__ = indgen(5) + 3
v__ = orderedhash()
v__['n2o_norm_grid'] = ma_var(n2o_norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['mission'] = ma_var(mission_names[mi__], ['mission'], '', 'ER-2 mission')
v__['mission_color'] = ma_var(colors[mi__], ['mission'], '', 'plot colour')
v__['age_archived'] = ma_var(reform(bins_co2_a0.n2o_norm_co2_age_m[*,0,mi__]), ['norm','mission'], 'years', 'archived CO2 mean age averaged in normalized N2O bins')
v__['age_updated'] = ma_var(reform(bins_co2_a0.n2o_norm_age_all_m[*,0,mi__]), ['norm','mission'], 'years', 'updated CO2 mean age (this study)')
v__['age_updated_uncert'] = ma_var(reform(bins_co2_a0.n2o_norm_age_all_m[*,1,mi__]), ['norm','mission'], 'years', 'uncertainty of the bin average')
ma_nc_write, fig_dir__+'supp_fig10.nc', v__, global=hash('title', 'Supplementary Figure 10: archived mean age bias vs N2O', $
  'note', 'The blue average line is the 7-point smoothed mean over missions of age_archived - age_updated (computed in supp_fig10.pro).')
