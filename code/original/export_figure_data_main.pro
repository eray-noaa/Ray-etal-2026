; Included at the end of age_time_series.pro with @export_figure_data_main.
;
; Writes the data plotted in each main-text and Extended Data figure to its
; own NetCDF file in <output>/figure_data/. The plotting routines in
; code/figures/ read only these files.

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
flask_dates__ = dates[i_flask]
flask__ = ma_flatten_profiles(bal, flask_dates__, $
  ['alt','n2o_norm','n2o_norm_q','ch4_norm','ch4_norm_q','f12_norm','f12_norm_q','co2_age','co2_age_q','sf6_age','sf6_age_q'], prefix='flask_')
flask__['flask_alt','units'] = 'km'
foreach k__, ['flask_co2_age','flask_co2_age_q','flask_sf6_age','flask_sf6_age_q'] do flask__[k__,'units'] = 'years'
flask__['flask_co2_age','long_name'] = 'cryo-flask CO2 mean age'
flask__['flask_sf6_age','long_name'] = 'cryo-flask SF6 mean age'
flask__['flask_co2_age_q','long_name'] = 'cryo-flask CO2 mean age after quality control (NaN = removed)'
flask__['flask_sf6_age_q','long_name'] = 'cryo-flask SF6 mean age after quality control (NaN = removed)'
foreach k__, ['n2o','ch4','f12'] do begin
  flask__['flask_'+k__+'_norm','long_name'] = 'cryo-flask normalized '+k__
  flask__['flask_'+k__+'_norm_q','long_name'] = 'cryo-flask normalized '+k__+' after quality control (NaN = removed)'
endforeach
stat5__ = ['weighted_mean','uncertainty_of_weighted_mean','mean','standard_deviation','n_measurements']

; Tropical Leaky Pipe model curves (base run fit to in situ data, and w+20%).
tlp__ = orderedhash()
foreach k__, ['n2o','ch4','f12'] do begin
  nend__ = (k__ eq 'ch4') ? [-18,-14] : [-20,-16]
  tlp__['tlp_base_mean_age_'+k__] = ma_var(reform(tlp_low[run_key_low[26]+'_mean_age',2,0:nend__[0]]), ['tlp_base_'+k__+'_level'], 'years', $
    'TLP base run NH mean age (curve vs '+k__+')')
  tlp__['tlp_base_'+k__+'_norm'] = ma_var(reform(tlp_low[run_key_low[26]+'_'+k__,2,0:nend__[0]]), ['tlp_base_'+k__+'_level'], '', 'TLP base run NH normalized '+k__)
  tlp__['tlp_w20_mean_age_'+k__] = ma_var(reform(tlp[run_key[47]+'_mean_age',2,0:nend__[1]]), ['tlp_w20_'+k__+'_level'], 'years', $
    'TLP run with 20% faster tropical upwelling: NH mean age (curve vs '+k__+')')
  tlp__['tlp_w20_'+k__+'_norm'] = ma_var(reform(tlp[run_key[47]+'_'+k__,2,0:nend__[1]]), ['tlp_w20_'+k__+'_level'], '', 'TLP w+20% run NH normalized '+k__)
endforeach

; ---------- Figure 1 and Extended Data Figure 1: mean age vs normalized tracers
foreach tr__, ['n2o','ch4','f12'] do begin
  v__ = orderedhash()
  v__['norm_grid_hist'] = ma_var(norm_grid2, ['norm_hist'], '', 'normalized tracer bin centres for the 1990s in situ distribution')
  v__['age_grid_hist'] = ma_var((tr__ eq 'n2o') ? age_grid_h : age_grid_hist, ['age_hist'], 'years', 'mean age bin centres for the 1990s in situ distribution')
  hist__ = (tr__ eq 'n2o') ? mean_age_all_on_n2o_hist : ((tr__ eq 'ch4') ? mean_age_all_on_ch4_hist : mean_age_all_on_f12_hist)
  v__['insitu_1990s_hist'] = ma_var(hist__, ['age_hist','norm_hist'], '', $
    'fraction of 1990s in situ measurements in each normalized-tracer bin falling in each mean age bin (shaded where >= 0.03)')
  v__['norm_grid2'] = ma_var(norm_grid2, ['norm2'], '', 'normalized tracer bin centres (interval 0.02)')
  a0__ = (tr__ eq 'n2o') ? mean_age_all_on_n2o_a0 : ((tr__ eq 'ch4') ? mean_age_all_on_ch4_a0 : mean_age_all_on_f12_a0)
  b0__ = (tr__ eq 'n2o') ? mean_age_all_on_n2o_b0 : ((tr__ eq 'ch4') ? mean_age_all_on_ch4_b0 : mean_age_all_on_f12_b0)
  v__['aircraft_1990s_mean_age'] = ma_var(a0__, ['norm2','stat5'], 'years', '1990s in situ aircraft mean age in each normalized-tracer bin', $
    statistics='columns: '+strjoin(stat5__,', ')+'; the figure plots column 2 with column 3 as error bar')
  v__['balloon_1990s_mean_age'] = ma_var(b0__, ['norm2','stat5'], 'years', '1990s in situ balloon mean age in each normalized-tracer bin', $
    statistics='columns: '+strjoin(stat5__,', '))
  v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized tracer bin centres (interval 0.05)')
  fl__ = (tr__ eq 'n2o') ? age_q_combo_n2o_bin_fl_avg : ((tr__ eq 'ch4') ? age_q_combo_ch4_bin_fl_avg : age_q_combo_f12_bin_fl_avg)
  v__['flask_avg_mean_age'] = ma_var(fl__, ['norm','stat5'], 'years', 'cryo-flask (1970s-2000s) mean age averaged in each normalized-tracer bin', $
    statistics='columns: '+strjoin(stat5__,', '))
  v__ += flask__
  foreach k__, tlp__.keys() do if strpos(k__, '_'+tr__) ge 0 then v__[k__] = tlp__[k__]
  fname__ = (tr__ eq 'n2o') ? 'fig1.nc' : ((tr__ eq 'ch4') ? 'ed_fig1a.nc' : 'ed_fig1b.nc')
  ma_nc_write, fig_dir__+fname__, v__, global=hash('title', 'Mean age vs normalized '+strupcase(tr__)+' (' + file_basename(fname__,'.nc') + ')')
endforeach

; ---------- Figure 2: NH midlatitude mean age and N2O trend profiles
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt_obs'], 'km', 'altitude bin centres for observed trends (1 km bins)')
v__['insitu_age_trend'] = ma_var(reform(age_nh_lat_adj_trends_coarse[*,1]), ['alt_obs'], 'years/year', 'in situ 30-60N mean age trend (latitude adjusted)')
v__['insitu_age_trend_sigma'] = ma_var(reform(age_nh_lat_adj_trends_coarse_sigma[*,1]), ['alt_obs'], 'years/year', '1-sigma uncertainty of insitu_age_trend')
v__['insitu_n2o_norm_trend'] = ma_var(reform(n2o_nh_lat_adj_trends_coarse[*,1]), ['alt_obs'], '1/year', 'in situ 30-60N normalized N2O trend (latitude adjusted)')
v__['insitu_n2o_norm_trend_sigma'] = ma_var(reform(n2o_nh_lat_adj_trends_coarse_sigma[*,1]), ['alt_obs'], '1/year', '1-sigma uncertainty of insitu_n2o_norm_trend')
v__['ace_alt'] = ma_var(ace_alts, ['alt_ace'], 'km', 'ACE-FTS altitude grid')
v__['ace_n2o_trend'] = ma_var(reform(ace_n2o_nh_trends[1,*]), ['alt_ace'], '1/year', 'ACE-FTS v4.1 30-60N normalized N2O trend 2004-2025')
v__['ace_n2o_trend_sigma'] = ma_var(ace_n2o_nh_trends_chi, ['alt_ace'], '1/year', '1-sigma uncertainty of ace_n2o_trend')
v__['ace_age_trend_alt'] = ma_var([15.5,18.5], ['ace_age'], 'km', 'altitudes of published ACE-FTS mean age trends (Saunders et al.; values entered by hand in age_time_series.pro)')
v__['ace_age_trend'] = ma_var([-0.014,-0.010], ['ace_age'], 'years/year', 'published ACE-FTS 40-50N mean age trends 2004-2021')
v__['ace_age_trend_sigma'] = ma_var([0.006,0.005], ['ace_age'], 'years/year', '1-sigma uncertainty of ace_age_trend')
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km', 'CCMI-2022 altitude (-7 ln(p/1000 hPa))')
foreach d__, ['d1','d2'] do begin
  v__['ccmi_ref'+d__+'_age_trend_avg'] = ma_var((d__ eq 'd1') ? aoa_trends_d1_avg : aoa_trends_d2_avg, ['alt_ccmi'], 'years/year', 'CCMI-2022 ref'+strupcase(d__)+' multi-model mean 30-60N mean age trend')
  v__['ccmi_ref'+d__+'_age_trend_range'] = ma_var((d__ eq 'd1') ? aoa_trends_d1_minmax : aoa_trends_d2_minmax, ['alt_ccmi','minmax'], 'years/year', 'CCMI-2022 ref'+strupcase(d__)+' min and max over ensemble members and periods')
  v__['ccmi_ref'+d__+'_n2o_norm_trend_avg'] = ma_var((d__ eq 'd1') ? n2o_norm_trends_d1_avg : n2o_norm_trends_d2_avg, ['alt_ccmi'], '1/year', 'CCMI-2022 ref'+strupcase(d__)+' multi-model mean 30-60N normalized N2O trend')
  v__['ccmi_ref'+d__+'_n2o_norm_trend_range'] = ma_var((d__ eq 'd1') ? n2o_norm_trends_d1_minmax : n2o_norm_trends_d2_minmax, ['alt_ccmi','minmax'], '1/year', 'CCMI-2022 ref'+strupcase(d__)+' min and max')
endforeach
ma_nc_write, fig_dir__+'fig2.nc', v__, global=hash('title', 'Figure 2: NH mean age and N2O trend profiles', 'note', 'Trends are per year; the figure plots them per decade (x10).')

; ---------- Figure 3
z3a__ = 12
v__ = orderedhash()
v__['years'] = ma_var(years, ['time'], 'year', 'decimal year (bins of 1/48 year)')
v__['n2o_norm_bin'] = ma_var([norm_grid[z3a__]-dn/2., norm_grid[z3a__]+dn/2.], ['bin_edge'], '', 'normalized N2O range of panel a')
v__['age_in_bin'] = ma_var(reform(age_on_n2o_adj_tseries[2,z3a__,*]), ['time'], 'years', 'mean age averaged in the panel a normalized N2O bin')
v__['age_in_bin_sd'] = ma_var(reform(age_on_n2o_adj_tseries[3,z3a__,*]), ['time'], 'years', 'standard deviation of age_in_bin')
v__['sample_type'] = ma_var(fix(reform(age_on_n2o_adj_tseries[5,z3a__,*])), ['time'], '', 'platform', flag_values='1 2 3', flag_meanings='aircraft in_situ_balloon aircore')
v__['n2o_from_ch4'] = ma_var(fix(n2o_from_ch4_yrs), ['time'], '', '1 where AirCore N2O was derived from a CH4-N2O relationship')
v__['bin_trend_fit'] = ma_var(reform(n2o_nofl_age_trends[z3a__,*]), ['coef'], '', 'linear fit (intercept years, slope years/year) to age_in_bin')
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres (interval 0.05)')
v__['insitu_age_trend_on_n2o'] = ma_var(reform(n2o_nofl_age_trends[*,1]), ['norm'], 'years/year', 'in situ NH extratropical mean age trend in normalized N2O bins')
v__['insitu_age_trend_on_n2o_sigma'] = ma_var(reform(n2o_nofl_age_trends_sigma[*,1]), ['norm'], 'years/year', '1-sigma uncertainty')
v__['ccmi_n2o_norm_grid'] = ma_var(n2o_norm_grid_model, ['norm_ccmi'], '', 'normalized N2O grid of CCMI trends')
v__['ccmi_refd1_age_trend_on_n2o_avg'] = ma_var(aoa_on_n2o_trends_d1_avg, ['norm_ccmi'], 'years/year', 'CCMI-2022 refD1 multi-model mean')
v__['ccmi_refd1_age_trend_on_n2o_range'] = ma_var(reform(aoa_on_n2o_trends_d1_minmax[*,*,1]), ['norm_ccmi','minmax'], 'years/year', 'CCMI-2022 refD1 min and max')
v__['ccmi_refd2_age_trend_on_n2o_avg'] = ma_var(aoa_on_n2o_trends_d2_avg, ['norm_ccmi'], 'years/year', 'CCMI-2022 refD2 multi-model mean')
v__['ccmi_refd2_age_trend_on_n2o_range'] = ma_var(aoa_on_n2o_trends_d2_minmax, ['norm_ccmi','minmax'], 'years/year', 'CCMI-2022 refD2 min and max')
v__['norm_grid2'] = ma_var(norm_grid2, ['norm2'], '', 'normalized N2O bin centres (interval 0.02)')
v__['age_on_n2o_1990s'] = ma_var(mean_age_on_n2o_90s, ['norm2'], 'years', 'NH extratropical mean age vs normalized N2O, 1990s')
v__['age_on_n2o_2020s'] = ma_var(mean_age_on_n2o_20s_from_trend, ['norm2'], 'years', 'NH extratropical mean age vs normalized N2O, 2020s (1990s plus trend)')
v__['arrow_alt'] = ma_var(alt_grid2[[11,13,14,15,17,19]], ['arrow'], 'km', 'altitudes of the constant-altitude change arrows in panel c')
v__['arrow_age'] = ma_var(transpose([[reform(age_alt_profile_avg[0,[11,13,14,15,17,19]])],[age_alt_profile_from_trend[[11,13,14,15,17,19]]]]), ['start_end','arrow'], 'years')
v__['arrow_n2o_norm'] = ma_var(transpose([[reform(n2o_norm_alt_profile_avg[0,[11,13,14,15,17,19]])],[n2o_norm_alt_profile_from_trend[[11,13,14,15,17,19]]]]), ['start_end','arrow'], '')
ma_nc_write, fig_dir__+'fig3.nc', v__, global=hash('title', 'Figure 3: mean age trends as a function of N2O')

; ---------- Figure 4 (computed in age_time_series.pro just above)
v__ = orderedhash()
v__['tlp_alt'] = ma_var(alt_tlp, ['alt_tlp'], 'km', 'TLP model altitude')
v__['tlp_w_trend_range'] = ma_var(diffs_minmax/3., ['alt_tlp','minmax'], '%/decade', 'range of TLP tropical upwelling trends within tolerance of the observed changes')
v__['tlp_w_trend_wmean'] = ma_var(wmean/3., ['alt_tlp'], '%/decade', 'weighted mean TLP tropical upwelling trend (figure applies a 3-point smooth)')
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km')
v__['ccmi_refd1_w_trend_avg'] = ma_var(10.*w_turn_trends_per_d1_avg, ['alt_ccmi'], '%/decade', 'CCMI-2022 refD1 tropical mean upwelling trend')
v__['ccmi_refd1_w_trend_range'] = ma_var(10.*w_turn_trends_per_d1_minmax, ['alt_ccmi','minmax'], '%/decade')
v__['ccmi_refd2_w_trend_avg'] = ma_var(10.*w_turn_trends_per_d2_avg, ['alt_ccmi'], '%/decade', 'CCMI-2022 refD2 tropical mean upwelling trend')
v__['ccmi_refd2_w_trend_range'] = ma_var(10.*w_turn_trends_per_d2_minmax, ['alt_ccmi','minmax'], '%/decade')
ma_nc_write, fig_dir__+'fig4.nc', v__, global=hash('title', 'Figure 4: tropical upwelling trends')

; ---------- Extended Data Figures 2 and 3: profiles and tracer correlations
v__ = orderedhash()
v__['alt_grid'] = ma_var(alt_grid, ['alt'], 'km', 'altitude grid of the OMS in situ balloon profiles')
pct__ = ', columns are percentiles; the figure uses column 3 (median) and columns 0 and 6 (range)'
foreach k__, ['n2o_norm','ch4_norm','f12_norm','age'] do begin
  case k__ of
    'n2o_norm': begin & tr__ = alt_n2o_norm_avg_tr_percent_oms & ml__ = alt_n2o_norm_avg_percent_oms & vx__ = alt_n2o_norm_vx_oms & end
    'ch4_norm': begin & tr__ = alt_ch4_norm_avg_tr_percent_oms & ml__ = alt_ch4_norm_avg_percent_oms & vx__ = alt_ch4_norm_vx_oms & end
    'f12_norm': begin & tr__ = alt_f12_norm_sf6_tr_percent_oms & ml__ = alt_f12_norm_sf6_percent_oms & vx__ = alt_f12_norm_vx_oms & end
    'age':      begin & tr__ = alt_age_avg_tr_percent_oms & ml__ = alt_age_avg_percent_oms & vx__ = alt_age_vx_oms & end
  endcase
  u__ = (k__ eq 'age') ? 'years' : ''
  v__['oms_tropics_'+k__] = ma_var(tr__, ['alt','pct'], u__, 'OMS in situ balloon tropical profile'+pct__)
  v__['oms_midlat_'+k__] = ma_var(ml__, ['alt','pct'], u__, 'OMS in situ balloon midlatitude profile'+pct__)
  v__['oms_vortex_'+k__] = ma_var(vx__, ['alt','pct'], u__, 'OMS in situ balloon vortex profile'+pct__)
endforeach
v__['alt_grid2'] = ma_var(alt_grid2, ['alt2'], 'km', '1 km altitude bins for flask averages')
v__['flask_avg_n2o_norm'] = ma_var(n2o_norm_alt_grid_fl_avg, ['alt2','stat5'], '', 'cryo-flask average normalized N2O (columns 0 mean, 1 sd)')
v__['flask_avg_ch4_norm'] = ma_var(ch4_norm_alt_grid_fl_avg, ['alt2','stat5'], '', 'cryo-flask average normalized CH4 (columns 0 mean, 1 sd)')
v__['flask_avg_f12_norm'] = ma_var(f12_norm_alt_grid_fl_avg, ['alt2','stat5'], '', 'cryo-flask average normalized CFC-12 (columns 0 mean, 1 sd)')
v__['flask_avg_co2_age'] = ma_var(co2_age_alt_grid_fl_avg, ['alt2','stat5'], 'years', 'cryo-flask CO2 mean age', statistics='columns: '+strjoin(stat5__,', '))
v__['flask_avg_sf6_age'] = ma_var(sf6_age_alt_grid_fl_avg, ['alt2','stat5'], 'years', 'cryo-flask SF6 mean age', statistics='columns: '+strjoin(stat5__,', '))
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['aircore_ch4_norm_on_n2o'] = ma_var(n2o_norm_ch4_norm_ac_sm_coarse, ['norm'], '', 'AirCore CH4-N2O relationship (Extended Data Fig. 3a green symbols: x = norm_grid, y = this)')
v__ += flask__
ma_nc_write, fig_dir__+'ed_fig2_3.nc', v__, global=hash('title', 'Extended Data Figures 2 and 3: tracer and mean age profiles and correlations')

; ---------- Extended Data Figure 4
v__ = orderedhash()
v__['years'] = ma_var(years, ['time'], 'year', 'decimal year (bins of 1/48 year)')
v__['balloon_max_alt'] = ma_var(max_alt_b_yrs, ['time'], 'km', 'maximum altitude of balloon profiles')
v__['balloon_sample_type'] = ma_var(samp_type_b_yrs, ['time'], '', 'balloon sampling type', flag_values='0 1 2 3', flag_meanings='flask in_situ_balloon aircore flask')
v__['aircraft_years'] = ma_var(n2o_years, ['time_aircraft'], 'year', 'decimal year')
v__['aircraft_max_alt'] = ma_var(max_alts_a, ['time_aircraft'], 'km', 'maximum altitude of aircraft flights')
v__['age_24_35km_with_flask'] = ma_var(bins_all_combo_fl_age_midstrat_nh_years[*,0:1], ['time','mean_unc'], 'years', '24-35 km weighted mean age including flask data (column 0 mean, 1 uncertainty)')
v__['age_24_35km_no_flask'] = ma_var(bins_all_combo_nofl_age_midstrat_nh_years[*,0:1], ['time','mean_unc'], 'years', '24-35 km weighted mean age without flask data')
v__['fit_with_flask'] = ma_var(age_fl_midstrat_linfit, ['coef'], '', 'linear fit (intercept, slope per year) through all data')
v__['fit_no_flask'] = ma_var(age_nofl_midstrat_linfit, ['coef'], '', 'linear fit through in situ and AirCore data')
ma_nc_write, fig_dir__+'ed_fig4.nc', v__, global=hash('title', 'Extended Data Figure 4: sampling altitudes and 24-35 km mean age')

; ---------- Extended Data Figure 5
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt'], 'km')
v__['age_trend_no_flask_lat_adj'] = ma_var(age_nh_lat_adj_trends_coarse[*,1], ['alt'], 'years/year', 'in situ without flasks, latitude adjusted')
v__['age_trend_no_flask_lat_adj_sigma'] = ma_var(age_nh_lat_adj_trends_coarse_sigma[*,1], ['alt'], 'years/year')
v__['age_trend_no_flask'] = ma_var(age_nh_trends_coarse[*,1], ['alt'], 'years/year', 'in situ without flasks')
v__['age_trend_no_flask_sigma'] = ma_var(age_nh_trends_coarse_sigma[*,1], ['alt'], 'years/year')
v__['age_trend_with_flask'] = ma_var(age_fl_alt_trends[*,1], ['alt'], 'years/year', 'in situ with flasks')
v__['age_trend_with_flask_sigma'] = ma_var(age_fl_alt_trends_sigma[*,1], ['alt'], 'years/year')
v__['ray2014_alt'] = ma_var(findgen(4)*5+17.5, ['ray2014'], 'km', 'Ray et al. (2014) trend altitudes (entered by hand)')
v__['ray2014_trend'] = ma_var([-0.07,0.09,0.2,0.22]/10., ['ray2014'], 'years/year')
v__['ray2014_trend_sigma'] = ma_var([0.037,0.075,0.075,0.082]/10., ['ray2014'], 'years/year')
v__['fritsch2020_alt'] = ma_var(28.5, '', 'km', 'Fritsch et al. (2020) trend altitude (entered by hand)')
v__['fritsch2020_trend'] = ma_var(0.0083, '', 'years/year')
v__['fritsch2020_trend_sigma'] = ma_var(0.0156, '', 'years/year')
v__['n_total'] = ma_var(tot_samp, ['alt'], '', 'total number of NH extratropical measurements 1975-2025')
v__['n_aircraft'] = ma_var(tot_samp_aircraft, ['alt'], '')
v__['n_insitu_balloon'] = ma_var(nsamp_oms, ['alt'], '')
v__['n_aircore'] = ma_var(tot_samp_aircore, ['alt'], '')
v__['n_flask'] = ma_var(nsamp_flask, ['alt'], '')
ma_nc_write, fig_dir__+'ed_fig5.nc', v__, global=hash('title', 'Extended Data Figure 5: mean age trend and measurement number profiles')

; ---------- Extended Data Figure 6
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt_obs'], 'km')
v__['insitu_ch4_norm_trend'] = ma_var(ch4_nh_trends_coarse[*,1], ['alt_obs'], '1/year', 'in situ 30-60N normalized CH4 trend')
v__['insitu_ch4_norm_trend_sigma'] = ma_var(ch4_nh_trends_coarse_sigma[*,1], ['alt_obs'], '1/year')
v__['ace_alt'] = ma_var(ace_alts, ['alt_ace'], 'km')
v__['ace_ch4_trend'] = ma_var(reform(ace_ch4_nh_trends[1,*]), ['alt_ace'], '1/year', 'ACE-FTS 30-60N normalized CH4 trend')
v__['ace_ch4_trend_sigma'] = ma_var(ace_ch4_nh_trends_chi, ['alt_ace'], '1/year')
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km')
v__['ccmi_refd1_ch4_norm_trend_avg'] = ma_var(ch4_norm_trends_d1_avg, ['alt_ccmi'], '1/year')
v__['ccmi_refd1_ch4_norm_trend_range'] = ma_var(ch4_norm_trends_d1_minmax, ['alt_ccmi','minmax'], '1/year')
v__['ccmi_refd2_ch4_norm_trend_avg'] = ma_var(ch4_norm_trends_d2_avg, ['alt_ccmi'], '1/year')
v__['ccmi_refd2_ch4_norm_trend_range'] = ma_var(ch4_norm_trends_d2_minmax, ['alt_ccmi','minmax'], '1/year')
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized CH4 bin centres')
v__['insitu_age_trend_on_ch4'] = ma_var(age_on_ch4_trends[*,1], ['norm'], 'years/year')
v__['insitu_age_trend_on_ch4_sigma'] = ma_var(age_on_ch4_trends_sigma[*,1], ['norm'], 'years/year')
v__['ccmi_norm_grid'] = ma_var(n2o_norm_grid_model, ['norm_ccmi'], '')
v__['ccmi_refd1_age_trend_on_ch4_avg'] = ma_var(aoa_on_ch4_trends_d1_avg, ['norm_ccmi'], 'years/year')
v__['ccmi_refd1_age_trend_on_ch4_range'] = ma_var(aoa_on_ch4_trends_d1_minmax, ['norm_ccmi','minmax'], 'years/year')
v__['ccmi_refd2_age_trend_on_ch4_avg'] = ma_var(aoa_on_ch4_trends_d2_avg, ['norm_ccmi'], 'years/year')
v__['ccmi_refd2_age_trend_on_ch4_range'] = ma_var(aoa_on_ch4_trends_d2_minmax, ['norm_ccmi','minmax'], 'years/year')
ma_nc_write, fig_dir__+'ed_fig6.nc', v__, global=hash('title', 'Extended Data Figure 6: CH4 and mean age on CH4 trend profiles')

print, 'Figure data written to ', fig_dir__
