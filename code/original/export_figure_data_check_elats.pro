; Included in check_elats.pro (after the trend profiles are saved) with @export_figure_data_check_elats.
;
; Writes the data plotted in Supplementary Figures 3 and 5 to NetCDF in
; <output>/figure_data/. The plotting code that drew these figures before the
; revision (as Supplement Figs. 1 and 3) is at the end of n2o_elats.pro.

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
m_ex__ = 7       ; Supp. Fig. 3: August
z_ex__ = 44      ; Supp. Fig. 3: altitude bin alt_grid[44] = 17 km
z2_ex__ = 28     ; Supp. Fig. 5d: alt_grid2[28] = 20 km
z3_ex__ = 14     ; Supp. Fig. 5a-c: alt_grid3[14] = 20 km

; ---------- Supplementary Figure 3: N2O equivalent latitude example (AirCore, August, one altitude bin)
gd__ = where(months_all[ico2] eq m_ex__+1 and finite(lat[ico2]) and alt[ico2] ge alt_grid[z_ex__]-dz/2. and alt[ico2] lt alt_grid[z_ex__]+dz/2., ngd__)
sel__ = ico2[gd__]
v__ = orderedhash()
v__['alt_bin'] = ma_var([alt_grid[z_ex__]-dz/2., alt_grid[z_ex__]+dz/2.], ['bin_edge'], 'km', 'altitude bin of the example')
v__['month'] = ma_var(m_ex__+1, '', '', 'calendar month of the example (all years after 2015)')
v__['swoosh_lat'] = ma_var(slat, ['lat'], 'degrees_north')
v__['swoosh_n2o_norm_seasonal'] = ma_var(reform(combo_n2o_seas_hires[*,z_ex__,m_ex__]), ['lat'], '', $
  'SWOOSH (with WACCM below 15 km) monthly mean normalized N2O latitudinal gradient used for the N2O equivalent latitude')
v__['n2o_norm'] = ma_var(n2o_norm[sel__], ['obs'], '', 'AirCore normalized N2O')
v__['lat'] = ma_var(lat[sel__], ['obs'], 'degrees_north', 'flight latitude')
v__['equiv_lat_merra2'] = ma_var(elat_m[sel__], ['obs'], 'degrees_north', 'equivalent latitude from MERRA-2 PV')
v__['equiv_lat_n2o'] = ma_var(elat_m_adj_b_late[sel__], ['obs'], 'degrees_north', 'N2O equivalent latitude (this study)')
ma_nc_write, fig_dir__+'supp_fig3.nc', v__, global=hash('title', 'Supplementary Figure 3: N2O equivalent latitude example')

; ---------- Supplementary Figure 5: NH midlatitude mean age, sampled latitude and latitude adjustment at 20 km
v__ = orderedhash()
v__['years'] = ma_var(years, ['time'], 'year', 'decimal year (bins of 1/48 year)')
v__['alt'] = ma_var(alt_grid3[z3_ex__], '', 'km', 'altitude of panels a-c')
v__['age'] = ma_var(reform(age_nh_coarse[2,z3_ex__,*]), ['time'], 'years', '30-60N average mean age, unadjusted')
v__['age_lat_adj'] = ma_var(reform(age_nh_coarse_lat_adj[2,z3_ex__,*]), ['time'], 'years', '30-60N average mean age, latitude adjusted')
v__['age_lat_adj_uncert'] = ma_var(reform(age_nh_coarse_lat_adj[3,z3_ex__,*]), ['time'], 'years', 'uncertainty of the weighted mean')
v__['sample_type'] = ma_var(fix(reform(age_nh_coarse[5,z3_ex__,*])), ['time'], '', 'platform', flag_values='1 2 3', flag_meanings='aircraft in_situ_balloon aircore')
v__['fit'] = ma_var(reform(age_nh_trends_coarse[z3_ex__,*]), ['coef'], '', 'linear fit to age (intercept years, slope years/year)')
v__['fit_lat_adj'] = ma_var(reform(age_nh_lat_adj_trends_coarse[z3_ex__,*]), ['coef'], '', 'linear fit to age_lat_adj')
v__['rmse'] = ma_var(rse_orig[z3_ex__], '', 'years', 'RMSE of fit')
v__['rmse_lat_adj'] = ma_var(rse_lat_adj[z3_ex__], '', 'years', 'RMSE of fit_lat_adj')
v__['avg_sampled_lat'] = ma_var(reform(lat_avg_nh_coarse[z3_ex__,*]), ['time'], 'degrees_north', 'average sampled latitude within 30-60N')
v__['clim_lat'] = ma_var(lat_grid, ['lat'], 'degrees_north')
v__['clim_alt'] = ma_var(alt_grid2[z2_ex__], '', 'km', 'altitude of panel d')
v__['clim_age_rel_45n'] = ma_var(reform(age_grid_seas_avg[0,*,z2_ex__,*]) - rebin(reform(age_grid_seas_avg[0,26,z2_ex__,*],1,4), n_elements(lat_grid), 4), $
  ['lat','season'], 'years', 'climatological seasonal mean age relative to 45N', seasons='DJF MAM JJA SON')
ma_nc_write, fig_dir__+'supp_fig5.nc', v__, global=hash('title', 'Supplementary Figure 5: NH midlatitude mean ages, latitude sampling and adjustments')
