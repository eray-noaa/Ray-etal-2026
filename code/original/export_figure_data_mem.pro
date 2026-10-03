; Included at the end of make_engel_merge_file.pro with @export_figure_data_mem.
;
; Writes the data plotted in Supplementary Figures 1, 2, 6 and 7 to NetCDF in
; <output>/figure_data/. The original plotting code is in the commented blocks
; above that save Mean_age_vs_N2O_norm_compare.png, Mean_age_vs_CH4_norm_compare.png,
; Mean_age_trend_abs_N2O_profiles_norm_compare.png, Mean_age_trend_CH4_profiles_norm_compare.png,
; Mean_age_trend_alt_profiles_compare_periods.png, Mean_age_trend_N2O_profiles_compare_periods.png
; and Mean_age_trend_alt_profiles_models_compare_periods.png.

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
norm_note__ = '"simple" = normalized by the surface mean at the measurement time (used in the paper); ' + $
  '"age_conv" = normalized by the surface time series convolved with the full age spectrum'

; ---------- Supplementary Figure 1: mean age vs N2O (a) and CH4 (b) for two normalizations
; Needs the age-convolution bins from Aircraft_mean_ages_early_bins_comb.sav (post-revision version).
if n_elements(bins_a0_n2o_norm_e_age_nh) gt 0 && n_elements(bins_a0_ch4_norm_e_age_nh) gt 0 then begin
  v__ = orderedhash()
  v__['norm_grid'] = ma_var(n2o_norm_grid, ['norm'], '', 'normalized tracer bin centres')
  v__['n2o_age_1990s_simple'] = ma_var(bins_a0_n2o_norm_age_nh, ['norm'], 'years', '1990s in situ mean age vs N2O, simple normalization')
  v__['n2o_age_1990s_age_conv'] = ma_var(bins_a0_n2o_norm_e_age_nh, ['norm'], 'years', '1990s in situ mean age vs N2O, age-convolution normalization')
  v__['n2o_age_1990s_uncert'] = ma_var(bins_a0_n2o_norm_age_sd_nh, ['norm'], 'years', 'uncertainty of the 1990s weighted means (used for both normalizations)')
  v__['n2o_age_2020s_simple'] = ma_var(bins_a2_n2o_norm_age_nh, ['norm'], 'years', '2020s in situ mean age vs N2O, simple normalization')
  v__['n2o_age_2020s_age_conv'] = ma_var(bins_a2_n2o_norm_e_age_nh, ['norm'], 'years', '2020s in situ mean age vs N2O, age-convolution normalization')
  v__['n2o_age_2020s_uncert'] = ma_var(bins_a2_n2o_norm_age_sd_nh, ['norm'], 'years', 'uncertainty of the 2020s weighted means')
  v__['ch4_age_1990s_simple'] = ma_var(bins_a0_ch4_norm_age_nh[*,0:1], ['norm','mean_unc'], 'years', '1990s mean age vs CH4, simple normalization (mean, uncertainty)')
  v__['ch4_age_1990s_age_conv'] = ma_var(bins_a0_ch4_norm_e_age_nh[*,0:1], ['norm','mean_unc'], 'years', '1990s mean age vs CH4, age-convolution normalization')
  v__['ch4_age_2020s_simple'] = ma_var(bins_a2_ch4_norm_age_nh[*,0:1], ['norm','mean_unc'], 'years', '2020s mean age vs CH4, simple normalization')
  v__['ch4_age_2020s_age_conv'] = ma_var(bins_a2_ch4_norm_e_age_nh[*,0:1], ['norm','mean_unc'], 'years', '2020s mean age vs CH4, age-convolution normalization')
  ma_nc_write, fig_dir__+'supp_fig1.nc', v__, global=hash('title', 'Supplementary Figure 1: mean age vs N2O and CH4 for two normalizations', 'normalization', norm_note__)
endif else print, 'Supplementary Figure 1 skipped: age-convolution bins not found in Aircraft_mean_ages_early_bins_comb.sav'

; ---------- Supplementary Figure 2: mean age trends on N2O (a) and CH4 (b) for two normalizations
v__ = orderedhash()
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized tracer bin centres')
v__['trend_on_n2o_simple'] = ma_var(n2o_nofl_age_trends[*,1], ['norm'], 'years/year')
v__['trend_on_n2o_simple_sigma'] = ma_var(n2o_nofl_age_trends_sigma[*,1], ['norm'], 'years/year')
v__['trend_on_n2o_age_conv'] = ma_var(n2o_e_nofl_age_trends[*,1], ['norm'], 'years/year')
v__['trend_on_n2o_age_conv_sigma'] = ma_var(n2o_e_nofl_age_trends_sigma[*,1], ['norm'], 'years/year')
v__['trend_on_ch4_simple'] = ma_var(age_on_ch4_trends[*,1], ['norm'], 'years/year')
v__['trend_on_ch4_simple_sigma'] = ma_var(age_on_ch4_trends_sigma[*,1], ['norm'], 'years/year')
v__['trend_on_ch4_age_conv'] = ma_var(age_on_ch4_e_trends[*,1], ['norm'], 'years/year')
v__['trend_on_ch4_age_conv_sigma'] = ma_var(age_on_ch4_e_trends_sigma[*,1], ['norm'], 'years/year')
ma_nc_write, fig_dir__+'supp_fig2.nc', v__, global=hash('title', 'Supplementary Figure 2: mean age trends for two normalizations', 'normalization', norm_note__)

; ---------- Supplementary Figure 6: trend profiles over the full period and the AirCore period
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt'], 'km')
v__['trend_alt_1993_2025'] = ma_var(age_nh_lat_adj_trends_coarse[*,1], ['alt'], 'years/year', '30-60N mean age trend vs altitude, 1993-2025')
v__['trend_alt_1993_2025_sigma'] = ma_var(age_nh_lat_adj_trends_coarse_sigma[*,1], ['alt'], 'years/year')
v__['trend_alt_2011_2025'] = ma_var(age_nh_lat_adj_trends_ac_coarse[*,1], ['alt'], 'years/year', '30-60N mean age trend vs altitude, 2011-2025 (AirCore period)')
v__['trend_alt_2011_2025_sigma'] = ma_var(age_nh_lat_adj_trends_ac_coarse_sigma[*,1], ['alt'], 'years/year')
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['trend_n2o_1993_2025'] = ma_var(n2o_nofl_age_trends[*,1], ['norm'], 'years/year', 'NH mean age trend vs normalized N2O, 1993-2025')
v__['trend_n2o_1993_2025_sigma'] = ma_var(n2o_nofl_age_trends_sigma[*,1], ['norm'], 'years/year')
v__['trend_n2o_2011_2025'] = ma_var(n2o_nofl_age_trends_ac[*,1], ['norm'], 'years/year', 'NH mean age trend vs normalized N2O, 2011-2025')
v__['trend_n2o_2011_2025_sigma'] = ma_var(n2o_nofl_age_trends_ac_sigma[*,1], ['norm'], 'years/year')
ma_nc_write, fig_dir__+'supp_fig6.nc', v__, global=hash('title', 'Supplementary Figure 6: NH midlatitude mean age trend profiles for two periods')

; ---------- Supplementary Figure 7: CCMI-2022 refD2 trend profiles for two periods
v__ = orderedhash()
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km')
v__['refd2_1990s_2020s_avg'] = ma_var(aoa_trends_d2_avg[*,0], ['alt_ccmi'], 'years/year', 'refD2 multi-model mean 30-60N mean age trend, 1990s-2020s')
v__['refd2_1990s_2020s_range'] = ma_var(aoa_trends_d2_minmax, ['alt_ccmi','minmax'], 'years/year')
v__['refd2_2010_2025_avg'] = ma_var(aoa_trends_d2_s_avg[*,0], ['alt_ccmi'], 'years/year', 'refD2 multi-model mean 30-60N mean age trend, 2010-2025')
v__['refd2_2010_2025_range'] = ma_var(aoa_trends_d2_s_minmax, ['alt_ccmi','minmax'], 'years/year')
ma_nc_write, fig_dir__+'supp_fig7.nc', v__, global=hash('title', 'Supplementary Figure 7: CCMI-2022 refD2 mean age trend profiles for two periods')

print, 'Supplementary figure data written to ', fig_dir__
