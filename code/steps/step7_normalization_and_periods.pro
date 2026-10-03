;+
;  Step 7: normalization and period comparisons (Supplementary Figs 1, 2, 6, 7).
;
;  Mean age trends for the two tracer normalizations and for the AirCore period,
;  and CCMI-2022 refD2 trends for two periods.
;
;  Reads:  step 3 outputs; Aircraft_mean_ages_early_bins_comb.sav; Andrews et al. N2O-age; CCMI-2022 trends
;  Writes: <output>/figure_data/supp_fig1.nc, supp_fig2.nc, supp_fig6.nc, supp_fig7.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step7_normalization_and_periods

dir = mean_age_data_dir()

nz2 = 30
dz2 = 1.0
alt_grid2 = findgen(nz2)*dz2+6.0
nn = 20
dn = 0.05
norm_grid = findgen(nn)*dn+0.025

nt = 2644
dt = 1./48.
years = findgen(nt)*dt+1975.
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_time_series.nc', $
  ['age_co2_b_grid_adj_tseries','age_co2_grid_adj_tseries','age_grid_adj_tseries','age_on_ch4_adj_tseries', $
  'age_on_ch4_e_adj_tseries','age_on_n2o_adj_tseries','age_on_n2o_e_adj_tseries','age_opt_b_grid_adj_tseries', $
  'age_opt_grid_adj_tseries','age_sf6_b_grid_adj_tseries','age_sf6_grid_adj_tseries','ch4_grid_adj_tseries', $
  'n2o_a_grid_adj_tseries','n2o_b_grid_adj_tseries','n2o_grid_adj_tseries','n2o_on_age_adj_tseries', $
  'n2o_on_age_both_adj_tseries','n2o_on_age_b_adj_tseries']

vv = replicate(!values.f_nan,nn,2)
n2o_nofl_age_trends = vv
n2o_nofl_age_trends_sigma = vv
age_on_ch4_trends = vv
age_on_ch4_trends_sigma = vv
n2o_e_nofl_age_trends = vv
n2o_e_nofl_age_trends_sigma = vv
age_on_ch4_e_trends = vv
age_on_ch4_e_trends_sigma = vv
n2o_nofl_age_trends_ac = vv
n2o_nofl_age_trends_ac_sigma = vv
for i = 0, nn-1 do begin
  gd = where(finite(age_on_n2o_adj_tseries[2,i,*]),ngd)
  if ngd ge 20 then begin
    n2o_nofl_age_trends[i,*] = linfit(years[gd],age_on_n2o_adj_tseries[2,i,gd],sigma=sigma)
    n2o_nofl_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_n2o_adj_tseries[2,i,*]) and (age_on_n2o_adj_tseries[5,i, $
    *] eq 3 or (age_on_n2o_adj_tseries[5,i,*] eq 1 and years ge 2015)),ngd)
  if ngd ge 10 then begin
    n2o_nofl_age_trends_ac[i,*] = linfit(years[gd],age_on_n2o_adj_tseries[2,i,gd],sigma=sigma)
    n2o_nofl_age_trends_ac_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_n2o_e_adj_tseries[2,i,*]),ngd)
  if ngd ge 20 then begin
    n2o_e_nofl_age_trends[i,*] = linfit(years[gd],age_on_n2o_e_adj_tseries[2,i,gd],sigma=sigma)
    n2o_e_nofl_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_ch4_adj_tseries[2,i,*]),ngd)
  if ngd ge 20 then begin
    age_on_ch4_trends[i,*] = linfit(years[gd],age_on_ch4_adj_tseries[2,i,gd],sigma=sigma)
    age_on_ch4_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_ch4_e_adj_tseries[2,i,*]),ngd)
  if ngd ge 20 then begin
    age_on_ch4_e_trends[i,*] = linfit(years[gd],age_on_ch4_e_adj_tseries[2,i,gd],sigma=sigma)
    age_on_ch4_e_trends_sigma[i,*] = sigma
  endif
endfor
restore,dir+'Aircraft/Airborne_Save_Files/Aircraft_mean_ages_early_bins_comb.sav'

ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'trend_profiles.nc', $
  ['age_nhe_adj_tseries','age_nhe_coarse','age_nh_adj_tseries','age_nh_adj_tseries_lat_adj','age_nh_coarse', $
  'age_nh_coarse_lat_adj','age_nh_lat_adj_trends_ac_coarse','age_nh_lat_adj_trends_ac_coarse_sigma', $
  'age_nh_lat_adj_trends_coarse','age_nh_lat_adj_trends_coarse_sigma','age_nh_lat_adj_trends_per_coarse', $
  'age_nh_lat_adj_trends_per_coarse_sigma','age_nh_trends_coarse','age_nh_trends_coarse_sigma','age_nh_trends_per_coarse', $
  'age_nh_trends_per_coarse_sigma','ch4_nh_trends_coarse','ch4_nh_trends_coarse_sigma','n2o_nh_lat_adj_trends_coarse', $
  'n2o_nh_lat_adj_trends_coarse_sigma','n2o_nh_lat_adj_trends_per_coarse','n2o_nh_lat_adj_trends_per_coarse_sigma']
restore,dir+'Models/CCMI/CCMI-2022_refd1_trends.sav'
restore,dir+'Models/CCMI/CCMI-2022_refd2_trends.sav'
nz_c = n_elements(pres_c)
alt_c = -7.*alog(pres_c/1e3)
zz = fltarr(nz_c,2)
zz2 = fltarr(nz_c,2)
aoa_trends_d2_minmax = zz
aoa_trends_d2_avg = zz2
aoa_trends_d2_s_minmax = zz
aoa_trends_d2_s_avg = zz2
for z = 0, nz_c-1 do begin
  aoa_trends_d2_minmax[z,0] = min(aoa_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_trends_d2_minmax[z,1] = max(aoa_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  for i = 0, 1 do aoa_trends_d2_avg[z,i] = mean(aoa_trends_avg_d2[i+1,2,z,10:15])
  aoa_trends_d2_s_minmax[z,0] = min(aoa_trends_d2_s[1,2,z,mi_d2_common,29:32],/nan)
  aoa_trends_d2_s_minmax[z,1] = max(aoa_trends_d2_s[1,2,z,mi_d2_common,29:32],/nan)
  for i = 0, 1 do aoa_trends_d2_s_avg[z,i] = mean(aoa_trends_avg_d2_s[i+1,2,z,29:32])
endfor

restore,dir+'Aircraft/Airborne_Save_Files/Andrews_N2O_age.sav'

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
norm_note__ = '"simple" = normalized by the surface mean at the measurement time (used in the paper); ' +  '"age_conv" = normalized by the surface time series convolved with the full age spectrum'

; Needs the age-convolution bins from Aircraft_mean_ages_early_bins_comb.sav (post-revision version).
if n_elements(bins_a0_n2o_norm_e_age_nh) gt 0 && n_elements(bins_a0_ch4_norm_e_age_nh) gt 0 then begin
  v__ = orderedhash()
  v__['norm_grid'] = ma_var(n2o_norm_grid, ['norm'], '', 'normalized tracer bin centres')
  v__['n2o_age_1990s_simple'] = ma_var(bins_a0_n2o_norm_age_nh, ['norm'], 'years', $
    '1990s in situ mean age vs N2O, simple normalization')
  v__['n2o_age_1990s_age_conv'] = ma_var(bins_a0_n2o_norm_e_age_nh, ['norm'], 'years', $
    '1990s in situ mean age vs N2O, age-convolution normalization')
  v__['n2o_age_1990s_uncert'] = ma_var(bins_a0_n2o_norm_age_sd_nh, ['norm'], 'years', $
    'uncertainty of the 1990s weighted means (used for both normalizations)')
  v__['n2o_age_2020s_simple'] = ma_var(bins_a2_n2o_norm_age_nh, ['norm'], 'years', $
    '2020s in situ mean age vs N2O, simple normalization')
  v__['n2o_age_2020s_age_conv'] = ma_var(bins_a2_n2o_norm_e_age_nh, ['norm'], 'years', $
    '2020s in situ mean age vs N2O, age-convolution normalization')
  v__['n2o_age_2020s_uncert'] = ma_var(bins_a2_n2o_norm_age_sd_nh, ['norm'], 'years', 'uncertainty of the 2020s weighted means')
  v__['ch4_age_1990s_simple'] = ma_var(bins_a0_ch4_norm_age_nh[*,0:1], ['norm','mean_unc'], 'years', $
    '1990s mean age vs CH4, simple normalization (mean, uncertainty)')
  v__['ch4_age_1990s_age_conv'] = ma_var(bins_a0_ch4_norm_e_age_nh[*,0:1], ['norm','mean_unc'], 'years', $
    '1990s mean age vs CH4, age-convolution normalization')
  v__['ch4_age_2020s_simple'] = ma_var(bins_a2_ch4_norm_age_nh[*,0:1], ['norm','mean_unc'], 'years', $
    '2020s mean age vs CH4, simple normalization')
  v__['ch4_age_2020s_age_conv'] = ma_var(bins_a2_ch4_norm_e_age_nh[*,0:1], ['norm','mean_unc'], 'years', $
    '2020s mean age vs CH4, age-convolution normalization')
  ma_nc_write, fig_dir__+'supp_fig1.nc', v__, $
    global=hash('title', 'Supplementary Figure 1: mean age vs N2O and CH4 for two normalizations', 'normalization', norm_note__)
endif

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
ma_nc_write, fig_dir__+'supp_fig2.nc', v__, $
  global=hash('title', 'Supplementary Figure 2: mean age trends for two normalizations', 'normalization', norm_note__)

; ---------- Supplementary Figure 6: trend profiles over the full period and the AirCore period
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt'], 'km')
v__['trend_alt_1993_2025'] = ma_var(age_nh_lat_adj_trends_coarse[*,1], ['alt'], 'years/year', $
  '30-60N mean age trend vs altitude, 1993-2025')
v__['trend_alt_1993_2025_sigma'] = ma_var(age_nh_lat_adj_trends_coarse_sigma[*,1], ['alt'], 'years/year')
v__['trend_alt_2011_2025'] = ma_var(age_nh_lat_adj_trends_ac_coarse[*,1], ['alt'], 'years/year', $
  '30-60N mean age trend vs altitude, 2011-2025 (AirCore period)')
v__['trend_alt_2011_2025_sigma'] = ma_var(age_nh_lat_adj_trends_ac_coarse_sigma[*,1], ['alt'], 'years/year')
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['trend_n2o_1993_2025'] = ma_var(n2o_nofl_age_trends[*,1], ['norm'], 'years/year', $
  'NH mean age trend vs normalized N2O, 1993-2025')
v__['trend_n2o_1993_2025_sigma'] = ma_var(n2o_nofl_age_trends_sigma[*,1], ['norm'], 'years/year')
v__['trend_n2o_2011_2025'] = ma_var(n2o_nofl_age_trends_ac[*,1], ['norm'], 'years/year', $
  'NH mean age trend vs normalized N2O, 2011-2025')
v__['trend_n2o_2011_2025_sigma'] = ma_var(n2o_nofl_age_trends_ac_sigma[*,1], ['norm'], 'years/year')
ma_nc_write, fig_dir__+'supp_fig6.nc', v__, $
  global=hash('title', 'Supplementary Figure 6: NH midlatitude mean age trend profiles for two periods')

; ---------- Supplementary Figure 7: CCMI-2022 refD2 trend profiles for two periods
v__ = orderedhash()
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km')
v__['refd2_1990s_2020s_avg'] = ma_var(aoa_trends_d2_avg[*,0], ['alt_ccmi'], 'years/year', $
  'refD2 multi-model mean 30-60N mean age trend, 1990s-2020s')
v__['refd2_1990s_2020s_range'] = ma_var(aoa_trends_d2_minmax, ['alt_ccmi','minmax'], 'years/year')
v__['refd2_2010_2025_avg'] = ma_var(aoa_trends_d2_s_avg[*,0], ['alt_ccmi'], 'years/year', $
  'refD2 multi-model mean 30-60N mean age trend, 2010-2025')
v__['refd2_2010_2025_range'] = ma_var(aoa_trends_d2_s_minmax, ['alt_ccmi','minmax'], 'years/year')
ma_nc_write, fig_dir__+'supp_fig7.nc', v__, $
  global=hash('title', 'Supplementary Figure 7: CCMI-2022 refD2 mean age trend profiles for two periods')
end

