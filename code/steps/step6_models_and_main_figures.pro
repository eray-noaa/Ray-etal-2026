;+
;  Step 6: model comparisons and main-text figure data.
;
;  ACE-FTS N2O and CH4 trends, CCMI-2022 refD1/refD2 trends on the observation
;  grids, and the Tropical Leaky Pipe (TLP) fits to the observed changes; writes the
;  data for Figs 1-4 and Extended Data Figs 1 and 6.
;
;  Reads:  step 2-5 outputs; ACE-FTS gridded data; CCMI-2022 means and trends; TLP model runs
;  Writes: <output>/figure_data/fig1.nc ... fig4.nc, ed_fig1a.nc, ed_fig1b.nc, ed_fig6.nc, supp_fig4.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step6_models_and_main_figures

dir = mean_age_data_dir()

;  Results from step4_balloon_profiles.
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'balloon_profiles.nc', $
  ['age_q_combo_ch4_bin_fl_avg','age_q_combo_f12_bin_fl_avg','age_q_combo_n2o_bin_fl_avg','alt_grid2','dn','n2o_from_ch4_yrs', $
  'norm_grid','nz2','years']
;  Results from step5_observed_trends.
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'observed_trends.nc', $
  ['age_alt_profile_avg','age_alt_profile_from_trend','age_grid_hist','age_on_ch4_trends','age_on_ch4_trends_sigma', $
  'alt_grid_a0','lat_grid_a','lat_tp','mean_age_on_n2o_20s_from_trend','mean_age_on_n2o_90s','n2o_alt_trends', $
  'n2o_nofl_age_trends','n2o_nofl_age_trends_sigma','n2o_norm_alt_profile_avg','n2o_norm_alt_profile_from_trend','nz', $
  'tp_alt_grid','tpause_alt','y','z']
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_time_series.nc', $
  ['age_co2_b_grid_adj_tseries','age_co2_grid_adj_tseries','age_grid_adj_tseries','age_on_ch4_adj_tseries', $
  'age_on_ch4_e_adj_tseries','age_on_n2o_adj_tseries','age_on_n2o_e_adj_tseries','age_opt_b_grid_adj_tseries', $
  'age_opt_grid_adj_tseries','age_sf6_b_grid_adj_tseries','age_sf6_grid_adj_tseries','ch4_grid_adj_tseries', $
  'n2o_a_grid_adj_tseries','n2o_b_grid_adj_tseries','n2o_grid_adj_tseries','n2o_on_age_adj_tseries', $
  'n2o_on_age_both_adj_tseries','n2o_on_age_b_adj_tseries']
restore,dir+'Balloon/Mean_age_relationships.sav'
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'mean_age_vs_tracers_1990s.nc', $
  ['age_grid_h','mean_age_all_on_ch4_a0','mean_age_all_on_ch4_b0','mean_age_all_on_ch4_hist','mean_age_all_on_ch4_hist_cent', $
  'mean_age_all_on_f12_a0','mean_age_all_on_f12_b0','mean_age_all_on_f12_hist','mean_age_all_on_f12_hist_cent', $
  'mean_age_all_on_n2o_0','mean_age_all_on_n2o_a0','mean_age_all_on_n2o_a2','mean_age_all_on_n2o_b0','mean_age_all_on_n2o_e_a0', $
  'mean_age_all_on_n2o_hist','mean_age_all_on_n2o_hist_cent','norm_grid2']
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'seasonal_grids_1990s.nc', $
  ['age_co2_old_grid_seas_sm','age_combo_grid_early','age_combo_grid_seas_sm_early','n2o_norm_grid_early', $
  'n2o_norm_grid_seas_both']
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'trend_profiles.nc', $
  ['age_nhe_adj_tseries','age_nhe_coarse','age_nh_adj_tseries','age_nh_adj_tseries_lat_adj','age_nh_coarse', $
  'age_nh_coarse_lat_adj','age_nh_lat_adj_trends_ac_coarse','age_nh_lat_adj_trends_ac_coarse_sigma', $
  'age_nh_lat_adj_trends_coarse','age_nh_lat_adj_trends_coarse_sigma','age_nh_lat_adj_trends_per_coarse', $
  'age_nh_lat_adj_trends_per_coarse_sigma','age_nh_trends_coarse','age_nh_trends_coarse_sigma','age_nh_trends_per_coarse', $
  'age_nh_trends_per_coarse_sigma','ch4_nh_trends_coarse','ch4_nh_trends_coarse_sigma','n2o_nh_lat_adj_trends_coarse', $
  'n2o_nh_lat_adj_trends_coarse_sigma','n2o_nh_lat_adj_trends_per_coarse','n2o_nh_lat_adj_trends_per_coarse_sigma']
;  Read in the gridded ACE data.
restore,dir+'/Satellite/ACE/ACE_gridding_5p3.sav'
yim = where(lat_grid gt 30 and lat_grid lt 60)
nza = n_elements(ace_alts)
nta = n_elements(time_grid)
ace_n2o_nh = replicate(!values.f_nan,nza,nta)
ace_n2o_nh_trends = replicate(!values.f_nan,2,nza)
ace_n2o_nh_trends_chi = replicate(!values.f_nan,nza)
ace_n2o_nh_avg = replicate(!values.f_nan,nza)
ace_ch4_nh = ace_n2o_nh
ace_ch4_nh_trends = ace_n2o_nh_trends
ace_ch4_nh_trends_chi = ace_n2o_nh_trends_chi
for t = 0, nta-1 do for zz = 6, nza-1 do begin
  ace_n2o_nh[zz,t] = mean(ace_tracers_lat_time_grid[yim,zz,t,2],/nan)
  ace_ch4_nh[zz,t] = mean(ace_tracers_lat_time_grid[yim,zz,t,8],/nan)
endfor
for zz = 6, nza-1 do begin
  chk = where(finite(ace_n2o_nh[zz,*]))
  ace_n2o_nh_avg[zz] = mean(ace_n2o_nh[zz,chk])
  ace_n2o_nh_trends[*,zz] = linfit(time_grid[chk],ace_n2o_nh[zz,chk],sigma=sigma)
  ace_n2o_nh_trends_chi[zz] = sigma[1]
  chk = where(finite(ace_ch4_nh[zz,*]))
  ace_ch4_nh_trends[*,zz] = linfit(time_grid[chk],ace_ch4_nh[zz,chk],sigma=sigma)
  ace_ch4_nh_trends_chi[zz] = sigma[1]
endfor

n2o_nh_alt_trends_combo = ace_n2o_nh_trends
n2o_nh_alt_trends_combo[0,*] = ace_n2o_nh_trends[1,*]
n2o_nh_alt_trends_combo[1,*] = ace_n2o_nh_trends_chi
zia = interpol(indgen(nz2),alt_grid2,ace_alts)
n2o_nh_alt_trends_combo[0,10:26] = interpolate(n2o_alt_trends[0:16,1],zia[10:26])
n2o_nh_alt_trends_combo[0,26] = 0.5 * (n2o_nh_alt_trends_combo[0,26] + ace_n2o_nh_trends[1,26])
n2o_nh_alt_trends_combo[0,27:28] += [5e-4,2e-4]
restore,dir+'Models/CCMI/CCMI-2022_refd1_means.sav'
restore,dir+'Models/CCMI/CCMI-2022_refd1_trends.sav'
restore,dir+'Models/CCMI/CCMI-2022_refd2_trends.sav'
ny_c = n_elements(lat_c)
nz_c = n_elements(pres_c)
nnm = n_elements(n2o_norm_grid_model)
alt_c = -7.*alog(pres_c/1e3)
zz = fltarr(nz_c,2)
zzz = fltarr(nz_c)
aoa_trends_d1_minmax = zz
aoa_trends_d1_avg = zzz
aoa_trends_d2_minmax = zz
aoa_trends_d2_avg = zzz
n2o_norm_trends_d1_minmax = zz
n2o_norm_trends_d1_avg = zzz
n2o_norm_trends_d2_minmax = zz
n2o_norm_trends_d2_avg = zzz
ch4_norm_trends_d1_minmax = zz
ch4_norm_trends_d1_avg = zzz
ch4_norm_trends_d2_minmax = zz
ch4_norm_trends_d2_avg = zzz
w_turn_trends_per_d2_minmax = zz
w_turn_trends_per_d2_avg = zzz
w_turn_trends_per_d1_minmax = zz
w_turn_trends_per_d1_avg = zzz
for z = 0, nz_c-1 do begin
  aoa_trends_d1_minmax[z,0] = min(aoa_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  aoa_trends_d1_minmax[z,1] = max(aoa_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  aoa_trends_d1_avg[z] = mean(aoa_trends_avg_d1[1,2,z,10:14])
  w_turn_trends_per_d1_minmax[z,0] = min(w_turn_trends_per_d1[1,1,z,mi_d1_common,10:14],/nan)
  w_turn_trends_per_d1_minmax[z,1] = max(w_turn_trends_per_d1[1,1,z,mi_d1_common,10:14],/nan)
  w_turn_trends_per_d1_avg[z] = mean(w_turn_trends_per_avg_d1[1,1,z,10:14])
  aoa_trends_d2_minmax[z,0] = min(aoa_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_trends_d2_minmax[z,1] = max(aoa_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_trends_d2_avg[z] = mean(aoa_trends_avg_d2[1,2,z,10:15])
  w_turn_trends_per_d2_minmax[z,0] = min(w_turn_trends_per_d2[1,1,z,mi_d2_common,10:15],/nan)
  w_turn_trends_per_d2_minmax[z,1] = max(w_turn_trends_per_d2[1,1,z,mi_d2_common,10:15],/nan)
  w_turn_trends_per_d2_avg[z] = mean(w_turn_trends_per_avg_d2[1,1,z,10:15])
  n2o_norm_trends_d1_minmax[z,0] = min(n2o_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  n2o_norm_trends_d1_minmax[z,1] = max(n2o_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  n2o_norm_trends_d1_avg[z] = mean(n2o_norm_trends_avg_d1[1,2,z,10:14])
  n2o_norm_trends_d2_minmax[z,0] = min(n2o_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  n2o_norm_trends_d2_minmax[z,1] = max(n2o_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  n2o_norm_trends_d2_avg[z] = mean(n2o_norm_trends_avg_d2[1,2,z,10:15])
  ch4_norm_trends_d1_minmax[z,0] = min(ch4_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  ch4_norm_trends_d1_minmax[z,1] = max(ch4_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  ch4_norm_trends_d1_avg[z] = mean(ch4_norm_trends_avg_d1[1,2,z,10:14])
  ch4_norm_trends_d2_minmax[z,0] = min(ch4_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  ch4_norm_trends_d2_minmax[z,1] = max(ch4_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  ch4_norm_trends_d2_avg[z] = mean(ch4_norm_trends_avg_d2[1,2,z,10:15])
endfor

yy = fltarr(nnm,2)
aoa_on_n2o_trends_d1_avg = fltarr(nnm)
aoa_on_n2o_trends_d2_avg = fltarr(nnm)
aoa_on_n2o_trends_d2_minmax = yy
aoa_on_ch4_trends_d1_minmax = yy
aoa_on_ch4_trends_d1_avg = fltarr(nnm)
aoa_on_ch4_trends_d2_minmax = yy
aoa_on_ch4_trends_d2_avg = fltarr(nnm)
ww = fltarr(nnm,2,3)
aoa_on_n2o_trends_d1_minmax = ww
for z = 0, nnm-1 do begin
  for j = 0, 2 do begin
    aoa_on_n2o_trends_d1_minmax[z,0,j] = min(aoa_on_n2o_trends_d1[j,2,z,mi_d1_common,10:14],/nan)
    aoa_on_n2o_trends_d1_minmax[z,1,j] = max(aoa_on_n2o_trends_d1[j,2,z,mi_d1_common,10:14],/nan)
  endfor
  aoa_on_n2o_trends_d1_avg[z] = mean(aoa_on_n2o_trends_avg_d1[1,2,z,10:14])
  aoa_on_n2o_trends_d2_minmax[z,0] = min(aoa_on_n2o_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_on_n2o_trends_d2_minmax[z,1] = max(aoa_on_n2o_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_on_n2o_trends_d2_avg[z] = mean(aoa_on_n2o_trends_avg_d2[1,2,z,10:15])
  aoa_on_ch4_trends_d1_minmax[z,0] = min(aoa_on_ch4_trends_d1[0,2,z,mi_d1_common,10:14],/nan)
  aoa_on_ch4_trends_d1_minmax[z,1] = max(aoa_on_ch4_trends_d1[0,2,z,mi_d1_common,10:14],/nan)
  aoa_on_ch4_trends_d1_avg[z] = mean(aoa_on_ch4_trends_avg_d1[0,2,z,10:14])
  aoa_on_ch4_trends_d2_minmax[z,0] = min(aoa_on_ch4_trends_d2[0,2,z,mi_d2_common,10:15],/nan)
  aoa_on_ch4_trends_d2_minmax[z,1] = max(aoa_on_ch4_trends_d2[0,2,z,mi_d2_common,10:15],/nan)
  aoa_on_ch4_trends_d2_avg[z] = mean(aoa_on_ch4_trends_avg_d2[0,2,z,10:15])
endfor

;  Make CCMI time averages.
ccc = replicate(!values.f_nan,ny_c,nz_c,nnm)
n2o_norm_c_1990s = ccc
age_c_1990s = ccc
ti1 = where(years_c ge 1990 and years_c lt 2000)
for m = 0, nnm-1 do for z = 0, nz_c-1 do for y = 0, ny_c-1 do begin
  n2o_norm_c_1990s[y,z,m] = mean(n2o_norm_d1[y,z,ti1,m])
  age_c_1990s[y,z,m] = mean(age_d1[y,z,ti1,m])
endfor

; Interpolate to obs vertical grid.
ddd = replicate(!values.f_nan,ny_c,nz,nnm)
n2o_norm_c_1990s_zobs = ddd
age_c_1990s_zobs = ddd
zic = interpol(indgen(nz_c),alt_c,alt_grid_a0)
for m = 0, nnm-1 do for y = 0, ny_c-1 do begin
  n2o_norm_c_1990s_zobs[y,*,m] = interpolate(n2o_norm_c_1990s[y,*,m],zic)
  age_c_1990s_zobs[y,*,m] = interpolate(age_c_1990s[y,*,m],zic)
endfor

;  Read in TLP model output

restore,dir+'Models/TLP/ideal_sweep_w_LowStrat_UpStrat_Budget.sav'
tlp_lu = tlp
run_key_lu = run_key

nz_tlp = n_elements(alt_tlp)
nnt = 21
dn2ot = 0.05
n2o_norm_grid_tlp = findgen(nnt)*dn2ot
ntlp = n_elements(run_key_lu)-1

age_nh_diff_alt_profile = 30.*age_nh_lat_adj_trends_coarse[*,1]
n2o_nh_diff_alt_profile = 30.*n2o_nh_lat_adj_trends_per_coarse[*,1]

n2o_nh_alt_trends_per_combo = n2o_nh_alt_trends_combo
for i = 0, 1 do n2o_nh_alt_trends_per_combo[i,*] = 1e2*reform(n2o_nh_alt_trends_combo[i,*])/ace_n2o_nh_avg

for z = 23, 24 do begin
  n2o_nh_diff_alt_profile[z] = 30 * 0.5 * total(n2o_nh_alt_trends_per_combo[0,z+5:z+6])
endfor

tlp_diffs = hash()

ni = interpol(indgen(nnt),n2o_norm_grid_tlp,n2o_norm_grid)

for ii = 1, ntlp do begin
  tlp_diffs[run_key_lu[ii]+'_mean_age'] = tlp_lu[run_key_lu[ii]+'_mean_age'] - tlp_lu[run_key_lu[0]+'_mean_age']
  tlp_diffs[run_key_lu[ii]+'_n2o'] = tlp_lu[run_key_lu[ii]+'_n2o'] - tlp_lu[run_key_lu[0]+'_n2o']
  tlp_diffs[run_key_lu[ii]+'_n2o_per'] = 1e2*tlp_diffs[run_key_lu[ii]+'_n2o'] / tlp_lu[run_key_lu[0]+'_n2o']

  n2o_mean_age_diffs = replicate(!values.f_nan,3,n_elements(n2o_norm_grid))
  for y = 0, 2 do n2o_mean_age_diffs[y,*] = interpolate(tlp_lu[run_key_lu[ii]+'_mean_age_n2o',y, $
    *]-tlp_lu[run_key_lu[0]+'_mean_age_n2o',y,*],ni)
  tlp_diffs[run_key_lu[ii]+'_n2o_mean_age'] = n2o_mean_age_diffs
endfor

aa = replicate(!values.f_nan,ntlp+1)
rmse_age_nh = aa
rmse_n2o_nh = aa
rmse_mid_age_nh = aa

;  Mean age profile differences from 16-17km.
nage = 2.
age_nh_i = indgen(nage)
;  Mean age profile differences from 22-27km.
nage_mid = 6.
age_mid_nh_i = indgen(nage_mid)+6
;  N2O profile differences from 20-27km.
nn2o = 10.
n2o_nh_i = indgen(nn2o)+5

for i = 1, ntlp do begin
  mean_age_tlp_nh_diffs = reform(tlp_diffs[run_key_lu[i]+'_mean_age',2,*])
  rmse_age_nh[i] = sqrt(total(mean_age_tlp_nh_diffs[age_nh_i] - age_nh_diff_alt_profile[age_nh_i+10])^2/nage)
  rmse_mid_age_nh[i] = sqrt(total(mean_age_tlp_nh_diffs[age_mid_nh_i] - age_nh_diff_alt_profile[age_mid_nh_i+10])^2/nage)

  n2o_tlp_nh_diffs = reform(tlp_diffs[run_key_lu[i]+'_n2o_per',2,*])
  rmse_n2o_nh[i] = sqrt(total(n2o_tlp_nh_diffs[n2o_nh_i] - n2o_nh_diff_alt_profile[n2o_nh_i+10])^2/nn2o)
endfor

rmse_age_thresh = 0.16
rmse_up_n2o_thresh = 15
best_age_nh_i = where(rmse_age_nh le rmse_age_thresh,nb_ra)

;  Best N2O within the age NH set of solutions.
best_age_n2o_nh_i = where(rmse_n2o_nh[best_age_nh_i] le rmse_up_n2o_thresh and rmse_mid_age_nh[best_age_nh_i] le 2.5*rmse_age_thresh,nbu_rn)

restore,dir+'Models/TLP/ideal_sweep_w_LowStrat.sav'
tlp_low = tlp
run_key_low = run_key
restore,dir+'Models/TLP/ideal_sweep_w_UpStrat_budget.sav'

wdiffs = fltarr(nbu_rn,nz_tlp)
for i = 0, nbu_rn-1 do wdiffs[i,*] = 1e2*(tlp_lu[run_key_lu[best_age_nh_i[best_age_n2o_nh_i[i]]]+'_w',1, $
  *]-tlp_lu[run_key_lu[0]+'_w',1,*])/tlp_lu[run_key_lu[0]+'_w',1,*]
diffs_minmax = fltarr(nz_tlp,2)
for z = 0, nz_tlp-1 do diffs_minmax[z,*] = [min(wdiffs[1:-1,z]),max(wdiffs[1:-1,z])]
diffs = 0.5 * (rmse_age_nh[best_age_nh_i[best_age_n2o_nh_i]] + rmse_n2o_nh[best_age_nh_i[best_age_n2o_nh_i]])
weights = 1. / diffs[1:-1]^2
wmean = fltarr(nz_tlp)
for z = 0, nz_tlp-1 do wmean[z] = total(weights * wdiffs[1:-1,z]) / total(weights)

inter_dir__ = ma_output_dir() + 'intermediate' + path_sep()
fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
stat5__ = ['weighted_mean','uncertainty_of_weighted_mean','mean','standard_deviation','n_measurements']
flask__ = ma_nc_read_full(inter_dir__+'flask_profiles.nc')
; Tropical Leaky Pipe model curves (base run fit to in situ data, and w+20%).
tlp__ = orderedhash()
foreach k__, ['n2o','ch4','f12'] do begin
  nend__ = (k__ eq 'ch4') ? [-18,-14] : [-20,-16]
  tlp__['tlp_base_mean_age_'+k__] = ma_var(reform(tlp_low[run_key_low[26]+'_mean_age',2,0:nend__[0]]), $
    ['tlp_base_'+k__+'_level'], 'years',  'TLP base run NH mean age (curve vs '+k__+')')
  tlp__['tlp_base_'+k__+'_norm'] = ma_var(reform(tlp_low[run_key_low[26]+'_'+k__,2,0:nend__[0]]), ['tlp_base_'+k__+'_level'], $
    '', 'TLP base run NH normalized '+k__)
  tlp__['tlp_w20_mean_age_'+k__] = ma_var(reform(tlp[run_key[47]+'_mean_age',2,0:nend__[1]]), ['tlp_w20_'+k__+'_level'], $
    'years',  'TLP run with 20% faster tropical upwelling: NH mean age (curve vs '+k__+')')
  tlp__['tlp_w20_'+k__+'_norm'] = ma_var(reform(tlp[run_key[47]+'_'+k__,2,0:nend__[1]]), ['tlp_w20_'+k__+'_level'], '', $
    'TLP w+20% run NH normalized '+k__)
endforeach
; ---------- Figure 1 and Extended Data Figure 1: mean age vs normalized tracers
foreach tr__, ['n2o','ch4','f12'] do begin
  v__ = orderedhash()
  v__['norm_grid_hist'] = ma_var(norm_grid2, ['norm_hist'], '', $
    'normalized tracer bin centres for the 1990s in situ distribution')
  v__['age_grid_hist'] = ma_var((tr__ eq 'n2o') ? age_grid_h : age_grid_hist, ['age_hist'], 'years', $
    'mean age bin centres for the 1990s in situ distribution')
  hist__ = (tr__ eq 'n2o') ? mean_age_all_on_n2o_hist : ((tr__ eq 'ch4') ? mean_age_all_on_ch4_hist : mean_age_all_on_f12_hist)
  v__['insitu_1990s_hist'] = ma_var(hist__, ['age_hist','norm_hist'], '', $
    'fraction of 1990s in situ measurements in each normalized-tracer bin falling in each mean age bin (shaded where >= 0.03)')
  v__['norm_grid2'] = ma_var(norm_grid2, ['norm2'], '', 'normalized tracer bin centres (interval 0.02)')
  a0__ = (tr__ eq 'n2o') ? mean_age_all_on_n2o_a0 : ((tr__ eq 'ch4') ? mean_age_all_on_ch4_a0 : mean_age_all_on_f12_a0)
  b0__ = (tr__ eq 'n2o') ? mean_age_all_on_n2o_b0 : ((tr__ eq 'ch4') ? mean_age_all_on_ch4_b0 : mean_age_all_on_f12_b0)
  v__['aircraft_1990s_mean_age'] = ma_var(a0__, ['norm2','stat5'], 'years', $
    '1990s in situ aircraft mean age in each normalized-tracer bin', $
    statistics='columns: '+strjoin(stat5__,', ')+'; the figure plots column 2 with column 3 as error bar')
  v__['balloon_1990s_mean_age'] = ma_var(b0__, ['norm2','stat5'], 'years', $
    '1990s in situ balloon mean age in each normalized-tracer bin',  statistics='columns: '+strjoin(stat5__,', '))
  v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized tracer bin centres (interval 0.05)')
  fl__ = (tr__ eq 'n2o') ? age_q_combo_n2o_bin_fl_avg : ((tr__ eq 'ch4') ? age_q_combo_ch4_bin_fl_avg : age_q_combo_f12_bin_fl_avg)
  v__['flask_avg_mean_age'] = ma_var(fl__, ['norm','stat5'], 'years', $
    'cryo-flask (1970s-2000s) mean age averaged in each normalized-tracer bin',  statistics='columns: '+strjoin(stat5__,', '))
  v__ += flask__
  foreach k__, tlp__.keys() do begin
    if strpos(k__, '_'+tr__) ge 0 then v__[k__] = tlp__[k__]
  endforeach
  fname__ = (tr__ eq 'n2o') ? 'fig1.nc' : ((tr__ eq 'ch4') ? 'ed_fig1a.nc' : 'ed_fig1b.nc')
  ma_nc_write, fig_dir__+fname__, v__, global=hash('title', $
    'Mean age vs normalized '+strupcase(tr__)+' (' + file_basename(fname__,'.nc') + ')')
endforeach

; ---------- Figure 2: NH midlatitude mean age and N2O trend profiles
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt_obs'], 'km', 'altitude bin centres for observed trends (1 km bins)')
v__['insitu_age_trend'] = ma_var(reform(age_nh_lat_adj_trends_coarse[*,1]), ['alt_obs'], 'years/year', $
  'in situ 30-60N mean age trend (latitude adjusted)')
v__['insitu_age_trend_sigma'] = ma_var(reform(age_nh_lat_adj_trends_coarse_sigma[*,1]), ['alt_obs'], 'years/year', $
  '1-sigma uncertainty of insitu_age_trend')
v__['insitu_n2o_norm_trend'] = ma_var(reform(n2o_nh_lat_adj_trends_coarse[*,1]), ['alt_obs'], '1/year', $
  'in situ 30-60N normalized N2O trend (latitude adjusted)')
v__['insitu_n2o_norm_trend_sigma'] = ma_var(reform(n2o_nh_lat_adj_trends_coarse_sigma[*,1]), ['alt_obs'], '1/year', $
  '1-sigma uncertainty of insitu_n2o_norm_trend')
v__['ace_alt'] = ma_var(ace_alts, ['alt_ace'], 'km', 'ACE-FTS altitude grid')
v__['ace_n2o_trend'] = ma_var(reform(ace_n2o_nh_trends[1,*]), ['alt_ace'], '1/year', $
  'ACE-FTS v4.1 30-60N normalized N2O trend 2004-2025')
v__['ace_n2o_trend_sigma'] = ma_var(ace_n2o_nh_trends_chi, ['alt_ace'], '1/year', '1-sigma uncertainty of ace_n2o_trend')
v__['ace_age_trend_alt'] = ma_var([15.5,18.5], ['ace_age'], 'km', $
  'altitudes of published ACE-FTS mean age trends (Saunders et al.; values entered by hand in age_time_series.pro)')
v__['ace_age_trend'] = ma_var([-0.014,-0.010], ['ace_age'], 'years/year', 'published ACE-FTS 40-50N mean age trends 2004-2021')
v__['ace_age_trend_sigma'] = ma_var([0.006,0.005], ['ace_age'], 'years/year', '1-sigma uncertainty of ace_age_trend')
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km', 'CCMI-2022 altitude (-7 ln(p/1000 hPa))')
foreach d__, ['d1','d2'] do begin
  v__['ccmi_ref'+d__+'_age_trend_avg'] = ma_var((d__ eq 'd1') ? aoa_trends_d1_avg : aoa_trends_d2_avg, ['alt_ccmi'], $
    'years/year', 'CCMI-2022 ref'+strupcase(d__)+' multi-model mean 30-60N mean age trend')
  v__['ccmi_ref'+d__+'_age_trend_range'] = ma_var((d__ eq 'd1') ? aoa_trends_d1_minmax : aoa_trends_d2_minmax, $
    ['alt_ccmi','minmax'], 'years/year', 'CCMI-2022 ref'+strupcase(d__)+' min and max over ensemble members and periods')
  v__['ccmi_ref'+d__+'_n2o_norm_trend_avg'] = ma_var((d__ eq 'd1') ? n2o_norm_trends_d1_avg : n2o_norm_trends_d2_avg, $
    ['alt_ccmi'], '1/year', 'CCMI-2022 ref'+strupcase(d__)+' multi-model mean 30-60N normalized N2O trend')
  v__['ccmi_ref'+d__+'_n2o_norm_trend_range'] = ma_var((d__ eq 'd1') ? n2o_norm_trends_d1_minmax : n2o_norm_trends_d2_minmax, $
    ['alt_ccmi','minmax'], '1/year', 'CCMI-2022 ref'+strupcase(d__)+' min and max')
endforeach
ma_nc_write, fig_dir__+'fig2.nc', v__, global=hash('title', 'Figure 2: NH mean age and N2O trend profiles', 'note', $
  'Trends are per year; the figure plots them per decade (x10).')

; ---------- Figure 3
z3a__ = 12
v__ = orderedhash()
v__['years'] = ma_var(years, ['time'], 'year', 'decimal year (bins of 1/48 year)')
v__['n2o_norm_bin'] = ma_var([norm_grid[z3a__]-dn/2., norm_grid[z3a__]+dn/2.], ['bin_edge'], '', $
  'normalized N2O range of panel a')
v__['age_in_bin'] = ma_var(reform(age_on_n2o_adj_tseries[2,z3a__,*]), ['time'], 'years', $
  'mean age averaged in the panel a normalized N2O bin')
v__['age_in_bin_sd'] = ma_var(reform(age_on_n2o_adj_tseries[3,z3a__,*]), ['time'], 'years', 'standard deviation of age_in_bin')
v__['sample_type'] = ma_var(fix(reform(age_on_n2o_adj_tseries[5,z3a__,*])), ['time'], '', 'platform', flag_values='1 2 3', $
  flag_meanings='aircraft in_situ_balloon aircore')
v__['n2o_from_ch4'] = ma_var(fix(n2o_from_ch4_yrs), ['time'], '', '1 where AirCore N2O was derived from a CH4-N2O relationship')
v__['bin_trend_fit'] = ma_var(reform(n2o_nofl_age_trends[z3a__,*]), ['coef'], '', $
  'linear fit (intercept years, slope years/year) to age_in_bin')
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres (interval 0.05)')
v__['insitu_age_trend_on_n2o'] = ma_var(reform(n2o_nofl_age_trends[*,1]), ['norm'], 'years/year', $
  'in situ NH extratropical mean age trend in normalized N2O bins')
v__['insitu_age_trend_on_n2o_sigma'] = ma_var(reform(n2o_nofl_age_trends_sigma[*,1]), ['norm'], 'years/year', $
  '1-sigma uncertainty')
v__['ccmi_n2o_norm_grid'] = ma_var(n2o_norm_grid_model, ['norm_ccmi'], '', 'normalized N2O grid of CCMI trends')
v__['ccmi_refd1_age_trend_on_n2o_avg'] = ma_var(aoa_on_n2o_trends_d1_avg, ['norm_ccmi'], 'years/year', $
  'CCMI-2022 refD1 multi-model mean')
v__['ccmi_refd1_age_trend_on_n2o_range'] = ma_var(reform(aoa_on_n2o_trends_d1_minmax[*,*,1]), ['norm_ccmi','minmax'], $
  'years/year', 'CCMI-2022 refD1 min and max')
v__['ccmi_refd2_age_trend_on_n2o_avg'] = ma_var(aoa_on_n2o_trends_d2_avg, ['norm_ccmi'], 'years/year', $
  'CCMI-2022 refD2 multi-model mean')
v__['ccmi_refd2_age_trend_on_n2o_range'] = ma_var(aoa_on_n2o_trends_d2_minmax, ['norm_ccmi','minmax'], 'years/year', $
  'CCMI-2022 refD2 min and max')
v__['norm_grid2'] = ma_var(norm_grid2, ['norm2'], '', 'normalized N2O bin centres (interval 0.02)')
v__['age_on_n2o_1990s'] = ma_var(mean_age_on_n2o_90s, ['norm2'], 'years', 'NH extratropical mean age vs normalized N2O, 1990s')
v__['age_on_n2o_2020s'] = ma_var(mean_age_on_n2o_20s_from_trend, ['norm2'], 'years', $
  'NH extratropical mean age vs normalized N2O, 2020s (1990s plus trend)')
v__['arrow_alt'] = ma_var(alt_grid2[[11,13,14,15,17,19]], ['arrow'], 'km', $
  'altitudes of the constant-altitude change arrows in panel c')
v__['arrow_age'] = ma_var(transpose([[reform(age_alt_profile_avg[0,[11,13,14,15,17,19]])], $
  [age_alt_profile_from_trend[[11,13,14,15,17,19]]]]), ['start_end','arrow'], 'years')
v__['arrow_n2o_norm'] = ma_var(transpose([[reform(n2o_norm_alt_profile_avg[0,[11,13,14,15,17,19]])], $
  [n2o_norm_alt_profile_from_trend[[11,13,14,15,17,19]]]]), ['start_end','arrow'], '')
ma_nc_write, fig_dir__+'fig3.nc', v__, global=hash('title', 'Figure 3: mean age trends as a function of N2O')

v__ = orderedhash()
v__['tlp_alt'] = ma_var(alt_tlp, ['alt_tlp'], 'km', 'TLP model altitude')
v__['tlp_w_trend_range'] = ma_var(diffs_minmax/3., ['alt_tlp','minmax'], '%/decade', $
  'range of TLP tropical upwelling trends within tolerance of the observed changes')
v__['tlp_w_trend_wmean'] = ma_var(wmean/3., ['alt_tlp'], '%/decade', $
  'weighted mean TLP tropical upwelling trend (figure applies a 3-point smooth)')
v__['ccmi_alt'] = ma_var(alt_c, ['alt_ccmi'], 'km')
v__['ccmi_refd1_w_trend_avg'] = ma_var(10.*w_turn_trends_per_d1_avg, ['alt_ccmi'], '%/decade', $
  'CCMI-2022 refD1 tropical mean upwelling trend')
v__['ccmi_refd1_w_trend_range'] = ma_var(10.*w_turn_trends_per_d1_minmax, ['alt_ccmi','minmax'], '%/decade')
v__['ccmi_refd2_w_trend_avg'] = ma_var(10.*w_turn_trends_per_d2_avg, ['alt_ccmi'], '%/decade', $
  'CCMI-2022 refD2 tropical mean upwelling trend')
v__['ccmi_refd2_w_trend_range'] = ma_var(10.*w_turn_trends_per_d2_minmax, ['alt_ccmi','minmax'], '%/decade')
ma_nc_write, fig_dir__+'fig4.nc', v__, global=hash('title', 'Figure 4: tropical upwelling trends')

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

; ---------- Supplementary Figure 4: 1990s gridded mean age and N2O, 20 km latitude distributions
z20__ = 28
v__ = orderedhash()
v__['equiv_lat'] = ma_var(lat_grid_a, ['lat'], 'degrees_north', 'N2O equivalent latitude bin centres')
v__['alt'] = ma_var(alt_grid, ['alt'], 'km', 'altitude bin centres')
v__['tropopause_alt_on_grid'] = ma_var(tp_alt_grid, ['lat'], 'km', $
  'climatological tropopause altitude on the equivalent latitude grid')
v__['age_1990s'] = ma_var(age_combo_grid_early, ['stat','lat','alt'], 'years', '1990s gridded mean age', $
  statistics='index 0 mean; index 2 number of seasons sampled (full-size symbols where >= 3, small where 1-2)')
v__['n2o_norm_1990s'] = ma_var(n2o_norm_grid_early, ['stat','lat','alt'], '', '1990s gridded normalized N2O', $
  statistics='index 0 mean; index 2 number of seasons sampled')
v__['merra2_lat'] = ma_var(lat_tp, ['lat_tp'], 'degrees_north')
v__['merra2_tropopause_alt'] = ma_var(tpause_alt, ['lat_tp'], 'km', 'MERRA-2 climatological tropopause altitude (magenta line)')
v__['alt_20km'] = ma_var(alt_grid[z20__], '', 'km', 'altitude of panels c and d')
v__['ccmi_lat'] = ma_var(lat_c, ['lat_ccmi'], 'degrees_north')
v__['ccmi_age_20km_1990s'] = ma_var(reform(age_c_1990s_zobs[*,z20__,*]), ['lat_ccmi','member'], 'years', $
  'CCMI-2022 refD1 1990s mean age at 20 km, each model member')
v__['ccmi_n2o_norm_20km_1990s'] = ma_var(reform(n2o_norm_c_1990s_zobs[*,z20__,*]), ['lat_ccmi','member'], '', $
  'CCMI-2022 refD1 1990s normalized N2O at 20 km, each model member (the figure omits members 0-2 and 15)')
ma_nc_write, fig_dir__+'supp_fig4.nc', v__, global=hash('title', 'Supplementary Figure 4: 1990s gridded mean age and N2O')
end

