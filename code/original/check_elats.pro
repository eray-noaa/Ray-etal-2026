pro check_elats, STAGE=stage

;  STAGE (string) selects what is rebuilt; run_all.pro runs the stages in order:
;    'elats'       - N2O equivalent latitudes and gridded fields
;                    (writes Elat_adj_all, OMS_grid_elat_adj, Aircore_grid_elat_adj), then returns
;    'time_series' - gridded time series, mean age vs N2O relationships and trends
;                    (writes Age_N2O_gridded_N2O_elat_adj_time_series, Mean_age_vs_N2O_1990s,
;                    Mean_age_and_N2O_trend_alt_profile_elat_adj and the Supp. Fig. 3 and 5 data), then returns
;  With no STAGE the cached files are used and the program runs to the end as before.
if n_elements(stage) eq 0 then stage = ''

dir = mean_age_data_dir()

;  Tropopause data.
ncdf_get,dir+'Reanalysis/MERRA2/tp.monmean.zm.nc',['time','tpp','tpt','tpz','lat'],tp,/quiet
ny_tp = tp['lat','dim_sizes']
lat_tp = tp['lat','value']
tp_jdays = julday(1,1,1900,0,0) + double(tp['time','value'])
caldat,tp_jdays,mon,day,yr
tp_yr = float(yr) + float(mon)/12. + 1./24.
tpause = fltarr(ny_tp) & tpause_alt = tpause & tpause_mon = fltarr(ny_tp,12) & tpause_mon_alt = tpause_mon & tpause_mon_temp = tpause_mon & tpause_seas = fltarr(ny_tp,4)
tpause_seas_alt = tpause_seas & tpause_seas_temp = tpause_seas & tpause_temp = tpause
for y = 0, ny_tp[0]-1 do begin
  tpause[y] = mean(tp['tpp','value',y,*])
  tpause_alt[y] = mean(tp['tpz','value',y,*])
  tpause_temp[y] = mean(tp['tpt','value',y,*])
  for t = 0, 11 do begin
    ti = where(mon eq t+1)
    tpause_mon[y,t] = mean(tp['tpp','value',y,ti],/nan)
    tpause_mon_alt[y,t] = mean(tp['tpz','value',y,ti],/nan)
    tpause_mon_temp[y,t] = mean(tp['tpt','value',y,ti],/nan)
  endfor
  for t = 0, 3 do begin
    ti = indgen(3)-1+t*3
    tpause_seas[y,t] = mean(tpause_mon[y,ti])
    tpause_seas_alt[y,t] = mean(tpause_mon_alt[y,ti])
    tpause_seas_temp[y,t] = mean(tpause_mon_temp[y,ti])
  endfor
endfor
tpause_theta = tpause_temp * (1e3/tpause)^(2./7.)
tpause_mon_theta = tpause_mon_temp * (1e3/tpause_mon)^(2./7.)
tpause_seas_theta = tpause_seas_temp * (1e3/tpause_seas)^(2./7.)

;  Read in the gridded ACE data.
restore,dir+'/Satellite/ACE/ACE_gridding_5p3.sav'
yim = where(lat_grid gt 30 and lat_grid lt 70)
nza = n_elements(ace_alts)
nta = n_elements(time_grid)
ace_n2o_nh = replicate(!values.f_nan,nza,nta) & ace_n2o_nh_trends = replicate(!values.f_nan,2,nza) & ace_n2o_nh_trends_chi = replicate(!values.f_nan,nza)
ace_n2o_nh_trends_per = ace_n2o_nh_trends_chi & ace_n2o_nh_trends_chi_per = ace_n2o_nh_trends_chi & ace_n2o_nh_avg = replicate(!values.f_nan,nza)
ace_ch4_nh = ace_n2o_nh & ace_ch4_nh_avg = ace_n2o_nh_avg & ace_ch4_nh_trends = ace_n2o_nh_trends & ace_ch4_nh_trends_chi = ace_n2o_nh_trends_chi
ace_ch4_nh_trends_per = ace_n2o_nh_trends_per & ace_ch4_nh_trends_chi_per = ace_n2o_nh_trends_chi_per
for t = 0, nta-1 do for zz = 6, nza-1 do begin
  ace_n2o_nh[zz,t] = mean(ace_tracers_lat_time_grid[yim,zz,t,2],/nan)
  ace_ch4_nh[zz,t] = mean(ace_tracers_lat_time_grid[yim,zz,t,8],/nan)
endfor
for zz = 6, nza-1 do begin
  chk = where(finite(ace_n2o_nh[zz,*]))
  ace_n2o_nh_avg[zz] = mean(ace_n2o_nh[zz,chk])
  ace_n2o_nh_trends[*,zz] = linfit(time_grid[chk],ace_n2o_nh[zz,chk],sigma=sigma)
  ace_n2o_nh_trends_chi[zz] = sigma[1]
  ace_n2o_nh_trends_per[zz] = 1e2*ace_n2o_nh_trends[1,zz] / ace_n2o_nh_avg[zz]
  ace_n2o_nh_trends_chi_per[zz] = 1e2*ace_n2o_nh_trends_chi[zz] / ace_n2o_nh_avg[zz]

  chk = where(finite(ace_ch4_nh[zz,*]))
  ace_ch4_nh_avg[zz] = mean(ace_ch4_nh[zz,chk])
  ace_ch4_nh_trends[*,zz] = linfit(time_grid[chk],ace_ch4_nh[zz,chk],sigma=sigma)
  ace_ch4_nh_trends_chi[zz] = sigma[1]
  ace_ch4_nh_trends_per[zz] = 1e2*ace_ch4_nh_trends[1,zz] / ace_ch4_nh_avg[zz]
  ace_ch4_nh_trends_chi_per[zz] = 1e2*ace_ch4_nh_trends_chi[zz] / ace_ch4_nh_avg[zz]
endfor
lat_grid_ace = lat_grid

data = read_ascii(dir+'Surface_Trace_Gas/GMD/n2o/combined/GML_global_N2O.txt',data_start=68)
n2o_yr = reform(data.field01[0,*])
n2o_mon = reform(data.field01[1,*])
n2o_date = n2o_yr + (n2o_mon-0.5)/12.
n2o_global = reform(data.field01[6,*])
data = read_ascii(dir+'Surface_Trace_Gas/GMD/cfc12/combined/HATS_global_F12.txt',data_start=68)
f12_yr = reform(data.field01[0,*])
f12_mon = reform(data.field01[1,*])
f12_date = f12_yr + (f12_mon-0.5)/12.
f12_global = reform(data.field01[6,*])


restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_early.sav'
caldat,jday_all,months_all,days_all,years_all
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
ti = interpol(indgen(n_elements(f12_date)),f12_date,times_all)
f12_surface = interpolate(f12_global,ti)
f12_norm = f12/f12_surface
co2_a0 = co2 & sf6_a0 = sf6 & elat_m_a0 = elat_m & lat_a0 = lat & alt_a0 = alt & n2o_a0 = n2o & n2o_norm_a0 = n2o_norm & months_all_a0 = months_all & times_all_a0 = times_all
ch4_all_norm_a0 = ch4_all_norm & f12_norm_a0 = f12_norm

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_mid.sav'
caldat,jday_all,months_all,days_all,years_all
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
ti = interpol(indgen(n_elements(f12_date)),f12_date,times_all)
f12_surface = interpolate(f12_global,ti)
f12_norm = f12/f12_surface
co2_a1 = co2 & sf6_a1 = sf6 & elat_m_a1 = elat_m & lat_a1 = lat & alt_a1 = alt & n2o_a1 = n2o & n2o_norm_a1 = n2o_norm & months_all_a1 = months_all & times_all_a1 = times_all
ch4_all_norm_a1 = ch4_all_norm & f12_norm_a1 = f12_norm

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_late.sav'
caldat,jday_all,months_all,days_all,years_all
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
ti = interpol(indgen(n_elements(f12_date)),f12_date,times_all)
f12_surface = interpolate(f12_global,ti)
f12_norm = f12/f12_surface
co2_a = co2 & sf6_a = sf6 & elat_m_a = elat_m & lat_a = lat & alt_a = alt & n2o_a = n2o & n2o_norm_a = n2o_norm & months_all_a = months_all & times_all_a = times_all
ch4_all_norm_a = ch4_all_norm & f12_norm_a = f12_norm

restore,dir+'Balloon/OMS_common_gc_merge.sav'
co2_b = co2 & sf6_b = sf6 & elat_m_b = elat_m & lat_b = lat & alt_b = alt & n2o_b = n2o & n2o_norm_b = n2o_norm & months_all_b = months_all & times_all_b = times_all
co2_all_b = co2_all & ch4_all_norm_b = ch4_all_norm & f12_norm_b = f12_norm
;  Fill in the N2O.
chk_n2o = where(finite(n2o),nchk_n2o)
n2o_norm_all = replicate(!values.f_nan,n_elements(n2o))
n2o_norm_all[chk_n2o] = n2o_norm[chk_n2o]
for i = 0, nchk_n2o-2 do begin
  ni = chk_n2o[i+1] - (chk_n2o[i]+1)
  if ni gt 0 and ni lt 100 then begin
    ii = findgen(ni) + chk_n2o[i]+1
    n2o_norm_all[ii] = n2o_norm_all[chk_n2o[i]] + (n2o_norm_all[chk_n2o[i+1]] - n2o_norm_all[chk_n2o[i]]) * (findgen(ni)+1.)/(ni+1.)
  endif
endfor
n2o_norm_all_b = n2o_norm_all


restore,dir+'Balloon/Aircore_common_merge.sav'

restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sweep_sf6_co2_early.sav'
dat_a0 = dat
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_early.sav'
dat_co2_a0 = dat_co2
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sf6_early.sav'
dat_sf6_a0 = dat_sf6
restore,dir+'Balloon/Balloon_mean_ages_sweep_sf6_co2_early.sav'
dat_b0 = dat
restore,dir+'Balloon/Balloon_mean_ages_co2_early.sav'
dat_co2_b0 = dat_co2
restore,dir+'Balloon/Balloon_mean_ages_sf6_early.sav'
dat_sf6_b0 = dat_sf6
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sweep_sf6_co2_mid.sav'
dat_a1 = dat
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_mid.sav'
dat_co2_a1 = dat_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_mid.sav'
dat_sf6_a1 = dat_sf6
restore,dir+'Balloon/Balloon_mean_ages_sweep_sf6_co2_late_v2025.sav'
dat_b2 = dat
restore,dir+'Balloon/Balloon_mean_ages_co2_late_v2025.sav'
dat_co2_b2 = dat_co2
restore,dir+'Balloon/Balloon_mean_ages_sf6_late.sav'
dat_sf6_b2 = dat_sf6
restore,dir+'Balloon/Balloon_mean_ages_sf6_late_v2025_UTC.sav'
help,utc_sf6,dat_sf6_b2.age_opt_all
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sweep_sf6_co2_late.sav'
dat_a2 = dat
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_late.sav'
dat_co2_a2 = dat_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late.sav'
dat_sf6_a2 = dat_sf6

in2o = where(finite(n2o_norm),nn2o)
ich4 = where(finite(ch4_all_norm),nch4)
ico2 = where(finite(co2),nco2)
isf6 = where(finite(sf6),nsf6)
iage = where(finite(co2) and finite(sf6),nage)
in2o_a0 = where(finite(n2o_norm_a0),nn2o_a0)
ich4_a0 = where(finite(ch4_all_norm_a0),nch4_a0)
iage_a0 = where(finite(co2_a0) and finite(sf6_a0),nage_a0)
ico2_a0 = where(finite(co2_a0),nco2_a0)
isf6_a0 = where(finite(sf6_a0),nsf6_a0)
in2o_a1 = where(finite(n2o_norm_a1),nn2o_a1)
ich4_a1 = where(finite(ch4_all_norm_a1),nch4_a1)
iage_a1 = where(finite(co2_a1) and finite(sf6_a1),nage_a1)
ico2_a1 = where(finite(co2_a1),nco2_a1)
isf6_a1 = where(finite(sf6_a1),nsf6_a1)
in2o_a = where(finite(n2o_norm_a),nn2o_a)
ich4_a = where(finite(ch4_all_norm_a),nch4_a)
iage_a = where(finite(co2_a) and finite(sf6_a),nage_a)
ico2_a = where(finite(co2_a),nco2_a)
isf6_a = where(finite(sf6_a),nsf6_a)
in2o_b = where(finite(n2o_norm_b),nn2o_b)
ich4_b = where(finite(ch4_all_norm_b),nch4_b)
ico2_b = where(finite(co2_b),nco2_b)
isf6_b = where(finite(sf6_b),nsf6_b)
iage_b = where(finite(co2_all_b) and finite(sf6_b),nage_b)


restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_early_bins.sav'
bins_a0 = bins
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_early_bins.sav'
bins_co2_a0 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_early_bins.sav'
bins_sf6_a0 = bins_sf6
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_optimum_sf6_co2_mid_bins.sav'
bins_a1 = bins
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_mid_bins.sav'
bins_co2_a1 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_mid_bins.sav'
bins_sf6_a1 = bins_sf6
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_late_bins.sav'
bins_a2 = bins
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_late_bins.sav'
bins_co2_a2 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late_bins.sav'
bins_sf6_a2 = bins_sf6

restore,dir+'Balloon/Balloon_mean_ages_co2_early_bins.sav'
bins_co2_b0 = bins_co2
restore,dir+'Balloon/Balloon_mean_ages_sf6_early_bins.sav'
bins_sf6_b0 = bins_sf6
restore,dir+'Balloon/Balloon_mean_ages_co2_late_v2025_bins.sav'
bins_co2_b2 = bins_co2
restore,dir+'Balloon/Balloon_mean_ages_sf6_late_bins.sav'
bins_sf6_b2 = bins_sf6


restore,dir+'Models/WACCM/FWSD_means.sav'
nyw = n_elements(wlat)
nzw = n_elements(walt)
print,years_w[0],years_w[-1]
restore,dir+'swoosh/n2o_merge.sav'
nys = n_elements(slat)
nzs = n_elements(level)
nts = n_elements(swoosh_year)
alt_s = -7.0*alog(level/1e3)
swoosh_n2o_seas = replicate(!values.f_nan,nys,nzs,12) & swoosh_n2o_seas_early = swoosh_n2o_seas
for t = 0, 11 do for z = 0, nzs-1 do for y = 0, nys-1 do begin
  ti = indgen(10)*12+t+360
  swoosh_n2o_seas[y,z,t] = mean(reform(combn2oq_sc_o3_norm[y,z,ti]),/nan)
  ti = indgen(6)*12+t+240
  swoosh_n2o_seas_early[y,z,t] = mean(reform(combn2oq_sc_o3_norm[y,z,ti]),/nan)
endfor
print,swoosh_year[0],swoosh_year[-1]

;  Add on WACCM seasonal cycle to bottom levels of swoosh.
yi = interpol(findgen(nyw),wlat,slat)
zi = interpol(findgen(nzw),walt,alt_s)
n2o_w = fltarr(nys,nzs,12) & n2o_w1 = fltarr(nyw,nzs,12) & n2o_w_early = n2o_w
for t = 0, 11 do begin
  for y = 0, nyw-1 do n2o_w1[y,*,t] = interpolate(n2o_norm_seas[0,y,*,t],zi)
  for z = 0, nzs-1 do n2o_w[*,z,t] = interpolate(n2o_w1[*,z,t],yi)
  for z = 0, nzs-1 do n2o_w_early[*,z,t] = interpolate(n2o_w1[*,z,t],yi)
endfor
n2o_ws = replicate(!values.f_nan,nys,nzs,nts) & n2o_ws1 = replicate(!values.f_nan,nyw,nzs,nts)
for t = 0, nts-38 do begin
  for y = 0, nyw-1 do n2o_ws1[y,*,t] = interpolate(n2o_norm_w[y,*,t+48],zi)
  for z = 0, nzs-1 do n2o_ws[*,z,t] = interpolate(n2o_ws1[*,z,t],yi)
endfor
for t = nts-12, nts-1 do begin
  n2o_ws[*,*,t] = n2o_w[*,*,t-(nts-12)]
endfor

combo_n2o_seas = swoosh_n2o_seas & combo_n2o_seas_early = swoosh_n2o_seas_early
combo_n2o_seas[*,0:5,*] = n2o_w[*,0:5,*]
combo_n2o_seas_early[*,0:5,*] = n2o_w[*,0:5,*]
for z = 4, 5 do begin
  tmp1 = reform(n2o_w[*,z,*]) & tmp2 = reform(swoosh_n2o_seas[*,z,*]) & tmp3 = reform(combo_n2o_seas[*,z,*])
  gd = where(finite(tmp2),ngd)
  if z eq 4 then tmp3[gd] = 0.67*tmp1[gd] + 0.33*tmp2[gd]
  if z eq 5 then tmp3[gd] = 0.33*tmp1[gd] + 0.67*tmp2[gd]
  combo_n2o_seas[*,z,*] = tmp3
  tmp1 = reform(n2o_w[*,z,*]) & tmp2 = reform(swoosh_n2o_seas_early[*,z,*]) & tmp3 = reform(combo_n2o_seas_early[*,z,*])
  gd = where(finite(tmp2),ngd)
  if z eq 4 then tmp3[gd] = 0.67*tmp1[gd] + 0.33*tmp2[gd]
  if z eq 5 then tmp3[gd] = 0.33*tmp1[gd] + 0.67*tmp2[gd]
  combo_n2o_seas_early[*,z,*] = tmp3
endfor
combo_n2o_s = combn2oq_sc_o3_norm
combo_n2o_s[*,0:5,*] = n2o_ws[*,0:5,*]
for z = 4, 5 do begin
  tmp1 = reform(n2o_ws[*,z,*]) & tmp2 = reform(combn2oq_sc_o3_norm[*,z,*]) & tmp3 = reform(combo_n2o_s[*,z,*])
  gd = where(finite(tmp2),ngd)
  if z eq 4 then tmp3[gd] = 0.67*tmp1[gd] + 0.33*tmp2[gd]
  if z eq 5 then tmp3[gd] = 0.33*tmp1[gd] + 0.67*tmp2[gd]
  combo_n2o_s[*,z,*] = tmp3
endfor


;  Fill in missing lats.
;chk = where(~finite(lat),nchk)
;diffs = chk - shift(chk,1)
;ni = where(diffs gt 1,nn)
;ni1 = where(diffs eq 1,nn)
;for i = 0, nn-1 do begin
;  bi = 
;endfor
;print,nn


nyg = 36
dy = 5.
lat_grid = findgen(nyg)*dy-85.
nz = 120
dz = 0.25
alt_grid = findgen(nz)*dz+6.
nz2 = 60
dz2 = 0.5
alt_grid2 = findgen(nz2)*dz2+6.
nz3 = 30
dz3 = 1.0
alt_grid3 = findgen(nz3)*dz3+6.0

do_elats = (stage eq 'elats')

if do_elats then begin

  elat_m_adj = elat_m & elat_m_adj_a = elat_m_a & elat_m_adj_b = elat_m_b & elat_m_adj_a0 = elat_m_a0 & elat_m_adj_a1 = elat_m_a1
  elat_n2o_fit = replicate(!values.f_nan,2,nz,12) & aa = replicate(!values.f_nan,7,nyg,nz2,12) & n2o_norm_grid_mon_avg = aa & age_co2_grid_mon_avg = aa & n2o_norm_sf6_grid_mon_avg = aa
  bb = replicate(!values.f_nan,3,nyg,nz,4) & combo_n2o_seas_hires = replicate(!values.f_nan,nys,nz,12) & age_co2_grid_seas = bb & n2o_norm_grid_seas = bb & age_sf6_grid_mon_avg = aa
  n2o_norm_sf6_grid_seas = bb & age_sf6_grid_seas = bb & n2o_norm_grid_mon_avg_a = aa & age_co2_grid_mon_avg_a = aa & n2o_norm_sf6_grid_mon_avg_a = aa & age_sf6_grid_mon_avg_a = aa
  dd = replicate(!values.f_nan,nyg,nz2,12) & age_combo_grid_mon_avg = aa & n2o_norm_combo_grid_mon_avg = aa & n2o_norm_combo_grid_seas = bb & age_combo_grid_seas = bb
  n2o_norm_grid_elat_orig_mon_avg = aa & age_co2_grid_elat_orig_mon_avg = aa & combo_n2o_s_hires = replicate(!values.f_nan,nys,nz,nts) & n2o_norm_grid_mon_avg_b = aa
  age_co2_grid_mon_avg_b = aa & combo_n2o_seas_early_hires = combo_n2o_seas_hires & n2o_norm_sf6_grid_mon_avg_b = aa & age_sf6_grid_mon_avg_b = aa & n2o_norm_both_grid_mon_avg_b = aa
  age_both_grid_mon_avg_b = aa & age_grid_mon_avg_b = aa & n2o_norm_opt_grid_mon_avg = aa & age_opt_grid_mon_avg = aa & n2o_norm_opt_grid_mon_avg_a = aa & age_opt_grid_mon_avg_a = aa
  n2o_norm_grid_mon_avg_a0 = aa & age_co2_grid_mon_avg_a0 = aa & n2o_norm_sf6_grid_mon_avg_a0 = aa & age_sf6_grid_mon_avg_a0 = aa & n2o_norm_grid_seas_a0 = bb & n2o_norm_opt_grid_mon_avg_a0 = aa
  age_opt_grid_mon_avg_a0 = aa
  n2o_norm_grid_mon_avg_a1 = aa & age_co2_grid_mon_avg_a1 = aa & n2o_norm_sf6_grid_mon_avg_a1 = aa & age_sf6_grid_mon_avg_a1 = aa & n2o_norm_grid_seas_a1 = bb
  n2o_norm_opt_grid_mon_avg_a1 = aa & age_opt_grid_mon_avg_a1 = aa
  for t = 0, 11 do for z = 10, nz-1 do begin
  
    zis = interpol(findgen(nzs),alt_s,alt_grid[z])
    for y = 0, nys-1 do combo_n2o_seas_early_hires[y,z,t] = interpolate(combo_n2o_seas_early[y,*,t],zis)
    zis = interpol(findgen(nzs),alt_s,alt_grid[z])
    for y = 0, nys-1 do combo_n2o_seas_hires[y,z,t] = interpolate(combo_n2o_seas[y,*,t],zis)
  
    gd = where(months_all[in2o] eq t+1 and elat_m[in2o] gt 0 and alt[in2o] ge alt_grid[z]-dz/2. and alt[in2o] lt alt_grid[z]+dz/2.,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(nys/2-1),combo_n2o_seas_hires[nys/2+1:-1,z,t],n2o_norm[in2o[gd]])
      elat_m_adj[in2o[gd]] = interpolate(slat[nys/2+1:-1],yi)

      chk = where(elat_m_adj[in2o[gd]] ge 85,nchk)
      if nchk gt 0 then elat_m_adj[in2o[gd[chk]]] = 85.

      chk = where(n2o_norm[in2o[gd]] gt max(combo_n2o_seas_hires[nys/2+1:-1,z,t]),nchk)
      if nchk gt 0 then elat_m_adj[in2o[gd[chk]]] = 0. ;min(elat_m_adj[in2o[gd]])

      chk = where(n2o_norm[in2o[gd]] lt min(combo_n2o_seas_hires[nys/2+1:-1,z,t]),nchk)
      if nchk gt 0 then elat_m_adj[in2o[gd[chk]]] = 85. ;max(elat_m_adj[in2o[gd]])
    endif

    gd = where(months_all_a[in2o_a] eq t+1 and elat_m_a[in2o_a] gt 0 and alt_a[in2o_a] ge alt_grid[z]-dz/2. and alt_a[in2o_a] lt alt_grid[z]+dz/2.,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(nys/2-1),combo_n2o_seas_hires[nys/2+1:-1,z,t],n2o_norm_a[in2o_a[gd]])
      elat_m_adj_a[in2o_a[gd]] = interpolate(slat[nys/2+1:-1],yi)

      chk = where(elat_m_adj_a[in2o_a[gd]] ge 85,nchk)
      if nchk gt 0 then elat_m_adj_a[in2o_a[gd[chk]]] = 85.

      chk = where(n2o_norm_a[in2o_a[gd]] gt max(combo_n2o_seas_hires[nys/2+1:-1,z,t]),nchk)
      if nchk gt 0 then elat_m_adj_a[in2o_a[gd[chk]]] = 0.

      chk = where(n2o_norm_a[in2o_a[gd]] lt min(combo_n2o_seas_hires[nys/2+1:-1,z,t]),nchk)
      if nchk gt 0 then elat_m_adj_a[in2o_a[gd[chk]]] = 85.
    endif

    gd = where(months_all_b[in2o_b] eq t+1 and finite(n2o_norm_b[in2o_b]) and finite(lat_b[in2o_b]) and alt_b[in2o_b] ge alt_grid[z]-dz/2. and alt_b[in2o_b] lt alt_grid[z]+dz/2.,ngd)
    if ngd gt 1 then begin

      chk_s = where(lat_b[in2o_b[gd]] lt 0,nchk_s)
      if nchk_s gt 0 then begin
        yi = interpol(findgen(nys/2+1),combo_n2o_seas_early_hires[0:nys/2,z,t],n2o_norm_b[in2o_b[gd[chk_s]]])
        elat_m_adj_b[in2o_b[gd[chk_s]]] = interpolate(slat[0:nys/2],yi)

        chk = where(elat_m_adj_b[in2o_b[gd[chk_s]]] le -85,nchk)
        if nchk gt 0 then elat_m_adj_b[in2o_b[gd[chk_s[chk]]]] = -85.

        chk = where(n2o_norm_b[in2o_b[gd[chk_s]]] gt max(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_b[in2o_b[gd[chk_s[chk]]]] = 0.

        chk = where(n2o_norm_b[in2o_b[gd[chk_s]]] lt min(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_b[in2o_b[gd[chk_s[chk]]]] = -85.
      endif

      chk_n = where(lat_b[in2o_b[gd]] ge 0,nchk_n)
      if nchk_n gt 0 then begin
        yi = interpol(findgen(nys/2),combo_n2o_seas_early_hires[nys/2:-1,z,t],n2o_norm_b[in2o_b[gd[chk_n]]])
        elat_m_adj_b[in2o_b[gd[chk_n]]] = interpolate(slat[nys/2:-1],yi)

        chk = where(elat_m_adj_b[in2o_b[gd[chk_n]]] ge 85,nchk)
        if nchk gt 0 then elat_m_adj_b[in2o_b[gd[chk_n[chk]]]] = 85.

        chk = where(n2o_norm_b[in2o_b[gd[chk_n]]] gt max(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_b[in2o_b[gd[chk_n[chk]]]] = 0.

        chk = where(n2o_norm_b[in2o_b[gd[chk_n]]] lt min(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_b[in2o_b[gd[chk_n[chk]]]] = 85.
      endif
    endif

;    gd = where(months_all_b[isf6_b] eq t+1 and finite(lat_b[isf6_b]) and alt_b[isf6_b] ge alt_grid[z]-dz/2. and alt_b[isf6_b] lt alt_grid[z]+dz/2.,ngd)
;    if ngd gt 1 then begin
;      if t eq 1 or t eq 10 then begin
;  
;        chk = where(lat_b[isf6_b[gd]] lt 0)
;        yi = interpol(findgen(nys/2+1),combo_n2o_seas_early_hires[0:nys/2,z,t],n2o_norm_b[isf6_b[gd]])
;        elat_m_adj_b[isf6_b[gd]] = interpolate(slat[0:nys/2],yi)
;  
;        chk = where(n2o_norm_b[isf6_b[gd]] gt max(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
;        if nchk gt 0 then elat_m_adj_b[isf6_b[gd[chk]]] = elat_m_b[isf6_b[gd[chk]]]
;  
;        chk1 = where(lat_b[isf6_b[gd]] gt 0,nchk1)
;        if nchk1 gt 0 then begin
;          yi = interpol(findgen(nys/2),combo_n2o_seas_early_hires[nys/2:-1,z,t],n2o_norm_b[isf6_b[gd[chk1]]])
;          elat_m_adj_b[isf6_b[gd[chk1]]] = interpolate(slat[nys/2:-1],yi)
;  
;          chk = where(n2o_norm_b[isf6_b[gd[chk1]]] gt max(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
;          if nchk gt 0 then elat_m_adj_b[isf6_b[gd[chk1[chk]]]] = min(elat_m_adj_b[isf6_b[gd]])
;        endif
;        ;      if z eq 60 then print,t,alt_grid[z],ngd,elat_m_adj_b[isf6_b[gd]]
;      endif else begin
;        yi = interpol(findgen(nys/2),combo_n2o_seas_early_hires[nys/2:-1,z,t],n2o_norm_b[isf6_b[gd]])
;        elat_m_adj_b[isf6_b[gd]] = interpolate(slat[nys/2:-1],yi)
;  
;        chk = where(n2o_norm_b[isf6_b[gd]] gt max(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
;        if nchk gt 0 then elat_m_adj_b[isf6_b[gd[chk]]] = min(elat_m_adj_b[isf6_b[gd]])
;      endelse
;      chk = where(elat_m_adj_b[isf6_b[gd]] ge 85,nchk)
;      if nchk gt 0 then elat_m_adj_b[isf6_b[gd[chk]]] = 85.
;  
;      chk = where(n2o_norm_b[isf6_b[gd]] lt min(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
;      if nchk gt 0 then elat_m_adj_b[isf6_b[gd[chk]]] = 85.
;    endif
  
    gd = where(months_all_a0[in2o_a0] eq t+1 and finite(n2o_norm_a0[in2o_a0]) and finite(lat_a0[in2o_a0]) and alt_a0[in2o_a0] ge alt_grid[z]-dz/2. and alt_a0[in2o_a0] lt alt_grid[z]+dz/2.,ngd)
    if ngd gt 1 then begin

      chk_s = where(lat_a0[in2o_a0[gd]] lt 0,nchk_s)
      if nchk_s gt 0 then begin
        yi = interpol(findgen(nys/2+1),combo_n2o_seas_early_hires[0:nys/2,z,t],n2o_norm_a0[in2o_a0[gd[chk_s]]])
        elat_m_adj_a0[in2o_a0[gd[chk_s]]] = interpolate(slat[0:nys/2],yi)

        chk = where(elat_m_adj_a0[in2o_a0[gd[chk_s]]] le -85,nchk)
        if nchk gt 0 then elat_m_adj_a0[in2o_a0[gd[chk_s[chk]]]] = -85.

        chk = where(n2o_norm_a0[in2o_a0[gd[chk_s]]] gt max(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a0[in2o_a0[gd[chk_s[chk]]]] = 0.

        chk = where(n2o_norm_a0[in2o_a0[gd[chk_s]]] lt min(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a0[in2o_a0[gd[chk_s[chk]]]] = -85.
;        if t eq 8 and z eq 32 then print,alt_grid[z],nchk_s,n2o_norm_a0[in2o_a0[gd[chk_s]]],lat_a0[in2o_a0[gd[chk_s]]],elat_m_adj_a0[in2o_a0[gd[chk_s]]]
      endif

      chk_n = where(lat_a0[in2o_a0[gd]] ge 0,nchk_n)
      if nchk_n gt 0 then begin
        yi = interpol(findgen(nys/2),combo_n2o_seas_early_hires[nys/2:-1,z,t],n2o_norm_a0[in2o_a0[gd[chk_n]]])
        elat_m_adj_a0[in2o_a0[gd[chk_n]]] = interpolate(slat[nys/2:-1],yi)

        chk = where(elat_m_adj_a0[in2o_a0[gd[chk_n]]] ge 85,nchk)
        if nchk gt 0 then elat_m_adj_a0[in2o_a0[gd[chk_n[chk]]]] = 85.

        chk = where(n2o_norm_a0[in2o_a0[gd[chk_n]]] gt max(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a0[in2o_a0[gd[chk_n[chk]]]] = 0.

        chk = where(n2o_norm_a0[in2o_a0[gd[chk_n]]] lt min(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a0[in2o_a0[gd[chk_n[chk]]]] = 85.
;        if t eq 8 and z eq 32 then begin
;          test = where(elat_m_adj_a0[in2o_a0[gd[chk_n]]] ge 60 and elat_m_adj_a0[in2o_a0[gd[chk_n]]] le 70)
;;          print,alt_grid[z],nchk_n,n2o_norm_a0[in2o_a0[gd[chk_n]]],elat_m_adj_a0[in2o_a0[gd[chk_n]]]
;          print,alt_grid[z],nchk_n,n2o_norm_a0[in2o_a0[gd[chk_n[test]]]],elat_m_adj_a0[in2o_a0[gd[chk_n[test]]]]
;        endif
      endif
    endif

    gd = where(months_all_a1[in2o_a1] eq t+1 and finite(n2o_norm_a1[in2o_a1]) and finite(lat_a1[in2o_a1]) and alt_a1[in2o_a1] ge alt_grid[z]-dz/2. and alt_a1[in2o_a1] lt alt_grid[z]+dz/2.,ngd)
    if ngd gt 1 then begin

      chk_s = where(lat_a1[in2o_a1[gd]] lt 0,nchk_s)
      if nchk_s gt 0 then begin
        yi = interpol(findgen(nys/2+1),combo_n2o_seas_early_hires[0:nys/2,z,t],n2o_norm_a1[in2o_a1[gd[chk_s]]])
        elat_m_adj_a1[in2o_a1[gd[chk_s]]] = interpolate(slat[0:nys/2],yi)

        chk = where(elat_m_adj_a1[in2o_a1[gd[chk_s]]] le -85,nchk)
        if nchk gt 0 then elat_m_adj_a1[in2o_a1[gd[chk_s[chk]]]] = -85.

        chk = where(n2o_norm_a1[in2o_a1[gd[chk_s]]] gt max(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a1[in2o_a1[gd[chk_s[chk]]]] = 0.

        chk = where(n2o_norm_a1[in2o_a1[gd[chk_s]]] lt min(combo_n2o_seas_early_hires[0:nys/2,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a1[in2o_a1[gd[chk_s[chk]]]] = -85.
      endif

      chk_n = where(lat_a1[in2o_a1[gd]] ge 0,nchk_n)
      if nchk_n gt 0 then begin
        yi = interpol(findgen(nys/2),combo_n2o_seas_early_hires[nys/2:-1,z,t],n2o_norm_a1[in2o_a1[gd[chk_n]]])
        elat_m_adj_a1[in2o_a1[gd[chk_n]]] = interpolate(slat[nys/2:-1],yi)

        chk = where(elat_m_adj_a1[in2o_a1[gd[chk_n]]] ge 85,nchk)
        if nchk gt 0 then elat_m_adj_a1[in2o_a1[gd[chk_n[chk]]]] = 85.

        chk = where(n2o_norm_a1[in2o_a1[gd[chk_n]]] gt max(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a1[in2o_a1[gd[chk_n[chk]]]] = 0.

        chk = where(n2o_norm_a1[in2o_a1[gd[chk_n]]] lt min(combo_n2o_seas_early_hires[nys/2:-1,z,t]),nchk)
        if nchk gt 0 then elat_m_adj_a1[in2o_a1[gd[chk_n[chk]]]] = 85.
      endif
    endif

  endfor
  
  for t = 360, nts-1 do begin
  ;  i_arr = ico2_a
  ;  gd = where(months_all_a[i_arr] eq ((t mod 12) + 1) and times_all_a[i_arr] ge swoosh_year[t]-1./24. and times_all_a[i_arr] le swoosh_year[t]+1./24.,ngd)
  ;  if ngd gt 1 then begin
  ;    for z = 10, nz-1 do begin
  ;      gd2 = where(elat_m_a[i_arr[gd]] gt 0 and alt_a[i_arr[gd]] ge alt_grid[z]-dz/2. and alt_a[i_arr[gd]] lt alt_grid[z]+dz/2.,ngd2)
  ;      if ngd2 gt 1 then begin
  ;        zis = interpol(findgen(nzs),alt_s,alt_grid[z])
  ;        for y = 0, nys-1 do combo_n2o_s_hires[y,z,t] = interpolate(combo_n2o_s[y,*,t],zis)
  ;        yi = interpol(findgen(nys/2-1),combo_n2o_s_hires[nys/2+1:-1,z,t],n2o_norm_a[i_arr[gd[gd2]]])
  ;        elat_m_adj_a[i_arr[gd[gd2]]] = interpolate(slat[nys/2+1:-1],yi)
  ;        chk_n = where(finite(n2o_norm_a[i_arr[gd[gd2]]]),nchk_n)
  ;        chk_c = where(finite(co2_a[i_arr[gd[gd2]]]),nchk_c)
  ;;        if z ge 47 and z le 49 then print,alt_grid[z],t,swoosh_year[t],ngd2,nchk_c,nchk_n
  ;    
  ;        chk = where(elat_m_adj_a[i_arr[gd[gd2]]] ge 85,nchk)
  ;        if nchk gt 0 then elat_m_adj_a[i_arr[gd[gd2[chk]]]] = elat_m_a[i_arr[gd[gd2[chk]]]]
  ;    
  ;        chk = where(n2o_norm_a[i_arr[gd[gd2]]] gt max(combo_n2o_s_hires[nys/2+1:-1,z,t]),nchk)
  ;        if nchk gt 0 then elat_m_adj_a[i_arr[gd[gd2[chk]]]] = min(elat_m_adj_a[i_arr[gd[gd2]]])
  ;    
  ;        chk = where(n2o_norm_a[i_arr[gd[gd2]]] lt min(combo_n2o_s_hires[nys/2+1:-1,z,t]),nchk)
  ;        if nchk gt 0 then elat_m_adj_a[i_arr[gd[gd2[chk]]]] = max(elat_m_adj_a[i_arr[gd[gd2]]])
  ;      endif
  ;    endfor
  ;  endif
  ;  gd = where(months_all_a[isf6_a] eq ((t mod 12) + 1) and times_all_a[isf6_a] ge swoosh_year[t]-1./24. and times_all_a[isf6_a] le swoosh_year[t]+1./24.,ngd)
  ;  if ngd gt 1 then begin
  ;    for z = 10, nz-1 do begin
  ;      gd2 = where(elat_m_a[isf6_a[gd]] gt 0 and alt_a[isf6_a[gd]] ge alt_grid[z]-dz/2. and alt_a[isf6_a[gd]] le alt_grid[z]+dz/2.,ngd2)
  ;      if ngd2 gt 1 then begin
  ;        zis = interpol(findgen(nzs),alt_s,alt_grid[z])
  ;        for y = 0, nys-1 do combo_n2o_s_hires[y,z,t] = interpolate(combo_n2o_s[y,*,t],zis)
  ;        yi = interpol(findgen(nys/2-1),combo_n2o_s_hires[nys/2+1:-1,z,t],n2o_norm_a[isf6_a[gd[gd2]]])
  ;        elat_m_adj_a[isf6_a[gd[gd2]]] = interpolate(slat[nys/2+1:-1],yi)
  ;    ;    if z eq 49 then print,alt_grid[z],t,swoosh_year[t],ngd
  ;    
  ;        chk = where(elat_m_adj_a[isf6_a[gd[gd2]]] ge 85,nchk)
  ;        if nchk gt 0 then elat_m_adj_a[isf6_a[gd[gd2[chk]]]] = elat_m_a[isf6_a[gd[gd2[chk]]]]
  ;    
  ;        chk = where(n2o_norm_a[isf6_a[gd[gd2]]] gt max(combo_n2o_s_hires[nys/2+1:-1,z,t],mi),nchk)
  ;        if nchk gt 0 then elat_m_adj_a[isf6_a[gd[gd2[chk]]]] = slat[mi+nys/2+1]
  ;    
  ;        chk = where(n2o_norm_a[isf6_a[gd[gd2]]] lt min(combo_n2o_s_hires[nys/2+1:-1,z,t]),nchk)
  ;        if nchk gt 0 then elat_m_adj_a[isf6_a[gd[gd2[chk]]]] = max(elat_m_adj_a[isf6_a[gd[gd2]]])
  ;      endif
  ;    endfor
  ;  endif
  endfor
  
  ;m = 1 & z = 60
  ;gd_b = where(months_all_b[ico2_b] eq m+1 and finite(lat_b[ico2_b]) and alt_b[ico2_b] ge alt_grid[z]-dz/2. and alt_b[ico2_b] lt alt_grid[z]+dz/2.,ngdb)
  ;print,ngdb,elat_m_adj_b[ico2_b[gd_b]]
  
  
  for t = 0, 11 do begin
    ii = in2o
    gd = where(months_all[ii] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj[ii[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt[ii[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm[ii[gd]]) and elat_m_adj[ii[gd]] gt 0 and lat[ii[gd]] gt 0 and alt[ii[gd]] ge alt_grid2[z]-dz2/2. and alt[ii[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm[ii[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
        ii_all = where(elat_m[ii[gd]] gt 0 and lat[ii[gd]] gt 0 and alt[ii[gd]] ge alt_grid2[z]-dz2/2. and alt[ii[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m[ii[gd]] ge lat_grid[y]-dy/2. and elat_m[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm[ii[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_grid_elat_orig_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all[ico2] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj[ico2[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt[ico2[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt[ico2[gd]] lt tp_elat,nchk)
;      elat_m_adj[ico2[gd[chk]]] = lat[ico2[gd[chk]]]
;      chk = where(elat_m_adj[ico2[gd]] gt 0 and elat_m_adj[ico2[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm[ico2[gd]]) and elat_m_adj[ico2[gd]] gt 0 and lat[ico2[gd]] gt 0 and alt[ico2[gd]] ge alt_grid2[z]-dz2/2. and alt[ico2[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj[ico2[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[ico2[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
;          stats = moment(n2o_norm[ico2[gd[ii_all]]],sd=sd,/nan)
;          n2o_norm_grid_mon_avg[*,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_co2_b2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_co2_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
        ii_all = where(elat_m[ico2[gd]] gt 0 and lat[ico2[gd]] gt 0 and alt[ico2[gd]] ge alt_grid2[z]-dz2/2. and alt[ico2[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m[ico2[gd]] ge lat_grid[y]-dy/2. and elat_m[ico2[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
;          stats = moment(n2o_norm[ico2[gd[ii_all]]],sd=sd,/nan)
;          n2o_norm_grid_elat_orig_mon_avg[*,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_co2_b2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_co2_grid_elat_orig_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all[isf6] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj[isf6[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt[isf6[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt[isf6[gd]] lt tp_elat,nchk)
;      elat_m_adj[isf6[gd[chk]]] = lat[isf6[gd[chk]]]
;      chk = where(elat_m_adj[isf6[gd]] gt 0 and elat_m_adj[isf6[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm[isf6[gd]]) and elat_m_adj[isf6[gd]] gt 0 and lat[isf6[gd]] gt 0 and alt[isf6[gd]] ge alt_grid2[z]-dz2/2. and alt[isf6[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj[isf6[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[isf6[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm[isf6[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_sf6_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_sf6_b2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_sf6_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all[iage] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj[iage[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt[iage[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt[iage[gd]] lt tp_elat,nchk)
;      elat_m_adj[iage[gd[chk]]] = lat[iage[gd[chk]]]
;      chk = where(elat_m_adj[iage[gd]] gt 0 and elat_m_adj[iage[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm[iage[gd]]) and elat_m_adj[iage[gd]] gt 0 and lat[iage[gd]] gt 0 and alt[iage[gd]] ge alt_grid2[z]-dz2/2. and alt[iage[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj[iage[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[iage[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm[iage[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_opt_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_b2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_opt_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif

    ;  Aircraft late.
    ii = in2o_a
    gd = where(months_all_a[ii] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a[ii[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a[ii[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm_a[ii[gd]]) and elat_m_adj_a[ii[gd]] gt 0 and lat_a[ii[gd]] gt 0 and alt_a[ii[gd]] ge alt_grid2[z]-dz2/2. and $
          alt_a[ii[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a[ii[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_a[ico2_a] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a[ico2_a[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a[ico2_a[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt_a[ico2_a[gd]] lt tp_elat,nchk)
;      elat_m_adj_a[ico2_a[gd[chk]]] = lat_a[ico2_a[gd[chk]]]
;      chk = where(elat_m_adj_a[ico2_a[gd]] gt 0 and elat_m_adj_a[ico2_a[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm_a[ico2_a[gd]]) and elat_m_adj_a[ico2_a[gd]] gt 0 and lat_a[ico2_a[gd]] gt 0 and alt_a[ico2_a[gd]] ge alt_grid2[z]-dz2/2. and $
          alt_a[ico2_a[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[ico2_a[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[ico2_a[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
;          stats = moment(n2o_norm_a[ico2_a[gd[ii_all]]],sd=sd,/nan)
;          n2o_norm_grid_mon_avg_a[*,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_co2_a2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_co2_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_a[iage_a] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a[iage_a[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a[iage_a[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt_a[iage_a[gd]] lt tp_elat,nchk)
;      elat_m_adj_a[iage_a[gd[chk]]] = lat_a[iage_a[gd[chk]]]
;      chk = where(elat_m_adj_a[iage_a[gd]] gt 0 and elat_m_adj_a[iage_a[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm_a[iage_a[gd]]) and elat_m_adj_a[iage_a[gd]] gt 0 and lat_a[iage_a[gd]] gt 0 and alt_a[iage_a[gd]] ge alt_grid2[z]-dz2/2. and $
          alt_a[iage_a[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[iage_a[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[iage_a[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a[iage_a[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_opt_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_a2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_opt_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_a[isf6_a] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a[isf6_a[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a[isf6_a[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt_a[isf6_a[gd]] lt tp_elat,nchk)
;      elat_m_adj_a[isf6_a[gd[chk]]] = lat_a[isf6_a[gd[chk]]]
;      chk = where(elat_m_adj_a[isf6_a[gd]] gt 0 and elat_m_adj_a[isf6_a[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(n2o_norm_a[isf6_a[gd]]) and elat_m_adj_a[isf6_a[gd]] gt 0 and lat_a[isf6_a[gd]] gt 0 and alt_a[isf6_a[gd]] ge alt_grid2[z]-dz2/2. and $
          alt_a[isf6_a[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[isf6_a[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[isf6_a[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a[isf6_a[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_sf6_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_sf6_a2.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_sf6_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
  ;        if t eq 2 and z eq 24 then print,elat_m_adj_a[isf6_a[gd[ii_all]]],n2o_norm_a[isf6_a[gd[ii_all]]]
        endif
      endfor
    endif
  
    ;  Balloon early.
    ii = in2o_b
    gd = where(months_all_b[ii] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_b[ii[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_b[ii[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_b[ii[gd]]) and alt_b[ii[gd]] ge alt_grid2[z]-dz2/2. and alt_b[ii[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_b[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_b[ii[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_b[ico2_b] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_b[ico2_b[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_b[ico2_b[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt_b[ico2_b[gd]] lt tp_elat,nchk)
;      elat_m_adj_b[ico2_b[gd[chk]]] = lat_b[ico2_b[gd[chk]]]
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_b[ico2_b[gd]]) and alt_b[ico2_b[gd]] ge alt_grid2[z]-dz2/2. and alt_b[ico2_b[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_b[ico2_b[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[ico2_b[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
;          stats = moment(n2o_norm_b[ico2_b[gd[ii_all]]],sd=sd,/nan)
;          n2o_norm_grid_mon_avg_b[*,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_co2_b0.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_co2_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif  
    gd = where(months_all_b[isf6_b] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_b[isf6_b[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_b[isf6_b[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt_b[isf6_b[gd]] lt tp_elat,nchk)
;      elat_m_adj_b[isf6_b[gd[chk]]] = lat_b[isf6_b[gd[chk]]]
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_b[isf6_b[gd]]) and alt_b[isf6_b[gd]] ge alt_grid2[z]-dz2/2. and alt_b[isf6_b[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_b[isf6_b[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[isf6_b[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_b[isf6_b[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_sf6_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_sf6_b0.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_sf6_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif  
    gd = where(months_all_b[iage_b] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_b[iage_b[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_b[iage_b[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
;      chk = where(alt_b[iage_b[gd]] lt tp_elat,nchk)
;      elat_m_adj_b[iage_b[gd[chk]]] = lat_b[iage_b[gd[chk]]]
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_b[iage_b[gd]]) and alt_b[iage_b[gd]] ge alt_grid2[z]-dz2/2. and alt_b[iage_b[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_b[iage_b[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[iage_b[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(dat_sf6_b0.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
  
    ;  Aircraft early.
    ii = in2o_a0
    gd = where(months_all_a0[ii] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a0[ii[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a0[ii[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_a0[ii[gd]]) and finite(n2o_norm_a0[ii[gd]]) and alt_a0[ii[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[ii[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a0[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a0[ii[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_grid_mon_avg_a0[0:2,y,z,t] = [stats[0],sd,tot_i]
          if t eq 3 and z eq 20 then print,alt_grid2[z],lat_grid[y],stats[0],sd,tot_i
        endif
      endfor
    endif
    gd = where(months_all_a0[ico2_a0] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a0[ico2_a0[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a0[ico2_a0[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_a0[ico2_a0[gd]]) and finite(n2o_norm_a0[ico2_a0[gd]]) and alt_a0[ico2_a0[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[ico2_a0[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a0[ico2_a0[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[ico2_a0[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          datsub = dat_co2_a0.age_opt_all[gd[ii_all]]
          datusub = dat_co2_a0.age_range_all[*,gd[ii_all]]
          age_err = 0.5 * (datusub[1,*] - datusub[0,*])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * datsub) / total(age_weights)
          ;  Weighted sample variance.
          wsv = total(age_weights * (datsub - wmean)^2) / (total(age_weights) - 1)
          ;  Standard error of the weighted mean.
          se = sqrt(wsv / total(age_weights))
          stats = moment(datsub,sd=sd,/nan)
          age_co2_grid_mon_avg_a0[*,y,z,t] = [stats[0],sd,tot_i,wmean,err,wsv,se]
          if t eq 3 and z eq 20 then print,alt_grid2[z],lat_grid[y],stats[0],sd,tot_i,wmean,err,wsv,se
        endif
      endfor
    endif  
    gd = where(months_all_a0[iage_a0] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a0[iage_a0[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a0[iage_a0[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      chk = where(elat_m_adj_a0[iage_a0[gd]] gt 0 and elat_m_adj_a0[iage_a0[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(lat_a0[iage_a0[gd]] gt 0 and finite(n2o_norm_a0[iage_a0[gd]]) and alt_a0[iage_a0[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[iage_a0[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a0[iage_a0[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[iage_a0[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a0[iage_a0[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_opt_grid_mon_avg_a0[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(datsub,sd=sd,/nan)
          age_opt_grid_mon_avg_a0[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif  
    gd = where(months_all_a0[isf6_a0] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a0[isf6_a0[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a0[isf6_a0[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_a0[isf6_a0[gd]]) and finite(n2o_norm_a0[isf6_a0[gd]]) and alt_a0[isf6_a0[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[isf6_a0[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a0[isf6_a0[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[isf6_a0[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a0[isf6_a0[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_sf6_grid_mon_avg_a0[0:2,y,z,t] = [stats[0],sd,tot_i]
          datsub = dat_sf6_a0.age_opt_all[gd[ii_all]]
          datusub = dat_sf6_a0.age_range_all[*,gd[ii_all]]
          age_err = 0.5 * (datusub[1,*] - datusub[0,*])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * datsub) / total(age_weights)
          wsv = total(age_weights * (datsub - wmean)^2) / (total(age_weights) - 1)
          se = sqrt(wsv / total(age_weights))
          stats = moment(datsub,sd=sd,/nan)
          age_sf6_grid_mon_avg_a0[*,y,z,t] = [stats[0],sd,tot_i,wmean,err,wsv,se]
        endif
      endfor
    endif
  
    ;  Aircraft mid.
    ii = in2o_a1
    gd = where(months_all_a1[ii] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a1[ii[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a1[ii[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_a1[ii[gd]]) and finite(n2o_norm_a1[ii[gd]]) and alt_a1[ii[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[ii[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a1[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a1[ii[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_a1[ico2_a1] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a1[ico2_a1[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a1[ico2_a1[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_a1[ico2_a1[gd]]) and finite(n2o_norm_a1[ico2_a1[gd]]) and alt_a1[ico2_a1[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[ico2_a1[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a1[ico2_a1[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[ico2_a1[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
;          stats = moment(n2o_norm_a1[ico2_a1[gd[ii_all]]],sd=sd,/nan)
;          n2o_norm_grid_mon_avg_a1[*,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_co2_a1.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_co2_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_a1[iage_a1] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a1[iage_a1[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a1[iage_a1[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      chk = where(elat_m_adj_a1[iage_a1[gd]] gt 0 and elat_m_adj_a1[iage_a1[gd]] lt 20,nchk)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(lat_a1[iage_a1[gd]] gt 0 and finite(n2o_norm_a1[iage_a1[gd]]) and alt_a1[iage_a1[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[iage_a1[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a1[iage_a1[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[iage_a1[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a1[iage_a1[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_opt_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_a1.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_opt_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif
    gd = where(months_all_a1[isf6_a1] eq t+1,ngd)
    if ngd gt 1 then begin
      yi = interpol(findgen(ny_tp),lat_tp,elat_m_adj_a1[isf6_a1[gd]])
      zi = interpol(findgen(ny_tp),tpause_mon_alt[*,t],alt_a1[isf6_a1[gd]])
      tp_elat = interpolate(tpause_mon_alt[*,t],yi)
      for z = 6, nz2-1 do for y = 0, nyg-1 do begin
        ii_all = where(finite(lat_a1[isf6_a1[gd]]) and finite(n2o_norm_a1[isf6_a1[gd]]) and alt_a1[isf6_a1[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[isf6_a1[gd]] lt alt_grid2[z]+dz2/2. and $
          elat_m_adj_a1[isf6_a1[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[isf6_a1[gd]] lt lat_grid[y]+dy/2.,tot_i)
        if tot_i gt 1 then begin
          stats = moment(n2o_norm_a1[isf6_a1[gd[ii_all]]],sd=sd,/nan)
          n2o_norm_sf6_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
          stats = moment(dat_sf6_a1.age_opt_all[gd[ii_all]],sd=sd,/nan)
          age_sf6_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
        endif
      endfor
    endif

    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      for j = 0, 1 do begin
        n2o_norm_combo_grid_mon_avg[j,y,z,t] = mean([n2o_norm_grid_mon_avg[j,y,z,t],n2o_norm_sf6_grid_mon_avg[j,y,z,t],n2o_norm_grid_mon_avg_a[j,y,z,t],n2o_norm_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
        age_combo_grid_mon_avg[j,y,z,t] = mean([age_co2_grid_mon_avg[j,y,z,t],age_sf6_grid_mon_avg[j,y,z,t],age_co2_grid_mon_avg_a[j,y,z,t],age_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
        n2o_norm_both_grid_mon_avg_b[j,y,z,t] = mean([n2o_norm_grid_mon_avg_b[j,y,z,t],n2o_norm_sf6_grid_mon_avg_b[j,y,z,t]],/nan)
        age_both_grid_mon_avg_b[j,y,z,t] = mean([age_co2_grid_mon_avg_b[j,y,z,t],age_grid_mon_avg_b[j,y,z,t]],/nan)
      endfor
      j = 2
      n2o_norm_combo_grid_mon_avg[j,y,z,t] = total([n2o_norm_grid_mon_avg[j,y,z,t],n2o_norm_sf6_grid_mon_avg[j,y,z,t],n2o_norm_grid_mon_avg_a[j,y,z,t],n2o_norm_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
      age_combo_grid_mon_avg[j,y,z,t] = total([age_co2_grid_mon_avg[j,y,z,t],age_sf6_grid_mon_avg[j,y,z,t],age_co2_grid_mon_avg_a[j,y,z,t],age_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
      n2o_norm_both_grid_mon_avg_b[j,y,z,t] = total([n2o_norm_grid_mon_avg_b[j,y,z,t],n2o_norm_sf6_grid_mon_avg_b[j,y,z,t]],/nan)
      age_both_grid_mon_avg_b[j,y,z,t] = total([age_co2_grid_mon_avg_b[j,y,z,t],age_grid_mon_avg_b[j,y,z,t]],/nan)
    endfor
  endfor
  
  for z = 0, nz2-1 do for y = 0, nyg-1 do for j = 0, 1 do begin
    n2o_norm_grid_seas[j,y,z,0] = mean([n2o_norm_grid_mon_avg[j,y,z,11],reform(n2o_norm_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas[j,y,z,s] = mean(n2o_norm_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    age_co2_grid_seas[j,y,z,0] = mean([age_co2_grid_mon_avg[j,y,z,11],reform(age_co2_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_co2_grid_seas[j,y,z,s] = mean(age_co2_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    n2o_norm_sf6_grid_seas[j,y,z,0] = mean([n2o_norm_sf6_grid_mon_avg[j,y,z,11],reform(n2o_norm_sf6_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_sf6_grid_seas[j,y,z,s] = mean(n2o_norm_sf6_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    age_sf6_grid_seas[j,y,z,0] = mean([age_sf6_grid_mon_avg[j,y,z,11],reform(age_sf6_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_sf6_grid_seas[j,y,z,s] = mean(age_sf6_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    n2o_norm_combo_grid_seas[j,y,z,0] = mean([n2o_norm_combo_grid_mon_avg[j,y,z,11],reform(n2o_norm_combo_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_combo_grid_seas[j,y,z,s] = mean(n2o_norm_combo_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    age_combo_grid_seas[j,y,z,0] = mean([age_combo_grid_mon_avg[j,y,z,11],reform(age_combo_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_combo_grid_seas[j,y,z,s] = mean(age_combo_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    n2o_norm_grid_seas_a0[j,y,z,0] = mean([n2o_norm_grid_mon_avg_a0[j,y,z,11],reform(n2o_norm_grid_mon_avg_a0[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas_a0[j,y,z,s] = mean(n2o_norm_grid_mon_avg_a0[j,y,z,indgen(3)+3*s-1],/nan)
    n2o_norm_grid_seas_a1[j,y,z,0] = mean([n2o_norm_grid_mon_avg_a1[j,y,z,11],reform(n2o_norm_grid_mon_avg_a1[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas_a1[j,y,z,s] = mean(n2o_norm_grid_mon_avg_a1[j,y,z,indgen(3)+3*s-1],/nan)
  endfor
  
  ;  Fill in missing locations.
  for s = 0, 3 do for z = 1, nz2-2 do for y = 1, nyg-2 do for j = 0, 1 do begin
    if ~finite(n2o_norm_combo_grid_seas[j,y,z,s]) and finite(n2o_norm_combo_grid_seas[j,y-1,z,s]) and finite(n2o_norm_combo_grid_seas[j,y+1,z,s]) then $
      n2o_norm_combo_grid_seas[j,y,z,s] = 0.5 * n2o_norm_combo_grid_seas[j,y-1,z,s] + 0.5 * n2o_norm_combo_grid_seas[j,y+1,z,s]
    if ~finite(n2o_norm_combo_grid_seas[j,y,z,s]) and finite(n2o_norm_combo_grid_seas[j,y,z-1,s]) and finite(n2o_norm_combo_grid_seas[j,y,z+1,s]) then $
      n2o_norm_combo_grid_seas[j,y,z,s] = 0.5 * n2o_norm_combo_grid_seas[j,y,z-1,s] + 0.5 * n2o_norm_combo_grid_seas[j,y,z+1,s]
    if ~finite(age_combo_grid_seas[j,y,z,s]) and finite(age_combo_grid_seas[j,y-1,z,s]) and finite(age_combo_grid_seas[j,y+1,z,s]) then $
      age_combo_grid_seas[j,y,z,s] = 0.5 * age_combo_grid_seas[j,y-1,z,s] + 0.5 * age_combo_grid_seas[j,y+1,z,s]
    if ~finite(age_combo_grid_seas[j,y,z,s]) and finite(age_combo_grid_seas[j,y,z-1,s]) and finite(age_combo_grid_seas[j,y,z+1,s]) then $
      age_combo_grid_seas[j,y,z,s] = 0.5 * age_combo_grid_seas[j,y,z-1,s] + 0.5 * age_combo_grid_seas[j,y,z+1,s]
  endfor
  
  n2o_norm_grid_seas_late = n2o_norm_grid_seas
  age_combo_grid_seas_late = age_combo_grid_seas
  n2o_norm_combo_grid_seas_late = n2o_norm_combo_grid_seas
  
  cc = replicate(!values.f_nan,3,nyg,nz2) & age_combo_grid_late = cc & age_combo_grid_minmax_late = cc & n2o_norm_combo_grid_late = cc & n2o_norm_combo_grid_minmax_late = cc
  for z = 1, nz2-2 do for y = 0, nyg-1 do begin
    gd = where(finite(age_combo_grid_seas_late[0,y,z,*]),ngd)
  ;  err = 1. / sqrt(total(cgrid[1,y,z,gd]^2))
  ;  wmean = total(1. / cgrid[1,y,z,gd]^2 * cgrid[0,y,z,gd]) / total(1. / cgrid[1,y,z,gd]^2)
    if ngd ge 1 then begin
      stats = moment(age_combo_grid_seas_late[0,y,z,gd],sdev=sdev)
      age_combo_grid_late[*,y,z] = [stats[0],sdev,ngd]
      age_combo_grid_minmax_late[0:1,y,z] = [min(age_combo_grid_seas_late[0,y,z,gd]-age_combo_grid_seas_late[1,y,z,gd]),max(age_combo_grid_seas_late[0,y,z,gd]+age_combo_grid_seas_late[1,y,z,gd])]
    endif
    gd = where(finite(n2o_norm_combo_grid_seas_late[0,y,z,*]),ngd)
    if ngd ge 1 then begin
      stats = moment(n2o_norm_combo_grid_seas_late[0,y,z,gd],sdev=sdev)
      n2o_norm_combo_grid_late[*,y,z] = [stats[0],sdev,ngd]
      n2o_norm_combo_grid_minmax_late[0:1,y,z] = [min(n2o_norm_combo_grid_seas_late[0,y,z,gd]-n2o_norm_combo_grid_seas_late[1,y,z,gd]), $
        max(n2o_norm_combo_grid_seas_late[0,y,z,gd]+n2o_norm_combo_grid_seas_late[1,y,z,gd])]
    endif
  endfor
  
  n2o_norm_opt_grid_mon_avg_b_late = n2o_norm_opt_grid_mon_avg
  n2o_norm_grid_mon_avg_b_late = n2o_norm_grid_mon_avg
  n2o_norm_sf6_grid_mon_avg_b_late = n2o_norm_sf6_grid_mon_avg
  n2o_norm_grid_mon_avg_a_late = n2o_norm_grid_mon_avg_a
  n2o_norm_sf6_grid_mon_avg_a_late = n2o_norm_sf6_grid_mon_avg_a
  n2o_norm_grid_mon_avg_a_early = n2o_norm_grid_mon_avg_a0
  n2o_norm_sf6_grid_mon_avg_a_early = n2o_norm_sf6_grid_mon_avg_a0
  n2o_norm_grid_mon_avg_a_mid = n2o_norm_grid_mon_avg_a1
  n2o_norm_sf6_grid_mon_avg_a_mid = n2o_norm_sf6_grid_mon_avg_a1
  n2o_norm_both_grid_mon_avg_b_early = n2o_norm_both_grid_mon_avg_b
  age_co2_grid_mon_avg_b_late = age_co2_grid_mon_avg
  age_sf6_grid_mon_avg_b_late = age_sf6_grid_mon_avg
  age_opt_grid_mon_avg_b_late = age_opt_grid_mon_avg
  age_co2_grid_mon_avg_a_late = age_co2_grid_mon_avg_a
  age_sf6_grid_mon_avg_a_late = age_sf6_grid_mon_avg_a
  age_opt_grid_mon_avg_a_late = age_opt_grid_mon_avg_a
  age_co2_grid_mon_avg_b_early = age_co2_grid_mon_avg_b
  age_grid_mon_avg_b_early = age_grid_mon_avg_b
  age_both_grid_mon_avg_b_early = age_both_grid_mon_avg_b
  age_co2_grid_mon_avg_a_early = age_co2_grid_mon_avg_a0
  age_sf6_grid_mon_avg_a_early = age_sf6_grid_mon_avg_a0
  age_opt_grid_mon_avg_a_early = age_opt_grid_mon_avg_a0
  age_co2_grid_mon_avg_a_mid = age_co2_grid_mon_avg_a1
  age_sf6_grid_mon_avg_a_mid = age_sf6_grid_mon_avg_a1
  age_opt_grid_mon_avg_a_mid = age_opt_grid_mon_avg_a1
  
  elat_m_adj_b_late = elat_m_adj
  elat_m_adj_a_late = elat_m_adj_a
  elat_m_adj_b_early = elat_m_adj_b
  elat_m_adj_a_early = elat_m_adj_a0
  elat_m_adj_a_mid = elat_m_adj_a1
  
  save,n2o_norm_grid_seas_late,n2o_norm_combo_grid_seas_late,age_co2_grid_seas,age_combo_grid_seas_late,age_combo_grid_late,n2o_norm_combo_grid_late,n2o_norm_combo_grid_minmax_late, $
    n2o_norm_combo_grid_mon_avg,n2o_norm_grid_mon_avg_b_late,age_co2_grid_mon_avg_b_late,n2o_norm_sf6_grid_mon_avg_b_late,age_combo_grid_mon_avg,age_sf6_grid_mon_avg_b_late, $
    age_opt_grid_mon_avg_b_late,age_co2_grid_mon_avg_a_late,n2o_norm_grid_mon_avg_a_late,age_sf6_grid_mon_avg_a_late,n2o_norm_opt_grid_mon_avg_b_late,age_opt_grid_mon_avg_a_late, $
    filename=dir+'Balloon/Aircore_grid_elat_adj.sav'
  save,n2o_norm_both_grid_mon_avg_b_early,age_co2_grid_mon_avg_b_early,age_grid_mon_avg_b_early,age_both_grid_mon_avg_b_early,n2o_norm_grid_mon_avg_a_early, $
    n2o_norm_sf6_grid_mon_avg_a_early,age_co2_grid_mon_avg_a_early,age_sf6_grid_mon_avg_a_early,age_opt_grid_mon_avg_a_early,combo_n2o_seas_hires, $
    combo_n2o_seas_early_hires,n2o_norm_grid_mon_avg_a_mid,n2o_norm_sf6_grid_mon_avg_a_mid,age_co2_grid_mon_avg_a_mid,age_sf6_grid_mon_avg_a_mid, $
    age_opt_grid_mon_avg_a_mid,n2o_norm_grid_seas_a1,filename=dir+'Balloon/OMS_grid_elat_adj.sav'
  save,elat_m_adj_a_early,elat_m_adj_b_early,elat_m_adj_a_late,elat_m_adj_b_late,elat_m_adj_a_mid,filename=dir+'Balloon/Elat_adj_all.sav'

endif
if stage eq 'elats' then return

restore,dir+'Balloon/Elat_adj_all.sav'
restore,dir+'Balloon/OMS_grid_elat_adj.sav'
restore,dir+'Balloon/Aircore_grid_elat_adj.sav'
restore,dir+'Balloon/balloon_mean_age_seas_grid_early.sav'

z = 28
;print,reform(age_co2_grid_mon_avg_a_early[0,*,z,3]),reform(age_co2_grid_mon_avg_a_early[3,*,z,3])

yi = where(lat_grid ge 30 and lat_grid le 60)
nsamp_oms = fltarr(nz3)
for z = 3, nz3-1 do nsamp_oms[z] = total(age_both_grid_mon_avg_b_early[2,yi,2*z-1:2*z,*],/nan)
;save,nsamp_oms,filename=dir+'Balloon/NSamp_OMS.sav'


;  Make mean age seasonal 'climatology' from early and late grids to use for lat gradients.
age_grid_seas_avg = age_combo_grid_seas_late
for s = 0, 3 do for z = 0, nz2-1 do for y = 0, nyg-1 do for j = 0, 1 do begin
  age_grid_seas_avg[j,y,z,s] = mean([age_combo_grid_seas_late[j,y,z,s],age_combo_grid_seas_sm_early[j,y,z,s]],/nan)
endfor
for s = 0, 3 do begin
  for z = 0, nz2-1 do begin
    gd = where(finite(age_grid_seas_avg[0,*,z,s]),ngd)
    if ngd ge 7 then age_grid_seas_avg[0,gd,z,s] = smooth(age_grid_seas_avg[0,gd,z,s],5,/edge_truncate)
  endfor
  for y = 0, nyg-1 do begin
    gd = where(finite(age_grid_seas_avg[0,y,*,s]),ngd)
    if ngd ge 7 then age_grid_seas_avg[0,y,gd,s] = smooth(age_grid_seas_avg[0,y,gd,s],5,/edge_truncate)
  endfor
endfor
grid_fill_in_situ_seas,nyg,nz2,age_grid_seas_avg
for s = 0, 3 do begin
  for z = 0, nz2-1 do begin
    gd = where(finite(age_grid_seas_avg[0,*,z,s]),ngd)
    if ngd ge 7 then age_grid_seas_avg[0,gd,z,s] = smooth(age_grid_seas_avg[0,gd,z,s],5,/edge_truncate)
  endfor
  for y = 0, nyg-1 do begin
    gd = where(finite(age_grid_seas_avg[0,y,*,s]),ngd)
    if ngd ge 7 then age_grid_seas_avg[0,y,gd,s] = smooth(age_grid_seas_avg[0,y,gd,s],5,/edge_truncate)
  endfor
endfor

bb = replicate(!values.f_nan,3,nyg,nz,4)

;  New time series in lat/alt grid using adjusted elats.
nt = 2644
dt = 1./48.
years = findgen(nt)*dt+1975.
nnorm = 20
dn = 0.05
norm_grid = findgen(nnorm)*dn+0.025
nn = 20
da1 = 0.5
age_grid1 = findgen(nn)*da1

fff = replicate(!values.f_nan,6,nyg,nz2,nt) & ggg = replicate(!values.f_nan,6,nz2,nt) & hhh = replicate(!values.f_nan,nz2,nt) & jjj = replicate(!values.f_nan,3,nn,nt)
mmm = replicate(!values.f_nan,6,nn,nt)

do_time_series = (stage eq 'time_series')

if do_time_series then begin
  n2o_norm_grid_adj_tseries = fff & age_grid_adj_tseries = fff & age_co2_grid_adj_tseries = fff & age_sf6_grid_adj_tseries = fff
  age_opt_grid_adj_tseries = fff & age_co2_b_grid_adj_tseries = fff & age_sf6_b_grid_adj_tseries = fff & age_opt_b_grid_adj_tseries = fff
  n2o_a_grid_adj_tseries = fff & n2o_b_grid_adj_tseries = fff & n2o_grid_adj_tseries = fff & n2o_on_age_adj_tseries = jjj & n2o_on_age_b_adj_tseries = jjj
  n2o_on_age_both_adj_tseries = jjj & age_on_n2o_adj_tseries = mmm & age_sf6_on_n2o_b_adj_tseries = mmm & age_sf6_on_n2o_adj_tseries = mmm
  age_co2_on_n2o_b_adj_tseries = mmm & age_co2_on_n2o_adj_tseries = mmm & ch4_a_grid_adj_tseries = fff & ch4_b_grid_adj_tseries = fff & ch4_grid_adj_tseries = fff
  age_sf6_on_ch4_b_adj_tseries = mmm & age_sf6_on_ch4_adj_tseries = mmm & age_co2_on_ch4_b_adj_tseries = mmm & age_co2_on_ch4_adj_tseries = mmm
  age_on_ch4_adj_tseries = mmm & age_on_n2o_e_adj_tseries = mmm & age_co2_on_n2o_e_adj_tseries = mmm & age_sf6_on_n2o_e_adj_tseries = mmm
  age_co2_on_n2o_e_b_adj_tseries = mmm & age_sf6_on_n2o_e_b_adj_tseries = mmm & age_co2_on_ch4_e_adj_tseries = mmm & age_sf6_on_ch4_e_adj_tseries = mmm
  age_co2_on_ch4_e_b_adj_tseries = mmm & age_sf6_on_ch4_e_b_adj_tseries = mmm & age_on_ch4_e_adj_tseries = mmm

  dy_int = 2*dy/3.
  dz_int = 2*dz2/3.

  for t = 0, nt-1 do begin
    
    ;  Early aircraft.   -------------------------------------------------------------------------------------------------------------------
    chk = where(times_all_a0[in2o_a0] ge years[t]-dt/2. and times_all_a0[in2o_a0] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a0[in2o_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[in2o_a0[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_early[in2o_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[in2o_a0[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_norm_a0[in2o_a0[chk[gd]]],sdev=sdev,/nan)
          n2o_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
  ;        if t eq 802 and z eq 24 and y eq 29 then print,years[t],alt_grid2[z],lat_grid[y],ngd,n2o_a_grid_adj_tseries[0,y,z,t],n2o_norm_a0[in2o_a0[chk[gd]]]
        endif
      endfor
    endif
    chk = where(times_all_a0[ich4_a0] ge years[t]-dt/2. and times_all_a0[ich4_a0] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a0[ich4_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[ich4_a0[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_early[ich4_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[ich4_a0[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_a0[ich4_a0[chk]] gt 0 and $
          ch4_all_norm_a0[ich4_a0[chk]] lt 1.1,ngd)
        if ngd ge 1 then begin
          stats = moment(ch4_all_norm_a0[ich4_a0[chk[gd]]],sdev=sdev,/nan)
          ch4_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a0[ico2_a0] ge years[t]-dt/2. and times_all_a0[ico2_a0] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a0[ico2_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[ico2_a0[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_early[ico2_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[ico2_a0[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a0.age_range_all[1,chk[gd]] - dat_co2_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_a0[ico2_a0[chk]] ge norm_grid[n]-dn/2. and n2o_norm_a0[ico2_a0[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a0[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a0.age_range_all[1,chk[gd]] - dat_co2_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_a0.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a0.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a0[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a0.age_range_all[1,chk[gd]] - dat_co2_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_a0[ico2_a0[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a0[ico2_a0[chk]] lt norm_grid[n]+dn/2. and $
          elat_m_adj_a_early[ico2_a0[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a0.age_range_all[1,chk[gd]] - dat_co2_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_a0.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a0.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and $
          elat_m_adj_a_early[ico2_a0[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a0.age_range_all[1,chk[gd]] - dat_co2_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_co2_a0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_a0.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
           elat_m_adj_a_early[ico2_a0[chk]] ge 30 and elat_m_adj_a_early[ico2_a0[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_a0[ico2_a0[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_adj_tseries[*,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a0[isf6_a0] ge years[t]-dt/2. and times_all_a0[isf6_a0] lt years[t]+dt/2. and times_all_a0[isf6_a0] lt 1999,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a0[isf6_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[isf6_a0[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_early[isf6_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[isf6_a0[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a0.age_range_all[1,chk[gd]] - dat_sf6_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_a0[isf6_a0[chk]] ge norm_grid[n]-dn/2. and n2o_norm_a0[isf6_a0[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a0[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] lt 1999 then begin
          age_err = 0.5 * (dat_sf6_a0.age_range_all[1,chk[gd]] - dat_sf6_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_a0.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a0.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a0[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] lt 1999 then begin
          age_err = 0.5 * (dat_sf6_a0.age_range_all[1,chk[gd]] - dat_sf6_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_a0[isf6_a0[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a0[isf6_a0[chk]] lt norm_grid[n]+dn/2. and $
          elat_m_adj_a_early[isf6_a0[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] lt 1999 then begin
          age_err = 0.5 * (dat_sf6_a0.age_range_all[1,chk[gd]] - dat_sf6_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_a0.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a0.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and $
          elat_m_adj_a_early[isf6_a0[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] lt 1999 then begin
          age_err = 0.5 * (dat_sf6_a0.age_range_all[1,chk[gd]] - dat_sf6_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_sf6_a0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_a0.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_a_early[isf6_a0[chk]] ge 30 and elat_m_adj_a_early[isf6_a0[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_a0[isf6_a0[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_adj_tseries[*,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a0[iage_a0] ge years[t]-dt/2. and times_all_a0[iage_a0] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a0[iage_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[iage_a0[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_early[iage_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[iage_a0[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_a0.age_range_all[1,chk[gd]] - dat_a0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_a0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_opt_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
    endif

    ;  Early balloon.   -------------------------------------------------------------------------------------------------------------------
    chk = where(times_all_b[in2o_b] ge years[t]-dt/2. and times_all_b[in2o_b] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 2
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_b[in2o_b[chk]] ge alt_grid2[z] - dz_int and alt_b[in2o_b[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_early[in2o_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[in2o_b[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_norm_b[in2o_b[chk[gd]]],sdev=sdev,/nan)
          n2o_b_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_b[ich4_b] ge years[t]-dt/2. and times_all_b[ich4_b] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 2
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_b[ich4_b[chk]] ge alt_grid2[z] - dz_int and alt_b[ich4_b[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_early[ich4_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[ich4_b[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_b[ich4_b[chk]] gt 0 and $
          ch4_all_norm_b[ich4_b[chk]] lt 1.1,ngd)
        if ngd ge 1 then begin
          stats = moment(ch4_all_norm_b[ich4_b[chk[gd]]],sdev=sdev,/nan)
          ch4_b_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_b[ico2_b] ge years[t]-dt/2. and times_all_b[ico2_b] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 2
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_b[ico2_b[chk]] ge alt_grid2[z] - dz_int and alt_b[ico2_b[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_early[ico2_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[ico2_b[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b0.age_range_all[1,chk[gd]] - dat_co2_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_b_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_b[ico2_b[chk]] ge norm_grid[n]-dn/2. and n2o_norm_b[ico2_b[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[ico2_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b0.age_range_all[1,chk[gd]] - dat_co2_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_b0.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_b0.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[ico2_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b0.age_range_all[1,chk[gd]] - dat_co2_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_b[ico2_b[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_b[ico2_b[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[ico2_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b0.age_range_all[1,chk[gd]] - dat_co2_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_b0.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_b0.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[ico2_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b0.age_range_all[1,chk[gd]] - dat_co2_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_co2_b0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_b0.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_b_early[ico2_b[chk]] ge 30 and elat_m_adj_b_early[ico2_b[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_b[ico2_b[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_b_adj_tseries[0:2,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_b[isf6_b] ge years[t]-dt/2. and times_all_b[isf6_b] lt years[t]+dt/2. and times_all_b[isf6_b] lt 2002,ni_all)
    if ni_all ge 1 then begin
      samp_type = 2
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_b[isf6_b[chk]] ge alt_grid2[z] - dz_int and alt_b[isf6_b[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_early[isf6_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[isf6_b[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b0.age_range_all[1,chk[gd]] - dat_sf6_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_b_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_b[isf6_b[chk]] ge norm_grid[n]-dn/2. and n2o_norm_b[isf6_b[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b0.age_range_all[1,chk[gd]] - dat_sf6_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_b0.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_b0.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b0.age_range_all[1,chk[gd]] - dat_sf6_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_b[isf6_b[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_b[isf6_b[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b0.age_range_all[1,chk[gd]] - dat_sf6_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_b0.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_b0.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6_b[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b0.age_range_all[1,chk[gd]] - dat_sf6_b0.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_sf6_b0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_b0.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_b_early[isf6_b[chk]] ge 30 and elat_m_adj_b_early[isf6_b[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_b[isf6_b[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_b_adj_tseries[0:2,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_b[iage_b] ge years[t]-dt/2. and times_all_b[iage_b] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 2
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_b[iage_b[chk]] ge alt_grid2[z] - dz_int and alt_b[iage_b[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_early[iage_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[iage_b[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
;          age_err = 0.5 * (dat_b0.age_range_all[1,chk[gd]] - dat_b0.age_range_all[0,chk[gd]])
          age_err = replicate(0.5,ngd)
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_b0.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_b0.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_opt_b_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
    endif

    ;  Mid aircraft.   -------------------------------------------------------------------------------------------------------------------
    chk = where(times_all_a1[in2o_a1] ge years[t]-dt/2. and times_all_a1[in2o_a1] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a1[in2o_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[in2o_a1[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_mid[in2o_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[in2o_a1[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_norm_a1[in2o_a1[chk[gd]]],sdev=sdev,/nan)
          n2o_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a1[ich4_a1] ge years[t]-dt/2. and times_all_a1[ich4_a1] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a1[ich4_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[ich4_a1[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_mid[ich4_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[ich4_a1[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_a1[ich4_a1[chk]] gt 0 and $
          ch4_all_norm_a1[ich4_a1[chk]] lt 1.1,ngd)
        if ngd ge 1 then begin
          stats = moment(ch4_all_norm_a1[ich4_a1[chk[gd]]],sdev=sdev,/nan)
          ch4_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a1[ico2_a1] ge years[t]-dt/2. and times_all_a1[ico2_a1] lt years[t]+dt/2. and times_all_a1[ico2_a1] ge 2014,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a1[ico2_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[ico2_a1[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_mid[ico2_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[ico2_a1[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a1.age_range_all[1,chk[gd]] - dat_co2_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_a1[ico2_a1[chk]] ge norm_grid[n]-dn/2. and n2o_norm_a1[ico2_a1[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a1[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a1.age_range_all[1,chk[gd]] - dat_co2_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_a1.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a1.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a1[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a1.age_range_all[1,chk[gd]] - dat_co2_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_a1[ico2_a1[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a1[ico2_a1[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a1[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a1.age_range_all[1,chk[gd]] - dat_co2_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_a1.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a1.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a1[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a1.age_range_all[1,chk[gd]] - dat_co2_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_co2_a1.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_a1.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_a_mid[ico2_a1[chk]] ge 30 and elat_m_adj_a_mid[ico2_a1[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_a1[ico2_a1[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_adj_tseries[*,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a1[isf6_a1] ge years[t]-dt/2. and times_all_a1[isf6_a1] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a1[isf6_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[isf6_a1[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_mid[isf6_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[isf6_a1[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a1.age_range_all[1,chk[gd]] - dat_sf6_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_a1[isf6_a1[chk]] ge norm_grid[n]-dn/2. and n2o_norm_a1[isf6_a1[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a1[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] ge 2011 then begin
          age_err = 0.5 * (dat_sf6_a1.age_range_all[1,chk[gd]] - dat_sf6_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_a1.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a1.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a1[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] ge 2011 then begin
          age_err = 0.5 * (dat_sf6_a1.age_range_all[1,chk[gd]] - dat_sf6_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_a1[isf6_a1[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a1[isf6_a1[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a1[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] ge 2011 then begin
          age_err = 0.5 * (dat_sf6_a1.age_range_all[1,chk[gd]] - dat_sf6_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_a1.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a1.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a1[chk]] ge 30,ngd)
        if ngd ge 1 and years[t] ge 2011 then begin
          age_err = 0.5 * (dat_sf6_a1.age_range_all[1,chk[gd]] - dat_sf6_a1.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_sf6_a1.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_a1.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_a_mid[isf6_a1[chk]] ge 30 and elat_m_adj_a_mid[isf6_a1[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_a1[isf6_a1[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_adj_tseries[*,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a1[iage_a1] ge years[t]-dt/2. and times_all_a1[iage_a1] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a1[iage_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[iage_a1[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_mid[iage_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[iage_a1[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
;          age_err = 0.5 * (dat_a1.age_range_all[1,chk[gd]] - dat_a1.age_range_all[0,chk[gd]])
          age_err = replicate(0.5,ngd)
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_a1.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_a1.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_opt_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
    endif

    ;  Late aircraft.   -------------------------------------------------------------------------------------------------------------------
    chk = where(times_all_a[in2o_a] ge years[t]-dt/2. and times_all_a[in2o_a] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a[in2o_a[chk]] ge alt_grid2[z] - dz_int and alt_a[in2o_a[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_late[in2o_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[in2o_a[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_norm_a[in2o_a[chk[gd]]],sdev=sdev,/nan)
          n2o_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a[ich4_a] ge years[t]-dt/2. and times_all_a[ich4_a] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a[ich4_a[chk]] ge alt_grid2[z] - dz_int and alt_a[ich4_a[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_late[ich4_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[ich4_a[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_a[ich4_a[chk]] gt 0 and $
          ch4_all_norm_a[ich4_a[chk]] lt 1.1,ngd)
        if ngd ge 1 then begin
          stats = moment(ch4_all_norm_a[ich4_a[chk[gd]]],sdev=sdev,/nan)
          ch4_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a[ico2_a] ge years[t]-dt/2. and times_all_a[ico2_a] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a[ico2_a[chk]] ge alt_grid2[z] - dz_int and alt_a[ico2_a[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_late[ico2_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[ico2_a[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a2.age_range_all[1,chk[gd]] - dat_co2_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_a[ico2_a[chk]] ge norm_grid[n]-dn/2. and n2o_norm_a[ico2_a[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a2.age_range_all[1,chk[gd]] - dat_co2_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_a2.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a2.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a2.age_range_all[1,chk[gd]] - dat_co2_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_a[ico2_a[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a[ico2_a[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a2.age_range_all[1,chk[gd]] - dat_co2_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_a2.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a2.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[ico2_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_a2.age_range_all[1,chk[gd]] - dat_co2_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_co2_a2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_a2.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_a_late[ico2_a[chk]] ge 30 and elat_m_adj_a_late[ico2_a[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_a[ico2_a[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_adj_tseries[*,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a[isf6_a] ge years[t]-dt/2. and times_all_a[isf6_a] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a[isf6_a[chk]] ge alt_grid2[z] - dz_int and alt_a[isf6_a[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_late[isf6_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[isf6_a[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a2.age_range_all[1,chk[gd]] - dat_sf6_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm_a[isf6_a[chk]] ge norm_grid[n]-dn/2. and n2o_norm_a[isf6_a[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a2.age_range_all[1,chk[gd]] - dat_sf6_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_a2.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a2.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a2.age_range_all[1,chk[gd]] - dat_sf6_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm_a[isf6_a[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a[isf6_a[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a2.age_range_all[1,chk[gd]] - dat_sf6_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_a2.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a2.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_a_early[isf6_a[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_a2.age_range_all[1,chk[gd]] - dat_sf6_a2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_e_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_sf6_a2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_a2.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_a_late[isf6_a[chk]] ge 30 and elat_m_adj_a_late[isf6_a[chk]] le 60,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_a[isf6_a[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_adj_tseries[*,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all_a[iage_a] ge years[t]-dt/2. and times_all_a[iage_a] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 1
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt_a[iage_a[chk]] ge alt_grid2[z] - dz_int and alt_a[iage_a[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_a_late[iage_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[iage_a[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
;          age_err = 0.5 * (dat_a2.age_range_all[1,chk[gd]] - dat_a2.age_range_all[0,chk[gd]])
          age_err = replicate(1,ngd)
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_a2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_a2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_opt_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
    endif
    
    ;  Late balloon.   -------------------------------------------------------------------------------------------------------------------
    chk = where(times_all[in2o] ge years[t]-dt/2. and times_all[in2o] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 3
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt[in2o[chk]] ge alt_grid2[z] - dz_int and alt[in2o[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_late[in2o[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[in2o[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o_norm[in2o[chk[gd]]],sdev=sdev,/nan)
          n2o_b_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all[ich4] ge years[t]-dt/2. and times_all[ich4] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 3
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt[ich4[chk]] ge alt_grid2[z] - dz_int and alt[ich4[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_late[ich4[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[ich4[chk]] lt lat_grid[y]+dy_int and ch4_all_norm[ich4[chk]] gt 0 and $
          ch4_all_norm[ich4[chk]] lt 1.1,ngd)
        if ngd ge 1 then begin
          stats = moment(ch4_all_norm[ich4[chk[gd]]],sdev=sdev,/nan)
          ch4_b_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all[ico2] ge years[t]-dt/2. and times_all[ico2] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 3
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt[ico2[chk]] ge alt_grid2[z] - dz_int and alt[ico2[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_late[ico2[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[ico2[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b2.age_range_all[1,chk[gd]] - dat_co2_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_b_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm[ico2[chk]] ge norm_grid[n]-dn/2. and n2o_norm[ico2[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_late[ico2[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b2.age_range_all[1,chk[gd]] - dat_co2_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_b2.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_b2.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_late[ico2[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b2.age_range_all[1,chk[gd]] - dat_co2_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_n2o_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm[ico2[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm[ico2[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_late[ico2[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b2.age_range_all[1,chk[gd]] - dat_co2_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_co2_b2.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_b2.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_late[ico2[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_co2_b2.age_range_all[1,chk[gd]] - dat_co2_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_co2_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_co2_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_co2_on_ch4_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_co2_b2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_b2.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_b_late[ico2[chk]] ge 30 and elat_m_adj_b_late[ico2[chk]] le 50,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o[ico2[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_b_adj_tseries[0:2,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    chk = where(times_all[isf6] ge years[t]-dt/2. and times_all[isf6] lt years[t]+dt/2.,ni_all)
    if ni_all ge 1 then begin
      samp_type = 3
      for z = 0, nz2-1 do for y = 0, nyg-1 do begin
        gd = where(alt[isf6[chk]] ge alt_grid2[z] -dz_int and alt[isf6[chk]] lt alt_grid2[z] + dz_int and $
          elat_m_adj_b_late[isf6[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[isf6[chk]] lt lat_grid[y]+dy_int,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b2.age_range_all[1,chk[gd]] - dat_sf6_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_b_grid_adj_tseries[*,y,z,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for n = 0, nn-1 do begin
        gd = where(n2o_norm[isf6[chk]] ge norm_grid[n]-dn/2. and n2o_norm[isf6[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b2.age_range_all[1,chk[gd]] - dat_sf6_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_b2.n2o_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_b2.n2o_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b2.age_range_all[1,chk[gd]] - dat_sf6_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_n2o_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(ch4_all_norm[isf6[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm[isf6[chk]] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b2.age_range_all[1,chk[gd]] - dat_sf6_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
        gd = where(bins_sf6_b2.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_b2.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and elat_m_adj_b_early[isf6[chk]] ge 30,ngd)
        if ngd ge 1 then begin
          age_err = 0.5 * (dat_sf6_b2.age_range_all[1,chk[gd]] - dat_sf6_b2.age_range_all[0,chk[gd]])
          age_weights = 1. / age_err^2
          err = 1. / sqrt(total(age_weights))
          wmean = total(age_weights * dat_sf6_b2.age_opt_all[chk[gd]]) / total(age_weights)
          stats = moment(dat_sf6_b2.age_opt_all[chk[gd]],sdev=sdev,/nan)
          age_sf6_on_ch4_e_b_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
        endif
      endfor
      for a = 0, nn-1 do begin
        gd = where(dat_sf6_b2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_b2.age_opt_all[chk] lt age_grid1[a]+da1/2. and $
          elat_m_adj_b_late[isf6[chk]] ge 30 and elat_m_adj_b_late[isf6[chk]] le 50,ngd)
        if ngd ge 1 then begin
          stats = moment(n2o[isf6[chk[gd]]],sdev=sdev,/nan)
          n2o_on_age_b_adj_tseries[0:2,a,t] = [stats[0],sdev,samp_type]
        endif
      endfor
    endif
    
    for z = 0, nz2-1 do for y = 0, nyg-1 do begin
      for j = 0, 2 do n2o_grid_adj_tseries[j,y,z,t] = mean([n2o_a_grid_adj_tseries[j,y,z,t],n2o_b_grid_adj_tseries[j,y,z,t]],/nan)
      for j = 0, 2 do ch4_grid_adj_tseries[j,y,z,t] = mean([ch4_a_grid_adj_tseries[j,y,z,t],ch4_b_grid_adj_tseries[j,y,z,t]],/nan)
      age1 = [age_co2_grid_adj_tseries[0,y,z,t],age_sf6_grid_adj_tseries[0,y,z,t],age_co2_b_grid_adj_tseries[0,y,z,t],age_sf6_b_grid_adj_tseries[0,y,z,t]]
      err1 = [age_co2_grid_adj_tseries[1,y,z,t],age_sf6_grid_adj_tseries[1,y,z,t],age_co2_b_grid_adj_tseries[1,y,z,t],age_sf6_b_grid_adj_tseries[1,y,z,t]]
      samps = [age_co2_grid_adj_tseries[5,y,z,t],age_sf6_grid_adj_tseries[5,y,z,t],age_co2_b_grid_adj_tseries[5,y,z,t],age_sf6_b_grid_adj_tseries[5,y,z,t]]
      gd = where(finite(age1),ngd)
      if ngd gt 0 then begin
        age_weights = 1. / err1[gd]^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * age1[gd]) / total(age_weights)
        age_grid_adj_tseries[0:1,y,z,t] = [wmean,err]
        age_grid_adj_tseries[5,y,z,t] = samps[gd[0]]
      endif
      for j = 2, 3 do age_grid_adj_tseries[j,y,z,t] = mean([age_co2_grid_adj_tseries[j,y,z,t],age_sf6_grid_adj_tseries[j,y,z,t], $
        age_co2_b_grid_adj_tseries[j,y,z,t],age_sf6_b_grid_adj_tseries[j,y,z,t]],/nan)
      age_grid_adj_tseries[4,y,z,t] = total([age_co2_grid_adj_tseries[4,y,z,t],age_sf6_grid_adj_tseries[4,y,z,t],age_co2_b_grid_adj_tseries[4,y,z,t], $
        age_sf6_b_grid_adj_tseries[4,y,z,t]],/nan)      
    endfor
    for z = 0, nn-1 do begin
      age1 = [age_co2_on_n2o_adj_tseries[0,z,t],age_sf6_on_n2o_adj_tseries[0,z,t],age_co2_on_n2o_b_adj_tseries[0,z,t],age_sf6_on_n2o_b_adj_tseries[0,z,t]]
      err1 = [age_co2_on_n2o_adj_tseries[1,z,t],age_sf6_on_n2o_adj_tseries[1,z,t],age_co2_on_n2o_b_adj_tseries[1,z,t],age_sf6_on_n2o_b_adj_tseries[1,z,t]]
      samps = [age_co2_on_n2o_adj_tseries[5,z,t],age_sf6_on_n2o_adj_tseries[5,z,t],age_co2_on_n2o_b_adj_tseries[5,z,t],age_sf6_on_n2o_b_adj_tseries[5,z,t]]
      gd = where(finite(age1),ngd)
      if ngd gt 0 then begin
        age_weights = 1. / err1[gd]^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * age1[gd]) / total(age_weights)
        age_on_n2o_adj_tseries[0:1,z,t] = [wmean,err]
        age_on_n2o_adj_tseries[5,z,t] = samps[gd[0]]
      endif
      for j = 2, 3 do age_on_n2o_adj_tseries[j,z,t] = mean([age_co2_on_n2o_adj_tseries[j,z,t],age_sf6_on_n2o_adj_tseries[j,z,t], $
        age_co2_on_n2o_b_adj_tseries[j,z,t],age_sf6_on_n2o_b_adj_tseries[j,z,t]],/nan)
      age_on_n2o_adj_tseries[4,z,t] = total([age_co2_on_n2o_adj_tseries[4,z,t],age_sf6_on_n2o_adj_tseries[4,z,t],age_co2_on_n2o_b_adj_tseries[4,z,t], $
        age_sf6_on_n2o_b_adj_tseries[4,z,t]],/nan)

      age1 = [age_co2_on_n2o_e_adj_tseries[0,z,t],age_sf6_on_n2o_e_adj_tseries[0,z,t],age_co2_on_n2o_e_b_adj_tseries[0,z,t],age_sf6_on_n2o_e_b_adj_tseries[0,z,t]]
      err1 = [age_co2_on_n2o_e_adj_tseries[1,z,t],age_sf6_on_n2o_e_adj_tseries[1,z,t],age_co2_on_n2o_e_b_adj_tseries[1,z,t],age_sf6_on_n2o_e_b_adj_tseries[1,z,t]]
      samps = [age_co2_on_n2o_e_adj_tseries[5,z,t],age_sf6_on_n2o_e_adj_tseries[5,z,t],age_co2_on_n2o_e_b_adj_tseries[5,z,t],age_sf6_on_n2o_e_b_adj_tseries[5,z,t]]
      gd = where(finite(age1),ngd)
      if ngd gt 0 then begin
        age_weights = 1. / err1[gd]^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * age1[gd]) / total(age_weights)
        age_on_n2o_e_adj_tseries[0:1,z,t] = [wmean,err]
        age_on_n2o_e_adj_tseries[5,z,t] = samps[gd[0]]
      endif
      for j = 2, 3 do age_on_n2o_e_adj_tseries[j,z,t] = mean([age_co2_on_n2o_e_adj_tseries[j,z,t],age_sf6_on_n2o_e_adj_tseries[j,z,t], $
        age_co2_on_n2o_e_b_adj_tseries[j,z,t],age_sf6_on_n2o_e_b_adj_tseries[j,z,t]],/nan)
      age_on_n2o_e_adj_tseries[4,z,t] = total([age_co2_on_n2o_e_adj_tseries[4,z,t],age_sf6_on_n2o_e_adj_tseries[4,z,t],age_co2_on_n2o_e_b_adj_tseries[4,z,t], $
        age_sf6_on_n2o_e_b_adj_tseries[4,z,t]],/nan)

      age1 = [age_co2_on_ch4_adj_tseries[0,z,t],age_sf6_on_ch4_adj_tseries[0,z,t],age_co2_on_ch4_b_adj_tseries[0,z,t],age_sf6_on_ch4_b_adj_tseries[0,z,t]]
      err1 = [age_co2_on_ch4_adj_tseries[1,z,t],age_sf6_on_ch4_adj_tseries[1,z,t],age_co2_on_ch4_b_adj_tseries[1,z,t],age_sf6_on_ch4_b_adj_tseries[1,z,t]]
      samps = [age_co2_on_ch4_adj_tseries[5,z,t],age_sf6_on_ch4_adj_tseries[5,z,t],age_co2_on_ch4_b_adj_tseries[5,z,t],age_sf6_on_ch4_b_adj_tseries[5,z,t]]
      gd = where(finite(age1),ngd)
      if ngd gt 0 then begin
        age_weights = 1. / err1[gd]^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * age1[gd]) / total(age_weights)
        age_on_ch4_adj_tseries[0:1,z,t] = [wmean,err]
        age_on_ch4_adj_tseries[5,z,t] = samps[gd[0]]
      endif
      for j = 2, 3 do age_on_ch4_adj_tseries[j,z,t] = mean([age_co2_on_ch4_adj_tseries[j,z,t],age_sf6_on_ch4_adj_tseries[j,z,t], $
        age_co2_on_ch4_b_adj_tseries[j,z,t],age_sf6_on_ch4_b_adj_tseries[j,z,t]],/nan)
      age_on_ch4_adj_tseries[4,z,t] = total([age_co2_on_ch4_adj_tseries[4,z,t],age_sf6_on_ch4_adj_tseries[4,z,t],age_co2_on_ch4_b_adj_tseries[4,z,t], $
        age_sf6_on_ch4_b_adj_tseries[4,z,t]],/nan)

      age1 = [age_co2_on_ch4_e_adj_tseries[0,z,t],age_sf6_on_ch4_e_adj_tseries[0,z,t],age_co2_on_ch4_e_b_adj_tseries[0,z,t],age_sf6_on_ch4_e_b_adj_tseries[0,z,t]]
      err1 = [age_co2_on_ch4_e_adj_tseries[1,z,t],age_sf6_on_ch4_e_adj_tseries[1,z,t],age_co2_on_ch4_e_b_adj_tseries[1,z,t],age_sf6_on_ch4_e_b_adj_tseries[1,z,t]]
      samps = [age_co2_on_ch4_e_adj_tseries[5,z,t],age_sf6_on_ch4_e_adj_tseries[5,z,t],age_co2_on_ch4_e_b_adj_tseries[5,z,t],age_sf6_on_ch4_e_b_adj_tseries[5,z,t]]
      gd = where(finite(age1),ngd)
      if ngd gt 0 then begin
        age_weights = 1. / err1[gd]^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * age1[gd]) / total(age_weights)
        age_on_ch4_e_adj_tseries[0:1,z,t] = [wmean,err]
        age_on_ch4_e_adj_tseries[5,z,t] = samps[gd[0]]
      endif
      for j = 2, 3 do age_on_ch4_e_adj_tseries[j,z,t] = mean([age_co2_on_ch4_e_adj_tseries[j,z,t],age_sf6_on_ch4_e_adj_tseries[j,z,t], $
        age_co2_on_ch4_e_b_adj_tseries[j,z,t],age_sf6_on_ch4_e_b_adj_tseries[j,z,t]],/nan)
      age_on_ch4_e_adj_tseries[4,z,t] = total([age_co2_on_ch4_e_adj_tseries[4,z,t],age_sf6_on_ch4_e_adj_tseries[4,z,t],age_co2_on_ch4_e_b_adj_tseries[4,z,t], $
        age_sf6_on_ch4_e_b_adj_tseries[4,z,t]],/nan)

    endfor
    for a = 0, nn-1 do for j = 0, 2 do n2o_on_age_both_adj_tseries[j,a,t] = mean([n2o_on_age_adj_tseries[j,a,t],n2o_on_age_b_adj_tseries[j,a,t]],/nan)
    
  endfor
  
  save,n2o_a_grid_adj_tseries,age_co2_grid_adj_tseries,age_sf6_grid_adj_tseries,age_opt_grid_adj_tseries,n2o_b_grid_adj_tseries,age_co2_b_grid_adj_tseries, $
    age_sf6_b_grid_adj_tseries,age_opt_b_grid_adj_tseries,age_grid_adj_tseries,n2o_grid_adj_tseries,n2o_on_age_adj_tseries,n2o_on_age_b_adj_tseries,n2o_on_age_both_adj_tseries, $
    age_on_n2o_adj_tseries,ch4_grid_adj_tseries,age_on_ch4_adj_tseries,age_on_n2o_e_adj_tseries,age_on_ch4_e_adj_tseries, $
    filename=dir+'Balloon/Age_N2O_gridded_N2O_elat_adj_time_series.sav'
  
endif

restore,dir+'Balloon/Age_N2O_gridded_N2O_elat_adj_time_series.sav'


;  Initial trends
ww = replicate(!values.f_nan,nyg,nz2,2) & xx = replicate(!values.f_nan,nz2,2) & xxx = replicate(!values.f_nan,nz3,2) & age_grid_trends = ww & age_grid_trends_sigma = ww
n2o_grid_trends = ww & n2o_grid_trends_sigma = ww & age_nh_trends = xx & age_nh_trends_sigma = xx & age_nh_lat_adj_trends = xx & age_nh_lat_adj_trends_sigma = xx
age_nh_lat_adj_trends_coarse = xxx & age_nh_lat_adj_trends_coarse_sigma = xxx & n2o_nh_trends = xx & n2o_nh_trends_sigma = xx & n2o_nh_lat_adj_trends = xx
n2o_nh_lat_adj_trends_sigma = xx & n2o_nh_lat_adj_trends_coarse = xxx & n2o_nh_lat_adj_trends_coarse_sigma = xxx & age_nh_trends_coarse = xxx
age_nh_trends_coarse_sigma = xxx & ch4_grid_trends = ww & ch4_grid_trends_sigma = ww & ch4_nh_trends = xx & ch4_nh_trends_sigma = xx
ch4_nh_trends_coarse = xxx & ch4_nh_trends_coarse_sigma = xxx & age_nh_lat_adj_trends_per_coarse = xxx & age_nh_trends_per_coarse = xxx
age_nh_trends_per_coarse_sigma = xxx & age_nh_lat_adj_trends_per_coarse_sigma = xxx & age_nh_trends_ac_coarse = xxx & age_nh_lat_adj_trends_ac_coarse_sigma = xxx
age_nh_trends_ac_coarse_sigma = xxx & age_nh_lat_adj_trends_ac_coarse = xxx &  & n2o_nh_lat_adj_trends_per_coarse = xxx & n2o_nh_lat_adj_trends_per_coarse_sigma = xxx
for z = 0, nz2-1 do for y = 0, nyg-1 do begin
  gd = where(finite(age_grid_adj_tseries[0,y,z,*]),ngd)
  if ngd ge 10 then begin
    chk1 = where(years[gd] lt 2001,nchk1)
    chk = where(years[gd] ge 2010,nchk)
    if years[gd[-1]]-years[gd[0]] ge 20 and nchk1 ge 2 and nchk ge 2 then begin
      age_grid_trends[y,z,*] = linfit(years[gd],age_grid_adj_tseries[0,y,z,gd],sigma=sigma)
      age_grid_trends_sigma[y,z,*] = sigma
    endif
  endif
  gd = where(finite(n2o_grid_adj_tseries[0,y,z,*]),ngd)  
  if ngd ge 10 then begin
    chk1 = where(years[gd] lt 2001,nchk1)
    chk = where(years[gd] ge 2010,nchk)
    if years[gd[-1]]-years[gd[0]] ge 20 and nchk1 ge 2 and nchk ge 2 then begin
      n2o_grid_trends[y,z,*] = linfit(years[gd],n2o_grid_adj_tseries[0,y,z,gd],sigma=sigma)
      n2o_grid_trends_sigma[y,z,*] = sigma
    endif
  endif
  gd = where(finite(ch4_grid_adj_tseries[0,y,z,*]) and years ge 1995,ngd)
  if ngd ge 10 then begin
    chk1 = where(years[gd] lt 2001,nchk1)
    chk = where(years[gd] ge 2010,nchk)
    if years[gd[-1]]-years[gd[0]] ge 20 and nchk1 ge 2 and nchk ge 2 then begin
      ch4_grid_trends[y,z,*] = linfit(years[gd],ch4_grid_adj_tseries[0,y,z,gd],sigma=sigma)
      ch4_grid_trends_sigma[y,z,*] = sigma
    endif
  endif
endfor


;  Seasonal cycles.
seas_frac = findgen(12)/12.+1./24.
mm = replicate(!values.f_nan,nyg,nz2,12) & mean_age_grid_seas = mm
seas = years-fix(years)
ti = interpol(findgen(12),seas_frac,seas)
for z = 0, nz2-1 do for y = 0, nyg-1 do begin
  seas_anom = age_grid_adj_tseries[0,y,z,*] - (age_grid_trends[y,z,0] + age_grid_trends[y,z,1]*years)
  for t = 0, 11 do begin
    chk = where(seas ge seas_frac[t]-1./24. and seas le seas_frac[t]+1./24.,nchk)
    if nchk gt 0 then begin
      mean_age_grid_seas[y,z,t] = mean(seas_anom[chk],/nan)
    endif
  endfor
  chk = where(finite(mean_age_grid_seas[y,z,*]),nchk)
;  if nchk eq 12 then mean_age_grid_seas[y,z,*] = lowpass_cfc(mean_age_grid_seas[y,z,*], BOX=3, EDGE_PFCAST=1)
endfor


yi = interpol(findgen(ny_tp),lat_tp,lat_grid)
tp_alt_grid = interpolate(tpause_alt,yi)


;  Make NH midlat avg time series.
lat_avg_n2o_nh_adj_tseries = hhh & n2o_nh_adj_tseries = ggg & nh_adj_tseries_mon = hhh & n2o_nh_adj_tseries_lat_adj = ggg & lat_avg_nh_adj_tseries = hhh
mmm = replicate(!values.f_nan,3,nn,nt) & age_nh_adj_tseries = ggg & nh_adj_tseries_seas = hhh & age_nh_adj_tseries_lat_adj = ggg & age_on_n2o_nh_adj_tseries = mmm
ooo = replicate(!values.f_nan,nn,nt) & ppp = replicate(!values.f_nan,nn,2) & age_on_n2o_nh_adj_tseries_seas = ooo & age_on_n2o_nh_lat_adj_trends = ppp
age_on_n2o_nh_lat_adj_trends_sigma = ppp & ch4_nh_adj_tseries = ggg & age_nhe_adj_tseries = ggg

yi = where(lat_grid ge 30 and lat_grid le 60)
yi2 = where(lat_grid ge 30)
yi3 = where(lat_grid ge 20)
yi_mid = where(lat_grid eq 45)

for t = 0, nt-1 do begin
  for z = 0, nz2-1 do begin
    gd = where(finite(age_grid_adj_tseries[0,yi,z,t]) and tp_alt_grid[yi] le alt_grid2[z],ngd)
    if ngd gt 0 then begin
      lat_avg_nh_adj_tseries[z,t] = mean(lat_grid[gd+yi[0]])
      age_weights = 1. / age_grid_adj_tseries[1,yi[gd],z,t]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age_grid_adj_tseries[0,yi[gd],z,t]) / total(age_weights)
      age_nh_adj_tseries[0:1,z,t] = [wmean,err]
      for j = 2, 3 do age_nh_adj_tseries[j,z,t] = mean(age_grid_adj_tseries[j,yi,z,t],/nan)
      age_nh_adj_tseries[4,z,t] = total(age_grid_adj_tseries[4,yi[gd],z,t],/nan)
      age_nh_adj_tseries[5,z,t] = age_grid_adj_tseries[5,yi[gd[0]],z,t]
      if t mod 12 eq 11 or t mod 12 eq 0 or t mod 12 eq 1 then nh_adj_tseries_seas[z,t] = 0
      if t mod 12 ge 2 and t mod 12 le 4 then nh_adj_tseries_seas[z,t] = 1
      if t mod 12 ge 5 and t mod 12 le 7 then nh_adj_tseries_seas[z,t] = 2
      if t mod 12 ge 8 and t mod 12 le 10 then nh_adj_tseries_seas[z,t] = 3
  
      ;  Make adjustments based on sampled latitudes and seasonal gradients.
      yii = interpol(indgen(nyg),lat_grid,lat_avg_nh_adj_tseries[z,t])
      age_clim_samp = interpolate(age_grid_seas_avg[0,*,z,nh_adj_tseries_seas[z,t]],yii)
      age_nh_adj_tseries_lat_adj[*,z,t] = age_nh_adj_tseries[*,z,t]
      age_nh_adj_tseries_lat_adj[0,z,t] -= (age_clim_samp - age_grid_seas_avg[0,yi_mid,z,nh_adj_tseries_seas[z,t]])
      age_nh_adj_tseries_lat_adj[2,z,t] -= (age_clim_samp - age_grid_seas_avg[0,yi_mid,z,nh_adj_tseries_seas[z,t]])
    endif
    gd = where(finite(age_grid_adj_tseries[0,yi3,z,t]) and tp_alt_grid[yi3] le alt_grid2[z],ngd)
    if ngd gt 0 then begin
      age_weights = 1. / age_grid_adj_tseries[1,yi3[gd],z,t]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age_grid_adj_tseries[0,yi3[gd],z,t]) / total(age_weights)
      age_nhe_adj_tseries[0:1,z,t] = [wmean,err]
      for j = 2, 3 do age_nhe_adj_tseries[j,z,t] = mean(age_grid_adj_tseries[j,yi3,z,t],/nan)
      age_nhe_adj_tseries[4,z,t] = total(age_grid_adj_tseries[4,yi3[gd],z,t],/nan)
      age_nhe_adj_tseries[5,z,t] = age_grid_adj_tseries[5,yi3[gd[0]],z,t]
    endif
    gd = where(finite(n2o_grid_adj_tseries[0,yi,z,t]) and tp_alt_grid[yi] le alt_grid2[z],ngd)
    if ngd gt 0 then begin
;      print,alt_grid2[z],lat_grid[yi[gd]]
      lat_avg_n2o_nh_adj_tseries[z,t] = mean(lat_grid[gd+yi[0]])
      for j = 0, 2 do n2o_nh_adj_tseries[j,z,t] = mean(n2o_grid_adj_tseries[j,yi[gd],z,t],/nan)
      nh_adj_tseries_mon[z,t] = t mod 12

      ;  Make adjustments based on sampled latitudes and seasonal gradients.
      yii = interpol(indgen(nys),slat,lat_avg_n2o_nh_adj_tseries[z,t])
      n2o_clim_samp = interpolate(combo_n2o_seas_hires[*,z,nh_adj_tseries_mon[z,t]],yii)
      n2o_nh_adj_tseries_lat_adj[*,z,t] = n2o_nh_adj_tseries[*,z,t]
      n2o_nh_adj_tseries_lat_adj[0,z,t] -= (n2o_clim_samp - mean(combo_n2o_seas_hires[13:14,z,nh_adj_tseries_mon[z,t]]))
    endif
    gd = where(finite(ch4_grid_adj_tseries[0,yi,z,t]) and tp_alt_grid[yi] le alt_grid2[z],ngd)
    if ngd gt 0 then begin
      for j = 0, 2 do ch4_nh_adj_tseries[j,z,t] = mean(ch4_grid_adj_tseries[j,yi[gd],z,t],/nan)
    endif
  endfor

  for n = 0, nn-1 do begin
    n2o_time_series_sub = reform(n2o_grid_adj_tseries[0,yi2,*,t])
    age_time_series_sub = reform(age_grid_adj_tseries[0,yi2,*,t])
    gd = where(finite(n2o_time_series_sub) and finite(age_time_series_sub) and n2o_time_series_sub ge norm_grid[n]-dn/2. and $
      n2o_time_series_sub le norm_grid[n]+dn/2.,ngd)
    if ngd gt 0 then begin
      stats = moment(age_time_series_sub[gd],sdev=sdev,/nan)
      age_on_n2o_nh_adj_tseries[0:1,n,t] = [stats[0],sdev]
      if t mod 12 eq 11 or t mod 12 eq 0 or t mod 12 eq 1 then age_on_n2o_nh_adj_tseries_seas[n,t] = 0
      if t mod 12 ge 2 and t mod 12 le 4 then age_on_n2o_nh_adj_tseries_seas[n,t] = 1
      if t mod 12 ge 5 and t mod 12 le 7 then age_on_n2o_nh_adj_tseries_seas[n,t] = 2
      if t mod 12 ge 8 and t mod 12 le 10 then age_on_n2o_nh_adj_tseries_seas[n,t] = 3
    endif
  endfor
endfor


;  NH midlat avg profiles.
aaa = replicate(!values.f_nan,6,nz2) & age_nh_avg_profile_early = aaa & n2o_nh_avg_profile_early = aaa & age_nh_avg_profile_late = aaa & n2o_nh_avg_profile_late = aaa
ch4_nh_avg_profile_early = aaa & ch4_nh_avg_profile_late = aaa
ti1 = where(years lt 2000)
ti2 = where(years ge 2015)
for z = 0, nz2-1 do begin
  for j = 0, 3 do begin
    n2o_nh_avg_profile_early[j,z] = mean(n2o_nh_adj_tseries[j,z,ti1],/nan)
    n2o_nh_avg_profile_late[j,z] = mean(n2o_nh_adj_tseries[j,z,ti2],/nan)
    ch4_nh_avg_profile_early[j,z] = mean(ch4_nh_adj_tseries[j,z,ti1],/nan)
    ch4_nh_avg_profile_late[j,z] = mean(ch4_nh_adj_tseries[j,z,ti2],/nan)
    age_nh_avg_profile_early[j,z] = mean(age_nh_adj_tseries[j,z,ti1],/nan)
    age_nh_avg_profile_late[j,z] = mean(age_nh_adj_tseries[j,z,ti2],/nan)
  endfor
  age_nh_avg_profile_early[4,z] = total(age_nh_adj_tseries_lat_adj[4,z,ti1[gd]],/nan)
  age_nh_avg_profile_late[4,z] = total(age_nh_adj_tseries_lat_adj[4,z,ti2[gd]],/nan)
endfor
gd = where(finite(n2o_nh_avg_profile_early[0,*]))
for j = 0, 2, 2 do n2o_nh_avg_profile_early[j,gd] = smooth(n2o_nh_avg_profile_early[j,gd],3,/edge_truncate)
gd = where(finite(n2o_nh_avg_profile_late[0,*]))
for j = 0, 2, 2 do n2o_nh_avg_profile_late[j,gd] = smooth(n2o_nh_avg_profile_late[j,gd],3,/edge_truncate)
gd = where(finite(ch4_nh_avg_profile_early[0,*]))
for j = 0, 2, 2 do ch4_nh_avg_profile_early[j,gd] = smooth(ch4_nh_avg_profile_early[j,gd],3,/edge_truncate)
gd = where(finite(ch4_nh_avg_profile_late[0,*]))
for j = 0, 2, 2 do ch4_nh_avg_profile_late[j,gd] = smooth(ch4_nh_avg_profile_late[j,gd],3,/edge_truncate)
gd = where(finite(age_nh_avg_profile_early[0,*]))
for j = 0, 2, 2 do age_nh_avg_profile_early[j,gd] = smooth(age_nh_avg_profile_early[j,gd],3,/edge_truncate)
gd = where(finite(age_nh_avg_profile_late[0,*]))
for j = 0, 2, 2 do age_nh_avg_profile_late[j,gd] = smooth(age_nh_avg_profile_late[j,gd],3,/edge_truncate)


;  Make NH extratropical avg mean age vs tracer relationships.
nn2 = 50
dn2 = 0.02
norm_grid2 = findgen(nn2)*dn2+0.01
na = 81
da = 0.1
age_grid = findgen(na)*da
cent_interval = 0.95
ppp = replicate(!values.f_nan,nn2,5) & mean_age_co2_on_n2o_a0 = ppp & mean_age_sf6_on_n2o_a0 = ppp & mean_age_co2_on_n2o_b0 = ppp & mean_age_sf6_on_n2o_b0 = ppp
mean_age_on_n2o_b0 = ppp & mean_age_sf6_opt_on_n2o_b0 = ppp & mean_age_all_on_n2o_a0 = ppp & mean_age_all_on_n2o_b0 = ppp & elat_co2_on_n2o_a0 = ppp
elat_sf6_on_n2o_a0 = ppp & elat_all_on_n2o_a0 = ppp & elat_all_on_n2o_b0 = ppp & mean_age_co2_on_n2o_b0_mlat = ppp & mean_age_co2_on_n2o_b0_hlat = ppp
mean_age_co2_on_ch4_a0 = ppp & mean_age_sf6_on_ch4_a0 = ppp & mean_age_co2_on_ch4_b0 = ppp & mean_age_sf6_on_ch4_b0 = ppp & mean_age_all_on_n2o_e_a0 = ppp
mean_age_on_ch4_b0 = ppp & mean_age_sf6_opt_on_ch4_b0 = ppp & mean_age_all_on_ch4_a0 = ppp & mean_age_all_on_ch4_b0 = ppp & elat_all_on_n2o_e_a0 = ppp
mean_age_co2_on_f12_a0 = ppp & mean_age_sf6_on_f12_a0 = ppp & mean_age_co2_on_f12_b0 = ppp & mean_age_sf6_on_f12_b0 = ppp
mean_age_on_f12_b0 = ppp & mean_age_sf6_opt_on_f12_b0 = ppp & mean_age_all_on_f12_a0 = ppp & mean_age_all_on_f12_b0 = ppp
mean_age_co2_on_n2o_a2 = ppp & mean_age_sf6_on_n2o_a2 = ppp & mean_age_all_on_n2o_a2 = ppp & elat_all_on_n2o_a2 = ppp
ss = replicate(!values.f_nan,na,nn2) & mean_age_co2_on_n2o_a0_hist = ss & mean_age_sf6_on_n2o_a0_hist = ss & mean_age_co2_on_n2o_b0_hist = ss
mean_age_sf6_on_n2o_b0_hist = ss & mean_age_on_n2o_b0_hist = ss & mean_age_sf6_opt_on_n2o_b0_hist = ss & mean_age_all_on_n2o_hist = ss
mean_age_co2_on_ch4_a0_hist = ss & mean_age_sf6_on_ch4_a0_hist = ss & mean_age_co2_on_ch4_b0_hist = ss
mean_age_sf6_on_ch4_b0_hist = ss & mean_age_on_ch4_b0_hist = ss & mean_age_sf6_opt_on_ch4_b0_hist = ss & mean_age_all_on_ch4_hist = ss
mean_age_co2_on_f12_a0_hist = ss & mean_age_sf6_on_f12_a0_hist = ss & mean_age_co2_on_f12_b0_hist = ss
mean_age_sf6_on_f12_b0_hist = ss & mean_age_on_f12_b0_hist = ss & mean_age_sf6_opt_on_f12_b0_hist = ss & mean_age_all_on_f12_hist = ss
mean_age_all_on_n2o_hist_cent = replicate(!values.f_nan,2,nn2) & mean_age_all_on_ch4_hist_cent = mean_age_all_on_n2o_hist_cent
mean_age_all_on_f12_hist_cent = mean_age_all_on_n2o_hist_cent
for n = 0, nn2-1 do begin
  gd = where(n2o_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_a0.age_opt_all[gd]
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_n2o_a0[n,*] = [stats[0],sdev,wmean,err,ngd]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_co2_on_n2o_a0_hist[*,n] = H / total(H)
    stats = moment(elat_m_adj_a_early[ico2_a0[gd]],sdev=sdev,/nan)
    elat_co2_on_n2o_a0[n,0:1] = [stats[0],sdev]
;    if n le 7 then print,'a0 co2 ',norm_grid2[n],stats[0],sdev,wmean,err,ngd;,dat_co2_a0.age_opt_all[gd],elat_m_adj_a_early[ico2_a0[gd]]
  endif
  gd = where(n2o_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_a0.age_opt_all[gd]
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    chk = where(finite(datsub))
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_on_n2o_a0[n,*] = [stats[0],sdev,wmean,err,ngd]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_sf6_on_n2o_a0_hist[*,n] = H / total(H)
;    if n ge nn2-2 then print,'a0 sf6 ',norm_grid2[n],stats[0],sdev,wmean,err,ngd,dat_sf6_a0.age_opt_all[gd];,elat_m_adj_a_early[ico2_a0[gd]]
    stats = moment(elat_m_adj_a_early[isf6_a0[gd]],sdev=sdev,/nan)
    elat_sf6_on_n2o_a0[n,0:1] = [stats[0],sdev]
  endif 
   
  datsub_b = 0. & age_err_b = 0. & elat_sub_b = 0. & ngd_b = 0
  datsub_b_mlat = 0. & age_err_b_mlat = 0. & elat_sub_b_mlat = 0. & ngd_b_mlat = 0
  datsub_b_hlat = 0. & age_err_b_hlat = 0. & elat_sub_b_hlat = 0. & ngd_b_hlat = 0
  gd = where(n2o_norm_all_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and $
    finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    elat_sub_b = [elat_sub_b,elat_m_adj_b_early[ico2_b[gd]]]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    chk1 = where(~finite(age_err),nchk1)
    if nchk1 gt 0 then age_err[chk1] = 1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_n2o_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_co2_on_n2o_b0_hist[*,n] = H / total(H)
;    if n eq 20 then print,'b0 co2 ',norm_grid2[n],stats[0],sdev,wmean,err,ngd,dat_co2_b0.age_opt_all[gd],elat_m_adj_b_early[ico2_b[gd]],times_all_b[ico2_b[gd]]
  endif
  gd = where(n2o_norm_all_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and $
     elat_m_adj_b_early[ico2_b] le 65 and finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b_mlat = [datsub_b_mlat,datsub]
    ngd_b_mlat += ngd
    elat_sub_b_mlat = [elat_sub_b_mlat,elat_m_adj_b_early[ico2_b[gd]]]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    chk1 = where(~finite(age_err),nchk1)
    if nchk1 gt 0 then age_err[chk1] = 1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_n2o_b0_mlat[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(n2o_norm_all_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] gt 75 and $
    finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b_hlat = [datsub_b_hlat,datsub]
    ngd_b_hlat += ngd
    elat_sub_b_hlat = [elat_sub_b_hlat,elat_m_adj_b_early[ico2_b[gd]]]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    chk1 = where(~finite(age_err),nchk1)
    if nchk1 gt 0 then age_err[chk1] = 1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_n2o_b0_hlat[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(n2o_norm_all_b[isf6_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_b0.age_opt_all[gd]
    datusub = dat_sf6_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_on_n2o_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_sf6_on_n2o_b0_hist[*,n] = H / total(H)
;    if n le 7 then print,'b0 sf6 ',norm_grid2[n],stats[0],sdev,wmean,err,ngd;,dat_co2_a0.age_opt_all[gd],elat_m_adj_a_early[ico2_a0[gd]]
  endif
  gd = where(n2o_norm_all_b[iage_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_b0.age_opt_all[gd]
    datusub = replicate(1.0,ngd)
    age_err = 0.5 * (datusub)
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_on_n2o_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_on_n2o_b0_hist[*,n] = H / total(H)
;    if n le 7 then print,'b0 opt ',norm_grid2[n],stats[0],sdev,wmean,err,ngd;,dat_co2_a0.age_opt_all[gd],elat_m_adj_a_early[ico2_a0[gd]]
  endif
  ;  Combine the balloon opt ages and SF6 ages with no coincident opt value.  
  gd = where(n2o_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and n2o_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  gd2 = where(n2o_norm_all_b[iage_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd2)
  if ngd gt 0 and ngd2 gt 0 then begin
    lat1 = elat_m_adj_b_early[isf6_b[gd]]
    lat2 = elat_m_adj_b_early[iage_b[gd2]]
    sf61 = sf6_b[isf6_b[gd]]
    sf62 = sf6_b[iage_b[gd2]]
    datsub1 = dat_sf6_b0.age_opt_all[gd]
    datsub2 = dat_b0.age_opt_all[gd2]
    datsub = datsub1
    elatsub = lat1
    for i = 0, ngd-1 do begin
      chk = where(lat2 eq lat1[i] and sf62 eq sf61[i],nchk)
      if nchk gt 0 then begin
        datsub[i] = datsub2[chk]
        elatsub[i] = lat2[chk]
      endif
    endfor
    datsub_b = [datsub_b,datsub]
    elat_sub_b = [elat_sub_b,elatsub]
    ngd_b += ngd
    datusub = replicate(1.0,ngd)
    age_err = 0.5 * (datusub)
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_opt_on_n2o_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_sf6_opt_on_n2o_b0_hist[*,n] = H / total(H)
  endif
  if n_elements(datsub_b) gt 1 then begin
    age_err_b = age_err_b[1:-1]
    datsub_b = datsub_b[1:-1]
    chk = where(finite(age_err_b) and finite(datsub_b))
    age_weights = 1. / age_err_b[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub_b[chk]) / total(age_weights)
    stats = moment(datsub_b[chk],sdev=sdev,/nan)
    mean_age_all_on_n2o_b0[n,*] = [stats[0],sdev,wmean,err,ngd_b]
;    if n eq 20 then print,norm_grid2[n],stats[0],sdev,wmean,err,ngd_b,datsub_b[chk],age_err_b[chk],age_weights
    stats = moment(elat_sub_b[1:-1],sdev=sdev,/nan)
    elat_all_on_n2o_b0[n,0:1] = [stats[0],sdev]
  endif
  
  ;  Combine all of the data.
  datsub = 0. & age_err = 0. & elat_sub = 0. & ngd_a = 0
  gd = where(n2o_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_early[ico2_a0[gd]]]
  endif
  gd = where(n2o_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_early[isf6_a0[gd]]]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    age_err = age_err[1:-1]
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub[chk],sdev=sdev,/nan)
    mean_age_all_on_n2o_a0[n,*] = [stats[0],sdev,wmean,err,ngd_a]
;    if n ge nn2-4 then print,norm_grid2[n],stats[0],sdev,wmean,err,ngd_a;,datsub[chk],age_err_b[chk],age_weights
    stats = moment(elat_sub[1:-1],sdev=sdev,/nan)
    elat_all_on_n2o_a0[n,0:1] = [stats[0],sdev]
  endif
  gd = where(n2o_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_b0.age_opt_all[gd]]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = [age_err,reform(0.5 * (datusub[1,*] - datusub[0,*]))]
  endif
  gd = where(n2o_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and n2o_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_b0.age_opt_all[gd]]
    datusub = dat_sf6_b0.age_range_all[*,gd]
    age_err = [age_err,reform(0.5 * (datusub[1,*] - datusub[0,*]))]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_all_on_n2o_hist[*,n] = H / total(H)
    csum = 0.
    for i = 0, na-1 do begin
      csum += mean_age_all_on_n2o_hist[i,n]
      if csum ge (1.-cent_interval) / 2. then begin
        mean_age_all_on_n2o_hist_cent[0,n] = age_grid[i]
        i = na-1
      endif
    endfor
    csum = 0.
    for i = 0, na-1 do begin
      csum += mean_age_all_on_n2o_hist[i,n]
      if csum ge cent_interval + (1.-cent_interval) / 2. then begin
        mean_age_all_on_n2o_hist_cent[1,n] = age_grid[i-1]
        i = na-1
      endif
    endfor
  endif

  ;  Use age adjusted entry values for N2O normalizing.
  datsub = 0. & age_err = 0. & elat_sub = 0. & ngd_a = 0
  gd = where(bins_co2_a0.n2o_norm_e ge norm_grid2[n]-dn2/2. and bins_co2_a0.n2o_norm_e le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_early[ico2_a0[gd]]]
  endif
  gd = where(bins_sf6_a0.n2o_norm_e ge norm_grid2[n]-dn2/2. and bins_sf6_a0.n2o_norm_e le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_early[isf6_a0[gd]]]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    age_err = age_err[1:-1]
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub[chk],sdev=sdev,/nan)
    mean_age_all_on_n2o_e_a0[n,*] = [stats[0],sdev,wmean,err,ngd_a]
    stats = moment(elat_sub[1:-1],sdev=sdev,/nan)
    elat_all_on_n2o_e_a0[n,0:1] = [stats[0],sdev]
  endif


  ;  Late time binning.
  datsub = 0. & age_err = 0. & elat_sub = 0. & ngd_a = 0
;  gd = where(n2o_norm_a[iage_a] ge norm_grid2[n]-dn2/2. and n2o_norm_a[iage_a] le norm_grid2[n]+dn2/2. and elat_m_adj_a_late[iage_a] ge 30,ngd)
;  if ngd gt 0 then begin
;    datsub = [datsub,dat_a2.age_opt_all[gd]]
;    ngd_a += ngd
;    datusub = dat_sf6_a2.age_range_all[*,gd]
;    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
;    chk = where(age_e lt 0.1,nchk)
;    if nchk gt 0 then age_e[chk] = 0.1
;    age_err = [age_err,age_e]
;    elat_sub = [elat_sub,elat_m_adj_a_late[iage_a[gd]]]
;  endif
  gd = where(n2o_norm_a[isf6_a] ge norm_grid2[n]-dn2/2. and n2o_norm_a[isf6_a] le norm_grid2[n]+dn2/2. and elat_m_adj_a_late[isf6_a] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_a2.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_sf6_a2.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_late[isf6_a0[gd]]]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    age_err = age_err[1:-1]
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub[chk],sdev=sdev,/nan)
    mean_age_all_on_n2o_a2[n,*] = [stats[0],sdev,wmean,err,ngd_a]
    stats = moment(elat_sub[1:-1],sdev=sdev,/nan)
    elat_all_on_n2o_a2[n,0:1] = [stats[0],sdev]
  endif


  ;  CH4 binning
  gd = where(ch4_all_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_a0.age_opt_all[gd]
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_ch4_a0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(ch4_all_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_a0.age_opt_all[gd]
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    chk = where(finite(datsub))
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_on_ch4_a0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  datsub_b = 0. & age_err_b = 0. & ngd_b = 0
  gd = where(ch4_all_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and $
    finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    chk1 = where(~finite(age_err),nchk1)
    if nchk1 gt 0 then age_err[chk1] = 1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_ch4_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(ch4_all_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_b0.age_opt_all[gd]
    datusub = dat_sf6_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_on_ch4_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(ch4_all_norm_b[iage_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_b0.age_opt_all[gd]
    datusub = replicate(1.0,ngd)
    age_err = 0.5 * (datusub)
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_on_ch4_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  ;  Combine the balloon opt ages and SF6 ages with no coincident opt value.
  gd = where(ch4_all_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  gd2 = where(ch4_all_norm_b[iage_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd2)
  if ngd gt 0 and ngd2 gt 0 then begin
    lat1 = elat_m_adj_b_early[isf6_b[gd]]
    lat2 = elat_m_adj_b_early[iage_b[gd2]]
    sf61 = sf6_b[isf6_b[gd]]
    sf62 = sf6_b[iage_b[gd2]]
    ch41 = ch4_all_norm_b[isf6_b[gd]]
    ch42 = ch4_all_norm_b[iage_b[gd2]]
    datsub1 = dat_sf6_b0.age_opt_all[gd]
    datsub2 = dat_b0.age_opt_all[gd2]
    datsub = datsub1
    for i = 0, ngd-1 do begin
      chk = where(lat2 eq lat1[i] and sf62 eq sf61[i] and ch42 eq ch41[i],nchk)
      if nchk gt 0 then datsub[i] = datsub2[chk]
    endfor
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    datusub = replicate(1.0,ngd)
    age_err = 0.5 * (datusub)
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_opt_on_ch4_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  if n_elements(datsub_b) gt 1 then begin
    age_err_b = age_err_b[1:-1]
    datsub_b = datsub_b[1:-1]
    chk = where(finite(age_err_b) and finite(datsub_b))
    age_weights = 1. / age_err_b[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub_b[chk]) / total(age_weights)
    stats = moment(datsub_b[chk],sdev=sdev,/nan)
    mean_age_all_on_ch4_b0[n,*] = [stats[0],sdev,wmean,err,ngd_b]
  endif

  ;  Combine all of the data.
  datsub = 0. & age_err = 0. & ngd_a = 0
  gd = where(ch4_all_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_early[ico2_a0[gd]]]
  endif
  gd = where(ch4_all_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    age_err = age_err[1:-1]
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub[chk],sdev=sdev,/nan)
    mean_age_all_on_ch4_a0[n,*] = [stats[0],sdev,wmean,err,ngd_a]
  endif
  gd = where(ch4_all_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_b0.age_opt_all[gd]]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = [age_err,reform(0.5 * (datusub[1,*] - datusub[0,*]))]
  endif
  gd = where(ch4_all_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_b0.age_opt_all[gd]]
    datusub = dat_sf6_b0.age_range_all[*,gd]
    age_err = [age_err,reform(0.5 * (datusub[1,*] - datusub[0,*]))]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_all_on_ch4_hist[*,n] = H / total(H)
    csum = 0.
    for i = 0, na-1 do begin
      csum += mean_age_all_on_ch4_hist[i,n]
      if csum ge (1.-cent_interval) / 2. then begin
        mean_age_all_on_ch4_hist_cent[0,n] = age_grid[i]
        i = na-1
      endif
    endfor
    csum = 0.
    for i = 0, na-1 do begin
      csum += mean_age_all_on_ch4_hist[i,n]
      if csum ge cent_interval + (1.-cent_interval) / 2. then begin
        mean_age_all_on_ch4_hist_cent[1,n] = age_grid[i-1]
        i = na-1
      endif
    endfor
  endif

  ;  F12 binning
  gd = where(f12_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_a0.age_opt_all[gd]
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_f12_a0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(f12_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_a0.age_opt_all[gd]
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    chk = where(finite(datsub))
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_on_f12_a0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  datsub_b = 0. & age_err_b = 0. & ngd_b = 0
  gd = where(f12_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and $
    finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    chk1 = where(~finite(age_err),nchk1)
    if nchk1 gt 0 then age_err[chk1] = 1
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_co2_on_f12_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(f12_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_b0.age_opt_all[gd]
    datusub = dat_sf6_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_on_f12_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  gd = where(f12_norm_b[iage_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_b0.age_opt_all[gd]
    datusub = replicate(1.0,ngd)
    age_err = 0.5 * (datusub)
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_on_f12_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  ;  Combine the balloon opt ages and SF6 ages with no coincident opt value.
  gd = where(f12_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  gd2 = where(f12_norm_b[iage_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd2)
  if ngd gt 0 and ngd2 gt 0 then begin
    lat1 = elat_m_adj_b_early[isf6_b[gd]]
    lat2 = elat_m_adj_b_early[iage_b[gd2]]
    sf61 = sf6_b[isf6_b[gd]]
    sf62 = sf6_b[iage_b[gd2]]
    f121 = f12_norm_b[isf6_b[gd]]
    f122 = f12_norm_b[iage_b[gd2]]
    datsub1 = dat_sf6_b0.age_opt_all[gd]
    datsub2 = dat_b0.age_opt_all[gd2]
    datsub = datsub1
    for i = 0, ngd-1 do begin
      chk = where(lat2 eq lat1[i] and sf62 eq sf61[i] and f122 eq f121[i],nchk)
      if nchk gt 0 then datsub[i] = datsub2[chk]
    endfor
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    datusub = replicate(1.0,ngd)
    age_err = 0.5 * (datusub)
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    age_weights = 1. / age_err^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub) / total(age_weights)
    stats = moment(datsub,sdev=sdev,/nan)
    mean_age_sf6_opt_on_f12_b0[n,*] = [stats[0],sdev,wmean,err,ngd]
  endif
  if n_elements(datsub_b) gt 1 then begin
    age_err_b = age_err_b[1:-1]
    datsub_b = datsub_b[1:-1]
    chk = where(finite(age_err_b) and finite(datsub_b))
    age_weights = 1. / age_err_b[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub_b[chk]) / total(age_weights)
    stats = moment(datsub_b[chk],sdev=sdev,/nan)
    mean_age_all_on_f12_b0[n,*] = [stats[0],sdev,wmean,err,ngd_b]
  endif

  ;  Combine all of the data.
  datsub = 0. & age_err = 0. & ngd_a = 0
  gd = where(f12_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
    elat_sub = [elat_sub,elat_m_adj_a_early[ico2_a0[gd]]]
  endif
  gd = where(f12_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and $
    times_all_a0[isf6_a0] lt 1999,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_sf6_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    age_err = age_err[1:-1]
    chk = where(finite(datsub),nchk)
    age_weights = 1. / age_err[chk]^2
    err = 1. / sqrt(total(age_weights))
    wmean = total(age_weights * datsub[chk]) / total(age_weights)
    stats = moment(datsub[chk],sdev=sdev,/nan)
    mean_age_all_on_f12_a0[n,*] = [stats[0],sdev,wmean,err,ngd_a]
  endif
  gd = where(f12_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_b0.age_opt_all[gd]]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = [age_err,reform(0.5 * (datusub[1,*] - datusub[0,*]))]
  endif
  gd = where(f12_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and $
    times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_b0.age_opt_all[gd]]
    datusub = dat_sf6_b0.age_range_all[*,gd]
    age_err = [age_err,reform(0.5 * (datusub[1,*] - datusub[0,*]))]
  endif
  if n_elements(datsub) gt 1 then begin
    datsub = datsub[1:-1]
    H = histogram(datsub,binsize=da,max=age_grid[-1],min=0)
    mean_age_all_on_f12_hist[*,n] = H / total(H)
    csum = 0.
    for i = 0, na-1 do begin
      csum += mean_age_all_on_f12_hist[i,n]
      if csum ge (1.-cent_interval) / 2. then begin
        mean_age_all_on_f12_hist_cent[0,n] = age_grid[i]
        i = na-1
      endif
    endfor
    csum = 0.
    for i = 0, na-1 do begin
      csum += mean_age_all_on_f12_hist[i,n]
      if csum ge cent_interval + (1.-cent_interval) / 2. then begin
        mean_age_all_on_f12_hist_cent[1,n] = age_grid[i-1]
        i = na-1
      endif
    endfor
  endif
endfor

;  Use SF6 values for mean ages less than 1 year.
mean_age_all_on_n2o_a0[-1,*] = mean_age_sf6_on_n2o_a0[-1,*]
mean_age_all_on_n2o_b0[-1,*] = mean_age_sf6_on_n2o_b0[-1,*]
mean_age_all_on_ch4_a0[-1,*] = mean_age_sf6_on_ch4_a0[-1,*]
mean_age_all_on_f12_a0[-2:-1,*] = mean_age_sf6_on_f12_a0[-2:-1,*]
mean_age_all_on_f12_b0[-1,*] = mean_age_sf6_on_f12_b0[-1,*]

mean_age_all_on_n2o_0 = mean_age_all_on_n2o_a0
for n = 0, nn2-1 do begin
  for j = 0, 3 do begin
    mean_age_all_on_n2o_0[n,j] = mean([mean_age_all_on_n2o_a0[n,j],mean_age_all_on_n2o_b0[n,j]],/nan)
  endfor
  mean_age_all_on_n2o_0[n,4] = total([mean_age_all_on_n2o_a0[n,4],mean_age_all_on_n2o_b0[n,4]],/nan)
endfor

age_grid_h = age_grid

save,norm_grid2,mean_age_all_on_n2o_a0,mean_age_all_on_n2o_b0,age_grid_h,mean_age_all_on_n2o_hist,mean_age_all_on_ch4_a0,mean_age_all_on_ch4_b0, $
  mean_age_all_on_ch4_hist,mean_age_all_on_f12_a0,mean_age_all_on_f12_b0,mean_age_all_on_f12_hist,mean_age_all_on_n2o_hist_cent,mean_age_all_on_ch4_hist_cent, $
  mean_age_all_on_f12_hist_cent,mean_age_all_on_n2o_0,mean_age_all_on_n2o_e_a0,mean_age_all_on_n2o_a2,filename=dir+'Balloon/Mean_age_vs_N2O_1990s.sav'


;  Trends
for z = 0, nz2-1 do for y = 0, nyg-1 do begin
  gd = where(finite(age_grid_adj_tseries[0,y,z,*]),ngd)
  if ngd ge 10 then begin
    if years[gd[-1]]-years[gd[0]] ge 20 then begin
      age_grid_trends[y,z,*] = linfit(years[gd],age_grid_adj_tseries[0,y,z,gd],sigma=sigma)
      age_grid_trends_sigma[y,z,*] = sigma
    endif
  endif
endfor
for z = 0, nz2-1 do begin
  gd = where(finite(n2o_nh_adj_tseries[0,z,*]),ngd)
  if ngd ge 10 then begin
    n2o_nh_trends[z,*] = linfit(years[gd],n2o_nh_adj_tseries[0,z,gd],sigma=sigma)
    n2o_nh_trends_sigma[z,*] = sigma
    n2o_nh_lat_adj_trends[z,*] = linfit(years[gd],n2o_nh_adj_tseries_lat_adj[0,z,gd],sigma=sigma)
    n2o_nh_lat_adj_trends_sigma[z,*] = sigma
  endif
  gd = where(finite(ch4_nh_adj_tseries[0,z,*]) and years ge 1995,ngd)
  if ngd ge 10 then begin
    ch4_nh_trends[z,*] = linfit(years[gd],ch4_nh_adj_tseries[0,z,gd],sigma=sigma)
    ch4_nh_trends_sigma[z,*] = sigma
  endif
  gd = where(finite(age_nh_adj_tseries[0,z,*]),ngd)
  if ngd ge 10 then begin
    age_nh_trends[z,*] = linfit(years[gd],age_nh_adj_tseries[0,z,gd],sigma=sigma)
    age_nh_trends_sigma[z,*] = sigma
    age_nh_lat_adj_trends[z,*] = linfit(years[gd],age_nh_adj_tseries_lat_adj[2,z,gd],sigma=sigma)
    age_nh_lat_adj_trends_sigma[z,*] = sigma
  endif
endfor
ddd = replicate(!values.f_nan,6,nz3,nt) & age_nh_coarse = ddd & n2o_nh_coarse = ddd & age_nh_coarse_lat_adj = ddd & ch4_nh_coarse = ddd & age_nhe_coarse = ddd
lat_avg_nh_coarse = replicate(!values.f_nan,nz3,nt) & rse_orig = replicate(!values.f_nan,nz3) & rse_lat_adj = rse_orig
age_nh_coarse_mean = replicate(!values.f_nan,nz3) & age_nh_coarse_lat_adj_mean = age_nh_coarse_mean & n2o_nh_coarse_mean = age_nh_coarse_mean
for z = 0, nz3-1 do begin
  zi = where(alt_grid2 ge alt_grid3[z]-0.5 and alt_grid2 le alt_grid3[z]+0.5,nzi)
  for t = 0, nt-1 do begin
    zii = where(finite(age_nh_adj_tseries_lat_adj[2,zi,t]),nzi)
    if nzi ge 1 then begin
      lat_avg_nh_coarse[z,t] = mean(lat_avg_nh_adj_tseries[zi[zii],t],/nan)
      age_weights = 1. / age_nh_adj_tseries[1,zi[zii],t]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age_nh_adj_tseries[0,zi[zii],t]) / total(age_weights)
      age_nh_coarse[0:1,z,t] = [wmean,err]
      age_weights = 1. / age_nh_adj_tseries_lat_adj[1,zi[zii],t]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age_nh_adj_tseries_lat_adj[0,zi[zii],t]) / total(age_weights)
      age_nh_coarse_lat_adj[0:1,z,t] = [wmean,err]
      for j = 2, 3 do begin
        age_nh_coarse[j,z,t] = mean(age_nh_adj_tseries[j,zi[zii],t],/nan)
        if j eq 3 and ~finite(age_nh_coarse[j,z,t]) then age_nh_coarse[j,z,t] = 0.1
        age_nh_coarse_lat_adj[j,z,t] = mean(age_nh_adj_tseries_lat_adj[j,zi[zii],t],/nan)        
      endfor
      age_nh_coarse[4,z,t] = total(age_nh_adj_tseries[4,zi[zii],t],/nan)
      age_nh_coarse[5,z,t] = age_nh_adj_tseries[5,zi[zii],t]
      age_nh_coarse_lat_adj[4,z,t] = total(age_nh_adj_tseries_lat_adj[4,zi[zii],t],/nan)
      age_nh_coarse_lat_adj[5,z,t] = age_nh_adj_tseries_lat_adj[5,zi[zii],t]
    endif
    zii = where(finite(age_nhe_adj_tseries[2,zi,t]),nzi)
    if nzi ge 1 then begin
      age_weights = 1. / age_nhe_adj_tseries[1,zi[zii],t]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age_nhe_adj_tseries[0,zi[zii],t]) / total(age_weights)
      age_nhe_coarse[0:1,z,t] = [wmean,err]
      for j = 2, 3 do begin
        age_nhe_coarse[j,z,t] = mean(age_nhe_adj_tseries[j,zi[zii],t],/nan)
        if j eq 3 and ~finite(age_nhe_coarse[j,z,t]) then age_nhe_coarse[j,z,t] = 0.1
      endfor
      age_nhe_coarse[4,z,t] = total(age_nhe_adj_tseries[4,zi[zii],t],/nan)
      age_nhe_coarse[5,z,t] = age_nhe_adj_tseries[5,zi[zii],t]
    endif
;    if alt_grid3[z] le 22 then begin
      zii = where(finite(n2o_nh_adj_tseries_lat_adj[0,zi,t]),nzi)
      if nzi ge 1 then n2o_nh_coarse[0,z,t] = mean(n2o_nh_adj_tseries_lat_adj[0,zi[zii],t],/nan)
      zii = where(finite(ch4_nh_adj_tseries[0,zi,t]),nzi)
      if nzi ge 1 then ch4_nh_coarse[0,z,t] = mean(ch4_nh_adj_tseries[0,zi[zii],t],/nan)
;    endif else begin
;      zii = where(finite(n2o_nh_adj_tseries[0,zi,t]),nzi)
;      if nzi ge 2 then n2o_nh_coarse[z,t] = mean(n2o_nh_adj_tseries[0,zi[zii],t],/nan)      
;    endelse
  endfor
  gd = where(finite(age_nh_coarse_lat_adj[2,z,*]),ngd)
  if ngd ge 10 then begin
    if years[gd[-1]]-years[gd[0]] ge 20 then begin
      age_nh_coarse_mean[z] = mean(age_nh_coarse[2,z,gd])
      age_nh_trends_coarse[z,*] = linfit(years[gd],age_nh_coarse[2,z,gd],sigma=sigma)
      age_nh_trends_coarse_sigma[z,*] = sigma
      age_nh_trends_per_coarse[z,*] = 100 * age_nh_trends_coarse[z,*] / age_nh_coarse_mean[z]
      age_nh_trends_per_coarse_sigma[z,*] = 100 * sigma / age_nh_coarse_mean[z]
      age_nh_coarse_lat_adj_mean[z] = mean(age_nh_coarse_lat_adj[2,z,gd])
      age_nh_lat_adj_trends_coarse[z,*] = linfit(years[gd],age_nh_coarse_lat_adj[2,z,gd],sigma=sigma)
      age_nh_lat_adj_trends_coarse_sigma[z,*] = sigma
      age_nh_lat_adj_trends_per_coarse[z,*] = 100 * age_nh_lat_adj_trends_coarse[z,*] / age_nh_coarse_lat_adj_mean[z]
      age_nh_lat_adj_trends_per_coarse_sigma[z,*] = 100 * sigma / age_nh_coarse_lat_adj_mean[z]
    endif
    gd = where(finite(age_nh_coarse_lat_adj[2,z,*]) and years ge 2011 and (age_nh_coarse[5,z,*] eq 3 or (age_nh_coarse[5,z,*] eq 1 and years ge 2015)),ngd)
    age_nh_trends_ac_coarse[z,*] = linfit(years[gd],age_nh_coarse[2,z,gd],sigma=sigma)
    age_nh_trends_ac_coarse_sigma[z,*] = sigma
    age_nh_lat_adj_trends_ac_coarse[z,*] = linfit(years[gd],age_nh_coarse_lat_adj[2,z,gd],sigma=sigma)
    age_nh_lat_adj_trends_ac_coarse_sigma[z,*] = sigma
  endif
  
  ;  Residual standard error of data around linear fit.
  gd = where(finite(age_nh_coarse[2,z,*]),ngd)
  rse_orig[z] = sqrt(total((age_nh_coarse[2,z,gd] - (age_nh_trends_coarse[z,1]*years[gd] + age_nh_trends_coarse[z,0]))^2)/(ngd-2))
  gd = where(finite(age_nh_coarse_lat_adj[2,z,*]),ngd)
  rse_lat_adj[z] = sqrt(total((age_nh_coarse_lat_adj[2,z,gd] - (age_nh_lat_adj_trends_coarse[z,1]*years[gd] + age_nh_lat_adj_trends_coarse[z,0]))^2)/(ngd-2))

  gd = where(finite(n2o_nh_coarse[0,z,*]),ngd)
  if ngd ge 10 then begin
    if years[gd[-1]]-years[gd[0]] ge 20 then begin
      n2o_nh_lat_adj_trends_coarse[z,*] = linfit(years[gd],n2o_nh_coarse[0,z,gd],sigma=sigma)
      n2o_nh_lat_adj_trends_coarse_sigma[z,*] = sigma
      n2o_nh_coarse_mean[z] = mean(n2o_nh_coarse[0,z,gd])
      n2o_nh_lat_adj_trends_per_coarse[z,*] = 100 * n2o_nh_lat_adj_trends_coarse[z,*] / n2o_nh_coarse_mean[z]
      n2o_nh_lat_adj_trends_per_coarse_sigma[z,*] = 100 * sigma / n2o_nh_coarse_mean[z]
    endif
  endif
  gd = where(finite(ch4_nh_coarse[0,z,*]) and years ge 1998,ngd)
  if ngd ge 10 then begin
    if years[gd[-1]]-years[gd[0]] ge 20 then begin
      ch4_nh_trends_coarse[z,*] = linfit(years[gd],ch4_nh_coarse[0,z,gd],sigma=sigma)
      ch4_nh_trends_coarse_sigma[z,*] = sigma
    endif
  endif
endfor

for n = 0, nn-1 do begin
  gd = where(finite(age_on_n2o_nh_adj_tseries[0,n,*]),ngd)
  if ngd ge 10 then begin
    age_on_n2o_nh_lat_adj_trends[n,*] = linfit(years[gd],age_on_n2o_nh_adj_tseries[0,n,gd],sigma=sigma)
    age_on_n2o_nh_lat_adj_trends_sigma[n,*] = sigma
  endif
endfor

save,age_nh_coarse,age_nh_coarse_lat_adj,age_nh_adj_tseries,age_nh_adj_tseries_lat_adj,age_nh_lat_adj_trends_coarse,age_nh_lat_adj_trends_coarse_sigma, $
  n2o_nh_lat_adj_trends_coarse,n2o_nh_lat_adj_trends_coarse_sigma,age_nh_trends_coarse,age_nh_trends_coarse_sigma,ch4_nh_trends_coarse,ch4_nh_trends_coarse_sigma, $
  age_nhe_adj_tseries,age_nhe_coarse,age_nh_trends_per_coarse,age_nh_lat_adj_trends_per_coarse,age_nh_trends_per_coarse_sigma,age_nh_lat_adj_trends_per_coarse_sigma, $
  age_nh_lat_adj_trends_ac_coarse,age_nh_lat_adj_trends_ac_coarse_sigma,n2o_nh_lat_adj_trends_per_coarse,n2o_nh_lat_adj_trends_per_coarse_sigma, $
  filename=dir+'Balloon/Mean_age_and_N2O_trend_alt_profile_elat_adj.sav'

;  Write the data plotted in Supplementary Figures 3 and 5 to NetCDF (code/figures/supp_fig3.pro, supp_fig5.pro).
@export_figure_data_check_elats
if stage eq 'time_series' then return

lat_grid_a = lat_grid



;  Write out mean age time series into an ascii file.
;z = 10
;openw,1,dir+'Balloon/Mean_age_time_series.txt'
;printf,1,'Year    Mean Age   Uncertainty'
;;for t = 0, nt-1 do printf,1,years[t],age_nh_coarse[2:3,z,t]
;for t = 0, nt-1 do printf,1,years[t],age_nh_coarse_lat_adj[2:3,z,t]
;close,1


restore,dir+'Balloon/Sample_type_time_series.sav'
restore,dir+'Aircraft/Airborne_Save_Files/N2O_vs_mean_age_coarse.sav'
ai = interpol(findgen(n_elements(norm_grid_obs)),bins_a0_n2o_norm_age_nh_coarse,3)
n2o_3yr_old = interpolate(norm_grid_obs,ai)
ni = interpol(findgen(n_elements(n2o_global)),n2o_date,1997)
n2o_surface = interpolate(n2o_global,ni)
n2o_3yr_old *= n2o_surface
ai = interpol(findgen(n_elements(norm_grid_obs)),bins_a2_n2o_norm_age_nh_trend,3)
n2o_3yr_new = interpolate(norm_grid_obs,ai)
ni = interpol(findgen(n_elements(n2o_global)),n2o_date,2022)
n2o_surface = interpolate(n2o_global,ni)
n2o_3yr_new *= n2o_surface

;restore,dir+'Models/CCMI/CCMI-2022_refd1_means.sav'
;alt_c = -7.*alog(pres_c/1e3)
;nyc = n_elements(lat_c)
;nzc = n_elements(alt_c)
;print,years_c[721]
;model_n2o_seas = replicate(!values.f_nan,nyc,nzc,12)
;for t = 0, 11 do for z = 0, nzc-1 do for y = 0, nyc-1 do begin
;  ti = indgen(10)*12+t+721
;  model_n2o_seas[y,z,t] = median(reform(n2o_norm_d1[y,z,ti,*]))
;endfor

restore,dir+'Models/WACCM/WACCM_FWSD_N2O_on_age.sav'
nyw = n_elements(n2o_norm_aoa_grid[*,0,0])

;  Read in the CLaMS output.
restore,dir+'Models/CLaMS/CLaMS_SD_UVT_ERA5_tracers.sav'
lat_c = lat
restore,dir+'Models/CLaMS/CLaMS_ERA5_N2O_on_mean_ages.sav'
nyc = n_elements(clams_lat)


colors = ['sky blue','Blue','dark blue','Purple','Medium Purple','Magenta','Lime Green','Green','gold','Orange','orange red','Red','Brown','Grey','Dark Grey','dodger blue','cornflower', $
  'Purple','Medium Purple','Magenta']

LOADCT, 39, RGB_TABLE = rgb
LOADCT, 70, RGB_TABLE = rgb_70

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb']
seas = ['DJF','MAM','JJA','SON']

a = 6
ai = 6
z = 32
z2 = 28
z3 = z2/2
m = 1
s = 3
tt = 802
ni = 12

;p = plot(indgen(2),/nodata)
;p = plot(lat[86033:87865],alt[86033:87865],symbol='o',sym_size=0.5,color='blue',/overplot)

;p = plot(indgen(2),/nodata)
;p = plot(n2o_norm[86033:87865],alt[86033:87865],symbol='o',sym_size=0.5,color='blue',linestyle=6,/overplot)

;p = plot(indgen(2),/nodata,yrange=[7,23],xtitle='N$_2$O',ytitle='Altitude (km)',font_size=11,title='N$_2$O '+months[m])
;p = plot(swoosh_n2o_seas[14,*,m],alt_s,symbol='o',/sym_filled,color='blue',sym_size=1,/overplot)
;p = plot(n2o_w[14,*,m],alt_s,symbol='o',/sym_filled,color='sky blue',sym_size=1,/overplot)
;p = plot(combo_n2o_seas[14,*,m],alt_s,symbol='o',/sym_filled,color='magenta',sym_size=1,/overplot)

;p = plot(indgen(2),/nodata,yrange=[13,32],xrange=[1,0.2],ytitle='Altitude (km)',xtitle='Normalized N$_2$O',font_size=11,title='Normalized N$_2$O Profiles', $
;  margin=[0.12,0.09,0.04,0.08],dimensions=[450,500])
;p = errorplot(n2o_nh_avg_profile_early[0,*],alt_grid2,n2o_nh_avg_profile_early[1,*],replicate(!values.f_nan,nz2),symbol='o',/sym_filled,thick=2,color='blue', $
;  errorbar_capsize=0,/overplot)
;p = errorplot(n2o_nh_avg_profile_late[0,*],alt_grid2,n2o_nh_avg_profile_late[1,*],replicate(!values.f_nan,nz2),symbol='o',/sym_filled,thick=2,color='red', $
;  errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,yrange=[13,32],xrange=[0,6.5],xtitle='Mean Age (years)',ytitle='Altitude (km)',font_size=11,title='Mean Age Profiles', $
;  margin=[0.12,0.09,0.04,0.08],dimensions=[450,500])
;p = errorplot(age_nh_avg_profile_early[0,*],alt_grid2,age_nh_avg_profile_early[1,*],replicate(!values.f_nan,nz2),symbol='o',/sym_filled,thick=2,color='sky blue', $
;  errorbar_capsize=0,/overplot)
;p = errorplot(age_nh_avg_profile_early[2,*],alt_grid2,age_nh_avg_profile_early[3,*],replicate(!values.f_nan,nz2),symbol='o',/sym_filled,thick=2,color='blue', $
;  errorbar_capsize=0,/overplot)
;p = errorplot(age_nh_avg_profile_late[0,*],alt_grid2,age_nh_avg_profile_late[1,*],replicate(!values.f_nan,nz2),symbol='o',/sym_filled,thick=2,color='red', $
;  errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1.05,0],yrange=[1.05,0],ytitle='Normalized CH$_4$',xtitle='Normalized N$_2$O',title='CH$_4$ vs. N$_2$O',font_size=11)
;p = plot([0,1],[0,1],linestyle=2,/overplot)
;ti = where(years lt 2000)
;p = plot(n2o_nh_coarse[0,*,ti],ch4_nh_coarse[0,*,ti],color='sky blue',symbol='o',/sym_filled,sym_size=0.5,linestyle=6,/overplot)
;ti = where(years ge 2000)
;p = plot(n2o_nh_coarse[0,*,ti],ch4_nh_coarse[0,*,ti],color='blue',symbol='o',/sym_filled,sym_size=0.5,linestyle=6,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Normalized N$_2$O',xtitle='Years', $
;title='N$_2$O at '+strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+'km Altitude',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = errorplot(years,n2o_nh_adj_tseries[0,z2,*],replicate(0,nt),n2o_nh_adj_tseries[1,z2,*],color='orange',symbol='o',/sym_filled,linestyle=6,errorbar_capsize=0,/overplot)
;;for y = 28, 28 do p = errorplot(years,n2o_b_grid_adj_tseries[0,y,z2,*],replicate(0,nt),n2o_b_grid_adj_tseries[1,y,z2,*],color='sky blue',symbol='s',linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;for y = 28, 28 do p = errorplot(years,n2o_a_grid_adj_tseries[0,y,z2,*],replicate(0,nt),n2o_a_grid_adj_tseries[1,y,z2,*],color=colors[y-21],symbol='s',linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;for y = 28, 28 do p = errorplot(years,n2o_grid_adj_tseries[0,y,z2,*],replicate(0,nt),n2o_grid_adj_tseries[1,y,z2,*],color=colors[y-21],symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = plot([2015,2025],replicate(n2o_norm_combo_grid_late[0,y,z2],2),thick=2,color='blue',/overplot)
;p = errorplot(years,n2o_nh_adj_tseries_lat_adj[0,z2,*],replicate(0,nt),n2o_nh_adj_tseries_lat_adj[1,z2,*],color='red',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = plot(years,n2o_nh_coarse[0,z3,*],color='purple',symbol='o',/sym_filled,linestyle=6,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Normalized CH$_4$',xtitle='Years', $
;title='CH$_4$ at '+strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+'km Altitude',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = errorplot(years,ch4_nh_adj_tseries[0,z2,*],replicate(0,nt),ch4_nh_adj_tseries[1,z2,*],color='orange',symbol='o',/sym_filled,linestyle=6,errorbar_capsize=0, $
;  /overplot)
;p = plot(years,ch4_nh_coarse[0,z3,*],color='purple',symbol='o',/sym_filled,linestyle=6,/overplot)
;chk = where(ch4_nh_coarse[0,z3,*] lt n2o_nh_coarse[0,z3,*],nchk)
;if nchk gt 0 then for rr = 0, nchk-1 do s = symbol(years[chk[rr]],ch4_nh_coarse[0,z3,chk[rr]],'s',sym_color='magenta',sym_size=1.5,/data)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Mean Age (years)',xtitle='Years',title='Mean Age at '+strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+ $
;  'km Altitude',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;for y = 24, 30 do p = errorplot(years,age_co2_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_co2_grid_adj_tseries[1,y,z2,*],color='blue',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for y = 24, 30 do p = errorplot(years,age_sf6_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_sf6_grid_adj_tseries[1,y,z2,*],color='sky blue',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for y = 24, 30 do p = errorplot(years,age_opt_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_opt_grid_adj_tseries[1,y,z2,*],color='orange',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for y = 24, 30 do p = errorplot(years,age_co2_b_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_co2_b_grid_adj_tseries[1,y,z2,*],color='purple',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for y = 24, 30 do p = errorplot(years,age_sf6_b_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_sf6_b_grid_adj_tseries[1,y,z2,*],color='pink',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for y = 24, 30 do p = errorplot(years,age_opt_b_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_opt_b_grid_adj_tseries[1,y,z2,*],color='magenta',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Mean Age (years)',xtitle='Years',title='Mean Age at '+strmid(strcompress(string(alt_grid2[z2]),/r),0,4)+ $
;  'km Altitude',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;;for y = 20, 32 do p = errorplot(years,age_co2_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_co2_grid_adj_tseries[1,y,z2,*],color=colors[y-20],symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;for y = 20, 32 do p = errorplot(years,age_sf6_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_sf6_grid_adj_tseries[1,y,z2,*],color=colors[y-20],symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;for y = 20, 32 do p = errorplot(years,age_co2_b_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_co2_b_grid_adj_tseries[1,y,z2,*],color=colors[y-20],symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;for y = 20, 32 do p = errorplot(years,age_sf6_b_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_sf6_b_grid_adj_tseries[1,y,z2,*],color=colors[y-20],symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;for y = 26, 34 do p = errorplot(years,age_grid_adj_tseries[0,y,z2,*],replicate(0,nt),age_grid_adj_tseries[1,y,z2,*],color=colors[y-20],symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for y = 26, 34 do p = plot([1990,2028],age_grid_trends[y,z2,0]+age_grid_trends[y,z2,1]*[1990,2028],color=colors[y-20],thick=2,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Mean Age (years)',xtitle='Years',title='30$\deg$N-60$\deg$N Average Mean Age at '+ $
;  strmid(strcompress(string(alt_grid3[z3]),/r),0,2)+'km Altitude',yrange=[2,5],font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;;p = errorplot(years,age_nh_adj_tseries[2,z2,*],replicate(0,nt),age_nh_adj_tseries[3,z2,*],color='orange',symbol='o',linestyle=6,errorbar_capsize=0,/overplot)
;;chk = where(age_nh_adj_tseries[5,z2,*] lt 2,nchk)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries[2,z2,chk],replicate(0,nchk),age_nh_adj_tseries[3,z2,chk],color='blue',symbol='o',linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;if nchk ge 1 then p = plot(years[chk],age_nh_adj_tseries[2,z2,chk],color='blue',symbol='o',linestyle=6,/overplot)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries_lat_adj[2,z2,chk],replicate(0,nchk),age_nh_adj_tseries_lat_adj[3,z2,chk],color='blue',symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;chk = where(age_nh_coarse[5,z3,*] eq 1,nchk)
;if nchk ge 1 then p = plot(years[chk],age_nh_coarse[2,z3,chk],color='blue',symbol='o',linestyle=6,/overplot)
;if nchk ge 1 then p = errorplot(years[chk],age_nh_coarse_lat_adj[2,z3,chk],replicate(0,nchk),age_nh_coarse_lat_adj[3,z3,chk],color='blue',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;;chk = where(age_nh_adj_tseries[5,z2,*] eq 2,nchk)
;chk = where(age_nh_coarse[5,z3,*] eq 2,nchk)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries[2,z2,chk],replicate(0,nchk),age_nh_adj_tseries[3,z2,chk],color='purple',symbol='o',linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;if nchk ge 1 then p = plot(years[chk],age_nh_adj_tseries[2,z2,chk],color='purple',symbol='o',linestyle=6,/overplot)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries_lat_adj[2,z2,chk],replicate(0,nchk),age_nh_adj_tseries_lat_adj[3,z2,chk],color='purple',symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;if nchk ge 1 then p = plot(years[chk],age_nh_coarse[2,z3,chk],color='purple',symbol='tu',linestyle=6,/overplot)
;if nchk ge 1 then p = errorplot(years[chk],age_nh_coarse_lat_adj[2,z3,chk],replicate(0,nchk),age_nh_coarse_lat_adj[3,z3,chk],color='purple',symbol='tu',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;;chk = where(age_nh_adj_tseries[5,z2,*] eq 3,nchk)
;chk = where(age_nh_coarse[5,z3,*] eq 3,nchk)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries[2,z2,chk],replicate(0,nchk),age_nh_adj_tseries[3,z2,chk],color='magenta',symbol='o',linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;if nchk ge 1 then p = plot(years[chk],age_nh_adj_tseries[2,z2,chk],color='magenta',symbol='o',linestyle=6,/overplot)
;;p = errorplot(years[chk],age_nh_adj_tseries_lat_adj[2,z2,chk],replicate(0,nchk),age_nh_adj_tseries_lat_adj[3,z2,chk],color='magenta',symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;if nchk ge 1 then p = plot(years[chk],age_nh_coarse[2,z3,chk],color='magenta',symbol='o',linestyle=6,/overplot)
;if nchk ge 1 then p = errorplot(years[chk],age_nh_coarse_lat_adj[2,z3,chk],replicate(0,nchk),age_nh_coarse_lat_adj[3,z3,chk],color='magenta',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;;chk = where(nh_adj_tseries_seas[z2,*] eq 0,nchk)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries[0,z2,chk],replicate(0,nchk),age_nh_adj_tseries[1,z2,chk],color=colors[0],symbol='o',/sym_filled, $
;;  linestyle=6,errorbar_capsize=0,/overplot)
;;chk = where(nh_adj_tseries_seas[z2,*] eq 1,nchk)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries[0,z2,chk],replicate(0,nchk),age_nh_adj_tseries[1,z2,chk],color=colors[1],symbol='o',/sym_filled, $
;;  linestyle=6,errorbar_capsize=0,/overplot)
;;chk = where(nh_adj_tseries_seas[z2,*] eq 2,nchk)
;;if nchk ge 1 then p = errorplot(years[chk],age_nh_adj_tseries[0,z2,chk],replicate(0,nchk),age_nh_adj_tseries[1,z2,chk],color=colors[2],symbol='o',/sym_filled, $
;;  linestyle=6,errorbar_capsize=0,/overplot)
;p = plot([1990,2026],age_nh_trends_coarse[z3,1]*[1990,2026] + age_nh_trends_coarse[z3,0],color='orange',thick=2,linestyle=1,/overplot)
;p = plot([1990,2026],age_nh_lat_adj_trends_coarse[z3,1]*[1990,2026] + age_nh_lat_adj_trends_coarse[z3,0],thick=2,color='orange',linestyle=2,/overplot)
;t = text(0.35,0.24,'In situ balloon',color='purple',font_size=11,/norm)
;t = text(0.35,0.2,'In situ aircraft',color='blue',font_size=11,/norm)
;t = text(0.35,0.16,'AirCore',color='magenta',font_size=11,/norm)
;s = symbol(0.35,0.835,'o',/norm)
;s = symbol(0.35,0.795,'o',/sym_filled,/norm)
;t = text(0.37,0.82,'= Orig',font_size=11,/norm)
;t = text(0.37,0.78,'= Lat Adj',font_size=11,/norm)
;p = plot([2007,2009],[4.7,4.7],color='orange',linestyle=1,thick=2,/overplot)
;p = plot([2007,2009],[4.55,4.55],color='orange',linestyle=2,thick=2,/overplot)
;t = text(0.54,0.82,' = Fit ('+strmid(strcompress(string(10.*age_nh_trends_coarse[z3,1]),/r),0,5)+'yr/dec, RMSE='+ $
;  strmid(strcompress(string(rse_orig[z3]),/r),0,4)+'yr)',font_size=11,/norm)
;t = text(0.54,0.78,' = Fit ('+strmid(strcompress(string(10.*age_nh_lat_adj_trends_coarse[z3,1]),/r),0,5)+'yr/dec, RMSE='+ $
;  strmid(strcompress(string(rse_lat_adj[z3]),/r),0,4)+'yr)',font_size=11,/norm)
;p.save,dir+'Plots/Balloon Mean Ages/Mean_age_NH_avg_lat_avg_'+strmid(strcompress(string(alt_grid3[z3]),/r),0,2)+'km.png'

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Mean Age (years)',xtitle='Years',title='Mean Age at '+strmid(strcompress(string(norm_grid[ni]),/r),0,4)+ $
;  ' N$_2$O',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;;p = errorplot(years,age_on_n2o_nh_adj_tseries[0,ni,*],replicate(0,nt),age_on_n2o_nh_adj_tseries[1,ni,*],color='red',symbol='o',/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;p = errorplot(years,age_on_n2o_adj_tseries[2,ni,*],replicate(0,nt),age_on_n2o_adj_tseries[3,ni,*],color='blue',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;chk = where(age_on_n2o_nh_adj_tseries_seas[ni,*] eq 0,nchk)
;if nchk ge 1 then p = errorplot(years[chk],age_on_n2o_nh_adj_tseries[0,ni,chk],replicate(0,nchk),age_on_n2o_nh_adj_tseries[1,ni,chk],color=colors[0],symbol='o',/sym_filled, $
;  linestyle=6,errorbar_capsize=0,/overplot)
;chk = where(age_on_n2o_nh_adj_tseries_seas[ni,*] eq 1,nchk)
;if nchk ge 1 then p = errorplot(years[chk],age_on_n2o_nh_adj_tseries[0,ni,chk],replicate(0,nchk),age_on_n2o_nh_adj_tseries[1,ni,chk],color=colors[1],symbol='o',/sym_filled, $
;  linestyle=6,errorbar_capsize=0,/overplot)
;chk = where(age_on_n2o_nh_adj_tseries_seas[ni,*] eq 2,nchk)
;if nchk ge 1 then p = errorplot(years[chk],age_on_n2o_nh_adj_tseries[0,ni,chk],replicate(0,nchk),age_on_n2o_nh_adj_tseries[1,ni,chk],color=colors[2],symbol='o',/sym_filled, $
;  linestyle=6,errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,yrange=[220,285],ytitle='N$_2$O',xtitle='Years',title='N$_2$O at '+strmid(strcompress(string(age_grid1[a]),/r),0,3)+ $
;  ' years mean age',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = plot(clams_yrs,clams_n2o_aoa_grid[3*nyc/4,ai,*],color='lime green',thick=2,/overplot)
;p = plot(years_waccm,n2o_surface_mean*n2o_norm_aoa_grid[3*nyw/4,ai,*],color='blue',thick=2,/overplot)
;p = plot(years,n2o_on_age_both_adj_tseries[0,a,*],color='purple',symbol='o',linestyle=6,/overplot)
;si = where(samp_type_b_yrs eq 2,nsi)
;p = plot(years[si],n2o_on_age_both_adj_tseries[0,a,si],color='magenta',symbol='o',linestyle=6,/overplot)
;si = where(samp_type_b_yrs eq 1,nsi)
;p = plot(years[si],n2o_on_age_both_adj_tseries[0,a,si],color='magenta',symbol='o',linestyle=6,/overplot)
;s = symbol(1999,240,'o',/sym_filled,sym_color='red',sym_size=2,/data)
;s = symbol(1997,n2o_3yr_old,'s',/sym_filled,sym_color='purple',sym_size=2,/data)
;s = symbol(2022,n2o_3yr_new,'s',/sym_filled,sym_color='purple',sym_size=2,/data)
;t = text(0.4,0.83,'In situ',color='purple',font_size=11,/norm)
;;t = text(0.4,0.79,'In situ aircraft',color='purple',font_size=11,/norm)
;t = text(0.4,0.79,'WACCM/MERRA2',color='blue',font_size=10,/norm)
;t = text(0.4,0.75,'CLaMS/ERA5',color='lime green',font_size=10,/norm)
;t = text(0.4,0.71,'Engel 2002',color='red',font_size=10,/norm)
;p.save,dir+'Plots/Balloon Mean Ages/N2O_on_3yr_mean_age_time_series.png'

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Number of samples',xtitle='Years',title='Number of measurements at '+ $
;  strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+'km Altitude',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = plot(years,age_nh_adj_tseries[4,z2,*],color='blue',symbol='o',/sym_filled,linestyle=6,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,ytitle='Mean Age Adjustment (years)',xtitle='Years',title='Mean Age Adjustment at '+ $
;  strmid(strcompress(string(alt_grid3[z3]),/r),0,2)+'km Altitude',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = plot([1990,2028],[0,0],linestyle=2,/overplot)
;;p = plot(years,age_nh_adj_tseries_lat_adj[2,z2,*]-age_nh_adj_tseries[2,z2,*],color='green',symbol='o',/sym_filled,linestyle=6,/overplot)
;p = plot(years,age_nh_coarse_lat_adj[2,z3,*]-age_nh_coarse[2,z3,*],color='green',symbol='o',/sym_filled,linestyle=6,/overplot)
;p.save,dir+'Plots/Balloon Mean Ages/Mean_age_NH_avg_age_adj_'+strmid(strcompress(string(alt_grid3[z3]),/r),0,2)+'km.png'

;p = plot(indgen(2),/nodata,xrange=[-0.35,0.15],yrange=[13.5,31],ytitle='Altitude (km)',xtitle='Years/Decade',title='NH Mean Age Trends',font_size=11,margin=[0.13,0.1,0.08,0.08],dimensions=[450,500])
;p = plot([0,0],[13,35],linestyle=2,/overplot)
;p = errorplot(10.*age_nh_trends[*,1],alt_grid2,10.*age_nh_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1,/sym_filled,color='sky blue', $
;  linestyle=6,errorbar_capsize=0,/overplot)
;p = errorplot(10.*age_nh_lat_adj_trends[*,1],alt_grid2,10.*age_nh_lat_adj_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1,/sym_filled, $
;  color='blue',linestyle=6,errorbar_capsize=0,/overplot)
;p = errorplot(10.*age_nh_lat_adj_trends_coarse[*,1],alt_grid3,10.*age_nh_lat_adj_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz3),symbol='o',sym_size=1,/sym_filled, $
;  color='purple',linestyle=6,errorbar_capsize=0,/overplot)
;p = errorplot([-0.105,-0.055],[15.5,18.5],[0.05,0.05],replicate(!values.f_nan,2),symbol='D',sym_size=1.5,/sym_filled,color='red',linestyle=6,errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,xrange=[-0.4,0.15],yrange=[1.02,0.2],ytitle='Normalized N$_2$O',xtitle='Years/Decade',title='NH Mean Age Trends on N$_2$O',font_size=11, $
;  margin=[0.13,0.1,0.08,0.08],dimensions=[450,500])
;p = plot([0,0],[1,0],linestyle=2,/overplot)
;p = errorplot(10.*age_on_n2o_nh_lat_adj_trends[*,1],norm_grid,10.*age_on_n2o_nh_lat_adj_trends_sigma[*,1],replicate(!values.f_nan,nn),symbol='o',sym_size=1.5,/sym_filled, $
;  color='blue',linestyle=6,errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,yrange=[13.5,31],xrange=[-0.04,0.04],ytitle='Altitude (km)',xtitle='Normalized Tracer Fraction/Decade',title='NH N$_2$O Trend Profiles',font_size=11, $
;  margin=[0.12,0.09,0.04,0.08],dimensions=[450,500])
;p = plot([0,0],[15,36],linestyle=2,/overplot)
;;poly = polygon([10*n2o_norm_trends_post2000_minmax[*,0],reverse(10*n2o_norm_trends_post2000_minmax[*,1])],[alt_c,reverse(alt_c)],/data,/fill_background,fill_color='lime green', $
;;  fill_transparency=70,color='white')
;;p = plot(10.*n2o_norm_trends_post2000_avg[0,2,*],alt_c,color='green',thick=2,/overplot)
;;poly = polygon([10*n2o_norm_trends_pre2000_minmax[*,0],reverse(10*n2o_norm_trends_pre2000_minmax[*,1])],[alt_c,reverse(alt_c)],/data,/fill_background,fill_color='violet', $
;;  fill_transparency=70,color='white')
;;;for m = 0, nm-1 do p = plot(10.*n2o_norm_trends_pre2000[0,2,*,m],alt_c,color='pink',thick=1,/overplot)
;;p = plot(10.*n2o_norm_trends_pre2000_avg[0,2,*],alt_c,color='medium orchid',thick=2,/overplot)
;p = errorplot(10.*ace_n2o_nh_trends[1,*],ace_alts,10.*ace_n2o_nh_trends_chi,replicate(!values.f_nan,nza),symbol='D',sym_size=1,/sym_filled,color='red',errorbar_capsize=0,linestyle=6,/overplot)
;;for z = 0, nz2-1 do s = symbol(10.*n2o_alt_trends[z,1],alt_grid2[z],'o',sym_size=n2o_nofl_alt_trends_count[z]/1.5e2+0.5,sym_color='blue',/sym_filled,/data)
;;p = errorplot(10.*n2o_alt_trends[*,1],alt_grid2,10.*n2o_alt_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=0.1,/sym_filled,color='blue',linestyle=6,errorbar_capsize=0, $
;;  /overplot)
;;p = errorplot(10.*n2o_nh_trends[*,1],alt_grid2,10.*n2o_nh_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1,/sym_filled,color='sky blue',linestyle=6,errorbar_capsize=0, $
;;  /overplot)
;;p = errorplot(10.*n2o_nh_lat_adj_trends[*,1],alt_grid2,10.*n2o_nh_lat_adj_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1,/sym_filled, $
;;  color='blue',linestyle=6,errorbar_capsize=0,/overplot)
;p = errorplot(10.*n2o_nh_lat_adj_trends_coarse[*,1],alt_grid3,10.*n2o_nh_lat_adj_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz3),symbol='o',sym_size=1, $
;  /sym_filled,color='purple',linestyle=6,errorbar_capsize=0,/overplot)
;p = errorplot(10.*n2o_nh_alt_trends_combo[0,*],ace_alts,10.*n2o_nh_alt_trends_combo[1,*],replicate(!values.f_nan,nza),symbol='o',sym_size=1,/sym_filled,color='sky blue',linestyle=6, $
;  errorbar_capsize=0,/overplot)
;for z = 0, 18 do s = symbol(10.*n2o_fl_alt_trends[z,1],alt_grid2[z],'o',sym_size=n2o_wfl_alt_trends_count[z]/1.5e2+0.5,sym_color='blue',/data)
;p = errorplot(10.*n2o_fl_alt_trends[0:18,1],alt_grid2[0:18],10.*n2o_fl_alt_trends_sigma[0:18,1],replicate(!values.f_nan,19),symbol='o',sym_size=0.1,color='blue',linestyle=6,errorbar_capsize=0, $
;  /overplot)
;t = text(0.17,0.85,'In situ',color='blue',font_size=11,/norm)
;;t = text(0.65,0.31,'w/flasks = ',color='blue',font_size=10,/norm)
;;s = symbol(0.83,0.32,'o',sym_size=1.5,sym_color='blue',/norm)
;;t = text(0.65,0.27,'w/o flasks = ',color='blue',font_size=10,/norm)
;;s = symbol(0.85,0.28,'o',sym_size=1.5,/sym_filled,sym_color='blue',/norm)
;t = text(0.17,0.81,'ACE (2004-25)',color='red',font_size=11,/norm)
;t = text(0.5,0.85,'CCMI-2022',color='green',font_size=11,/norm)
;t = text(0.65,0.52,'CCMI (1995-2019)',color='green',font_size=10,/norm)
;t = text(0.65,0.48,'CCMI (1975-1995)',color='medium orchid',font_size=10,/norm)

;p = plot(indgen(2),/nodata,yrange=[13.5,31],xrange=[-0.04,0.04],ytitle='Altitude (km)',xtitle='Normalized Tracer Fraction/Decade',title='NH CH$_4$ Trend Profiles', $
;  font_size=11,margin=[0.12,0.09,0.04,0.08],dimensions=[450,500])
;p = plot([0,0],[15,36],linestyle=2,/overplot)
;p = errorplot(10.*ch4_nh_trends[*,1],alt_grid2,10.*ch4_nh_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1,/sym_filled,color='sky blue',linestyle=6,errorbar_capsize=0, $
;  /overplot)
;p = errorplot(10.*ch4_nh_trends_coarse[*,1],alt_grid3,10.*ch4_nh_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz3),symbol='o',sym_size=1, $
;  /sym_filled,color='purple',linestyle=6,errorbar_capsize=0,/overplot)

;p = plot(indgen(2),/nodata,xrange=[1990,2028],xtickinterval=5,yrange=[27,63],ytitle='Avg sampled latitude',xtitle='Years',title='Time series of average sampled latitudes '+ $
;  strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+'km',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = plot([1990,2028],[45,45],linestyle=2,/overplot)
;;p = plot(years,lat_avg_nh_adj_tseries[z2,*],color='green',symbol='o',/sym_filled,linestyle=6,/overplot)
;p = plot(years,lat_avg_nh_coarse[z3,*],color='green',symbol='o',/sym_filled,linestyle=6,/overplot)
;;p = plot(years,lat_avg_n2o_nh_adj_tseries[z2,*],color='blue',symbol='o',/sym_filled,linestyle=6,/overplot)
;p.save,dir+'Plots/Balloon Mean Ages/Mean_age_NH_avg_lat_samp_'+strmid(strcompress(string(alt_grid3[z3]),/r),0,2)+'km.png'

;p = plot(indgen(2),/nodata,ytitle='years',xrange=[20,70],xminor=1,yrange=[-1.2,1.2],xtitle='Latitude',title='Climatological Seasonal Mean Age Relative to 45$\deg$N  Alt='+ $
;  strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+'km',font_size=11,margin=[0.1,0.1,0.03,0.08],dimensions=[700,400])
;p = plot([20,70],[0,0],linestyle=2,/overplot)
;p = plot([45,45],[-1.5,1.5],linestyle=2,/overplot)
;p = plot(lat_grid_a,age_grid_seas_avg[0,*,z2,0]-age_grid_seas_avg[0,26,z2,0],thick=3,color='blue',/overplot)
;p = plot(lat_grid_a,age_grid_seas_avg[0,*,z2,1]-age_grid_seas_avg[0,26,z2,1],thick=3,color='green',/overplot)
;p = plot(lat_grid_a,age_grid_seas_avg[0,*,z2,2]-age_grid_seas_avg[0,26,z2,2],thick=3,color='orange',/overplot)
;p = plot(lat_grid_a,age_grid_seas_avg[0,*,z2,3]-age_grid_seas_avg[0,26,z2,3],thick=3,color='red',/overplot)
;t = text(0.17,0.81,'DJF',color='blue',font_size=11,/norm)
;t = text(0.17,0.77,'MAM',color='green',font_size=11,/norm)
;t = text(0.17,0.73,'JJA',color='orange',font_size=11,/norm)
;t = text(0.17,0.69,'SON',color='red',font_size=11,/norm)
;p.save,dir+'Plots/Balloon Mean Ages/Mean_age_climatology_lat_gradient_20km.png'

;p = plot(indgen(2),/nodata,xrange=[0,1],ytitle='Mean Age (years)',title='NH Mean Age Seasonal Cycle '+strmid(strcompress(string(alt_grid2[z2]),/r),0,2)+ $
;  'km Altitude',font_size=11,margin=[0.1,0.1,0.04,0.08],dimensions=[700,400])
;p = plot([0,1],[0,0],linestyle=2,/overplot)
;for y = 26, 32 do p = plot(seas_frac,mean_age_grid_seas[y,z2,*],thick=2,symbol='o',/sym_filled,color='purple',/overplot)

gd = where((months_all[ico2] eq m or months_all[ico2] eq m+1 or months_all[ico2] eq m+2) and finite(lat[ico2]) and alt[ico2] ge alt_grid[z]-dz and alt[ico2] lt alt_grid[z]+dz)
gd_a0 = where((months_all_a0[ico2_a0] eq m or months_all_a0[ico2_a0] eq m+1 or months_all_a0[ico2_a0] eq m+2) and finite(lat_a0[ico2_a0]) and $
  alt_a0[ico2_a0] ge alt_grid[z]-dz and alt_a0[ico2_a0] lt alt_grid[z]+dz)

;p = plot(indgen(2),/nodata,yrange=[1.03,0.75],ytitle='Normalized N$_2$O',xrange=[0,90],xtitle='Latitude',title='AirCore N$_2$O  Alt='+ $
;  strmid(strcompress(string(alt_grid[z]),/r),0,5)+'km',font_size=11,margin=[0.11,0.1,0.04,0.08],dimensions=[700,400])
;;p = plot(lat[ico2[gd]],n2o_norm[ico2[gd]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m[ico2[gd]],n2o_norm[ico2[gd]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;p = plot(lat_a0[in2o_a0[gd_a0]],n2o_norm_a0[ico2_a0[gd_a0]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;p = plot(elat_m_a0[ico2_a0[gd_a0]],n2o_norm_a0[ico2_a0[gd_a0]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;p = plot(elat_m_adj_a0[ico2_a0[gd_a0]],n2o_norm_a[ico2_a[gd_a]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;;for t = 0, 9 do p = plot(slat,combn2oq_norm[*,7,t*12+m+360],color='lime green',thick=1,/overplot)
;p = plot(slat,combo_n2o_seas_hires[*,z,m],color='green',thick=2,/overplot)
;t = text(0.15,0.83,'Lat',font_size=11,color='sky blue',/norm)
;t = text(0.15,0.78,'Elat from MERRA2',font_size=11,color='blue',/norm)
;t = text(0.15,0.73,'SWOOSH seas avg',font_size=11,color='green',/norm)
;p.save,dir+'Plots/Balloon Mean Ages/AirCore_elat_example1.png'
;p.save,dir+'Plots/Balloon Mean Ages/AirCore_elat_example2.png'

;p = plot(indgen(2),/nodata,yrange=[1.03,0.6],ytitle='Normalized N$_2$O',xrange=[0,90],xtitle='Latitude',title='AirCore N$_2$O  Alt='+ $
;  strmid(strcompress(string(alt_grid[z]),/r),0,5)+'km',font_size=11,margin=[0.11,0.1,0.04,0.08],dimensions=[700,400])
;p = plot(slat,combo_n2o_seas_hires[*,z,m],color='green',thick=2,/overplot)
;p = errorplot(lat_grid_a,n2o_grid_adj_tseries[0,*,z/2,tt],replicate(0,nyg),n2o_grid_adj_tseries[1,*,z/2,tt],color='blue',symbol='o',/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)

gd = where(months_all[ico2] eq m+1 and finite(lat[ico2]) and alt[ico2] ge alt_grid[z]-dz/2. and alt[ico2] lt alt_grid[z]+dz/2.)
gd_a = where(months_all_a[ico2_a] eq m+1 and finite(lat_a[ico2_a]) and alt_a[ico2_a] ge alt_grid[z]-dz/2. and alt_a[ico2_a] lt alt_grid[z]+dz/2.)
gd_a0 = where(months_all_a0[ico2_a0] eq m+1 and finite(lat_a0[ico2_a0]) and alt_a0[ico2_a0] ge alt_grid[z]-dz/2. and alt_a0[ico2_a0] lt alt_grid[z]+dz/2.)
gd_an0 = where(months_all_a0[in2o_a0] eq m+1 and finite(n2o_norm_a0[in2o_a0]) and finite(lat_a0[in2o_a0]) and alt_a0[in2o_a0] ge alt_grid[z]-dz/2. and $
  alt_a0[in2o_a0] lt alt_grid[z]+dz/2.)
gd_as = where(months_all_a[isf6_a] eq m+1 and finite(lat_a[isf6_a]) and alt_a[isf6_a] ge alt_grid[z]-dz/2. and alt_a[isf6_a] lt alt_grid[z]+dz/2.)
gd_b = where(months_all_b[ico2_b] eq m+1 and finite(lat_b[ico2_b]) and alt_b[ico2_b] ge alt_grid[z]-dz/2. and alt_b[ico2_b] lt alt_grid[z]+dz/2.,ngdb)
gd_b1 = where(months_all_b[ico2_b] eq m+1)
;print,ngdb,times_all_b[ico2_b[gd_b[0:2]]],elat_m_adj_b[ico2_b[gd_b]]
;print,months_all_b[0:10]

;p = plot(indgen(2),/nodata,xrange=[1.05,0],yrange=[9,32],xtitle='N$_2$O',ytitle='Altitude (km)',font_size=11,title='N$_2$O '+months[m])
;p = plot(n2o_norm_b[ico2_b[gd_b1]],alt_b[ico2_b[gd_b1]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)

;p = plot(indgen(2),/nodata,yrange=[1.03,0.75],ytitle='Normalized N$_2$O',xrange=[-90,90],xtitle='Latitude',title='N$_2$O  Alt='+strmid(strcompress(string(alt_grid[z]),/r),0,5)+'km', $
;  font_size=11,margin=[0.11,0.1,0.04,0.08],dimensions=[700,400])
;;;for i = 0, 25 do p = plot(lat_c[45:-1],n2o_norm_d1[45:-1,15,721+m,i],symbol='o',color='purple',sym_size=0.5,linestyle=6,/overplot)
;;;p = plot(lat_c[45:-1],model_n2o_seas[45:-1,15,m],color='purple',thick=2,/overplot)
;;;p = plot(wlat[45:-1],n2o_norm_seas[0,45:-1,23,m],color='magenta',thick=2,/overplot)
;;for t = 8, 8 do p = plot(slat,combn2oq_norm[*,7,t*12+m+360],color='lime green',thick=1,/overplot)
;;p = plot(slat,combo_n2o_s_hires[*,z,462],color='orange',thick=1,/overplot)
;;p = plot(slat,swoosh_n2o_seas[*,6,m],color='green',thick=2,/overplot)
;;p = plot(slat,swoosh_n2o_seas_early[*,6,m],color='dodger blue',thick=2,/overplot)
;;;p = plot(slat,combo_n2o_seas[*,6,m],color='green',thick=2,/overplot)
;;p = plot(slat,combo_n2o_seas_hires[*,z,m],color='green',thick=2,/overplot)
;p = plot(slat,combo_n2o_seas_early_hires[*,z,m],color='green',thick=2,/overplot)
;;p = plot(slat,combo_n2o_seas_early_hires[*,z-1,m],color='lime green',thick=2,/overplot)
;;for mi = 0, 11 do p = plot(slat,combo_n2o_seas_early_hires[*,z,mi],color=colors[mi],thick=1,/overplot)
;p = plot(elat_m[ico2[gd]],n2o_ch4_norm[ico2[gd]],symbol='o',/sym_filled,color='sky blue',sym_size=0.5,linestyle=6,/overplot)
;p = plot(lat[ico2[gd]],n2o_norm[ico2[gd]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;p = plot(elat_m_adj[ico2[gd]],n2o_norm[ico2[gd]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;;;p = plot(elat_m_adj_a[isf6_a[gd_as]],n2o_norm_a[isf6_a[gd_as]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_a[ico2_a[gd_a]],n2o_norm_a[ico2_a[gd_a]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;;;p = plot(elat_m_a[isf6_a[gd_as]],n2o_norm_a[isf6_a[gd_as]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(lat_a[ico2_a[gd_a]],n2o_norm_a[ico2_a[gd_a]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;;;p = plot(lat_a[isf6_a[gd_as]],n2o_norm_a[isf6_a[gd_as]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_adj_a[ico2_a[gd_a]],n2o_norm_a[ico2_a[gd_a]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(lat_b[ico2_b[gd_b]],n2o_norm_b[ico2_b[gd_b]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_b[ico2_b[gd_b]],n2o_norm_b[ico2_b[gd_b]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_adj_b[ico2_b[gd_b]],n2o_norm_b[ico2_b[gd_b]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(lat_a0[ico2_a0[gd_a0]],n2o_norm_a0[ico2_a0[gd_a0]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_a0[ico2_a0[gd_a0]],n2o_norm_a0[ico2_a0[gd_a0]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_adj_a0[ico2_a0[gd_a0]],n2o_norm_a0[ico2_a0[gd_a0]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(lat_a0[in2o_a0[gd_an0]],n2o_norm_a0[in2o_a0[gd_an0]],symbol='o',/sym_filled,color='sky blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_a0[in2o_a0[gd_an0]],n2o_norm_a0[in2o_a0[gd_an0]],symbol='o',/sym_filled,color='blue',sym_size=0.75,linestyle=6,/overplot)
;;p = plot(elat_m_adj_a_early[in2o_a0[gd_an0]],n2o_norm_a0[in2o_a0[gd_an0]],symbol='o',/sym_filled,color='red',sym_size=0.75,linestyle=6,/overplot)
;t = text(0.15,0.83,'Lat',font_size=11,color='sky blue',/norm)
;t = text(0.15,0.78,'Elat from MERRA2',font_size=11,color='blue',/norm)
;t = text(0.15,0.73,'SWOOSH seas avg',font_size=11,color='green',/norm)
;t = text(0.15,0.68,'Elat from SWOOSH',font_size=11,color='red',/norm)
;p.save,dir+'Plots/Balloon Mean Ages/AirCore_elat_example3.png'

;p = plot(indgen(2),/nodata,yrange=[0,4],ytitle='years',xrange=[0,90],xtitle='Equivalent Latitude',title='Equivalent latitude vs. Mean Age  Alt='+ $
;  strmid(strcompress(string(alt_grid[z]),/r),0,4)+'km',font_size=11,margin=[0.11,0.1,0.04,0.08],dimensions=[700,400])
;p = plot(elat_m[ico2[gd]],dat_co2_b2.age_opt_all[gd],symbol='o',/sym_filled,color='blue',sym_size=0.5,linestyle=6,/overplot)
;p = plot(elat_m_adj[ico2[gd]],dat_co2_b2.age_opt_all[gd],symbol='o',/sym_filled,color='sky blue',sym_size=0.5,linestyle=6,/overplot)
;
;p = plot(indgen(2),/nodata,yrange=[1,0.65],ytitle='Normalized N$_2$O',xrange=[0,5],xtitle='years',title='Mean Age vs. N$_2$O  Alt='+ $
;  strmid(strcompress(string(alt_grid[z]),/r),0,4)+'km',font_size=11,margin=[0.11,0.1,0.04,0.08],dimensions=[700,400])
;p = plot(dat_co2_b2.age_opt_all[gd],n2o_norm[ico2[gd]],symbol='o',/sym_filled,color='blue',sym_size=0.5,linestyle=6,/overplot)

;p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[0,7.5],ytitle='Normalized N$_2$O',xtitle='Mean Age (years)',title='Mean Age vs. N$_2$O',font_size=11, $
;  position=[0.1,0.1,0.97,0.9],dimensions=[700,500])
;;c = contour(mean_age_co2_on_n2o_a0_hist,age_grid,norm_grid2,c_value=[0.1,0.1],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;c = contour(mean_age_sf6_on_n2o_a0_hist,age_grid,norm_grid2,c_value=[0.1,0.1],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;c = contour(mean_age_co2_on_n2o_b0_hist,age_grid,norm_grid2,c_value=[0.1,0.1],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;c = contour(mean_age_sf6_on_n2o_b0_hist,age_grid,norm_grid2,c_value=[0.1,0.1],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;c = contour(mean_age_on_n2o_b0_hist,age_grid,norm_grid2,c_value=[0.1,0.1],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;c = contour(mean_age_sf6_opt_on_n2o_b0_hist,age_grid,norm_grid2,c_value=[0.1,0.1],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;c = contour(mean_age_all_on_n2o_hist,age_grid,norm_grid2,c_value=[0.03,0.03],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;p = errorplot(mean_age_co2_on_n2o_a0[*,0],norm_grid2,mean_age_co2_on_n2o_a0[*,1],replicate(0,nn2),symbol='s',color='blue',sym_size=1,/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_on_n2o_a0[*,0],norm_grid2,mean_age_sf6_on_n2o_a0[*,1],replicate(0,nn2),symbol='s',color='sky blue',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_n2o_0[*,2],norm_grid2,mean_age_all_on_n2o_0[*,3],replicate(0,nn2),symbol='s',color='magenta',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_n2o_a0[*,2],norm_grid2,mean_age_all_on_n2o_a0[*,3],replicate(0,nn2),symbol='s',color='blue',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_co2_on_n2o_b0[*,0],norm_grid2,mean_age_co2_on_n2o_b0[*,1],replicate(0,nn2),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_on_n2o_b0[*,2],norm_grid2,mean_age_sf6_on_n2o_b0[*,3],replicate(0,nn2),symbol='s',color='purple',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_on_n2o_b0[*,0],norm_grid2,mean_age_on_n2o_b0[*,1],replicate(0,nn2),symbol='s',color='magenta',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_opt_on_n2o_b0[*,2],norm_grid2,mean_age_sf6_opt_on_n2o_b0[*,3],replicate(0,nn2),symbol='s',color='lime green',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_n2o_b0[*,2],norm_grid2,mean_age_all_on_n2o_b0[*,3],replicate(0,nn2),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_co2_on_n2o_b0_mlat[*,2],norm_grid2,mean_age_co2_on_n2o_b0_mlat[*,3],replicate(0,nn2),symbol='o',color='purple',sym_size=1,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_co2_on_n2o_b0_hlat[*,2],norm_grid2,mean_age_co2_on_n2o_b0_hlat[*,3],replicate(0,nn2),symbol='o',color='magenta',sym_size=1,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = plot(mean_age_all_on_n2o_hist_cent[0,*],norm_grid2,color='magenta',thick=2,/overplot)
;p = plot(mean_age_all_on_n2o_hist_cent[1,*],norm_grid2,color='magenta',thick=2,/overplot)

;p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[0,7.5],ytitle='Normalized CH$_4$',xtitle='Mean Age (years)',title='Mean Age vs. CH$_4$',font_size=11, $
;  position=[0.1,0.1,0.97,0.9],dimensions=[700,500])
;c = contour(mean_age_all_on_ch4_hist,age_grid,norm_grid2,c_value=[0.03,0.03],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;p = errorplot(mean_age_co2_on_ch4_a0[*,0],norm_grid2,mean_age_co2_on_ch4_a0[*,1],replicate(0,nn2),symbol='s',color='blue',sym_size=1,/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_on_ch4_a0[*,0],norm_grid2,mean_age_sf6_on_ch4_a0[*,1],replicate(0,nn2),symbol='s',color='sky blue',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_co2_on_ch4_b0[*,0],norm_grid2,mean_age_co2_on_ch4_b0[*,1],replicate(0,nn2),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_on_ch4_b0[*,2],norm_grid2,mean_age_sf6_on_ch4_b0[*,3],replicate(0,nn2),symbol='s',color='purple',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_on_ch4_b0[*,0],norm_grid2,mean_age_on_ch4_b0[*,1],replicate(0,nn2),symbol='s',color='magenta',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_ch4_a0[*,2],norm_grid2,mean_age_all_on_ch4_a0[*,3],replicate(0,nn2),symbol='s',color='blue',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_ch4_b0[*,2],norm_grid2,mean_age_all_on_ch4_b0[*,3],replicate(0,nn2),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = plot(mean_age_all_on_ch4_hist_cent[0,*],norm_grid2,color='magenta',thick=2,/overplot)
;p = plot(mean_age_all_on_ch4_hist_cent[1,*],norm_grid2,color='magenta',thick=2,/overplot)

;p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[0,7.5],ytitle='Normalized CFC-12',xtitle='Mean Age (years)',title='Mean Age vs. CFC-12',font_size=11, $
;  position=[0.1,0.1,0.97,0.9],dimensions=[700,500])
;c = contour(mean_age_all_on_f12_hist,age_grid,norm_grid2,c_value=[0.03,0.03],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;;p = errorplot(mean_age_co2_on_f12_a0[*,0],norm_grid2,mean_age_co2_on_f12_a0[*,1],replicate(0,nn2),symbol='s',color='blue',sym_size=1,/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_on_f12_a0[*,0],norm_grid2,mean_age_sf6_on_f12_a0[*,1],replicate(0,nn2),symbol='s',color='sky blue',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_co2_on_f12_b0[*,0],norm_grid2,mean_age_co2_on_f12_b0[*,1],replicate(0,nn2),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(mean_age_sf6_on_f12_b0[*,2],norm_grid2,mean_age_sf6_on_f12_b0[*,3],replicate(0,nn2),symbol='s',color='purple',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_f12_a0[*,2],norm_grid2,mean_age_all_on_f12_a0[*,3],replicate(0,nn2),symbol='s',color='blue',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = errorplot(mean_age_all_on_f12_b0[*,2],norm_grid2,mean_age_all_on_f12_b0[*,3],replicate(0,nn2),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;  errorbar_capsize=0,/overplot)
;p = plot(mean_age_all_on_f12_hist_cent[0,*],norm_grid2,color='magenta',thick=2,/overplot)
;p = plot(mean_age_all_on_f12_hist_cent[1,*],norm_grid2,color='magenta',thick=2,/overplot)

;p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[25,90],ytitle='Normalized N$_2$O',xtitle='Equivalent Latitude',title='Equivalent Latitude vs. N$_2$O',font_size=11, $
;  position=[0.1,0.1,0.97,0.9],dimensions=[700,500])
;;p = errorplot(elat_co2_on_n2o_a0[*,0],norm_grid2,elat_co2_on_n2o_a0[*,1],replicate(0,nn2),symbol='s',color='blue',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;;p = errorplot(elat_sf6_on_n2o_a0[*,0],norm_grid2,elat_sf6_on_n2o_a0[*,1],replicate(0,nn2),symbol='s',color='sky blue',sym_size=1,linestyle=6, $
;;  errorbar_capsize=0,/overplot)
;p = errorplot(elat_all_on_n2o_a0[*,0],norm_grid2,elat_all_on_n2o_a0[*,1],replicate(0,nn2),symbol='s',color='blue',sym_size=1,linestyle=6,/sym_filled, $
;  errorbar_capsize=0,/overplot)
;p = errorplot(elat_all_on_n2o_b0[*,0],norm_grid2,elat_all_on_n2o_b0[*,1],replicate(0,nn2),symbol='s',color='purple',sym_size=1,linestyle=6,/sym_filled, $
;  errorbar_capsize=0,/overplot)

;  Lat grid averages
nn2 = 51
dn2 = 0.02
norm_grid2 = findgen(nn2)*dn2+0.02
sss = replicate(!values.f_nan,nz2) & age_grid_late_tr = sss & n2o_grid_late_tr = sss
ttt = replicate(!values.f_nan,2,nn2) & age_grid_late_nh = sss & n2o_grid_late_nh = sss & age_n2o_late_tr = ttt & age_n2o_late_nh = ttt
age_grid_late_vx = sss & n2o_grid_late_vx = sss & age_n2o_late_sh = ttt & age_n2o_late_vx = ttt
for z = 0, nz2-1 do begin
  age_grid_late_tr[z] = mean(age_combo_grid_late[0,15:19,z],/nan)
  n2o_grid_late_tr[z] = mean(n2o_norm_combo_grid_late[0,15:19,z],/nan)
  age_grid_late_nh[z] = mean(age_combo_grid_late[0,24:28,z],/nan)
  n2o_grid_late_nh[z] = mean(n2o_norm_combo_grid_late[0,24:28,z],/nan)
  age_grid_late_vx[z] = mean(age_combo_grid_late[0,nyg-3:nyg-1,z],/nan)
  n2o_grid_late_vx[z] = mean(n2o_norm_combo_grid_late[0,nyg-3:nyg-1,z],/nan)
endfor
for n = 0, nn2-1 do begin
  gd = where(finite(age_grid_late_tr) and n2o_grid_late_tr ge norm_grid2[n]-dn2/2. and n2o_grid_late_tr lt norm_grid2[n]+dn2/2.,ngd)
  if ngd ge 1 then begin
    stats = moment(age_grid_late_tr[gd],sdev=sdev)
    age_n2o_late_tr[*,n] = [stats[0],sdev]
  endif
  gd = where(finite(age_grid_late_nh) and n2o_grid_late_nh ge norm_grid2[n]-dn2/2. and n2o_grid_late_nh lt norm_grid2[n]+dn2/2.,ngd)
  if ngd ge 1 then begin
    stats = moment(age_grid_late_nh[gd],sdev=sdev)
    age_n2o_late_nh[*,n] = [stats[0],sdev]
  endif
  gd = where(finite(age_grid_late_vx) and n2o_grid_late_vx ge norm_grid2[n]-dn2/2. and n2o_grid_late_vx lt norm_grid2[n]+dn2/2.,ngd)
  if ngd ge 1 then begin
    stats = moment(age_grid_late_vx[gd],sdev=sdev)
    age_n2o_late_vx[*,n] = [stats[0],sdev]
  endif
endfor

;p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[0,7],ytitle='Normalized N$_2$O',xtitle='Mean Age (years)',title='In Situ Mean Age vs. N$_2$O  2020s',font_size=11)
;;for y = 0, 14 do p = plot(age_combo_grid_late[0,y,*],n2o_norm_combo_grid_late[0,y,*],symbol='o',color='orange',sym_size=1,/sym_filled,linestyle=6,/overplot)
;;for y = 20, nyg-1 do p = plot(age_combo_grid_late[0,y,*],n2o_norm_combo_grid_late[0,y,*],symbol='o',color='blue',sym_size=1,/sym_filled,linestyle=6,/overplot)
;;for y = 15, 19 do p = plot(age_combo_grid_late[0,y,*],n2o_norm_combo_grid_late[0,y,*],symbol='o',color='red',sym_size=1,/sym_filled,linestyle=6,/overplot)
;;p = plot(age_grid_late_vx,n2o_grid_late_vx,symbol='o',color='purple',sym_size=1,linestyle=6,/overplot)
;;p = plot(age_grid_late_nh,n2o_grid_late_nh,symbol='o',color='blue',sym_size=1,/sym_filled,linestyle=6,/overplot)
;;p = plot(age_grid_late_tr,n2o_grid_late_tr,symbol='o',color='red',sym_size=1,linestyle=6,/overplot)
;p = errorplot(age_n2o_late_vx[0,*],norm_grid2,age_n2o_late_vx[1,*],replicate(0,nn2),symbol='o',color='lime green',sym_size=1,/sym_filled,linestyle=6,errorbar_capsize=0, $
;  /overplot)
;p = errorplot(age_n2o_late_nh[0,*],norm_grid2,age_n2o_late_nh[1,*],replicate(0,nn2),symbol='o',color='blue',sym_size=1,/sym_filled,linestyle=6,errorbar_capsize=0, $
;  /overplot)
;p = errorplot(age_n2o_late_tr[0,*],norm_grid2,age_n2o_late_tr[1,*],replicate(0,nn2),symbol='o',color='red',sym_size=1,/sym_filled,linestyle=6,errorbar_capsize=0, $
;  /overplot)
;t = text(0.18,0.8,'Tropics (15S-15N)',color='red',font_size=11,/norm)
;t = text(0.18,0.76,'NH Midlat (40-60N)',color='blue',font_size=11,/norm)
;t = text(0.18,0.72,'NH High Lat (80-90N)',color='lime green',font_size=11,/norm)
;p.save,dir+'Plots/Balloon Mean Ages/Mean_age_vs_N2O_NH_lat_bins_2020s.png'


tickname = strmid(strcompress(string(findgen(9)*0.1+0.2),/r),0,3)
tickname1 = strmid(strcompress(string(findgen(6)),/r),0,1)
tickname2 = strmid(strcompress(string(findgen(7)),/r),0,1)
tickname4 = strmid(strcompress(string(findgen(11)*0.01-0.05),/r),0,5)
tickname4[5] = '0'
tickname5 = strmid(strcompress(string(findgen(11)*0.1-0.5),/r),0,4)

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='AirCore N$_2$O '+months[m])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_grid_mon_avg[0,*,z,m]),nii)
;;  ii = where(finite(n2o_norm_grid_elat_orig_mon_avg[0,*,z,m]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_grid_mon_avg[0,ii,z,m]-0.2)*300.
;;    vcols = (n2o_norm_grid_elat_orig_mon_avg[0,ii,z,m]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;p.save,dir+'Plots/Balloon Mean Ages/AirCore_N2O_'+months[m]+'.png'
;p.save,dir+'Plots/Balloon Mean Ages/AirCore_N2O_elat_orig_'+months[m]+'.png'
;
;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Aircraft N$_2$O '+months[m])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_grid_mon_avg_a[0,*,z,m]),nii)
;;  ii = where(finite(n2o_norm_sf6_grid_mon_avg_a[0,*,z,m]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_grid_mon_avg_a[0,ii,z,m]-0.2)*300.
;;    vcols = (n2o_norm_sf6_grid_mon_avg_a[0,ii,z,m]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;
;p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Aircraft Early N$_2$O '+months[m])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_grid_mon_avg_a0[0,*,z,m]),nii)
;;  ii = where(finite(n2o_norm_sf6_grid_mon_avg_a0[0,*,z,m]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_grid_mon_avg_a0[0,ii,z,m]-0.2)*300.
;;    vcols = (n2o_norm_sf6_grid_mon_avg_a0[0,ii,z,m]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;
;p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Aircraft Early Mean Ages '+months[m])
;for z = 0, nz2-1 do begin
;;  ii = where(finite(age_co2_grid_mon_avg_a0[0,*,z,m]),nii)
;;  ii = where(finite(age_opt_grid_mon_avg_a0[0,*,z,m]),nii)
;  ii = where(finite(age_sf6_grid_mon_avg_a0[0,*,z,m]),nii)
;  if nii gt 0 then begin
;;    vcols = age_co2_grid_mon_avg_a0[0,ii,z,m]*40.
;;    vcols = age_opt_grid_mon_avg_a0[0,ii,z,m]*40.
;    vcols = age_sf6_grid_mon_avg_a0[0,ii,z,m]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')

;p = plot(indgen(2),/nodata,xrange=[-40,90],xtickinterval=20,xminor=1,yrange=[9,33],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Balloon OMS N$_2$O '+months[m])
;for z = 0, nz2-1 do begin
;;  ii = where(finite(n2o_norm_grid_mon_avg_b[0,*,z,m]),nii)
;;  ii = where(finite(n2o_norm_sf6_grid_mon_avg_b[0,*,z,m]),nii)
;  ii = where(finite(n2o_norm_both_grid_mon_avg_b[0,*,z,m]),nii)
;  if nii gt 0 then begin
;;    vcols = (n2o_norm_grid_mon_avg_b[0,ii,z,m]-0.2)*300.
;;    vcols = (n2o_norm_sf6_grid_mon_avg_b[0,ii,z,m]-0.2)*300.
;    vcols = (n2o_norm_both_grid_mon_avg_b[0,ii,z,m]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.2,/overplot)
;    p = plot(lat_grid[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.2,/overplot)
;    p = plot(lat_grid[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')

;p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Aircraft Mid N$_2$O '+months[m])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_grid_mon_avg_a_mid[0,*,z,m]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_grid_mon_avg_a_mid[0,ii,z,m]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Combined N$_2$O '+months[m])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_combo_grid_mon_avg[0,*,z,m]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_combo_grid_mon_avg[0,ii,z,m]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;
;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Mean Ages '+months[m])
;for z = 0, nz2-1 do begin
;;  ii = where(finite(age_co2_grid_mon_avg[0,*,z,m]),nii)
;  ii = where(finite(age_sf6_grid_mon_avg_a[0,*,z,m]),nii)
;  if nii gt 0 then begin
;;    vcols = age_co2_grid_mon_avg[0,ii,z,m]*40.
;    vcols = age_sf6_grid_mon_avg_a[0,ii,z,m]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')

;p = plot(indgen(2),/nodata,xrange=[-40,90],xtickinterval=20,xminor=1,yrange=[9,33],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='OMS Balloon Mean Ages '+months[m])
;for z = 0, nz2-1 do begin
;;  ii = where(finite(age_co2_grid_mon_avg_b[0,*,z,m]),nii)
;;  ii = where(finite(age_sf6_grid_mon_avg_b[0,*,z,m]),nii)
;  ii = where(finite(age_both_grid_mon_avg_b_early[0,*,z,m]),nii)
;  if nii gt 0 then begin
;;    vcols = age_co2_grid_mon_avg_b[0,ii,z,m]*40.
;;    vcols = age_sf6_grid_mon_avg_b[0,ii,z,m]*40.
;    vcols = age_both_grid_mon_avg_b_early[0,ii,z,m]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.2,/overplot)
;    p = plot(lat_grid[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.2,/overplot)
;    p = plot(lat_grid[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Combined Mean Ages '+months[m])
;for z = 0, nz2-1 do begin
;  ii = where(finite(age_combo_grid_mon_avg[*,z,m]),nii)
;  if nii gt 0 then begin
;    vcols = age_combo_grid_mon_avg[ii,z,m]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_mon_alt[*,m],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='AirCore N$_2$O '+seas[s])
;for z = 0, nz2-1 do begin
;;  ii = where(finite(n2o_norm_grid_seas[0,*,z,s]),nii)
;  ii = where(finite(n2o_norm_sf6_grid_seas[0,*,z,s]),nii)
;  if nii gt 0 then begin
;;    vcols = (n2o_norm_grid_seas[0,ii,z,s]-0.2)*300.
;    vcols = (n2o_norm_sf6_grid_seas[0,ii,z,s]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Combined N$_2$O '+seas[s])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_combo_grid_seas_late[0,*,z,s]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_combo_grid_seas_late[0,ii,z,s]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.7,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_N2O_'+seas[s]+'.png'

;p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Aircraft Early N$_2$O '+seas[s])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_grid_seas_a0[0,*,z,s]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_grid_seas_a0[0,ii,z,s]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')

;p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,23],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Aircraft Mid N$_2$O '+seas[s])
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_grid_seas_a1[0,*,z,s]),nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_grid_seas_a1[0,ii,z,s]-0.4)*416.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+0.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(7)*41.6,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Mean Ages '+seas[s])
;for z = 0, nz2-1 do begin
;;  ii = where(finite(age_co2_grid_seas[0,*,z,s]),nii)
;  ii = where(finite(age_sf6_grid_seas[0,*,z,s]),nii)
;  if nii gt 0 then begin
;;    vcols = age_co2_grid_seas[0,ii,z,s]*40.
;    vcols = age_sf6_grid_seas[0,ii,z,s]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Combined Mean Ages '+seas[s])
;for z = 0, nz2-1 do begin
;  ii = where(finite(age_combo_grid_seas[0,*,z,s]),nii)
;  if nii gt 0 then begin
;    vcols = age_combo_grid_seas[0,ii,z,s]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_mean_age_late_'+seas[s]+'.png'

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,27],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Mean Ages '+seas[s])
;for z = 0, nz2-1 do begin
;  ii = where(finite(age_grid_seas_avg[0,*,z,s]),nii)
;  if nii gt 0 then begin
;    vcols = age_grid_seas_avg[0,ii,z,s]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_seas_alt[*,s],thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname2,tickvalues=findgen(7)*40.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,27],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Normalized N$_2$O 2015-2025')
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_norm_combo_grid_late[0,*,z]) and n2o_norm_combo_grid_late[2,*,z] ge 3 and tp_alt_grid-0.5 le alt_grid2[z],nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_combo_grid_late[0,ii,z]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;  ii = where(finite(n2o_norm_combo_grid_late[0,*,z]) and n2o_norm_combo_grid_late[2,*,z] ge 1 and n2o_norm_combo_grid_late[2,*,z] le 2 and tp_alt_grid-0.5 le alt_grid2[z],nii)
;  if nii gt 0 then begin
;    vcols = (n2o_norm_combo_grid_late[0,ii,z]-0.2)*300.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1,/overplot)
;    p = plot(lat_grid[ii]-1,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1,/overplot)
;    p = plot(lat_grid[ii]+1,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(9)*30.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_N2O_late.png'

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='Mean Ages 2015-2025')
;for z = 0, nz2-1 do begin
;  ii = where(finite(age_combo_grid_late[0,*,z]) and age_combo_grid_late[2,*,z] ge 3,nii)
;  if nii gt 0 then begin
;    vcols = age_combo_grid_late[0,ii,z]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]-1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;    p = plot(lat_grid[ii]+1.4,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=2,/overplot)
;  endif
;  ii = where(finite(age_combo_grid_late[0,*,z]) and age_combo_grid_late[2,*,z] ge 1 and age_combo_grid_late[2,*,z] le 2,nii)
;  if nii gt 0 then begin
;    vcols = age_combo_grid_late[0,ii,z]*40.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1,/overplot)
;    p = plot(lat_grid[ii]-1,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1,/overplot)
;    p = plot(lat_grid[ii]+1,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(6)*50.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_mean_ages_late.png'

;p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xtitle='Latitude',yrange=[0,6],ytitle='Mean Age (years)',xticklen=0.02,yticklen=0.02,title='Mean age at 20 km', $
;  margin=[0.08,0.1,0.04,0.08],font_size=11,dimensions=[700,400])
;for t = 0, 3 do p = errorplot(lat_grid,age_combo_grid_seas[0,*,z2,t],replicate(0,nyg),age_combo_grid_seas[1,*,z2,t],symbol='s',/sym_filled,color=colors[t],linestyle=6,errorbar_capsize=0, $
;  /overplot)
;p = plot(lat_grid,age_combo_grid_late[0,*,z2],color='blue',thick=3,/overplot)

;p = plot(indgen(2),/nodata,xrange=[5,85],xtickinterval=10,xminor=1,yrange=[9,25],position=[0.08,0.1,0.87,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;  yticklen=0.03,dimensions=[700,450],title='In Situ Mean Age Trends')
;for z = 0, nz2-1 do begin
;  ii = where(finite(age_grid_trends[0:-3,z,1]) and tp_alt_grid[0:-3]-0.5 le alt_grid2[z],nii)
;  if nii gt 0 then begin
;;    vcols = (-10.*age_grid_trends[ii,z,1] + 0.8)*159.
;    vcols = (0.5 - 10.*age_grid_trends[ii,z,1])*255.
;    test = where(vcols ge 250.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 250.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.88,0.2,0.9,0.75],tickname=tickname5,tickvalues=findgen(11)*25.5,RGB_TABLE=reverse(rgb_70),BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years/decade')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_mean_ages_trends_grid_elat_adj.png'

;p = plot(indgen(2),/nodata,xrange=[5,85],xtickinterval=10,xminor=1,yrange=[9,25],position=[0.08,0.1,0.87,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)', $
;  font_size=11,xticklen=0.04,yticklen=0.03,dimensions=[700,450],title='In Situ Normalized N$_2$O Trends')
;for z = 0, nz2-1 do begin
;  ii = where(finite(n2o_grid_trends[0:-3,z,1]) and tp_alt_grid[0:-3]-0.5 le alt_grid2[z],nii)
;  if nii gt 0 then begin
;    vcols = (-10.*n2o_grid_trends[ii,z,1] + 0.05)*2550.
;    test = where(vcols ge 250.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 250.)] = 250.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;  endif
;endfor
;;c = contour(n2o_grid_trends[*,*,1],lat_grid,alt_grid2,c_value=[0,0],c_color=['blue','blue'],c_label_show=0,/overplot)
;p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.88,0.2,0.9,0.75],tickname=tickname4,tickvalues=findgen(11)*25.5,RGB_TABLE=reverse(rgb_70),BORDER=1, ORIENTATION=1,/textpos,/tickdir, $
;  font_size=10,title='Normalized N$_2$O/decade')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_n2o_trends_grid_elat_adj.png'

;p = plot(indgen(2),/nodata,xrange=[5,85],xtickinterval=10,xminor=1,yrange=[9,25],position=[0.08,0.1,0.87,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)', $
;  font_size=11,xticklen=0.04,yticklen=0.03,dimensions=[700,450],title='In Situ Normalized CH$_4$ Trends')
;for z = 0, nz2-1 do begin
;  ii = where(finite(ch4_grid_trends[0:-3,z,1]) and tp_alt_grid[0:-3]-0.5 le alt_grid2[z],nii)
;  if nii gt 0 then begin
;    vcols = (-10.*ch4_grid_trends[ii,z,1] + 0.05)*2550.
;    test = where(vcols ge 250.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 250.)] = 250.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;  endif
;endfor
;;c = contour(n2o_grid_trends[*,*,1],lat_grid,alt_grid2,c_value=[0,0],c_color=['blue','blue'],c_label_show=0,/overplot)
;p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.88,0.2,0.9,0.75],tickname=tickname4,tickvalues=findgen(11)*25.5,RGB_TABLE=reverse(rgb_70),BORDER=1, ORIENTATION=1,/textpos,/tickdir, $
;  font_size=10,title='Normalized CH$_4$/decade')
;p.save,dir+'Plots/Balloon Mean Ages/In_situ_ch4_trends_grid_elat_adj.png'

;p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,26],position=[0.09,0.1,0.87,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)', $
;  font_size=14,xticklen=0.04,yticklen=0.03,dimensions=[700,450],title='Mean Age Changes 1990s to 2015-25 '+seas[s])
;for z = 0, nz2-1 do begin
;  ii = where(finite(age_combo_grid_seas_late[0,*,z,s]) and finite(age_combo_grid_seas_sm_early[0,*,z,s]) and tp_alt_grid[0:-3]-0.5 le alt_grid[z],nii)
;  if nii gt 0 then begin
;    vcols = (-(age_combo_grid_seas_late[0,ii,z,s]-age_combo_grid_seas_sm_early[0,ii,z,s]) + 0.8)*156.
;    test = where(vcols ge 240.,ntest)
;    if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;    test = where(vcols lt 0.,ntest)
;    if ntest gt 0 then vcols[test] = 0.
;    p = plot(lat_grid_a[ii],replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid_a[ii]-1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    p = plot(lat_grid_a[ii]+1.6,replicate(alt_grid2[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;  endif
;endfor
;p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;cb = COLORBAR(POSITION=[0.88,0.2,0.9,0.75],tickname=tickname5,tickvalues=findgen(9)*31,RGB_TABLE=reverse(rgb_70),BORDER=1, ORIENTATION=1,/textpos,/tickdir, $
;  font_size=10,title='years')

end