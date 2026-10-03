;+
;  Step 2: seasonal 1990s grids.
;
;  Combines the gridded 1990s aircraft and balloon mean ages and N2O into seasonal
;  equivalent latitude-altitude grids, used for the latitude adjustment of the
;  time series (step 3) and for Supplementary Fig. 4.
;
;  Reads:  <output>/intermediate/gridded_1990s_balloon_aircraft.nc, gridded_aircore_aircraft_late.nc;
;          Aircraft/.../aircraft_mean_age_grid_early.sav
;  Writes: <output>/intermediate/seasonal_grids_1990s.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step2_seasonal_grids

dir = mean_age_data_dir()

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_age_grid_early.sav'
grid_a0 = grid
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_1990s_balloon_aircraft.nc', $
  ['age_both_grid_mon_avg_b_early','age_co2_grid_mon_avg_a_early','age_co2_grid_mon_avg_a_mid','age_co2_grid_mon_avg_b_early', $
  'age_grid_mon_avg_b_early','age_opt_grid_mon_avg_a_early','age_opt_grid_mon_avg_a_mid','age_sf6_grid_mon_avg_a_early', $
  'age_sf6_grid_mon_avg_a_mid','combo_n2o_seas_early_hires','combo_n2o_seas_hires','n2o_norm_both_grid_mon_avg_b_early', $
  'n2o_norm_grid_mon_avg_a_early','n2o_norm_grid_mon_avg_a_mid','n2o_norm_grid_seas_a1','n2o_norm_sf6_grid_mon_avg_a_early', $
  'n2o_norm_sf6_grid_mon_avg_a_mid']
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_aircore_aircraft_late.nc', $
  ['age_co2_grid_mon_avg_a_late','age_co2_grid_mon_avg_b_late','age_co2_grid_seas','age_combo_grid_late', $
  'age_combo_grid_mon_avg','age_combo_grid_seas_late','age_opt_grid_mon_avg_a_late','age_opt_grid_mon_avg_b_late', $
  'age_sf6_grid_mon_avg_a_late','age_sf6_grid_mon_avg_b_late','n2o_norm_combo_grid_late','n2o_norm_combo_grid_minmax_late', $
  'n2o_norm_combo_grid_mon_avg','n2o_norm_combo_grid_seas_late','n2o_norm_grid_mon_avg_a_late','n2o_norm_grid_mon_avg_b_late', $
  'n2o_norm_grid_seas_late','n2o_norm_opt_grid_mon_avg_b_late','n2o_norm_sf6_grid_mon_avg_b_late']
nyg = 36
nz = 60
rrr = replicate(!values.f_nan,3,nyg,nz,4)
age_co2_opt_grid_mon_avg_sm = age_co2_grid_mon_avg_a_early
age_opt_grid_mon_avg_sm = age_opt_grid_mon_avg_a_early
age_sf6_opt_grid_mon_avg_sm = age_sf6_grid_mon_avg_a_early
age_co2_opt_grid_mon_avg_sm_late = age_co2_grid_mon_avg_a_late
age_co2_opt_grid_mon_avg_sm_late_b = age_co2_grid_mon_avg_b_late
age_co2_opt_grid_mon_avg_sm_b = age_co2_grid_mon_avg_b_early
age_sf6_opt_grid_mon_avg_sm_late = age_sf6_grid_mon_avg_a_late
age_opt_grid_mon_avg_sm_late = age_opt_grid_mon_avg_a_late
age_sf6_opt_grid_mon_avg_sm_late_b = age_sf6_grid_mon_avg_b_late
age_opt_grid_mon_avg_sm_late_b = age_opt_grid_mon_avg_b_late
age_co2_old_grid_mon_avg_sm = grid_a0.age_co2_grid_mon_avg
age_co2_old_grid_seas_sm = rrr
age_opt_grid_mon_avg_sm_b = age_grid_mon_avg_b_early
n2o_norm_grid_mon_avg = n2o_norm_grid_mon_avg_a_early
age_combo_grid_seas_sm_both = rrr
n2o_norm_grid_mon_avg_b = n2o_norm_both_grid_mon_avg_b_early
n2o_norm_grid_seas_both = rrr
age_co2_opt_grid_mon_avg_sm_mid = age_co2_grid_mon_avg_a_mid
age_sf6_opt_grid_mon_avg_sm_mid = age_sf6_grid_mon_avg_a_mid
age_opt_grid_mon_avg_sm_mid = age_opt_grid_mon_avg_a_mid
n2o_norm_grid_mon_avg_mid = n2o_norm_grid_mon_avg_a_mid

grid_fill_in_situ,nyg,nz,age_co2_opt_grid_mon_avg_sm
grid_fill_in_situ,nyg,nz,age_sf6_opt_grid_mon_avg_sm
grid_fill_in_situ,nyg,nz,age_opt_grid_mon_avg_sm
grid_fill_in_situ,nyg,nz,age_co2_opt_grid_mon_avg_sm_b
grid_fill_in_situ,nyg,nz,age_opt_grid_mon_avg_sm_b
grid_fill_in_situ,nyg,nz,age_co2_opt_grid_mon_avg_sm_mid
grid_fill_in_situ,nyg,nz,age_sf6_opt_grid_mon_avg_sm_mid
grid_fill_in_situ,nyg,nz,age_opt_grid_mon_avg_sm_mid
grid_fill_in_situ,nyg,nz,age_co2_opt_grid_mon_avg_sm_late
grid_fill_in_situ,nyg,nz,age_sf6_opt_grid_mon_avg_sm_late
grid_fill_in_situ,nyg,nz,age_opt_grid_mon_avg_sm_late
grid_fill_in_situ,nyg,nz,age_co2_opt_grid_mon_avg_sm_late_b
grid_fill_in_situ,nyg,nz,age_sf6_opt_grid_mon_avg_sm_late_b
grid_fill_in_situ,nyg,nz,age_opt_grid_mon_avg_sm_late_b
grid_fill_in_situ,nyg-2,nz,age_co2_old_grid_mon_avg_sm
grid_fill_in_situ,nyg,nz,n2o_norm_grid_mon_avg
grid_fill_in_situ,nyg,nz,n2o_norm_grid_mon_avg_mid
grid_fill_in_situ,nyg,nz,n2o_norm_grid_mon_avg_b
grid_fill_in_situ,nyg,nz,age_both_grid_mon_avg_b_early

;  Aircraft early
age_combo_grid_mon_avg_sm = age_co2_opt_grid_mon_avg_sm

;  Balloon early
age_combo_grid_mon_avg_sm_b = age_both_grid_mon_avg_b_early

age_combo_grid_mon_avg_sm_both = age_combo_grid_mon_avg_sm
n2o_norm_grid_mon_avg_both = n2o_norm_grid_mon_avg
for m = 0, 11 do for z = 0, nz-1 do for y = 0, nyg-1 do begin
  for j = 0, 1 do begin
    n2o_norm_grid_mon_avg_both[j,y,z,m] = mean([n2o_norm_grid_mon_avg[j,y,z,m],n2o_norm_grid_mon_avg_b[j,y,z,m]],/nan)
    age_combo_grid_mon_avg_sm_both[j,y,z,m] = mean([age_combo_grid_mon_avg_sm[j,y,z,m],age_combo_grid_mon_avg_sm_b[j,y,z,m]],/nan)
  endfor
  j = 2
  n2o_norm_grid_mon_avg_both[j,y,z,m] = total([n2o_norm_grid_mon_avg[j,y,z,m],n2o_norm_grid_mon_avg_b[j,y,z,m]],/nan)
  age_combo_grid_mon_avg_sm_both[j,y,z,m] = total([age_combo_grid_mon_avg_sm[j,y,z,m],age_combo_grid_mon_avg_sm_b[j,y,z,m]],/nan)
endfor

for z = 1, nz-2 do for y = 0, nyg-1 do begin
  for j = 0, 1 do begin
    n2o_norm_grid_seas_both[j,y,z,0] = mean([n2o_norm_grid_mon_avg_both[j,y,z,11], $
      reform(n2o_norm_grid_mon_avg_both[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas_both[j,y,z,s] = mean(n2o_norm_grid_mon_avg_both[j,y,z,indgen(3)+3*s-1],/nan)
    age_combo_grid_seas_sm_both[j,y,z,0] = mean([age_combo_grid_mon_avg_sm_both[j,y,z,11], $
      reform(age_combo_grid_mon_avg_sm_both[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_combo_grid_seas_sm_both[j,y,z,s] = mean(age_combo_grid_mon_avg_sm_both[j,y,z,indgen(3)+3*s-1],/nan)
    if y ge 1 and y le nyg-2 then begin
      age_co2_old_grid_seas_sm[j,y,z,0] = mean([reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,11]), $
        reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1,0:1]),  reform(age_co2_old_grid_mon_avg_sm[j,y-1,z,0:1]), $
        reform(age_co2_old_grid_mon_avg_sm[j,y-1,z+1,0:1])],/nan)
      for s = 1, 3 do age_co2_old_grid_seas_sm[j,y,z,s] = mean(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,indgen(3)+3*s-1],/nan)
    endif
  endfor
  j = 2
  n2o_norm_grid_seas_both[j,y,z,0] = total([n2o_norm_grid_mon_avg_both[j,y,z,11],reform(n2o_norm_grid_mon_avg_both[j,y,z,0:1])], $
    /nan)
  for s = 1, 3 do n2o_norm_grid_seas_both[j,y,z,s] = total(n2o_norm_grid_mon_avg_both[j,y,z,indgen(3)+3*s-1],/nan)
  age_combo_grid_seas_sm_both[j,y,z,0] = total([age_combo_grid_mon_avg_sm_both[j,y,z,11], $
    reform(age_combo_grid_mon_avg_sm_both[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_combo_grid_seas_sm_both[j,y,z,s] = total(age_combo_grid_mon_avg_sm_both[j,y,z,indgen(3)+3*s-1],/nan)
  if y ge 1 and y le nyg-2 then begin
    age_co2_old_grid_seas_sm[j,y,z,0] = total([reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,11]), $
      reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1,0:1]),  reform(age_co2_old_grid_mon_avg_sm[j,y-1,z,0:1]), $
      reform(age_co2_old_grid_mon_avg_sm[j,y-1,z+1,0:1])],/nan)
    for s = 1, 3 do age_co2_old_grid_seas_sm[j,y,z,s] = total(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,indgen(3)+3*s-1],/nan)
  endif
endfor

age_co2_old_grid_seas_sm[*,0:10,28,3] = !values.f_nan

;  Fill in missing locations.
for s = 0, 3 do for z = 1, nz-2 do for y = 1, nyg-2 do for j = 0, 1 do begin
  if ~finite(age_combo_grid_seas_sm_both[j,y,z,s]) and finite(age_combo_grid_seas_sm_both[j,y-1,z, $
    s]) and finite(age_combo_grid_seas_sm_both[j,y+1,z,s]) then age_combo_grid_seas_sm_both[j,y,z, $
    s] = 0.5 * age_combo_grid_seas_sm_both[j,y-1,z,s] + 0.5 * age_combo_grid_seas_sm_both[j,y+1,z,s]
  if ~finite(age_combo_grid_seas_sm_both[j,y,z,s]) and finite(age_combo_grid_seas_sm_both[j,y,z-1, $
    s]) and finite(age_combo_grid_seas_sm_both[j,y,z+1,s]) then age_combo_grid_seas_sm_both[j,y,z, $
    s] = 0.5 * age_combo_grid_seas_sm_both[j,y,z-1,s] + 0.5 * age_combo_grid_seas_sm_both[j,y,z+1,s]
endfor

age_combo_grid_seas_sm_early = age_combo_grid_seas_sm_both

cc = replicate(!values.f_nan,4,nyg,nz)
age_combo_grid_early = cc
n2o_norm_grid_early = cc
cgrid = age_combo_grid_seas_sm_early
for z = 1, nz-2 do for y = 0, nyg-1 do begin
  gd = where(finite(cgrid[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(cgrid[0,y,z,gd],sdev=sdev)
    tot = total(cgrid[2,y,z,gd])
    age_combo_grid_early[*,y,z] = [stats[0],sdev,ngd,tot]
  endif
  gd = where(finite(n2o_norm_grid_seas_both[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(n2o_norm_grid_seas_both[0,y,z,gd],sdev=sdev)
    tot = total(n2o_norm_grid_seas_both[2,y,z,gd])
    n2o_norm_grid_early[*,y,z] = [stats[0],sdev,ngd,tot]
  endif
endfor
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'seasonal_grids_1990s.nc', $
  ['age_combo_grid_seas_sm_early','age_combo_grid_early','n2o_norm_grid_early','age_co2_old_grid_seas_sm', $
  'n2o_norm_grid_seas_both']
end

