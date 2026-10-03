pro make_balloon_seas_grid

;  Builds the seasonal equivalent latitude-altitude mean age and N2O grids for
;  the 1990s (early) period from the gridded aircraft and balloon data and the
;  equivalent-latitude-adjusted grids written by check_elats (stage 'elats').
;  Writes Balloon/balloon_mean_age_seas_grid_early.sav, which check_elats
;  (stage 'time_series') and age_time_series read.
;
;  This block was moved here unchanged from age_time_series.pro so that the
;  programs can run in a single forward order (see run_all.pro).

dir = mean_age_data_dir()

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_early.sav'
dat_co2_a0 = dat_co2

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_age_grid_early.sav'
grid_a0 = grid
restore,dir+'Balloon/balloon_mean_age_grid_early.sav'
grid_b0 = grid
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_age_grid_late.sav'
grid_a2 = grid
restore,dir+'Balloon/balloon_mean_age_grid_late_v2025.sav'
grid_b2 = grid
restore,dir+'Balloon/OMS_grid_elat_adj.sav'
restore,dir+'Balloon/Aircore_grid_elat_adj.sav'
restore,dir+'Balloon/Elat_adj_all.sav'

nyg = 36
dy = 5.
nz = 60
dz = 0.5
lat_grid_a = findgen(nyg)*dy-85.
nza = 14
dza = 1.0
alt_grid_a0 = findgen(nz)*dz+6.
alt_grid_a02 = findgen(nza)*dza+12.


age_co2_opt_grid_avg = replicate(!values.f_nan,3,nyg,nza) & rrr = replicate(!values.f_nan,3,nyg,nz,4) & age_co2_opt_grid_mon_avg_sm = age_co2_grid_mon_avg_a_early & age_co2_opt_grid_seas_sm = rrr
age_opt_grid_mon_avg_sm = age_opt_grid_mon_avg_a_early & age_opt_grid_seas_sm = rrr & age_sf6_opt_grid_mon_avg_sm = age_sf6_grid_mon_avg_a_early
vvv = replicate(!values.f_nan,2,nyg,nz) & age_combo_grid_sm_early = vvv & age_combo_grid_minmax_early = vvv & age_co2_opt_grid_mon_avg_sm_late = age_co2_grid_mon_avg_a_late
age_co2_opt_grid_mon_avg_sm_late_b = age_co2_grid_mon_avg_b_late & age_co2_opt_grid_mon_avg_sm_b = age_co2_grid_mon_avg_b_early & age_combo_grid_seas_sm = rrr
age_sf6_opt_grid_mon_avg_sm_late = age_sf6_grid_mon_avg_a_late & age_opt_grid_mon_avg_sm_late = age_opt_grid_mon_avg_a_late & age_combo_grid_seas_sm_late = rrr
age_sf6_opt_grid_mon_avg_sm_late_b = age_sf6_grid_mon_avg_b_late & age_opt_grid_mon_avg_sm_late_b = age_opt_grid_mon_avg_b_late & age_co2_old_grid_sm_early = vvv
age_co2_old_grid_mon_avg_sm = grid_a0.age_co2_grid_mon_avg & age_co2_old_grid_seas_sm = rrr & n2o_norm_grid_seas = rrr & age_combo_grid_seas_sm_mid = rrr
age_sf6_opt_grid_mon_avg_sm_b = grid_b0.age_opt_sf6_grid_mon_avg & age_opt_grid_mon_avg_sm_b = age_grid_mon_avg_b_early & n2o_norm_grid_mon_avg = n2o_norm_grid_mon_avg_a_early
age_combo_grid_seas_sm_both = rrr & n2o_norm_grid_mon_avg_b = n2o_norm_both_grid_mon_avg_b_early & n2o_norm_grid_seas_both = rrr & n2o_norm_grid_seas_mid = rrr
age_co2_opt_grid_mon_avg_sm_mid = age_co2_grid_mon_avg_a_mid & age_sf6_opt_grid_mon_avg_sm_mid = age_sf6_grid_mon_avg_a_mid
age_opt_grid_mon_avg_sm_mid = age_opt_grid_mon_avg_a_mid & n2o_norm_grid_mon_avg_mid = n2o_norm_grid_mon_avg_a_mid
for z = 0, nza-1 do for y = 0, nyg-1 do begin
  ii_all = where(dat_co2_a0.alt ge alt_grid_a02[z]-dza/2. and dat_co2_a0.alt lt alt_grid_a02[z]+dza/2. and dat_co2_a0.elat ge lat_grid_a[y]-dy/2. and dat_co2_a0.elat lt lat_grid_a[y]+dy/2.,tot_i)
  if tot_i gt 2 then begin
    stats = moment(dat_co2_a0.age_opt_all[ii_all],sd=sd,/nan)
    age_co2_opt_grid_avg[*,y,z] = [stats[0],sd,tot_i]
  endif
endfor

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
tmp1 = reform(age_combo_grid_mon_avg_sm[0,*,*,*]) & tmp11 = reform(age_combo_grid_mon_avg_sm[1,*,*,*]) & tmp2 = reform(age_opt_grid_mon_avg_sm[0,*,*,*])
tmp3 = reform(age_opt_grid_mon_avg_sm[1,*,*,*])
chk = where(tmp1 lt 1.5 and finite(tmp2))
tmp1[chk] = tmp2[chk]
tmp11[chk] = tmp3[chk]

;  Balloon early
age_combo_grid_mon_avg_sm_b = age_both_grid_mon_avg_b_early

;  Aircraft mid
age_combo_grid_mon_avg_sm_mid = age_co2_opt_grid_mon_avg_sm_mid
tmp1 = reform(age_combo_grid_mon_avg_sm_mid[0,*,*,*]) & tmp11 = reform(age_combo_grid_mon_avg_sm_mid[1,*,*,*]) & tmp2 = reform(age_sf6_opt_grid_mon_avg_sm_mid[0,*,*,*])
tmp3 = reform(age_sf6_opt_grid_mon_avg_sm_mid[1,*,*,*])
chk = where(finite(tmp2))
tmp1[chk] = tmp2[chk]
tmp11[chk] = tmp3[chk]
age_combo_grid_mon_avg_sm_mid[0,*,*,*] = tmp1
age_combo_grid_mon_avg_sm_mid[1,*,*,*] = tmp11
chk = where(finite(age_opt_grid_mon_avg_sm_mid))
age_combo_grid_mon_avg_sm_mid[chk] = age_opt_grid_mon_avg_sm_mid[chk]

age_combo_grid_mon_avg_sm_late = age_co2_opt_grid_mon_avg_sm_late
tmp1 = reform(age_combo_grid_mon_avg_sm_late[0,*,*,*]) & tmp11 = reform(age_combo_grid_mon_avg_sm_late[1,*,*,*]) & tmp2 = reform(age_sf6_opt_grid_mon_avg_sm_late[0,*,*,*])
tmp3 = reform(age_sf6_opt_grid_mon_avg_sm_late[1,*,*,*])
chk = where(finite(tmp2))
tmp1[chk] = tmp2[chk]
tmp11[chk] = tmp3[chk]
age_combo_grid_mon_avg_sm_late[0,*,*,*] = tmp1
age_combo_grid_mon_avg_sm_late[1,*,*,*] = tmp11
chk = where(finite(age_opt_grid_mon_avg_sm_late))
age_combo_grid_mon_avg_sm_late[chk] = age_opt_grid_mon_avg_sm_late[chk]

age_combo_grid_mon_avg_sm_late_b = age_co2_opt_grid_mon_avg_sm_late_b
tmp1 = reform(age_combo_grid_mon_avg_sm_late_b[0,*,*,*]) & tmp11 = reform(age_combo_grid_mon_avg_sm_late_b[1,*,*,*]) & tmp2 = reform(age_sf6_opt_grid_mon_avg_sm_late_b[0,*,*,*])
tmp3 = reform(age_sf6_opt_grid_mon_avg_sm_late_b[1,*,*,*])
chk = where(finite(tmp2))
tmp1[chk] = tmp2[chk]
tmp11[chk] = tmp3[chk]
age_combo_grid_mon_avg_sm_late_b[0,*,*,*] = tmp1
age_combo_grid_mon_avg_sm_late_b[1,*,*,*] = tmp11
chk = where(finite(age_opt_grid_mon_avg_sm_late_b))
age_combo_grid_mon_avg_sm_late_b[chk] = age_opt_grid_mon_avg_sm_late_b[chk]

age_combo_grid_mon_avg_sm_both = age_combo_grid_mon_avg_sm & n2o_norm_grid_mon_avg_both = n2o_norm_grid_mon_avg & age_combo_grid_mon_avg_sm_late_both = age_combo_grid_mon_avg_sm_late
for m = 0, 11 do for z = 0, nz-1 do for y = 0, nyg-1 do begin
  for j = 0, 1 do begin
    n2o_norm_grid_mon_avg_both[j,y,z,m] = mean([n2o_norm_grid_mon_avg[j,y,z,m],n2o_norm_grid_mon_avg_b[j,y,z,m]],/nan)
    age_combo_grid_mon_avg_sm_both[j,y,z,m] = mean([age_combo_grid_mon_avg_sm[j,y,z,m],age_combo_grid_mon_avg_sm_b[j,y,z,m]],/nan)
    age_combo_grid_mon_avg_sm_late_both[j,y,z,m] = mean([age_combo_grid_mon_avg_sm_late[j,y,z,m],age_combo_grid_mon_avg_sm_late_b[j,y,z,m]],/nan)
  endfor
  j = 2
  n2o_norm_grid_mon_avg_both[j,y,z,m] = total([n2o_norm_grid_mon_avg[j,y,z,m],n2o_norm_grid_mon_avg_b[j,y,z,m]],/nan)
  age_combo_grid_mon_avg_sm_both[j,y,z,m] = total([age_combo_grid_mon_avg_sm[j,y,z,m],age_combo_grid_mon_avg_sm_b[j,y,z,m]],/nan)
  age_combo_grid_mon_avg_sm_late_both[j,y,z,m] = total([age_combo_grid_mon_avg_sm_late[j,y,z,m],age_combo_grid_mon_avg_sm_late_b[j,y,z,m]],/nan)
endfor

for z = 1, nz-2 do for y = 0, nyg-1 do begin
  for j = 0, 1 do begin
    n2o_norm_grid_seas[j,y,z,0] = mean([n2o_norm_grid_mon_avg[j,y,z,11],reform(n2o_norm_grid_mon_avg[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas[j,y,z,s] = mean(n2o_norm_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
    n2o_norm_grid_seas_mid[j,y,z,0] = mean([n2o_norm_grid_mon_avg_mid[j,y,z,11],reform(n2o_norm_grid_mon_avg_mid[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas_mid[j,y,z,s] = mean(n2o_norm_grid_mon_avg_mid[j,y,z,indgen(3)+3*s-1],/nan)
    n2o_norm_grid_seas_both[j,y,z,0] = mean([n2o_norm_grid_mon_avg_both[j,y,z,11],reform(n2o_norm_grid_mon_avg_both[j,y,z,0:1])],/nan)
    for s = 1, 3 do n2o_norm_grid_seas_both[j,y,z,s] = mean(n2o_norm_grid_mon_avg_both[j,y,z,indgen(3)+3*s-1],/nan)
    age_co2_opt_grid_seas_sm[j,y,z,0] = mean([age_co2_opt_grid_mon_avg_sm[j,y,z,11],reform(age_co2_opt_grid_mon_avg_sm[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_co2_opt_grid_seas_sm[j,y,z,s] = mean(age_co2_opt_grid_mon_avg_sm[j,y,z,indgen(3)+3*s-1],/nan)
    age_opt_grid_seas_sm[j,y,z,0] = mean([age_opt_grid_mon_avg_sm[j,y,z,11],reform(age_opt_grid_mon_avg_sm[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_opt_grid_seas_sm[j,y,z,s] = mean(age_opt_grid_mon_avg_sm[j,y,z,indgen(3)+3*s-1],/nan)
    age_combo_grid_seas_sm[j,y,z,0] = mean([age_combo_grid_mon_avg_sm[j,y,z,11],reform(age_combo_grid_mon_avg_sm[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_combo_grid_seas_sm[j,y,z,s] = mean(age_combo_grid_mon_avg_sm[j,y,z,indgen(3)+3*s-1],/nan)
    age_combo_grid_seas_sm_both[j,y,z,0] = mean([age_combo_grid_mon_avg_sm_both[j,y,z,11],reform(age_combo_grid_mon_avg_sm_both[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_combo_grid_seas_sm_both[j,y,z,s] = mean(age_combo_grid_mon_avg_sm_both[j,y,z,indgen(3)+3*s-1],/nan)
    age_combo_grid_seas_sm_mid[j,y,z,0] = mean([age_combo_grid_mon_avg_sm_mid[j,y,z,11],reform(age_combo_grid_mon_avg_sm_mid[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_combo_grid_seas_sm_mid[j,y,z,s] = mean(age_combo_grid_mon_avg_sm_mid[j,y,z,indgen(3)+3*s-1],/nan)
    age_combo_grid_seas_sm_late[j,y,z,0] = mean([age_combo_grid_mon_avg_sm_late_both[j,y,z,11],reform(age_combo_grid_mon_avg_sm_late_both[j,y,z,0:1])],/nan)
    for s = 1, 3 do age_combo_grid_seas_sm_late[j,y,z,s] = mean(age_combo_grid_mon_avg_sm_late_both[j,y,z,indgen(3)+3*s-1],/nan)
    if y ge 1 and y le nyg-2 then begin
      age_co2_old_grid_seas_sm[j,y,z,0] = mean([reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,11]),reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1,0:1]), $
        reform(age_co2_old_grid_mon_avg_sm[j,y-1,z,0:1]),reform(age_co2_old_grid_mon_avg_sm[j,y-1,z+1,0:1])],/nan)
      for s = 1, 3 do age_co2_old_grid_seas_sm[j,y,z,s] = mean(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,indgen(3)+3*s-1],/nan)
    endif
  endfor
  j = 2
  n2o_norm_grid_seas[j,y,z,0] = total([n2o_norm_grid_mon_avg[j,y,z,11],reform(n2o_norm_grid_mon_avg[j,y,z,0:1])],/nan)
  for s = 1, 3 do n2o_norm_grid_seas[j,y,z,s] = total(n2o_norm_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
  n2o_norm_grid_seas_mid[j,y,z,0] = total([n2o_norm_grid_mon_avg_mid[j,y,z,11],reform(n2o_norm_grid_mon_avg_mid[j,y,z,0:1])],/nan)
  for s = 1, 3 do n2o_norm_grid_seas_mid[j,y,z,s] = total(n2o_norm_grid_mon_avg_mid[j,y,z,indgen(3)+3*s-1],/nan)
  n2o_norm_grid_seas_both[j,y,z,0] = total([n2o_norm_grid_mon_avg_both[j,y,z,11],reform(n2o_norm_grid_mon_avg_both[j,y,z,0:1])],/nan)
  for s = 1, 3 do n2o_norm_grid_seas_both[j,y,z,s] = total(n2o_norm_grid_mon_avg_both[j,y,z,indgen(3)+3*s-1],/nan)
  age_co2_opt_grid_seas_sm[j,y,z,0] = total([age_co2_opt_grid_mon_avg_sm[j,y,z,11],reform(age_co2_opt_grid_mon_avg_sm[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_co2_opt_grid_seas_sm[j,y,z,s] = total(age_co2_opt_grid_mon_avg_sm[j,y,z,indgen(3)+3*s-1],/nan)
  age_opt_grid_seas_sm[j,y,z,0] = total([age_opt_grid_mon_avg_sm[j,y,z,11],reform(age_opt_grid_mon_avg_sm[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_opt_grid_seas_sm[j,y,z,s] = total(age_opt_grid_mon_avg_sm[j,y,z,indgen(3)+3*s-1],/nan)
  age_combo_grid_seas_sm[j,y,z,0] = total([age_combo_grid_mon_avg_sm[j,y,z,11],reform(age_combo_grid_mon_avg_sm[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_combo_grid_seas_sm[j,y,z,s] = total(age_combo_grid_mon_avg_sm[j,y,z,indgen(3)+3*s-1],/nan)
  age_combo_grid_seas_sm_both[j,y,z,0] = total([age_combo_grid_mon_avg_sm_both[j,y,z,11],reform(age_combo_grid_mon_avg_sm_both[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_combo_grid_seas_sm_both[j,y,z,s] = total(age_combo_grid_mon_avg_sm_both[j,y,z,indgen(3)+3*s-1],/nan)
  age_combo_grid_seas_sm_mid[j,y,z,0] = total([age_combo_grid_mon_avg_sm_mid[j,y,z,11],reform(age_combo_grid_mon_avg_sm_mid[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_combo_grid_seas_sm_mid[j,y,z,s] = total(age_combo_grid_mon_avg_sm_mid[j,y,z,indgen(3)+3*s-1],/nan)
  age_combo_grid_seas_sm_late[j,y,z,0] = total([age_combo_grid_mon_avg_sm_late_both[j,y,z,11],reform(age_combo_grid_mon_avg_sm_late_both[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_combo_grid_seas_sm_late[j,y,z,s] = total(age_combo_grid_mon_avg_sm_late_both[j,y,z,indgen(3)+3*s-1],/nan)
  if y ge 1 and y le nyg-2 then begin
    age_co2_old_grid_seas_sm[j,y,z,0] = total([reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,11]),reform(age_co2_old_grid_mon_avg_sm[j,y-1,z-1,0:1]), $
      reform(age_co2_old_grid_mon_avg_sm[j,y-1,z,0:1]),reform(age_co2_old_grid_mon_avg_sm[j,y-1,z+1,0:1])],/nan)
    for s = 1, 3 do age_co2_old_grid_seas_sm[j,y,z,s] = total(age_co2_old_grid_mon_avg_sm[j,y-1,z-1:z+1,indgen(3)+3*s-1],/nan)
  endif
endfor

age_co2_old_grid_seas_sm[*,0:10,28,3] = !values.f_nan

;  Fill in missing locations.
for s = 0, 3 do for z = 1, nz-2 do for y = 1, nyg-2 do for j = 0, 1 do begin
  if ~finite(age_combo_grid_seas_sm[j,y,z,s]) and finite(age_combo_grid_seas_sm[j,y-1,z,s]) and finite(age_combo_grid_seas_sm[j,y+1,z,s]) then $
    age_combo_grid_seas_sm[j,y,z,s] = 0.5 * age_combo_grid_seas_sm[j,y-1,z,s] + 0.5 * age_combo_grid_seas_sm[j,y+1,z,s]
  if ~finite(age_combo_grid_seas_sm[j,y,z,s]) and finite(age_combo_grid_seas_sm[j,y,z-1,s]) and finite(age_combo_grid_seas_sm[j,y,z+1,s]) then $
    age_combo_grid_seas_sm[j,y,z,s] = 0.5 * age_combo_grid_seas_sm[j,y,z-1,s] + 0.5 * age_combo_grid_seas_sm[j,y,z+1,s]
  if ~finite(age_combo_grid_seas_sm_both[j,y,z,s]) and finite(age_combo_grid_seas_sm_both[j,y-1,z,s]) and finite(age_combo_grid_seas_sm_both[j,y+1,z,s]) then $
    age_combo_grid_seas_sm_both[j,y,z,s] = 0.5 * age_combo_grid_seas_sm_both[j,y-1,z,s] + 0.5 * age_combo_grid_seas_sm_both[j,y+1,z,s]
  if ~finite(age_combo_grid_seas_sm_both[j,y,z,s]) and finite(age_combo_grid_seas_sm_both[j,y,z-1,s]) and finite(age_combo_grid_seas_sm_both[j,y,z+1,s]) then $
    age_combo_grid_seas_sm_both[j,y,z,s] = 0.5 * age_combo_grid_seas_sm_both[j,y,z-1,s] + 0.5 * age_combo_grid_seas_sm_both[j,y,z+1,s]
  if ~finite(age_combo_grid_seas_sm_late[j,y,z,s]) and finite(age_combo_grid_seas_sm_late[j,y-1,z,s]) and finite(age_combo_grid_seas_sm_late[j,y+1,z,s]) then $
    age_combo_grid_seas_sm_late[j,y,z,s] = 0.5 * age_combo_grid_seas_sm_late[j,y-1,z,s] + 0.5 * age_combo_grid_seas_sm_late[j,y+1,z,s]
  if ~finite(age_combo_grid_seas_sm_late[j,y,z,s]) and finite(age_combo_grid_seas_sm_late[j,y,z-1,s]) and finite(age_combo_grid_seas_sm_late[j,y,z+1,s]) then $
    age_combo_grid_seas_sm_late[j,y,z,s] = 0.5 * age_combo_grid_seas_sm_late[j,y,z-1,s] + 0.5 * age_combo_grid_seas_sm_late[j,y,z+1,s]
endfor

age_combo_grid_seas_sm_early = age_combo_grid_seas_sm_both
tmp1 = reform(age_combo_grid_seas_sm_early[0,*,*,*]) & tmp2 = reform(age_opt_grid_seas_sm[0,*,*,*])
chk = where(tmp1 lt 1.5 and finite(tmp2))
tmp1[chk] = tmp2[chk]

cc = replicate(!values.f_nan,4,nyg,nz) & age_combo_grid_early = cc & age_combo_grid_minmax_early = cc & n2o_norm_grid_early = cc & n2o_norm_grid_minmax_early = cc & age_co2_old_grid_sm_early = cc
cgrid = age_combo_grid_seas_sm_early & age_combo_grid_all = cc & n2o_norm_grid_mid = cc & n2o_norm_grid_minmax_mid = cc
dd = replicate(!values.f_nan,4,nz) & age_nh_avg_profile_early = dd
for z = 1, nz-2 do for y = 0, nyg-1 do begin
  gd = where(finite(cgrid[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    err = 1. / sqrt(total(cgrid[1,y,z,gd]^2))
    wmean = total(1. / cgrid[1,y,z,gd]^2 * cgrid[0,y,z,gd]) / total(1. / cgrid[1,y,z,gd]^2)
    stats = moment(cgrid[0,y,z,gd],sdev=sdev)
    tot = total(cgrid[2,y,z,gd])
    age_combo_grid_early[*,y,z] = [stats[0],sdev,ngd,tot]
    age_combo_grid_minmax_early[0:2,y,z] = [min(cgrid[0,y,z,gd]-cgrid[1,y,z,gd]),max(cgrid[0,y,z,gd]+cgrid[1,y,z,gd]),ngd]
  endif
  gd = where(finite(cgrid[0,y,z,*]) or finite(age_combo_grid_seas_sm_late[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment([reform(cgrid[0,y,z,*]),reform(age_combo_grid_seas_sm_late[0,y,z,*])],sdev=sdev,/nan)
    tot = total([reform(cgrid[0,y,z,*]),reform(age_combo_grid_seas_sm_late[0,y,z,*])],/nan)
    age_combo_grid_all[*,y,z] = [stats[0],sdev,ngd,tot]
  endif
  gd = where(finite(age_co2_old_grid_seas_sm[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(age_co2_old_grid_seas_sm[0,y,z,gd],sdev=sdev)
    tot = total(age_co2_old_grid_seas_sm[2,y,z,gd])
    age_co2_old_grid_sm_early[*,y,z] = [stats[0],sdev,ngd,tot]
  endif
  gd = where(finite(n2o_norm_grid_seas_both[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(n2o_norm_grid_seas_both[0,y,z,gd],sdev=sdev)
    tot = total(n2o_norm_grid_seas_both[2,y,z,gd])
    n2o_norm_grid_early[*,y,z] = [stats[0],sdev,ngd,tot]
    n2o_norm_grid_minmax_early[0:1,y,z] = [min(n2o_norm_grid_seas_both[0,y,z,gd]-n2o_norm_grid_seas_both[1,y,z,gd]), $
      max(n2o_norm_grid_seas_both[0,y,z,gd]+n2o_norm_grid_seas_both[1,y,z,gd])]
  endif
  gd = where(finite(n2o_norm_grid_seas_mid[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(n2o_norm_grid_seas_mid[0,y,z,gd],sdev=sdev)
    tot = total(n2o_norm_grid_seas_mid[2,y,z,gd])
    n2o_norm_grid_mid[*,y,z] = [stats[0],sdev,ngd,tot]
    n2o_norm_grid_minmax_mid[0:1,y,z] = [min(n2o_norm_grid_seas_mid[0,y,z,gd]-n2o_norm_grid_seas_mid[1,y,z,gd]), $
      max(n2o_norm_grid_seas_mid[0,y,z,gd]+n2o_norm_grid_seas_mid[1,y,z,gd])]
  endif
endfor
yi = where(lat_grid_a ge 30 and lat_grid_a le 60)
for z = 1, nz-2 do begin
  for j = 0, 1 do begin
    age_nh_avg_profile_early[j,z] = mean(age_combo_grid_early[j,yi,z],/nan)
  endfor
  age_nh_avg_profile_early[2,z] = total(age_combo_grid_early[2,yi,z],/nan)
endfor
for z = 1, nz-2 do begin
  gd = where(finite(age_combo_grid_all[0,*,z]),ngd)
  if ngd ge 10 then age_combo_grid_all[0,gd,z] = smooth(age_combo_grid_all[0,gd,z],3,/edge_truncate)
endfor
for y = 1, nyg-2 do begin
  gd = where(finite(age_combo_grid_all[0,y,*]),ngd)
  if ngd ge 10 then age_combo_grid_all[0,y,gd] = smooth(age_combo_grid_all[0,y,gd],3,/edge_truncate)
endfor

save,age_combo_grid_seas_sm_early,age_combo_grid_early,n2o_norm_grid_early,age_co2_old_grid_seas_sm,n2o_norm_grid_seas_both, $
  filename=dir+'Balloon/balloon_mean_age_seas_grid_early.sav'

end
