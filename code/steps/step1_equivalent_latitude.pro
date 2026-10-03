;+
;  Step 1: N2O equivalent latitude.
;
;  Assigns an N2O-based equivalent latitude to every in situ aircraft, balloon and
;  AirCore measurement by matching its normalized N2O to the SWOOSH monthly
;  latitude-altitude N2O climatology (WACCM below 15 km), and grids the 1990s
;  (aircraft and OMS balloon) and 2010s-2020s (aircraft and AirCore) mean ages and
;  N2O by equivalent latitude, altitude and month (Supplementary Note 2).
;
;  Reads:  merged aircraft, balloon and AirCore mean age files; swoosh/n2o_merge.sav;
;          Models/WACCM/FWSD_means.sav; NOAA surface N2O
;  Writes: <output>/intermediate/equivalent_latitude_adjustments.nc,
;          gridded_1990s_balloon_aircraft.nc, gridded_aircore_aircraft_late.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step1_equivalent_latitude

dir = mean_age_data_dir()

data = read_ascii(dir+'Surface_Trace_Gas/GMD/n2o/combined/GML_global_N2O.txt',data_start=68)
n2o_yr = reform(data.field01[0,*])
n2o_mon = reform(data.field01[1,*])
n2o_date = n2o_yr + (n2o_mon-0.5)/12.
n2o_global = reform(data.field01[6,*])

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_early.sav'
caldat,jday_all,months_all,days_all,years_all
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
co2_a0 = co2
sf6_a0 = sf6
elat_m_a0 = elat_m
lat_a0 = lat
alt_a0 = alt
n2o_norm_a0 = n2o_norm
months_all_a0 = months_all

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_mid.sav'
caldat,jday_all,months_all,days_all,years_all
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
co2_a1 = co2
sf6_a1 = sf6
elat_m_a1 = elat_m
lat_a1 = lat
alt_a1 = alt
n2o_norm_a1 = n2o_norm
months_all_a1 = months_all

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_late.sav'
caldat,jday_all,months_all,days_all,years_all
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
co2_a = co2
sf6_a = sf6
elat_m_a = elat_m
lat_a = lat
alt_a = alt
n2o_norm_a = n2o_norm
months_all_a = months_all

restore,dir+'Balloon/OMS_common_gc_merge.sav'
co2_b = co2
sf6_b = sf6
elat_m_b = elat_m
lat_b = lat
alt_b = alt
n2o_norm_b = n2o_norm
months_all_b = months_all
co2_all_b = co2_all

restore,dir+'Balloon/Aircore_common_merge.sav'
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_early.sav'
dat_co2_a0 = dat_co2
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sf6_early.sav'
dat_sf6_a0 = dat_sf6
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
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sweep_sf6_co2_late.sav'
dat_a2 = dat
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_late.sav'
dat_co2_a2 = dat_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late.sav'
dat_sf6_a2 = dat_sf6

in2o = where(finite(n2o_norm),nn2o)
ico2 = where(finite(co2),nco2)
isf6 = where(finite(sf6),nsf6)
iage = where(finite(co2) and finite(sf6),nage)
in2o_a0 = where(finite(n2o_norm_a0),nn2o_a0)
iage_a0 = where(finite(co2_a0) and finite(sf6_a0),nage_a0)
ico2_a0 = where(finite(co2_a0),nco2_a0)
isf6_a0 = where(finite(sf6_a0),nsf6_a0)
in2o_a1 = where(finite(n2o_norm_a1),nn2o_a1)
iage_a1 = where(finite(co2_a1) and finite(sf6_a1),nage_a1)
ico2_a1 = where(finite(co2_a1),nco2_a1)
isf6_a1 = where(finite(sf6_a1),nsf6_a1)
in2o_a = where(finite(n2o_norm_a),nn2o_a)
iage_a = where(finite(co2_a) and finite(sf6_a),nage_a)
ico2_a = where(finite(co2_a),nco2_a)
isf6_a = where(finite(sf6_a),nsf6_a)
in2o_b = where(finite(n2o_norm_b),nn2o_b)
ico2_b = where(finite(co2_b),nco2_b)
isf6_b = where(finite(sf6_b),nsf6_b)
iage_b = where(finite(co2_all_b) and finite(sf6_b),nage_b)

restore,dir+'Models/WACCM/FWSD_means.sav'
nyw = n_elements(wlat)
nzw = n_elements(walt)
restore,dir+'swoosh/n2o_merge.sav'
nys = n_elements(slat)
nzs = n_elements(level)
alt_s = -7.0*alog(level/1e3)
swoosh_n2o_seas = replicate(!values.f_nan,nys,nzs,12)
swoosh_n2o_seas_early = swoosh_n2o_seas
for t = 0, 11 do for z = 0, nzs-1 do for y = 0, nys-1 do begin
  ti = indgen(10)*12+t+360
  swoosh_n2o_seas[y,z,t] = mean(reform(combn2oq_sc_o3_norm[y,z,ti]),/nan)
  ti = indgen(6)*12+t+240
  swoosh_n2o_seas_early[y,z,t] = mean(reform(combn2oq_sc_o3_norm[y,z,ti]),/nan)
endfor

;  Add on WACCM seasonal cycle to bottom levels of swoosh.
yi = interpol(findgen(nyw),wlat,slat)
zi = interpol(findgen(nzw),walt,alt_s)
n2o_w = fltarr(nys,nzs,12)
n2o_w1 = fltarr(nyw,nzs,12)
for t = 0, 11 do begin
  for y = 0, nyw-1 do n2o_w1[y,*,t] = interpolate(n2o_norm_seas[0,y,*,t],zi)
  for z = 0, nzs-1 do n2o_w[*,z,t] = interpolate(n2o_w1[*,z,t],yi)
endfor

combo_n2o_seas = swoosh_n2o_seas
combo_n2o_seas_early = swoosh_n2o_seas_early
combo_n2o_seas[*,0:5,*] = n2o_w[*,0:5,*]
combo_n2o_seas_early[*,0:5,*] = n2o_w[*,0:5,*]
for z = 4, 5 do begin
  tmp1 = reform(n2o_w[*,z,*])
  tmp2 = reform(swoosh_n2o_seas[*,z,*])
  tmp3 = reform(combo_n2o_seas[*,z,*])
  gd = where(finite(tmp2),ngd)
  if z eq 4 then tmp3[gd] = 0.67*tmp1[gd] + 0.33*tmp2[gd]
  if z eq 5 then tmp3[gd] = 0.33*tmp1[gd] + 0.67*tmp2[gd]
  combo_n2o_seas[*,z,*] = tmp3
  tmp1 = reform(n2o_w[*,z,*])
  tmp2 = reform(swoosh_n2o_seas_early[*,z,*])
  tmp3 = reform(combo_n2o_seas_early[*,z,*])
  gd = where(finite(tmp2),ngd)
  if z eq 4 then tmp3[gd] = 0.67*tmp1[gd] + 0.33*tmp2[gd]
  if z eq 5 then tmp3[gd] = 0.33*tmp1[gd] + 0.67*tmp2[gd]
  combo_n2o_seas_early[*,z,*] = tmp3
endfor

;  Fill in missing lats.

nyg = 36
dy = 5.
lat_grid = findgen(nyg)*dy-85.
nz = 120
dz = 0.25
alt_grid = findgen(nz)*dz+6.
nz2 = 60
dz2 = 0.5
alt_grid2 = findgen(nz2)*dz2+6.

elat_m_adj = elat_m
elat_m_adj_a = elat_m_a
elat_m_adj_b = elat_m_b
elat_m_adj_a0 = elat_m_a0
elat_m_adj_a1 = elat_m_a1
aa = replicate(!values.f_nan,7,nyg,nz2,12)
n2o_norm_grid_mon_avg = aa
age_co2_grid_mon_avg = aa
n2o_norm_sf6_grid_mon_avg = aa
bb = replicate(!values.f_nan,3,nyg,nz,4)
combo_n2o_seas_hires = replicate(!values.f_nan,nys,nz,12)
age_co2_grid_seas = bb
n2o_norm_grid_seas = bb
age_sf6_grid_mon_avg = aa
n2o_norm_grid_mon_avg_a = aa
age_co2_grid_mon_avg_a = aa
n2o_norm_sf6_grid_mon_avg_a = aa
age_sf6_grid_mon_avg_a = aa
age_combo_grid_mon_avg = aa
n2o_norm_combo_grid_mon_avg = aa
n2o_norm_combo_grid_seas = bb
age_combo_grid_seas = bb
n2o_norm_grid_mon_avg_b = aa
age_co2_grid_mon_avg_b = aa
combo_n2o_seas_early_hires = combo_n2o_seas_hires
n2o_norm_sf6_grid_mon_avg_b = aa
n2o_norm_both_grid_mon_avg_b = aa
age_both_grid_mon_avg_b = aa
age_grid_mon_avg_b = aa
n2o_norm_opt_grid_mon_avg = aa
age_opt_grid_mon_avg = aa
age_opt_grid_mon_avg_a = aa
n2o_norm_grid_mon_avg_a0 = aa
age_co2_grid_mon_avg_a0 = aa
n2o_norm_sf6_grid_mon_avg_a0 = aa
age_sf6_grid_mon_avg_a0 = aa
age_opt_grid_mon_avg_a0 = aa
n2o_norm_grid_mon_avg_a1 = aa
age_co2_grid_mon_avg_a1 = aa
n2o_norm_sf6_grid_mon_avg_a1 = aa
age_sf6_grid_mon_avg_a1 = aa
n2o_norm_grid_seas_a1 = bb
age_opt_grid_mon_avg_a1 = aa
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
    if nchk gt 0 then elat_m_adj[in2o[gd[chk]]] = 0.

    chk = where(n2o_norm[in2o[gd]] lt min(combo_n2o_seas_hires[nys/2+1:-1,z,t]),nchk)
    if nchk gt 0 then elat_m_adj[in2o[gd[chk]]] = 85.
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

for t = 0, 11 do begin
  ii = in2o
  gd = where(months_all[ii] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm[ii[gd]]) and elat_m_adj[ii[gd]] gt 0 and lat[ii[gd]] gt 0 and alt[ii[gd]] ge alt_grid2[z]-dz2/2. and alt[ii[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm[ii[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
      ii_all = where(elat_m[ii[gd]] gt 0 and lat[ii[gd]] gt 0 and alt[ii[gd]] ge alt_grid2[z]-dz2/2. and alt[ii[gd]] lt alt_grid2[z]+dz2/2. and  elat_m[ii[gd]] ge lat_grid[y]-dy/2. and elat_m[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm[ii[gd[ii_all]]],sd=sd,/nan)
      endif
    endfor
  endif
  gd = where(months_all[ico2] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm[ico2[gd]]) and elat_m_adj[ico2[gd]] gt 0 and lat[ico2[gd]] gt 0 and alt[ico2[gd]] ge alt_grid2[z]-dz2/2. and alt[ico2[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj[ico2[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[ico2[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(dat_co2_b2.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_co2_grid_mon_avg[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
      ii_all = where(elat_m[ico2[gd]] gt 0 and lat[ico2[gd]] gt 0 and alt[ico2[gd]] ge alt_grid2[z]-dz2/2. and alt[ico2[gd]] lt alt_grid2[z]+dz2/2. and  elat_m[ico2[gd]] ge lat_grid[y]-dy/2. and elat_m[ico2[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(dat_co2_b2.age_opt_all[gd[ii_all]],sd=sd,/nan)
      endif
    endfor
  endif
  gd = where(months_all[isf6] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm[isf6[gd]]) and elat_m_adj[isf6[gd]] gt 0 and lat[isf6[gd]] gt 0 and alt[isf6[gd]] ge alt_grid2[z]-dz2/2. and alt[isf6[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj[isf6[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[isf6[gd]] lt lat_grid[y]+dy/2.,tot_i)
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
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm[iage[gd]]) and elat_m_adj[iage[gd]] gt 0 and lat[iage[gd]] gt 0 and alt[iage[gd]] ge alt_grid2[z]-dz2/2. and alt[iage[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj[iage[gd]] ge lat_grid[y]-dy/2. and elat_m_adj[iage[gd]] lt lat_grid[y]+dy/2.,tot_i)
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
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm_a[ii[gd]]) and elat_m_adj_a[ii[gd]] gt 0 and lat_a[ii[gd]] gt 0 and alt_a[ii[gd]] ge alt_grid2[z]-dz2/2. and  alt_a[ii[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a[ii[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a[ico2_a] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm_a[ico2_a[gd]]) and elat_m_adj_a[ico2_a[gd]] gt 0 and lat_a[ico2_a[gd]] gt 0 and alt_a[ico2_a[gd]] ge alt_grid2[z]-dz2/2. and  alt_a[ico2_a[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[ico2_a[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[ico2_a[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(dat_co2_a2.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_co2_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a[iage_a] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm_a[iage_a[gd]]) and elat_m_adj_a[iage_a[gd]] gt 0 and lat_a[iage_a[gd]] gt 0 and alt_a[iage_a[gd]] ge alt_grid2[z]-dz2/2. and  alt_a[iage_a[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[iage_a[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[iage_a[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a[iage_a[gd[ii_all]]],sd=sd,/nan)
        stats = moment(dat_a2.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_opt_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a[isf6_a] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(n2o_norm_a[isf6_a[gd]]) and elat_m_adj_a[isf6_a[gd]] gt 0 and lat_a[isf6_a[gd]] gt 0 and alt_a[isf6_a[gd]] ge alt_grid2[z]-dz2/2. and  alt_a[isf6_a[gd]] lt alt_grid2[z]+dz2/2. and elat_m_adj_a[isf6_a[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a[isf6_a[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a[isf6_a[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_sf6_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
        stats = moment(dat_sf6_a2.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_sf6_grid_mon_avg_a[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif

  ;  Balloon early.
  ii = in2o_b
  gd = where(months_all_b[ii] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_b[ii[gd]]) and alt_b[ii[gd]] ge alt_grid2[z]-dz2/2. and alt_b[ii[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_b[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_b[ii[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_b[ico2_b] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_b[ico2_b[gd]]) and alt_b[ico2_b[gd]] ge alt_grid2[z]-dz2/2. and alt_b[ico2_b[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_b[ico2_b[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[ico2_b[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(dat_co2_b0.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_co2_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_b[isf6_b] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_b[isf6_b[gd]]) and alt_b[isf6_b[gd]] ge alt_grid2[z]-dz2/2. and alt_b[isf6_b[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_b[isf6_b[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[isf6_b[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_b[isf6_b[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_sf6_grid_mon_avg_b[0:2,y,z,t] = [stats[0],sd,tot_i]
        stats = moment(dat_sf6_b0.age_opt_all[gd[ii_all]],sd=sd,/nan)
      endif
    endfor
  endif
  gd = where(months_all_b[iage_b] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_b[iage_b[gd]]) and alt_b[iage_b[gd]] ge alt_grid2[z]-dz2/2. and alt_b[iage_b[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_b[iage_b[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_b[iage_b[gd]] lt lat_grid[y]+dy/2.,tot_i)
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
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_a0[ii[gd]]) and finite(n2o_norm_a0[ii[gd]]) and alt_a0[ii[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[ii[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a0[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a0[ii[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_grid_mon_avg_a0[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a0[ico2_a0] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_a0[ico2_a0[gd]]) and finite(n2o_norm_a0[ico2_a0[gd]]) and alt_a0[ico2_a0[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[ico2_a0[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a0[ico2_a0[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[ico2_a0[gd]] lt lat_grid[y]+dy/2.,tot_i)
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
      endif
    endfor
  endif
  gd = where(months_all_a0[iage_a0] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(lat_a0[iage_a0[gd]] gt 0 and finite(n2o_norm_a0[iage_a0[gd]]) and alt_a0[iage_a0[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[iage_a0[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a0[iage_a0[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[iage_a0[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a0[iage_a0[gd[ii_all]]],sd=sd,/nan)
        stats = moment(datsub,sd=sd,/nan)
        age_opt_grid_mon_avg_a0[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a0[isf6_a0] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_a0[isf6_a0[gd]]) and finite(n2o_norm_a0[isf6_a0[gd]]) and alt_a0[isf6_a0[gd]] ge alt_grid2[z]-dz2/2. and alt_a0[isf6_a0[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a0[isf6_a0[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a0[isf6_a0[gd]] lt lat_grid[y]+dy/2.,tot_i)
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
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_a1[ii[gd]]) and finite(n2o_norm_a1[ii[gd]]) and alt_a1[ii[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[ii[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a1[ii[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[ii[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a1[ii[gd[ii_all]]],sd=sd,/nan)
        n2o_norm_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a1[ico2_a1] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_a1[ico2_a1[gd]]) and finite(n2o_norm_a1[ico2_a1[gd]]) and alt_a1[ico2_a1[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[ico2_a1[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a1[ico2_a1[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[ico2_a1[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(dat_co2_a1.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_co2_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a1[iage_a1] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(lat_a1[iage_a1[gd]] gt 0 and finite(n2o_norm_a1[iage_a1[gd]]) and alt_a1[iage_a1[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[iage_a1[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a1[iage_a1[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[iage_a1[gd]] lt lat_grid[y]+dy/2.,tot_i)
      if tot_i gt 1 then begin
        stats = moment(n2o_norm_a1[iage_a1[gd[ii_all]]],sd=sd,/nan)
        stats = moment(dat_a1.age_opt_all[gd[ii_all]],sd=sd,/nan)
        age_opt_grid_mon_avg_a1[0:2,y,z,t] = [stats[0],sd,tot_i]
      endif
    endfor
  endif
  gd = where(months_all_a1[isf6_a1] eq t+1,ngd)
  if ngd gt 1 then begin
    for z = 6, nz2-1 do for y = 0, nyg-1 do begin
      ii_all = where(finite(lat_a1[isf6_a1[gd]]) and finite(n2o_norm_a1[isf6_a1[gd]]) and alt_a1[isf6_a1[gd]] ge alt_grid2[z]-dz2/2. and alt_a1[isf6_a1[gd]] lt alt_grid2[z]+dz2/2. and  elat_m_adj_a1[isf6_a1[gd]] ge lat_grid[y]-dy/2. and elat_m_adj_a1[isf6_a1[gd]] lt lat_grid[y]+dy/2.,tot_i)
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
      n2o_norm_combo_grid_mon_avg[j,y,z,t] = mean([n2o_norm_grid_mon_avg[j,y,z,t],n2o_norm_sf6_grid_mon_avg[j,y,z,t], $
        n2o_norm_grid_mon_avg_a[j,y,z,t],n2o_norm_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
      age_combo_grid_mon_avg[j,y,z,t] = mean([age_co2_grid_mon_avg[j,y,z,t],age_sf6_grid_mon_avg[j,y,z,t], $
        age_co2_grid_mon_avg_a[j,y,z,t],age_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
      n2o_norm_both_grid_mon_avg_b[j,y,z,t] = mean([n2o_norm_grid_mon_avg_b[j,y,z,t],n2o_norm_sf6_grid_mon_avg_b[j,y,z,t]],/nan)
      age_both_grid_mon_avg_b[j,y,z,t] = mean([age_co2_grid_mon_avg_b[j,y,z,t],age_grid_mon_avg_b[j,y,z,t]],/nan)
    endfor
    j = 2
    n2o_norm_combo_grid_mon_avg[j,y,z,t] = total([n2o_norm_grid_mon_avg[j,y,z,t],n2o_norm_sf6_grid_mon_avg[j,y,z,t], $
      n2o_norm_grid_mon_avg_a[j,y,z,t],n2o_norm_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
    age_combo_grid_mon_avg[j,y,z,t] = total([age_co2_grid_mon_avg[j,y,z,t],age_sf6_grid_mon_avg[j,y,z,t], $
      age_co2_grid_mon_avg_a[j,y,z,t],age_sf6_grid_mon_avg_a[j,y,z,t]],/nan)
    n2o_norm_both_grid_mon_avg_b[j,y,z,t] = total([n2o_norm_grid_mon_avg_b[j,y,z,t],n2o_norm_sf6_grid_mon_avg_b[j,y,z,t]],/nan)
    age_both_grid_mon_avg_b[j,y,z,t] = total([age_co2_grid_mon_avg_b[j,y,z,t],age_grid_mon_avg_b[j,y,z,t]],/nan)
  endfor
endfor

for z = 0, nz2-1 do for y = 0, nyg-1 do for j = 0, 1 do begin
  n2o_norm_grid_seas[j,y,z,0] = mean([n2o_norm_grid_mon_avg[j,y,z,11],reform(n2o_norm_grid_mon_avg[j,y,z,0:1])],/nan)
  for s = 1, 3 do n2o_norm_grid_seas[j,y,z,s] = mean(n2o_norm_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
  age_co2_grid_seas[j,y,z,0] = mean([age_co2_grid_mon_avg[j,y,z,11],reform(age_co2_grid_mon_avg[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_co2_grid_seas[j,y,z,s] = mean(age_co2_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
  n2o_norm_combo_grid_seas[j,y,z,0] = mean([n2o_norm_combo_grid_mon_avg[j,y,z,11], $
    reform(n2o_norm_combo_grid_mon_avg[j,y,z,0:1])],/nan)
  for s = 1, 3 do n2o_norm_combo_grid_seas[j,y,z,s] = mean(n2o_norm_combo_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
  age_combo_grid_seas[j,y,z,0] = mean([age_combo_grid_mon_avg[j,y,z,11],reform(age_combo_grid_mon_avg[j,y,z,0:1])],/nan)
  for s = 1, 3 do age_combo_grid_seas[j,y,z,s] = mean(age_combo_grid_mon_avg[j,y,z,indgen(3)+3*s-1],/nan)
  n2o_norm_grid_seas_a1[j,y,z,0] = mean([n2o_norm_grid_mon_avg_a1[j,y,z,11],reform(n2o_norm_grid_mon_avg_a1[j,y,z,0:1])],/nan)
  for s = 1, 3 do n2o_norm_grid_seas_a1[j,y,z,s] = mean(n2o_norm_grid_mon_avg_a1[j,y,z,indgen(3)+3*s-1],/nan)
endfor

;  Fill in missing locations.
for s = 0, 3 do for z = 1, nz2-2 do for y = 1, nyg-2 do for j = 0, 1 do begin
  if ~finite(n2o_norm_combo_grid_seas[j,y,z,s]) and finite(n2o_norm_combo_grid_seas[j,y-1,z, $
    s]) and finite(n2o_norm_combo_grid_seas[j,y+1,z,s]) then n2o_norm_combo_grid_seas[j,y,z, $
    s] = 0.5 * n2o_norm_combo_grid_seas[j,y-1,z,s] + 0.5 * n2o_norm_combo_grid_seas[j,y+1,z,s]
  if ~finite(n2o_norm_combo_grid_seas[j,y,z,s]) and finite(n2o_norm_combo_grid_seas[j,y,z-1, $
    s]) and finite(n2o_norm_combo_grid_seas[j,y,z+1,s]) then n2o_norm_combo_grid_seas[j,y,z, $
    s] = 0.5 * n2o_norm_combo_grid_seas[j,y,z-1,s] + 0.5 * n2o_norm_combo_grid_seas[j,y,z+1,s]
  if ~finite(age_combo_grid_seas[j,y,z,s]) and finite(age_combo_grid_seas[j,y-1,z,s]) and finite(age_combo_grid_seas[j,y+1,z, $
    s]) then age_combo_grid_seas[j,y,z,s] = 0.5 * age_combo_grid_seas[j,y-1,z,s] + 0.5 * age_combo_grid_seas[j,y+1,z,s]
  if ~finite(age_combo_grid_seas[j,y,z,s]) and finite(age_combo_grid_seas[j,y,z-1,s]) and finite(age_combo_grid_seas[j,y,z+1, $
    s]) then age_combo_grid_seas[j,y,z,s] = 0.5 * age_combo_grid_seas[j,y,z-1,s] + 0.5 * age_combo_grid_seas[j,y,z+1,s]
endfor

n2o_norm_grid_seas_late = n2o_norm_grid_seas
age_combo_grid_seas_late = age_combo_grid_seas
n2o_norm_combo_grid_seas_late = n2o_norm_combo_grid_seas

cc = replicate(!values.f_nan,3,nyg,nz2)
age_combo_grid_late = cc
n2o_norm_combo_grid_late = cc
n2o_norm_combo_grid_minmax_late = cc
for z = 1, nz2-2 do for y = 0, nyg-1 do begin
  gd = where(finite(age_combo_grid_seas_late[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(age_combo_grid_seas_late[0,y,z,gd],sdev=sdev)
    age_combo_grid_late[*,y,z] = [stats[0],sdev,ngd]
  endif
  gd = where(finite(n2o_norm_combo_grid_seas_late[0,y,z,*]),ngd)
  if ngd ge 1 then begin
    stats = moment(n2o_norm_combo_grid_seas_late[0,y,z,gd],sdev=sdev)
    n2o_norm_combo_grid_late[*,y,z] = [stats[0],sdev,ngd]
    n2o_norm_combo_grid_minmax_late[0:1,y, $
      z] = [min(n2o_norm_combo_grid_seas_late[0,y,z,gd]-n2o_norm_combo_grid_seas_late[1,y,z,gd]), $
      max(n2o_norm_combo_grid_seas_late[0,y,z,gd]+n2o_norm_combo_grid_seas_late[1,y,z,gd])]
  endif
endfor

n2o_norm_opt_grid_mon_avg_b_late = n2o_norm_opt_grid_mon_avg
n2o_norm_grid_mon_avg_b_late = n2o_norm_grid_mon_avg
n2o_norm_sf6_grid_mon_avg_b_late = n2o_norm_sf6_grid_mon_avg
n2o_norm_grid_mon_avg_a_late = n2o_norm_grid_mon_avg_a
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
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_aircore_aircraft_late.nc', $
  ['n2o_norm_grid_seas_late','n2o_norm_combo_grid_seas_late','age_co2_grid_seas','age_combo_grid_seas_late', $
  'age_combo_grid_late','n2o_norm_combo_grid_late','n2o_norm_combo_grid_minmax_late','n2o_norm_combo_grid_mon_avg', $
  'n2o_norm_grid_mon_avg_b_late','age_co2_grid_mon_avg_b_late','n2o_norm_sf6_grid_mon_avg_b_late','age_combo_grid_mon_avg', $
  'age_sf6_grid_mon_avg_b_late','age_opt_grid_mon_avg_b_late','age_co2_grid_mon_avg_a_late','n2o_norm_grid_mon_avg_a_late', $
  'age_sf6_grid_mon_avg_a_late','n2o_norm_opt_grid_mon_avg_b_late','age_opt_grid_mon_avg_a_late']
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_1990s_balloon_aircraft.nc', $
  ['n2o_norm_both_grid_mon_avg_b_early','age_co2_grid_mon_avg_b_early','age_grid_mon_avg_b_early', $
  'age_both_grid_mon_avg_b_early','n2o_norm_grid_mon_avg_a_early','n2o_norm_sf6_grid_mon_avg_a_early', $
  'age_co2_grid_mon_avg_a_early','age_sf6_grid_mon_avg_a_early','age_opt_grid_mon_avg_a_early','combo_n2o_seas_hires', $
  'combo_n2o_seas_early_hires','n2o_norm_grid_mon_avg_a_mid','n2o_norm_sf6_grid_mon_avg_a_mid','age_co2_grid_mon_avg_a_mid', $
  'age_sf6_grid_mon_avg_a_mid','age_opt_grid_mon_avg_a_mid','n2o_norm_grid_seas_a1']
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'equivalent_latitude_adjustments.nc', $
  ['elat_m_adj_a_early','elat_m_adj_b_early','elat_m_adj_a_late','elat_m_adj_b_late','elat_m_adj_a_mid']
end

