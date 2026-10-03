;+
;  Step 3: gridded time series, mean age vs tracer relationships and trend profiles.
;
;  Builds latitude-adjusted NH mean age and N2O time series in altitude and
;  normalized N2O/CH4 bins (Supplementary Note 3), the 1990s mean age vs N2O,
;  CH4 and CFC-12 relationships, and the altitude trend profiles.
;
;  Reads:  step 1 and 2 outputs; merged aircraft, balloon and AirCore files and their
;          normalized-tracer bins; SWOOSH N2O; MERRA-2 tropopause; NOAA surface N2O and CFC-12
;  Writes: <output>/intermediate/gridded_time_series.nc, mean_age_vs_tracers_1990s.nc,
;          trend_profiles.nc; <output>/figure_data/supp_fig3.nc, supp_fig5.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step3_time_series_and_trends

if n_elements(stage) eq 0 then stage = ''

dir = mean_age_data_dir()

;  Tropopause data.
ncdf_get,dir+'Reanalysis/MERRA2/tp.monmean.zm.nc',['time','tpp','tpt','tpz','lat'],tp,/quiet
ny_tp = tp['lat','dim_sizes']
lat_tp = tp['lat','value']
tpause = fltarr(ny_tp)
tpause_alt = tpause
for y = 0, ny_tp[0]-1 do begin
  tpause_alt[y] = mean(tp['tpz','value',y,*])
endfor

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
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
ti = interpol(indgen(n_elements(f12_date)),f12_date,times_all)
f12_surface = interpolate(f12_global,ti)
f12_norm = f12/f12_surface
co2_a0 = co2
sf6_a0 = sf6
alt_a0 = alt
n2o_a0 = n2o
n2o_norm_a0 = n2o_norm
times_all_a0 = times_all
ch4_all_norm_a0 = ch4_all_norm
f12_norm_a0 = f12_norm

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_mid.sav'
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
co2_a1 = co2
sf6_a1 = sf6
alt_a1 = alt
n2o_a1 = n2o
n2o_norm_a1 = n2o_norm
times_all_a1 = times_all
ch4_all_norm_a1 = ch4_all_norm

restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_late.sav'
ti = interpol(indgen(n_elements(n2o_date)),n2o_date,times_all)
n2o_surface = interpolate(n2o_global,ti)
n2o_norm = n2o/n2o_surface
co2_a = co2
sf6_a = sf6
alt_a = alt
n2o_a = n2o
n2o_norm_a = n2o_norm
times_all_a = times_all
ch4_all_norm_a = ch4_all_norm

restore,dir+'Balloon/OMS_common_gc_merge.sav'
co2_b = co2
sf6_b = sf6
alt_b = alt
n2o_b = n2o
n2o_norm_b = n2o_norm
times_all_b = times_all
co2_all_b = co2_all
ch4_all_norm_b = ch4_all_norm
f12_norm_b = f12_norm
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
ich4 = where(finite(ch4_all_norm),nch4)
ico2 = where(finite(co2),nco2)
isf6 = where(finite(sf6),nsf6)
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
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_early_bins.sav'
bins_co2_a0 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_early_bins.sav'
bins_sf6_a0 = bins_sf6
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_mid_bins.sav'
bins_co2_a1 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_mid_bins.sav'
bins_sf6_a1 = bins_sf6
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
restore,dir+'swoosh/n2o_merge.sav'
nys = n_elements(slat)

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
nz3 = 30
dz3 = 1.0
alt_grid3 = findgen(nz3)*dz3+6.0
if stage eq 'elats' then return

ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'equivalent_latitude_adjustments.nc', $
  ['elat_m_adj_a_early','elat_m_adj_a_late','elat_m_adj_a_mid','elat_m_adj_b_early','elat_m_adj_b_late']
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
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'seasonal_grids_1990s.nc', $
  ['age_co2_old_grid_seas_sm','age_combo_grid_early','age_combo_grid_seas_sm_early','n2o_norm_grid_early', $
  'n2o_norm_grid_seas_both']
z = 28

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

fff = replicate(!values.f_nan,6,nyg,nz2,nt)
ggg = replicate(!values.f_nan,6,nz2,nt)
hhh = replicate(!values.f_nan,nz2,nt)
jjj = replicate(!values.f_nan,3,nn,nt)
mmm = replicate(!values.f_nan,6,nn,nt)
age_grid_adj_tseries = fff
age_co2_grid_adj_tseries = fff
age_sf6_grid_adj_tseries = fff
age_opt_grid_adj_tseries = fff
age_co2_b_grid_adj_tseries = fff
age_sf6_b_grid_adj_tseries = fff
age_opt_b_grid_adj_tseries = fff
n2o_a_grid_adj_tseries = fff
n2o_b_grid_adj_tseries = fff
n2o_grid_adj_tseries = fff
n2o_on_age_adj_tseries = jjj
n2o_on_age_b_adj_tseries = jjj
n2o_on_age_both_adj_tseries = jjj
age_on_n2o_adj_tseries = mmm
age_sf6_on_n2o_b_adj_tseries = mmm
age_sf6_on_n2o_adj_tseries = mmm
age_co2_on_n2o_b_adj_tseries = mmm
age_co2_on_n2o_adj_tseries = mmm
ch4_a_grid_adj_tseries = fff
ch4_b_grid_adj_tseries = fff
ch4_grid_adj_tseries = fff
age_sf6_on_ch4_b_adj_tseries = mmm
age_sf6_on_ch4_adj_tseries = mmm
age_co2_on_ch4_b_adj_tseries = mmm
age_co2_on_ch4_adj_tseries = mmm
age_on_ch4_adj_tseries = mmm
age_on_n2o_e_adj_tseries = mmm
age_co2_on_n2o_e_adj_tseries = mmm
age_sf6_on_n2o_e_adj_tseries = mmm
age_co2_on_n2o_e_b_adj_tseries = mmm
age_sf6_on_n2o_e_b_adj_tseries = mmm
age_co2_on_ch4_e_adj_tseries = mmm
age_sf6_on_ch4_e_adj_tseries = mmm
age_co2_on_ch4_e_b_adj_tseries = mmm
age_sf6_on_ch4_e_b_adj_tseries = mmm
age_on_ch4_e_adj_tseries = mmm

dy_int = 2*dy/3.
dz_int = 2*dz2/3.

for t = 0, nt-1 do begin

  ;  Early aircraft.   -------------------------------------------------------------------------------------------------------------------
  chk = where(times_all_a0[in2o_a0] ge years[t]-dt/2. and times_all_a0[in2o_a0] lt years[t]+dt/2.,ni_all)
  if ni_all ge 1 then begin
    samp_type = 1
    for z = 0, nz2-1 do for y = 0, nyg-1 do begin
      gd = where(alt_a0[in2o_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[in2o_a0[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_early[in2o_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[in2o_a0[chk]] lt lat_grid[y]+dy_int,ngd)
      if ngd ge 1 then begin
        stats = moment(n2o_norm_a0[in2o_a0[chk[gd]]],sdev=sdev,/nan)
        n2o_a_grid_adj_tseries[0:2,y,z,t] = [stats[0],sdev,samp_type]
      endif
    endfor
  endif
  chk = where(times_all_a0[ich4_a0] ge years[t]-dt/2. and times_all_a0[ich4_a0] lt years[t]+dt/2.,ni_all)
  if ni_all ge 1 then begin
    samp_type = 1
    for z = 0, nz2-1 do for y = 0, nyg-1 do begin
      gd = where(alt_a0[ich4_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[ich4_a0[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_early[ich4_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[ich4_a0[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_a0[ich4_a0[chk]] gt 0 and  ch4_all_norm_a0[ich4_a0[chk]] lt 1.1,ngd)
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
      gd = where(alt_a0[ico2_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[ico2_a0[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_early[ico2_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[ico2_a0[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(ch4_all_norm_a0[ico2_a0[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a0[ico2_a0[chk]] lt norm_grid[n]+dn/2. and  elat_m_adj_a_early[ico2_a0[chk]] ge 30,ngd)
      if ngd ge 1 then begin
        age_err = 0.5 * (dat_co2_a0.age_range_all[1,chk[gd]] - dat_co2_a0.age_range_all[0,chk[gd]])
        age_weights = 1. / age_err^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * dat_co2_a0.age_opt_all[chk[gd]]) / total(age_weights)
        stats = moment(dat_co2_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
        age_co2_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
      endif
      gd = where(bins_co2_a0.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_co2_a0.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and  elat_m_adj_a_early[ico2_a0[chk]] ge 30,ngd)
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
      gd = where(dat_co2_a0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_a0.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_a_early[ico2_a0[chk]] ge 30 and elat_m_adj_a_early[ico2_a0[chk]] le 60,ngd)
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
      gd = where(alt_a0[isf6_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[isf6_a0[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_early[isf6_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[isf6_a0[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(ch4_all_norm_a0[isf6_a0[chk]] ge norm_grid[n]-dn/2. and ch4_all_norm_a0[isf6_a0[chk]] lt norm_grid[n]+dn/2. and  elat_m_adj_a_early[isf6_a0[chk]] ge 30,ngd)
      if ngd ge 1 and years[t] lt 1999 then begin
        age_err = 0.5 * (dat_sf6_a0.age_range_all[1,chk[gd]] - dat_sf6_a0.age_range_all[0,chk[gd]])
        age_weights = 1. / age_err^2
        err = 1. / sqrt(total(age_weights))
        wmean = total(age_weights * dat_sf6_a0.age_opt_all[chk[gd]]) / total(age_weights)
        stats = moment(dat_sf6_a0.age_opt_all[chk[gd]],sdev=sdev,/nan)
        age_sf6_on_ch4_adj_tseries[*,n,t] = [wmean,err,stats[0],sdev,ngd,samp_type]
      endif
      gd = where(bins_sf6_a0.ch4_norm_e[chk] ge norm_grid[n]-dn/2. and bins_sf6_a0.ch4_norm_e[chk] lt norm_grid[n]+dn/2. and  elat_m_adj_a_early[isf6_a0[chk]] ge 30,ngd)
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
      gd = where(dat_sf6_a0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_a0.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_a_early[isf6_a0[chk]] ge 30 and elat_m_adj_a_early[isf6_a0[chk]] le 60,ngd)
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
      gd = where(alt_a0[iage_a0[chk]] ge alt_grid2[z] - dz_int and alt_a0[iage_a0[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_early[iage_a0[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_early[iage_a0[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(alt_b[in2o_b[chk]] ge alt_grid2[z] - dz_int and alt_b[in2o_b[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_early[in2o_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[in2o_b[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(alt_b[ich4_b[chk]] ge alt_grid2[z] - dz_int and alt_b[ich4_b[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_early[ich4_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[ich4_b[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_b[ich4_b[chk]] gt 0 and  ch4_all_norm_b[ich4_b[chk]] lt 1.1,ngd)
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
      gd = where(alt_b[ico2_b[chk]] ge alt_grid2[z] - dz_int and alt_b[ico2_b[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_early[ico2_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[ico2_b[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_co2_b0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_b0.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_b_early[ico2_b[chk]] ge 30 and elat_m_adj_b_early[ico2_b[chk]] le 60,ngd)
      if ngd ge 1 then begin
        stats = moment(n2o_b[ico2_b[chk[gd]]],sdev=sdev,/nan)
        n2o_on_age_b_adj_tseries[0:2,a,t] = [stats[0],sdev,samp_type]
      endif
    endfor
  endif
  chk = where(times_all_b[isf6_b] ge years[t]-dt/2. and times_all_b[isf6_b] lt years[t]+dt/2. and times_all_b[isf6_b] lt 2002, $
    ni_all)
  if ni_all ge 1 then begin
    samp_type = 2
    for z = 0, nz2-1 do for y = 0, nyg-1 do begin
      gd = where(alt_b[isf6_b[chk]] ge alt_grid2[z] - dz_int and alt_b[isf6_b[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_early[isf6_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[isf6_b[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_sf6_b0.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_b0.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_b_early[isf6_b[chk]] ge 30 and elat_m_adj_b_early[isf6_b[chk]] le 60,ngd)
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
      gd = where(alt_b[iage_b[chk]] ge alt_grid2[z] - dz_int and alt_b[iage_b[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_early[iage_b[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_early[iage_b[chk]] lt lat_grid[y]+dy_int,ngd)
      if ngd ge 1 then begin
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
      gd = where(alt_a1[in2o_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[in2o_a1[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_mid[in2o_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[in2o_a1[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(alt_a1[ich4_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[ich4_a1[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_mid[ich4_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[ich4_a1[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_a1[ich4_a1[chk]] gt 0 and  ch4_all_norm_a1[ich4_a1[chk]] lt 1.1,ngd)
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
      gd = where(alt_a1[ico2_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[ico2_a1[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_mid[ico2_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[ico2_a1[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_co2_a1.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_a1.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_a_mid[ico2_a1[chk]] ge 30 and elat_m_adj_a_mid[ico2_a1[chk]] le 60,ngd)
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
      gd = where(alt_a1[isf6_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[isf6_a1[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_mid[isf6_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[isf6_a1[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_sf6_a1.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_a1.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_a_mid[isf6_a1[chk]] ge 30 and elat_m_adj_a_mid[isf6_a1[chk]] le 60,ngd)
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
      gd = where(alt_a1[iage_a1[chk]] ge alt_grid2[z] - dz_int and alt_a1[iage_a1[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_mid[iage_a1[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_mid[iage_a1[chk]] lt lat_grid[y]+dy_int,ngd)
      if ngd ge 1 then begin
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
      gd = where(alt_a[in2o_a[chk]] ge alt_grid2[z] - dz_int and alt_a[in2o_a[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_late[in2o_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[in2o_a[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(alt_a[ich4_a[chk]] ge alt_grid2[z] - dz_int and alt_a[ich4_a[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_late[ich4_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[ich4_a[chk]] lt lat_grid[y]+dy_int and ch4_all_norm_a[ich4_a[chk]] gt 0 and  ch4_all_norm_a[ich4_a[chk]] lt 1.1,ngd)
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
      gd = where(alt_a[ico2_a[chk]] ge alt_grid2[z] - dz_int and alt_a[ico2_a[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_late[ico2_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[ico2_a[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_co2_a2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_a2.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_a_late[ico2_a[chk]] ge 30 and elat_m_adj_a_late[ico2_a[chk]] le 60,ngd)
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
      gd = where(alt_a[isf6_a[chk]] ge alt_grid2[z] - dz_int and alt_a[isf6_a[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_late[isf6_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[isf6_a[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_sf6_a2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_a2.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_a_late[isf6_a[chk]] ge 30 and elat_m_adj_a_late[isf6_a[chk]] le 60,ngd)
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
      gd = where(alt_a[iage_a[chk]] ge alt_grid2[z] - dz_int and alt_a[iage_a[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_a_late[iage_a[chk]] ge lat_grid[y]-dy_int and elat_m_adj_a_late[iage_a[chk]] lt lat_grid[y]+dy_int,ngd)
      if ngd ge 1 then begin
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
      gd = where(alt[in2o[chk]] ge alt_grid2[z] - dz_int and alt[in2o[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_late[in2o[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[in2o[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(alt[ich4[chk]] ge alt_grid2[z] - dz_int and alt[ich4[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_late[ich4[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[ich4[chk]] lt lat_grid[y]+dy_int and ch4_all_norm[ich4[chk]] gt 0 and  ch4_all_norm[ich4[chk]] lt 1.1,ngd)
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
      gd = where(alt[ico2[chk]] ge alt_grid2[z] - dz_int and alt[ico2[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_late[ico2[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[ico2[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_co2_b2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_co2_b2.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_b_late[ico2[chk]] ge 30 and elat_m_adj_b_late[ico2[chk]] le 50,ngd)
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
      gd = where(alt[isf6[chk]] ge alt_grid2[z] -dz_int and alt[isf6[chk]] lt alt_grid2[z] + dz_int and  elat_m_adj_b_late[isf6[chk]] ge lat_grid[y]-dy_int and elat_m_adj_b_late[isf6[chk]] lt lat_grid[y]+dy_int,ngd)
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
      gd = where(dat_sf6_b2.age_opt_all[chk] ge age_grid1[a]-da1/2. and dat_sf6_b2.age_opt_all[chk] lt age_grid1[a]+da1/2. and  elat_m_adj_b_late[isf6[chk]] ge 30 and elat_m_adj_b_late[isf6[chk]] le 50,ngd)
      if ngd ge 1 then begin
        stats = moment(n2o[isf6[chk[gd]]],sdev=sdev,/nan)
        n2o_on_age_b_adj_tseries[0:2,a,t] = [stats[0],sdev,samp_type]
      endif
    endfor
  endif

  for z = 0, nz2-1 do for y = 0, nyg-1 do begin
    for j = 0, 2 do n2o_grid_adj_tseries[j,y,z,t] = mean([n2o_a_grid_adj_tseries[j,y,z,t],n2o_b_grid_adj_tseries[j,y,z,t]],/nan)
    for j = 0, 2 do ch4_grid_adj_tseries[j,y,z,t] = mean([ch4_a_grid_adj_tseries[j,y,z,t],ch4_b_grid_adj_tseries[j,y,z,t]],/nan)
    age1 = [age_co2_grid_adj_tseries[0,y,z,t],age_sf6_grid_adj_tseries[0,y,z,t],age_co2_b_grid_adj_tseries[0,y,z,t], $
      age_sf6_b_grid_adj_tseries[0,y,z,t]]
    err1 = [age_co2_grid_adj_tseries[1,y,z,t],age_sf6_grid_adj_tseries[1,y,z,t],age_co2_b_grid_adj_tseries[1,y,z,t], $
      age_sf6_b_grid_adj_tseries[1,y,z,t]]
    samps = [age_co2_grid_adj_tseries[5,y,z,t],age_sf6_grid_adj_tseries[5,y,z,t],age_co2_b_grid_adj_tseries[5,y,z,t], $
      age_sf6_b_grid_adj_tseries[5,y,z,t]]
    gd = where(finite(age1),ngd)
    if ngd gt 0 then begin
      age_weights = 1. / err1[gd]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age1[gd]) / total(age_weights)
      age_grid_adj_tseries[0:1,y,z,t] = [wmean,err]
      age_grid_adj_tseries[5,y,z,t] = samps[gd[0]]
    endif
    for j = 2, 3 do age_grid_adj_tseries[j,y,z, $
      t] = mean([age_co2_grid_adj_tseries[j,y,z,t],age_sf6_grid_adj_tseries[j,y,z,t],  age_co2_b_grid_adj_tseries[j,y,z,t], $
      age_sf6_b_grid_adj_tseries[j,y,z,t]],/nan)
    age_grid_adj_tseries[4,y,z,t] = total([age_co2_grid_adj_tseries[4,y,z,t],age_sf6_grid_adj_tseries[4,y,z,t], $
      age_co2_b_grid_adj_tseries[4,y,z,t],  age_sf6_b_grid_adj_tseries[4,y,z,t]],/nan)
  endfor
  for z = 0, nn-1 do begin
    age1 = [age_co2_on_n2o_adj_tseries[0,z,t],age_sf6_on_n2o_adj_tseries[0,z,t],age_co2_on_n2o_b_adj_tseries[0,z,t], $
      age_sf6_on_n2o_b_adj_tseries[0,z,t]]
    err1 = [age_co2_on_n2o_adj_tseries[1,z,t],age_sf6_on_n2o_adj_tseries[1,z,t],age_co2_on_n2o_b_adj_tseries[1,z,t], $
      age_sf6_on_n2o_b_adj_tseries[1,z,t]]
    samps = [age_co2_on_n2o_adj_tseries[5,z,t],age_sf6_on_n2o_adj_tseries[5,z,t],age_co2_on_n2o_b_adj_tseries[5,z,t], $
      age_sf6_on_n2o_b_adj_tseries[5,z,t]]
    gd = where(finite(age1),ngd)
    if ngd gt 0 then begin
      age_weights = 1. / err1[gd]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age1[gd]) / total(age_weights)
      age_on_n2o_adj_tseries[0:1,z,t] = [wmean,err]
      age_on_n2o_adj_tseries[5,z,t] = samps[gd[0]]
    endif
    for j = 2, 3 do age_on_n2o_adj_tseries[j,z, $
      t] = mean([age_co2_on_n2o_adj_tseries[j,z,t],age_sf6_on_n2o_adj_tseries[j,z,t],  age_co2_on_n2o_b_adj_tseries[j,z,t], $
      age_sf6_on_n2o_b_adj_tseries[j,z,t]],/nan)
    age_on_n2o_adj_tseries[4,z,t] = total([age_co2_on_n2o_adj_tseries[4,z,t],age_sf6_on_n2o_adj_tseries[4,z,t], $
      age_co2_on_n2o_b_adj_tseries[4,z,t],  age_sf6_on_n2o_b_adj_tseries[4,z,t]],/nan)

    age1 = [age_co2_on_n2o_e_adj_tseries[0,z,t],age_sf6_on_n2o_e_adj_tseries[0,z,t],age_co2_on_n2o_e_b_adj_tseries[0,z,t], $
      age_sf6_on_n2o_e_b_adj_tseries[0,z,t]]
    err1 = [age_co2_on_n2o_e_adj_tseries[1,z,t],age_sf6_on_n2o_e_adj_tseries[1,z,t],age_co2_on_n2o_e_b_adj_tseries[1,z,t], $
      age_sf6_on_n2o_e_b_adj_tseries[1,z,t]]
    samps = [age_co2_on_n2o_e_adj_tseries[5,z,t],age_sf6_on_n2o_e_adj_tseries[5,z,t],age_co2_on_n2o_e_b_adj_tseries[5,z,t], $
      age_sf6_on_n2o_e_b_adj_tseries[5,z,t]]
    gd = where(finite(age1),ngd)
    if ngd gt 0 then begin
      age_weights = 1. / err1[gd]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age1[gd]) / total(age_weights)
      age_on_n2o_e_adj_tseries[0:1,z,t] = [wmean,err]
      age_on_n2o_e_adj_tseries[5,z,t] = samps[gd[0]]
    endif
    for j = 2, 3 do age_on_n2o_e_adj_tseries[j,z, $
      t] = mean([age_co2_on_n2o_e_adj_tseries[j,z,t],age_sf6_on_n2o_e_adj_tseries[j,z,t], $
      age_co2_on_n2o_e_b_adj_tseries[j,z,t],age_sf6_on_n2o_e_b_adj_tseries[j,z,t]],/nan)
    age_on_n2o_e_adj_tseries[4,z,t] = total([age_co2_on_n2o_e_adj_tseries[4,z,t],age_sf6_on_n2o_e_adj_tseries[4,z,t], $
      age_co2_on_n2o_e_b_adj_tseries[4,z,t],  age_sf6_on_n2o_e_b_adj_tseries[4,z,t]],/nan)

    age1 = [age_co2_on_ch4_adj_tseries[0,z,t],age_sf6_on_ch4_adj_tseries[0,z,t],age_co2_on_ch4_b_adj_tseries[0,z,t], $
      age_sf6_on_ch4_b_adj_tseries[0,z,t]]
    err1 = [age_co2_on_ch4_adj_tseries[1,z,t],age_sf6_on_ch4_adj_tseries[1,z,t],age_co2_on_ch4_b_adj_tseries[1,z,t], $
      age_sf6_on_ch4_b_adj_tseries[1,z,t]]
    samps = [age_co2_on_ch4_adj_tseries[5,z,t],age_sf6_on_ch4_adj_tseries[5,z,t],age_co2_on_ch4_b_adj_tseries[5,z,t], $
      age_sf6_on_ch4_b_adj_tseries[5,z,t]]
    gd = where(finite(age1),ngd)
    if ngd gt 0 then begin
      age_weights = 1. / err1[gd]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age1[gd]) / total(age_weights)
      age_on_ch4_adj_tseries[0:1,z,t] = [wmean,err]
      age_on_ch4_adj_tseries[5,z,t] = samps[gd[0]]
    endif
    for j = 2, 3 do age_on_ch4_adj_tseries[j,z, $
      t] = mean([age_co2_on_ch4_adj_tseries[j,z,t],age_sf6_on_ch4_adj_tseries[j,z,t],  age_co2_on_ch4_b_adj_tseries[j,z,t], $
      age_sf6_on_ch4_b_adj_tseries[j,z,t]],/nan)
    age_on_ch4_adj_tseries[4,z,t] = total([age_co2_on_ch4_adj_tseries[4,z,t],age_sf6_on_ch4_adj_tseries[4,z,t], $
      age_co2_on_ch4_b_adj_tseries[4,z,t],  age_sf6_on_ch4_b_adj_tseries[4,z,t]],/nan)

    age1 = [age_co2_on_ch4_e_adj_tseries[0,z,t],age_sf6_on_ch4_e_adj_tseries[0,z,t],age_co2_on_ch4_e_b_adj_tseries[0,z,t], $
      age_sf6_on_ch4_e_b_adj_tseries[0,z,t]]
    err1 = [age_co2_on_ch4_e_adj_tseries[1,z,t],age_sf6_on_ch4_e_adj_tseries[1,z,t],age_co2_on_ch4_e_b_adj_tseries[1,z,t], $
      age_sf6_on_ch4_e_b_adj_tseries[1,z,t]]
    samps = [age_co2_on_ch4_e_adj_tseries[5,z,t],age_sf6_on_ch4_e_adj_tseries[5,z,t],age_co2_on_ch4_e_b_adj_tseries[5,z,t], $
      age_sf6_on_ch4_e_b_adj_tseries[5,z,t]]
    gd = where(finite(age1),ngd)
    if ngd gt 0 then begin
      age_weights = 1. / err1[gd]^2
      err = 1. / sqrt(total(age_weights))
      wmean = total(age_weights * age1[gd]) / total(age_weights)
      age_on_ch4_e_adj_tseries[0:1,z,t] = [wmean,err]
      age_on_ch4_e_adj_tseries[5,z,t] = samps[gd[0]]
    endif
    for j = 2, 3 do age_on_ch4_e_adj_tseries[j,z, $
      t] = mean([age_co2_on_ch4_e_adj_tseries[j,z,t],age_sf6_on_ch4_e_adj_tseries[j,z,t], $
      age_co2_on_ch4_e_b_adj_tseries[j,z,t],age_sf6_on_ch4_e_b_adj_tseries[j,z,t]],/nan)
    age_on_ch4_e_adj_tseries[4,z,t] = total([age_co2_on_ch4_e_adj_tseries[4,z,t],age_sf6_on_ch4_e_adj_tseries[4,z,t], $
      age_co2_on_ch4_e_b_adj_tseries[4,z,t],  age_sf6_on_ch4_e_b_adj_tseries[4,z,t]],/nan)
  endfor
  for a = 0, nn-1 do begin
    for j = 0, 2 do n2o_on_age_both_adj_tseries[j,a,t] = mean([n2o_on_age_adj_tseries[j,a,t],n2o_on_age_b_adj_tseries[j,a,t]], $
      /nan)
  endfor
endfor
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_time_series.nc', $
  ['n2o_a_grid_adj_tseries','age_co2_grid_adj_tseries','age_sf6_grid_adj_tseries','age_opt_grid_adj_tseries', $
  'n2o_b_grid_adj_tseries','age_co2_b_grid_adj_tseries','age_sf6_b_grid_adj_tseries','age_opt_b_grid_adj_tseries', $
  'age_grid_adj_tseries','n2o_grid_adj_tseries','n2o_on_age_adj_tseries','n2o_on_age_b_adj_tseries', $
  'n2o_on_age_both_adj_tseries','age_on_n2o_adj_tseries','ch4_grid_adj_tseries','age_on_ch4_adj_tseries', $
  'age_on_n2o_e_adj_tseries','age_on_ch4_e_adj_tseries']
xxx = replicate(!values.f_nan,nz3,2)
age_nh_lat_adj_trends_coarse = xxx
age_nh_lat_adj_trends_coarse_sigma = xxx
n2o_nh_lat_adj_trends_coarse = xxx
n2o_nh_lat_adj_trends_coarse_sigma = xxx
age_nh_trends_coarse = xxx
age_nh_trends_coarse_sigma = xxx
ch4_nh_trends_coarse = xxx
ch4_nh_trends_coarse_sigma = xxx
age_nh_lat_adj_trends_per_coarse = xxx
age_nh_trends_per_coarse = xxx
age_nh_trends_per_coarse_sigma = xxx
age_nh_lat_adj_trends_per_coarse_sigma = xxx
age_nh_lat_adj_trends_ac_coarse_sigma = xxx
age_nh_lat_adj_trends_ac_coarse = xxx
n2o_nh_lat_adj_trends_per_coarse = xxx
n2o_nh_lat_adj_trends_per_coarse_sigma = xxx

yi = interpol(findgen(ny_tp),lat_tp,lat_grid)
tp_alt_grid = interpolate(tpause_alt,yi)

;  Make NH midlat avg time series.
lat_avg_n2o_nh_adj_tseries = hhh
n2o_nh_adj_tseries = ggg
nh_adj_tseries_mon = hhh
n2o_nh_adj_tseries_lat_adj = ggg
lat_avg_nh_adj_tseries = hhh
age_nh_adj_tseries = ggg
nh_adj_tseries_seas = hhh
age_nh_adj_tseries_lat_adj = ggg
ch4_nh_adj_tseries = ggg
age_nhe_adj_tseries = ggg

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
    gd = where(finite(n2o_time_series_sub) and finite(age_time_series_sub) and n2o_time_series_sub ge norm_grid[n]-dn/2. and  n2o_time_series_sub le norm_grid[n]+dn/2.,ngd)
    if ngd gt 0 then begin
      stats = moment(age_time_series_sub[gd],sdev=sdev,/nan)
    endif
  endfor
endfor

;  Make NH extratropical avg mean age vs tracer relationships.
nn2 = 50
dn2 = 0.02
norm_grid2 = findgen(nn2)*dn2+0.01
na = 81
da = 0.1
age_grid = findgen(na)*da
cent_interval = 0.95
ppp = replicate(!values.f_nan,nn2,5)
mean_age_sf6_on_n2o_a0 = ppp
mean_age_sf6_on_n2o_b0 = ppp
mean_age_all_on_n2o_a0 = ppp
mean_age_all_on_n2o_b0 = ppp
mean_age_sf6_on_ch4_a0 = ppp
mean_age_all_on_n2o_e_a0 = ppp
mean_age_all_on_ch4_a0 = ppp
mean_age_all_on_ch4_b0 = ppp
mean_age_sf6_on_f12_a0 = ppp
mean_age_sf6_on_f12_b0 = ppp
mean_age_all_on_f12_a0 = ppp
mean_age_all_on_f12_b0 = ppp
mean_age_all_on_n2o_a2 = ppp
ss = replicate(!values.f_nan,na,nn2)
mean_age_all_on_n2o_hist = ss
mean_age_all_on_ch4_hist = ss
mean_age_all_on_f12_hist = ss
mean_age_all_on_n2o_hist_cent = replicate(!values.f_nan,2,nn2)
mean_age_all_on_ch4_hist_cent = mean_age_all_on_n2o_hist_cent
mean_age_all_on_f12_hist_cent = mean_age_all_on_n2o_hist_cent
for n = 0, nn2-1 do begin
  gd = where(n2o_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_a0.age_opt_all[gd]
    stats = moment(datsub,sdev=sdev,/nan)
    stats = moment(elat_m_adj_a_early[ico2_a0[gd]],sdev=sdev,/nan)
  endif
  gd = where(n2o_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
    stats = moment(elat_m_adj_a_early[isf6_a0[gd]],sdev=sdev,/nan)
  endif

  datsub_b = 0.
  age_err_b = 0.
  elat_sub_b = 0.
  ngd_b = 0
  gd = where(n2o_norm_all_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and  finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
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
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(n2o_norm_all_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and  elat_m_adj_b_early[ico2_b] le 65 and finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(n2o_norm_all_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] gt 75 and  finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(n2o_norm_all_b[isf6_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
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
  endif
  gd = where(n2o_norm_all_b[iage_b] ge norm_grid2[n]-dn2/2. and n2o_norm_all_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_b0.age_opt_all[gd]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  ;  Combine the balloon opt ages and SF6 ages with no coincident opt value.
  gd = where(n2o_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and n2o_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
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
    stats = moment(datsub,sdev=sdev,/nan)
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
    stats = moment(elat_sub_b[1:-1],sdev=sdev,/nan)
  endif

  ;  Combine all of the data.
  datsub = 0.
  age_err = 0.
  elat_sub = 0.
  ngd_a = 0
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
  gd = where(n2o_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and n2o_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
    stats = moment(elat_sub[1:-1],sdev=sdev,/nan)
  endif
  gd = where(n2o_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and n2o_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_b0.age_opt_all[gd]]
  endif
  gd = where(n2o_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and n2o_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_b0.age_opt_all[gd]]
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
  datsub = 0.
  age_err = 0.
  elat_sub = 0.
  ngd_a = 0
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
  gd = where(bins_sf6_a0.n2o_norm_e ge norm_grid2[n]-dn2/2. and bins_sf6_a0.n2o_norm_e le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
  endif

  ;  Late time binning.
  datsub = 0.
  age_err = 0.
  elat_sub = 0.
  ngd_a = 0
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
  endif

  ;  CH4 binning
  gd = where(ch4_all_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_a0.age_opt_all[gd]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(ch4_all_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
  datsub_b = 0.
  age_err_b = 0.
  ngd_b = 0
  gd = where(ch4_all_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and  finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(ch4_all_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_sf6_b0.age_opt_all[gd]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(ch4_all_norm_b[iage_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[iage_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[iage_b] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = dat_b0.age_opt_all[gd]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  ;  Combine the balloon opt ages and SF6 ages with no coincident opt value.
  gd = where(ch4_all_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
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
    stats = moment(datsub,sdev=sdev,/nan)
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
  datsub = 0.
  age_err = 0.
  ngd_a = 0
  gd = where(ch4_all_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
  endif
  gd = where(ch4_all_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and ch4_all_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
  endif
  gd = where(ch4_all_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and ch4_all_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_b0.age_opt_all[gd]]
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
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(f12_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
  datsub_b = 0.
  age_err_b = 0.
  ngd_b = 0
  gd = where(f12_norm_b[ico2_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[ico2_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[ico2_b] ge 30 and  finite(dat_co2_b0.age_opt_all) and times_all_b[ico2_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = dat_co2_b0.age_opt_all[gd]
    datsub_b = [datsub_b,datsub]
    ngd_b += ngd
    datusub = dat_co2_b0.age_range_all[*,gd]
    age_err = 0.5 * (datusub[1,*] - datusub[0,*])
    chk = where(age_err lt 0.1,nchk)
    if nchk gt 0 then age_err[chk] = 0.1
    age_err_b = [age_err_b,reform(age_err)]
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  gd = where(f12_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
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
    stats = moment(datsub,sdev=sdev,/nan)
  endif
  ;  Combine the balloon opt ages and SF6 ages with no coincident opt value.
  gd = where(f12_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
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
    stats = moment(datsub,sdev=sdev,/nan)
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
  datsub = 0.
  age_err = 0.
  ngd_a = 0
  gd = where(f12_norm_a0[ico2_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[ico2_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[ico2_a0] ge 30,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_co2_a0.age_opt_all[gd]]
    ngd_a += ngd
    datusub = dat_co2_a0.age_range_all[*,gd]
    age_e = reform(0.5 * (datusub[1,*] - datusub[0,*]))
    chk = where(age_e lt 0.1,nchk)
    if nchk gt 0 then age_e[chk] = 0.1
    age_err = [age_err,age_e]
  endif
  gd = where(f12_norm_a0[isf6_a0] ge norm_grid2[n]-dn2/2. and f12_norm_a0[isf6_a0] le norm_grid2[n]+dn2/2. and elat_m_adj_a_early[isf6_a0] ge 30 and  times_all_a0[isf6_a0] lt 1999,ngd)
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
  endif
  gd = where(f12_norm_b[isf6_b] ge norm_grid2[n]-dn2/2. and f12_norm_b[isf6_b] le norm_grid2[n]+dn2/2. and elat_m_adj_b_early[isf6_b] ge 30 and  times_all_b[isf6_b] lt 2001,ngd)
  if ngd gt 0 then begin
    datsub = [datsub,dat_sf6_b0.age_opt_all[gd]]
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
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'mean_age_vs_tracers_1990s.nc', $
  ['norm_grid2','mean_age_all_on_n2o_a0','mean_age_all_on_n2o_b0','age_grid_h','mean_age_all_on_n2o_hist', $
  'mean_age_all_on_ch4_a0','mean_age_all_on_ch4_b0','mean_age_all_on_ch4_hist','mean_age_all_on_f12_a0', $
  'mean_age_all_on_f12_b0','mean_age_all_on_f12_hist','mean_age_all_on_n2o_hist_cent','mean_age_all_on_ch4_hist_cent', $
  'mean_age_all_on_f12_hist_cent','mean_age_all_on_n2o_0','mean_age_all_on_n2o_e_a0','mean_age_all_on_n2o_a2']
ddd = replicate(!values.f_nan,6,nz3,nt)
age_nh_coarse = ddd
n2o_nh_coarse = ddd
age_nh_coarse_lat_adj = ddd
ch4_nh_coarse = ddd
age_nhe_coarse = ddd
lat_avg_nh_coarse = replicate(!values.f_nan,nz3,nt)
rse_orig = replicate(!values.f_nan,nz3)
rse_lat_adj = rse_orig
age_nh_coarse_mean = replicate(!values.f_nan,nz3)
age_nh_coarse_lat_adj_mean = age_nh_coarse_mean
n2o_nh_coarse_mean = age_nh_coarse_mean
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
    zii = where(finite(n2o_nh_adj_tseries_lat_adj[0,zi,t]),nzi)
    if nzi ge 1 then n2o_nh_coarse[0,z,t] = mean(n2o_nh_adj_tseries_lat_adj[0,zi[zii],t],/nan)
    zii = where(finite(ch4_nh_adj_tseries[0,zi,t]),nzi)
    if nzi ge 1 then ch4_nh_coarse[0,z,t] = mean(ch4_nh_adj_tseries[0,zi[zii],t],/nan)
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
    gd = where(finite(age_nh_coarse_lat_adj[2,z,*]) and years ge 2011 and (age_nh_coarse[5,z, $
      *] eq 3 or (age_nh_coarse[5,z,*] eq 1 and years ge 2015)),ngd)
    age_nh_lat_adj_trends_ac_coarse[z,*] = linfit(years[gd],age_nh_coarse_lat_adj[2,z,gd],sigma=sigma)
    age_nh_lat_adj_trends_ac_coarse_sigma[z,*] = sigma
  endif

  ;  Residual standard error of data around linear fit.
  gd = where(finite(age_nh_coarse[2,z,*]),ngd)
  rse_orig[z] = sqrt(total((age_nh_coarse[2,z,gd] - (age_nh_trends_coarse[z,1]*years[gd] + age_nh_trends_coarse[z,0]))^2)/(ngd-2))
  gd = where(finite(age_nh_coarse_lat_adj[2,z,*]),ngd)
  rse_lat_adj[z] = sqrt(total((age_nh_coarse_lat_adj[2,z, $
    gd] - (age_nh_lat_adj_trends_coarse[z,1]*years[gd] + age_nh_lat_adj_trends_coarse[z,0]))^2)/(ngd-2))

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
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'trend_profiles.nc', $
  ['age_nh_coarse','age_nh_coarse_lat_adj','age_nh_adj_tseries','age_nh_adj_tseries_lat_adj','age_nh_lat_adj_trends_coarse', $
  'age_nh_lat_adj_trends_coarse_sigma','n2o_nh_lat_adj_trends_coarse','n2o_nh_lat_adj_trends_coarse_sigma', $
  'age_nh_trends_coarse','age_nh_trends_coarse_sigma','ch4_nh_trends_coarse','ch4_nh_trends_coarse_sigma','age_nhe_adj_tseries', $
  'age_nhe_coarse','age_nh_trends_per_coarse','age_nh_lat_adj_trends_per_coarse','age_nh_trends_per_coarse_sigma', $
  'age_nh_lat_adj_trends_per_coarse_sigma','age_nh_lat_adj_trends_ac_coarse','age_nh_lat_adj_trends_ac_coarse_sigma', $
  'n2o_nh_lat_adj_trends_per_coarse','n2o_nh_lat_adj_trends_per_coarse_sigma']

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
m_ex__ = 7
z_ex__ = 44
z2_ex__ = 28
z3_ex__ = 14

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
v__['age_lat_adj'] = ma_var(reform(age_nh_coarse_lat_adj[2,z3_ex__,*]), ['time'], 'years', $
  '30-60N average mean age, latitude adjusted')
v__['age_lat_adj_uncert'] = ma_var(reform(age_nh_coarse_lat_adj[3,z3_ex__,*]), ['time'], 'years', $
  'uncertainty of the weighted mean')
v__['sample_type'] = ma_var(fix(reform(age_nh_coarse[5,z3_ex__,*])), ['time'], '', 'platform', flag_values='1 2 3', $
  flag_meanings='aircraft in_situ_balloon aircore')
v__['fit'] = ma_var(reform(age_nh_trends_coarse[z3_ex__,*]), ['coef'], '', $
  'linear fit to age (intercept years, slope years/year)')
v__['fit_lat_adj'] = ma_var(reform(age_nh_lat_adj_trends_coarse[z3_ex__,*]), ['coef'], '', 'linear fit to age_lat_adj')
v__['rmse'] = ma_var(rse_orig[z3_ex__], '', 'years', 'RMSE of fit')
v__['rmse_lat_adj'] = ma_var(rse_lat_adj[z3_ex__], '', 'years', 'RMSE of fit_lat_adj')
v__['avg_sampled_lat'] = ma_var(reform(lat_avg_nh_coarse[z3_ex__,*]), ['time'], 'degrees_north', $
  'average sampled latitude within 30-60N')
v__['clim_lat'] = ma_var(lat_grid, ['lat'], 'degrees_north')
v__['clim_alt'] = ma_var(alt_grid2[z2_ex__], '', 'km', 'altitude of panel d')
v__['clim_age_rel_45n'] = ma_var(reform(age_grid_seas_avg[0,*,z2_ex__,*]) - rebin(reform(age_grid_seas_avg[0,26,z2_ex__,*],1,4), $
  n_elements(lat_grid), 4),  ['lat','season'], 'years', 'climatological seasonal mean age relative to 45N', $
  seasons='DJF MAM JJA SON')
ma_nc_write, fig_dir__+'supp_fig5.nc', v__, $
  global=hash('title', 'Supplementary Figure 5: NH midlatitude mean ages, latitude sampling and adjustments')
end

