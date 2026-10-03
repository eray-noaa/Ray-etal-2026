;+
;  Step 5: observed time series and trends.
;
;  Combines the aircraft and balloon time series, computes the mean age trends in
;  altitude and normalized N2O/CH4 bins with and without flask data, the 24-35 km
;  mean age series, and the 1990s and 2020s mean age vs N2O relationships.
;
;  Reads:  step 3 and 4 outputs; aircraft normalized-tracer bins; Andrews et al. N2O-age;
;          MERRA-2 tropopause; NOAA CO2 (Mauna Loa, Samoa); ASHOE flight data
;  Writes: <output>/intermediate/observed_trends.nc;
;          <output>/figure_data/ed_fig4.nc, ed_fig5.nc, supp_fig8.nc, supp_fig9.nc, supp_fig10.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step5_observed_trends

dir = mean_age_data_dir()

;  Results from step4_balloon_profiles.
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'balloon_profiles.nc', $
  ['age_combo_alt_grid_fl_yrs','age_combo_alt_grid_nofl_yrs','alt_grid2','alt_grid3','dn','dz2','dz3','max_alt_b_yrs','mm', $
  'n2o_norm_alt_grid_yrs','nn','nnorm','norm_grid','nsamp_flask','nt','nz2','nz3','samp_type_b_yrs','sdev','years','z']
restore,dir+'Balloon/Mean_age_relationships.sav'

age_grid_hist = age_grid
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_early_bins.sav'
bins_a0 = bins
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_early_bins.sav'
bins_co2_a0 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_mid_bins.sav'
bins_co2_a1 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_mid_bins.sav'
bins_sf6_a1 = bins_sf6
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_late_bins.sav'
bins_a2 = bins
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late_bins.sav'
bins_sf6_a2 = bins_sf6
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_late_bins.sav'
bins_co2_a2 = bins_co2

;  Rebin the aircraft data into coarser balloon bins.
nz = 60
dz = 0.5
alt_grid_a = findgen(nz)*dz+6.
nt_n2o = n_elements(n2o_years)
uu = replicate(!values.f_nan,nz2,5,nt_n2o)
bins_a0_co2_age_alt_nh_years = uu
bins_a2_co2_age_alt_nh_years = uu
bins_a_combo_age_alt_nh_years = uu
bins_a0_co2_n2o_alt_nh_years = uu
bins_a2_co2_n2o_alt_nh_years = uu
bins_a_combo_n2o_alt_nh_years = uu
bins_a2_sf6_age_alt_nh_years = uu
bins_a2_sf6_n2o_alt_nh_years = uu
ww = replicate(!values.f_nan,nt_n2o)
max_alts_a = ww
for t = 0, nt_n2o-1 do begin

  for i = 0, nn-1 do begin

    ii = where(n2o_norm_grid ge norm_grid[i]-dn/2. and n2o_norm_grid le norm_grid[i]+dn/2.)
    stats = moment(bins_a0.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a0.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a1.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_sf6_a1.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_a2.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_sf6_a2.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a2.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)

    ii = where(n2o_norm_grid ge norm_grid[i]-dn/2. and n2o_norm_grid le norm_grid[i]+dn/2.)
    stats = moment(bins_a0.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a0.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a1.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    if years[t] lt 2010 or years[t] gt 2011 then begin
      stats = moment(bins_sf6_a1.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    endif
    stats = moment(bins_a2.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_sf6_a2.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a2.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
  endfor

  for z = 0, nz2-1 do begin
    zi = where(alt_grid_a ge alt_grid2[z]-dz2/2. and alt_grid_a le alt_grid2[z]+dz2/2.)
    err = 1. / sqrt(total(1./bins_co2_a0.age_uncert_alt_nh_years[zi,t]^2))
    npts = total(bins_co2_a0.age_npts_alt_nh_years[zi,t])
    stats = moment(bins_co2_a0.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a0_co2_age_alt_nh_years[z,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a2.age_uncert_alt_nh_years[zi,t]^2))
    npts = total(bins_co2_a2.age_npts_alt_nh_years[zi,t])
    stats = moment(bins_co2_a2.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_co2_age_alt_nh_years[z,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_sf6_a2.age_uncert_alt_nh_years[zi,t]^2))
    npts = total(bins_sf6_a2.age_npts_alt_nh_years[zi,t])
    stats = moment(bins_sf6_a2.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_sf6_age_alt_nh_years[z,*,t] = [stats[0],err,stats[0],sdev,npts]
    if finite(bins_a0_co2_age_alt_nh_years[z,0,t]) and n2o_years[t] lt 2000 then bins_a_combo_age_alt_nh_years[z,*, $
      t] = bins_a0_co2_age_alt_nh_years[z,*,t]
    if finite(bins_a2_co2_age_alt_nh_years[z,0,t]) then bins_a_combo_age_alt_nh_years[z,*,t] = bins_a2_co2_age_alt_nh_years[z,*,t]
    if finite(bins_a2_sf6_age_alt_nh_years[z,0,t]) then bins_a_combo_age_alt_nh_years[z,*,t] = bins_a2_sf6_age_alt_nh_years[z,*,t]
    stats = moment(bins_co2_a0.n2o_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a0_co2_n2o_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    stats = moment(bins_co2_a2.n2o_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_co2_n2o_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    stats = moment(bins_sf6_a2.n2o_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_sf6_n2o_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    if finite(bins_a0_co2_n2o_alt_nh_years[z,0,t]) and n2o_years[t] lt 2000 then bins_a_combo_n2o_alt_nh_years[z,*, $
      t] = bins_a0_co2_n2o_alt_nh_years[z,*,t]
    if finite(bins_a2_co2_n2o_alt_nh_years[z,0,t]) then bins_a_combo_n2o_alt_nh_years[z,*,t] = bins_a2_co2_n2o_alt_nh_years[z,*,t]
    if finite(bins_a2_sf6_n2o_alt_nh_years[z,0,t]) then bins_a_combo_n2o_alt_nh_years[z,*,t] = bins_a2_sf6_n2o_alt_nh_years[z,*,t]
    stats = moment(bins_co2_a0.ch4_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a2.ch4_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    stats = moment(bins_sf6_a2.ch4_alt_nh_years[zi,0,t],sdev=sdev,/nan)
  endfor

  for z = 0, nz3-1 do begin
    zi = where(alt_grid_a ge alt_grid3[z]-dz3/2. and alt_grid_a le alt_grid3[z]+dz3/2.)
    stats = moment(bins_co2_a0.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    stats = moment(bins_co2_a2.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    stats = moment(bins_sf6_a2.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
  endfor

  chk = where(finite(bins_co2_a0.age_alt_nh_years[*,0,t]),nchk)
  if nchk gt 0 then max_alts_a[t] = alt_grid_a[chk[-1]]
  chk = where(finite(bins_co2_a1.age_alt_nh_years[*,0,t]),nchk)
  if nchk gt 0 then max_alts_a[t] = alt_grid_a[chk[-1]]
  chk = where(finite(bins_co2_a2.age_alt_nh_years[*,0,t]),nchk)
  if nchk gt 0 then max_alts_a[t] = alt_grid_a[chk[-1]]
  chk = where(finite(bins_sf6_a1.age_alt_nh_years[*,0,t]),nchk)
  if nchk gt 0 then max_alts_a[t] = alt_grid_a[chk[-1]]
  chk = where(finite(bins_sf6_a2.age_alt_nh_years[*,0,t]),nchk)
  if nchk gt 0 then max_alts_a[t] = alt_grid_a[chk[-1]]
endfor
bins_all_combo_age_alt_nh_years = mm
bins_all_combo_n2o_norm_alt_nh_years = mm
bins_all_combo_fl_age_alt_nh_years = mm
for t = 0, nt-1 do begin
  for z = 0, nz2-1 do begin
    for j = 0, 3 do begin
      bins_all_combo_age_alt_nh_years[t,z,j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]], $
        /nan)
      bins_all_combo_fl_age_alt_nh_years[t,z, $
        j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      if finite(age_combo_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z, $
        j] = mean([age_combo_alt_grid_fl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      if finite(age_combo_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z, $
        j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      bins_all_combo_n2o_norm_alt_nh_years[t,z,j] = mean([n2o_norm_alt_grid_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
    endfor
    j = 4
    bins_all_combo_age_alt_nh_years[t,z,j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_fl_age_alt_nh_years[t,z,j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]], $
      /nan)
    if finite(age_combo_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z, $
      j] = total([age_combo_alt_grid_fl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    if finite(age_combo_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z, $
      j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_n2o_norm_alt_nh_years[t,z,j] = total([n2o_norm_alt_grid_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
  endfor
endfor

ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'trend_profiles.nc', $
  ['age_nhe_adj_tseries','age_nhe_coarse','age_nh_adj_tseries','age_nh_adj_tseries_lat_adj','age_nh_coarse', $
  'age_nh_coarse_lat_adj','age_nh_lat_adj_trends_ac_coarse','age_nh_lat_adj_trends_ac_coarse_sigma', $
  'age_nh_lat_adj_trends_coarse','age_nh_lat_adj_trends_coarse_sigma','age_nh_lat_adj_trends_per_coarse', $
  'age_nh_lat_adj_trends_per_coarse_sigma','age_nh_trends_coarse','age_nh_trends_coarse_sigma','age_nh_trends_per_coarse', $
  'age_nh_trends_per_coarse_sigma','ch4_nh_trends_coarse','ch4_nh_trends_coarse_sigma','n2o_nh_lat_adj_trends_coarse', $
  'n2o_nh_lat_adj_trends_coarse_sigma','n2o_nh_lat_adj_trends_per_coarse','n2o_nh_lat_adj_trends_per_coarse_sigma']
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'gridded_time_series.nc', $
  ['age_co2_b_grid_adj_tseries','age_co2_grid_adj_tseries','age_grid_adj_tseries','age_on_ch4_adj_tseries', $
  'age_on_ch4_e_adj_tseries','age_on_n2o_adj_tseries','age_on_n2o_e_adj_tseries','age_opt_b_grid_adj_tseries', $
  'age_opt_grid_adj_tseries','age_sf6_b_grid_adj_tseries','age_sf6_grid_adj_tseries','ch4_grid_adj_tseries', $
  'n2o_a_grid_adj_tseries','n2o_b_grid_adj_tseries','n2o_grid_adj_tseries','n2o_on_age_adj_tseries', $
  'n2o_on_age_both_adj_tseries','n2o_on_age_b_adj_tseries']
; ------------------------  Trends  ---------------------------------------------------------------------------------------

zi = where(alt_grid2 ge 24)
bins_all_combo_fl_age_midstrat_nh_years = replicate(!values.f_nan,nt,5)
bins_all_combo_nofl_age_midstrat_nh_years = bins_all_combo_fl_age_midstrat_nh_years
for t = 0, nt-1 do begin
  gd = where(finite(bins_all_combo_fl_age_alt_nh_years[t,zi,0]),ngd)
  if ngd gt 0 then begin
    err_fl = 1. / sqrt(total(1./bins_all_combo_fl_age_alt_nh_years[t,zi[gd],1]^2))
    wmean = total(1. / bins_all_combo_fl_age_alt_nh_years[t,zi[gd],1]^2 * bins_all_combo_fl_age_alt_nh_years[t,zi[gd], $
      0]) /  total(1. / bins_all_combo_fl_age_alt_nh_years[t,zi[gd],1]^2)
    bins_all_combo_fl_age_midstrat_nh_years[t,0:1] = [wmean,err_fl]
  endif
  gd = where(finite(bins_all_combo_age_alt_nh_years[t,zi,1]),ngd)
  if ngd gt 0 then begin
    err_nofl = 1. / sqrt(total(1./bins_all_combo_age_alt_nh_years[t,zi[gd],1]^2))
    wmean = total(1. / bins_all_combo_age_alt_nh_years[t,zi[gd],1]^2 * bins_all_combo_age_alt_nh_years[t,zi[gd], $
      0]) /  total(1. / bins_all_combo_age_alt_nh_years[t,zi[gd],1]^2)
    bins_all_combo_nofl_age_midstrat_nh_years[t,0:1] = [wmean,err_nofl]
  endif
  for j = 2, 3 do begin
    bins_all_combo_fl_age_midstrat_nh_years[t,j] = mean(bins_all_combo_fl_age_alt_nh_years[t,zi,j],/nan)
    bins_all_combo_nofl_age_midstrat_nh_years[t,j] = mean(bins_all_combo_age_alt_nh_years[t,zi,j],/nan)
  endfor
  bins_all_combo_fl_age_midstrat_nh_years[t,4] = total(bins_all_combo_fl_age_alt_nh_years[t,zi,4],/nan)
  bins_all_combo_nofl_age_midstrat_nh_years[t,4] = total(bins_all_combo_age_alt_nh_years[t,zi,4],/nan)
endfor
gd = where(finite(bins_all_combo_fl_age_midstrat_nh_years[*,0]))
age_fl_midstrat_linfit = linfit(years[gd],bins_all_combo_fl_age_midstrat_nh_years[gd,0],sigma=sigma)
gd = where(finite(bins_all_combo_nofl_age_midstrat_nh_years[*,0]))
age_nofl_midstrat_linfit = linfit(years[gd],bins_all_combo_nofl_age_midstrat_nh_years[gd,0],sigma=sigma)

vv = replicate(!values.f_nan,nn,2)
ww = replicate(!values.f_nan,nz2,2)
age_alt_trends = ww
n2o_alt_trends = ww
age_fl_alt_trends = ww
age_fl_alt_trends_sigma = ww
n2o_nofl_age_trends = vv
n2o_nofl_age_trends_sigma = vv
age_on_ch4_trends = vv
age_on_ch4_trends_sigma = vv
for i = 0, nn-1 do begin
  gd = where(finite(age_on_n2o_adj_tseries[2,i,*]),ngd)
  if ngd ge 20 then begin
    n2o_nofl_age_trends[i,*] = linfit(years[gd],age_on_n2o_adj_tseries[2,i,gd],sigma=sigma)
    n2o_nofl_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_ch4_adj_tseries[2,i,*]),ngd)
  if ngd ge 20 then begin
    age_on_ch4_trends[i,*] = linfit(years[gd],age_on_ch4_adj_tseries[2,i,gd],sigma=sigma)
    age_on_ch4_trends_sigma[i,*] = sigma
  endif
endfor

age_nh_coarse_wflask = age_nh_coarse

for z = 0, nz2-1 do begin
  gd = where(finite(bins_all_combo_age_alt_nh_years[*,z,0]),ngd)
  if ngd ge 15 then begin
    chk = where(~finite(bins_all_combo_age_alt_nh_years[gd,z,1]),nchk)
    if nchk gt 0 then bins_all_combo_age_alt_nh_years[gd[chk],z,1] = 0.25
    age_alt_trends[z,*] = linfit(years[gd],bins_all_combo_age_alt_nh_years[gd,z,0],sigma=sigma)
  endif

  ;  Trends including flask data.  Add in the elat adjusted non-flask data to the time series.
  si = where(samp_type_b_yrs eq 0,nsi)
  for j = 0, 1 do age_nh_coarse_wflask[j,z,si] = bins_all_combo_fl_age_alt_nh_years[si,z,j]
  for j = 2, 3 do age_nh_coarse_wflask[j,z,si] = bins_all_combo_fl_age_alt_nh_years[si,z,j-2]
  age_nh_coarse_wflask[4,z,si] = bins_all_combo_fl_age_alt_nh_years[si,z,4]
  chk = where(samp_type_b_yrs eq 0 and age_nh_coarse[5,z,*] gt 0,nchk)
  if nchk eq 0 then age_nh_coarse_wflask[5,z,si] = 0
  gd = where(finite(age_nh_coarse_wflask[2,z,*]),ngd)
  if ngd ge 15 then begin
    age_fl_alt_trends[z,*] = linfit(years[gd],age_nh_coarse_wflask[2,z,gd],sigma=sigma)
    age_fl_alt_trends_sigma[z,*] = sigma
  endif
  gd = where(finite(bins_all_combo_n2o_norm_alt_nh_years[*,z,0]) and years ge 1990,ngd)
  if ngd ge 15 then begin
    n2o_alt_trends[z,*] = linfit(years[gd],bins_all_combo_n2o_norm_alt_nh_years[gd,z,0],sigma=sigma)
  endif
endfor
restore,dir+'Aircraft/Airborne_Save_Files/Aircraft_mean_ages_early_bins_comb.sav'
nyg = 36
dy = 5.
nz = 60
dz = 0.5
lat_grid_a = findgen(nyg)*dy-85.
alt_grid_a0 = findgen(nz)*dz+6.

;  Read in mean age vs N2O from elat adj program.
ma_restore_nc, ma_output_dir() + 'intermediate' + path_sep() + 'mean_age_vs_tracers_1990s.nc', $
  ['age_grid_h','mean_age_all_on_ch4_a0','mean_age_all_on_ch4_b0','mean_age_all_on_ch4_hist','mean_age_all_on_ch4_hist_cent', $
  'mean_age_all_on_f12_a0','mean_age_all_on_f12_b0','mean_age_all_on_f12_hist','mean_age_all_on_f12_hist_cent', $
  'mean_age_all_on_n2o_0','mean_age_all_on_n2o_a0','mean_age_all_on_n2o_a2','mean_age_all_on_n2o_b0','mean_age_all_on_n2o_e_a0', $
  'mean_age_all_on_n2o_hist','mean_age_all_on_n2o_hist_cent','norm_grid2']
;  Calculate 2020s mean age vs N2O from 1990s avg and trend.
nic = interpol(indgen(nnorm),n2o_norm_grid,norm_grid)
bins_a0_n2o_norm_age_nh_coarse = interpolate(bins_a0_n2o_norm_age_nh,nic)
bins_a2_n2o_norm_age_nh_trend = bins_a0_n2o_norm_age_nh_coarse + 25.*n2o_nofl_age_trends[*,1]

mean_age_on_n2o_90s = mean_age_all_on_n2o_0[*,2]
mean_age_on_n2o_90s[0:-4] = smooth(mean_age_on_n2o_90s[0:-4],5,/edge_truncate)
mean_age_on_n2o_90s[0:-2] = smooth(mean_age_on_n2o_90s[0:-2],3,/edge_truncate)
ni = interpol(indgen(n_elements(norm_grid)),norm_grid,norm_grid2)
n2o_age_trend_interp = interpolate(n2o_nofl_age_trends[*,1],ni)
n2o_age_trend_interp[-3:-1] = [-0.012,-0.01,-0.005]
mean_age_on_n2o_20s_from_trend = mean_age_on_n2o_90s + 25.*n2o_age_trend_interp
mean_age_on_n2o_20s_from_trend[19:-3] = smooth(mean_age_on_n2o_20s_from_trend[19:-3],3,/edge_truncate)
bbb = replicate(!values.f_nan,2,nz2)
n2o_norm_alt_profile_avg = bbb
age_alt_profile_avg = bbb
ti = where(years lt 2000)
for z = 0, nz2-1 do begin
  gd = where(finite(bins_all_combo_n2o_norm_alt_nh_years[ti,z,0]),ngd)
  if ngd gt 1 then for j = 0, 1 do n2o_norm_alt_profile_avg[j,z] = mean(bins_all_combo_n2o_norm_alt_nh_years[ti[gd],z,j])
  gd = where(finite(bins_all_combo_age_alt_nh_years[ti,z,0]),ngd)
  if ngd gt 1 then for j = 0, 1 do age_alt_profile_avg[j,z] = mean(bins_all_combo_age_alt_nh_years[ti[gd],z,j])
endfor
gd = where(finite(n2o_norm_alt_profile_avg))
n2o_norm_alt_profile_avg[gd] = smooth(n2o_norm_alt_profile_avg[gd],5,/edge_truncate)
gd = where(finite(age_alt_profile_avg[0,*]))
age_alt_profile_avg[0,gd] = smooth(age_alt_profile_avg[0,gd],5,/edge_truncate)

age_alt_profile_from_trend = 25.*age_alt_trends[*,1] + reform(age_alt_profile_avg[0,*])
n2o_alt_trends_combo = n2o_nh_lat_adj_trends_coarse[*,1]
n2o_alt_trends_combo[22] = 0.
n2o_norm_alt_profile_from_trend = 25.*n2o_alt_trends_combo + reform(n2o_norm_alt_profile_avg[0,*])

;  Tropopause data.
ncdf_get,dir+'Reanalysis/MERRA2/tp.monmean.zm.nc',['time','tpp','tpt','tpz','lat'],tp,/quiet
ny_tp = tp['lat','dim_sizes']
lat_tp = tp['lat','value']
tpause = fltarr(ny_tp)
tpause_alt = tpause
for y = 0, ny_tp[0]-1 do begin
  tpause_alt[y] = mean(tp['tpz','value',y,*])
endfor

yi = interpol(findgen(ny_tp),lat_tp,lat_grid_a)
tp_alt_grid = interpolate(tpause_alt,yi)

restore,dir+'Aircraft/Airborne_Save_Files/Andrews_N2O_age.sav'

mission_names = ['aaoe','aase','aase2','ashoe','spade','strat','polaris','solve','accent','crystalf','pre_ave','ave_0506', $
  'cr_ave','tc4','glopac','attrex',  'seac4rs','ATom1','ATom2','ATom3','ATom4','dcotts1','dcotts2','sabre']

colors = ['sky blue','Blue','dark blue','Purple','Medium Purple','Magenta','Lime Green','Green','gold','Orange','orange red', $
  'Red','Brown','Grey','Dark Grey','dodger blue','cornflower',  'Purple','Medium Purple','Magenta']

restore,dir+'Balloon/NSamp_OMS.sav'
tot_samp = fltarr(nz2)
tot_samp_aircraft = tot_samp
tot_samp_aircore = tot_samp
for z = 0, nz2-1 do begin
  si = where(age_nh_coarse[5,z,*] eq 1)
  tot_samp_aircraft[z] = total(age_nh_coarse[4,z,si],/nan)
  si = where(age_nh_coarse[5,z,*] eq 3)
  tot_samp_aircore[z] = total(age_nh_coarse[4,z,si],/nan)
  tot_samp[z] = tot_samp_aircraft[z] + tot_samp_aircore[z] + nsamp_flask[z] + nsamp_oms[z]
endfor

;  Surface CO2 and the ASHOE flight for Supplementary Fig. 9
co2_m = read_ascii(dir+'Surface_Trace_Gas/CO2/CO2_mlo_monthly.txt',data_start=70)
yrs_m = reform(co2_m.field1[1,*])
mons_m = reform(co2_m.field1[2,*])
years_mlo = float(yrs_m) + float((mons_m-1.))/12. + 1./24.
co2_mlo = reform(co2_m.field1[3,*])

co2_s = read_ascii(dir+'Surface_Trace_Gas/CO2/CO2_smo_monthly.txt',data_start=70)
yrs_s = reform(co2_s.field1[1,*])
mons_s = reform(co2_s.field1[2,*])
co2_smo = reform(co2_s.field1[3,*])
years_smo = float(yrs_s) + float((mons_s-1.))/12. + 1./24.

ti = where(years_smo ge years_mlo[0],nt_sm)
years_sm = years_mlo
co2_sm = 0.5 * (co2_mlo + co2_smo[ti])
co2_sm_smooth = lowpass_cfc(co2_sm, BOX=12, EDGE_PFCAST=1)

co2_steady_growth_time = [1980,2000]
co2_steady_growth = [337.2,366.5]

restore,dir+'Aircraft/Airborne_Save_Files/ashoe.sav'

fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
; ---------- Extended Data Figure 4
v__ = orderedhash()
v__['years'] = ma_var(years, ['time'], 'year', 'decimal year (bins of 1/48 year)')
v__['balloon_max_alt'] = ma_var(max_alt_b_yrs, ['time'], 'km', 'maximum altitude of balloon profiles')
v__['balloon_sample_type'] = ma_var(samp_type_b_yrs, ['time'], '', 'balloon sampling type', flag_values='0 1 2 3', $
  flag_meanings='flask in_situ_balloon aircore flask')
v__['aircraft_years'] = ma_var(n2o_years, ['time_aircraft'], 'year', 'decimal year')
v__['aircraft_max_alt'] = ma_var(max_alts_a, ['time_aircraft'], 'km', 'maximum altitude of aircraft flights')
v__['age_24_35km_with_flask'] = ma_var(bins_all_combo_fl_age_midstrat_nh_years[*,0:1], ['time','mean_unc'], 'years', $
  '24-35 km weighted mean age including flask data (column 0 mean, 1 uncertainty)')
v__['age_24_35km_no_flask'] = ma_var(bins_all_combo_nofl_age_midstrat_nh_years[*,0:1], ['time','mean_unc'], 'years', $
  '24-35 km weighted mean age without flask data')
v__['fit_with_flask'] = ma_var(age_fl_midstrat_linfit, ['coef'], '', 'linear fit (intercept, slope per year) through all data')
v__['fit_no_flask'] = ma_var(age_nofl_midstrat_linfit, ['coef'], '', 'linear fit through in situ and AirCore data')
ma_nc_write, fig_dir__+'ed_fig4.nc', v__, global=hash('title', 'Extended Data Figure 4: sampling altitudes and 24-35 km mean age')

; ---------- Extended Data Figure 5
v__ = orderedhash()
v__['alt_grid2'] = ma_var(alt_grid2, ['alt'], 'km')
v__['age_trend_no_flask_lat_adj'] = ma_var(age_nh_lat_adj_trends_coarse[*,1], ['alt'], 'years/year', $
  'in situ without flasks, latitude adjusted')
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
ma_nc_write, fig_dir__+'ed_fig5.nc', v__, global=hash('title', $
  'Extended Data Figure 5: mean age trend and measurement number profiles')

; ---------- Supplementary Figure 8: mean age vs normalized N2O, this study and earlier relationships
n2o_ppb_engel__ = 313. * n2o_norm_grid
v__ = orderedhash()
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['age_1990s'] = ma_var(bins_a0_n2o_norm_age_nh_coarse, ['norm'], 'years', $
  'NH mean age vs normalized N2O, 1990s in situ (this study)')
v__['age_2020s'] = ma_var(bins_a2_n2o_norm_age_nh_trend, ['norm'], 'years', $
  'NH mean age vs normalized N2O, 2020s in situ (this study)')
v__['n2o_norm_grid_fine'] = ma_var(n2o_norm_grid, ['norm_fine'], '', 'normalized N2O grid of the published relationships')
v__['age_andrews2001'] = ma_var(n2o_age_a, ['norm_fine'], 'years', $
  'Andrews et al. (2001) ER-2 relationship, converted to normalized N2O')
v__['age_engel2002'] = ma_var(6.03 - 0.0136*n2o_ppb_engel__ + 8.5892e-5*n2o_ppb_engel__^2 - 3.38e-7*n2o_ppb_engel__^3, $
  ['norm_fine'], 'years',  'Engel et al. (2002) polynomial mean age vs N2O, evaluated at N2O = 313 ppb x normalized N2O')
ma_nc_write, fig_dir__+'supp_fig8.nc', v__, global=hash('title', 'Supplementary Figure 8: mean age vs normalized N2O')

; ---------- Supplementary Figure 9: surface CO2 and the 6 Aug 1994 ASHOE flight
i_ft__ = 29
v__ = orderedhash()
v__['years_surface'] = ma_var(years_sm+2./12., ['time'], 'year', 'decimal year, shifted by the 2-month delay used in the figure')
v__['co2_mlo_smo'] = ma_var(co2_sm, ['time'], 'ppmv', 'Mauna Loa and Samoa average CO2')
v__['co2_mlo_smo_smooth'] = ma_var(co2_sm_smooth, ['time'], 'ppmv', 'smoothed Mauna Loa and Samoa average CO2')
v__['steady_growth_time'] = ma_var(float(co2_steady_growth_time), ['endpoint'], 'year', $
  'constant growth rate CO2 used in a previous study')
v__['steady_growth_co2'] = ma_var(co2_steady_growth, ['endpoint'], 'ppmv')
v__['flight_date'] = ma_var(flight_track[i_ft__].date, '', '', 'ER-2 flight date (YYYYMMDD)')
v__['flight_date_float'] = ma_var(flight_track[i_ft__].date_float, '', 'year')
v__['flight_co2_adj'] = ma_var(*flight_track[i_ft__].co2_adj, ['obs'], 'ppmv', 'flight CO2 adjusted for CH4 oxidation')
v__['flight_co2_archive_adj'] = ma_var(*flight_track[i_ft__].co2 - (1.725-*flight_track[i_ft__].ch4/1e3), ['obs'], 'ppmv', $
  'flight CO2 with the archived CH4 adjustment (CO2 - (1.725 - CH4/1000))')
v__['flight_co2_age_archive'] = ma_var(*flight_track[i_ft__].co2_age, ['obs'], 'years', 'archived mean age')
v__['flight_co2_lag_age'] = ma_var(*flight_track[i_ft__].co2_lag_age, ['obs'], 'years', $
  'simple lag-technique mean age (this study)')
ma_nc_write, fig_dir__+'supp_fig9.nc', v__, global=hash('title', 'Supplementary Figure 9: CO2 time series and mean age example')

; ---------- Supplementary Figure 10: archived minus updated CO2 mean age vs N2O, by 1990s mission
mi__ = indgen(5) + 3
v__ = orderedhash()
v__['n2o_norm_grid'] = ma_var(n2o_norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['mission'] = ma_var(mission_names[mi__], ['mission'], '', 'ER-2 mission')
v__['mission_color'] = ma_var(colors[mi__], ['mission'], '', 'plot colour')
v__['age_archived'] = ma_var(reform(bins_co2_a0.n2o_norm_co2_age_m[*,0,mi__]), ['norm','mission'], 'years', $
  'archived CO2 mean age averaged in normalized N2O bins')
v__['age_updated'] = ma_var(reform(bins_co2_a0.n2o_norm_age_all_m[*,0,mi__]), ['norm','mission'], 'years', $
  'updated CO2 mean age (this study)')
v__['age_updated_uncert'] = ma_var(reform(bins_co2_a0.n2o_norm_age_all_m[*,1,mi__]), ['norm','mission'], 'years', $
  'uncertainty of the bin average')
ma_nc_write, fig_dir__+'supp_fig10.nc', v__, $
  global=hash('title', 'Supplementary Figure 10: archived mean age bias vs N2O',  'note', $
  'The blue average line is the 7-point smoothed mean over missions of age_archived - age_updated (computed in supp_fig10.pro).')

;  Hand the results needed by later steps to them.
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'observed_trends.nc', $
  ['age_alt_profile_avg','age_alt_profile_from_trend','age_grid_hist','age_on_ch4_trends','age_on_ch4_trends_sigma', $
  'alt_grid_a0','lat_grid_a','lat_tp','mean_age_on_n2o_20s_from_trend','mean_age_on_n2o_90s','n2o_alt_trends', $
  'n2o_nofl_age_trends','n2o_nofl_age_trends_sigma','n2o_norm_alt_profile_avg','n2o_norm_alt_profile_from_trend','nz', $
  'tp_alt_grid','tpause_alt','y','z']
end

