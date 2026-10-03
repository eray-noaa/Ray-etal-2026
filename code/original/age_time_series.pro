pro age_time_series, BUILD_PROFILES=build_profiles

;  BUILD_PROFILES: rebuild Balloon/Engel_profiles_v2025.sav (flask, OMS, AirCore and WAS
;  balloon profiles with mean ages) from the merged inputs. Without it the cached file is used.

dir = mean_age_data_dir()

data = read_ascii(dir+'Surface_Trace_Gas/GMD/n2o/combined/GML_global_N2O.txt',data_start=68)
n2o_yr = reform(data.field01[0,*])
n2o_mon = reform(data.field01[1,*])
n2o_date = n2o_yr + (n2o_mon-0.5)/12.
n2o_global = reform(data.field01[6,*])

;  Read in GMD marine boundary layer data for CH4.
dat = read_ascii(dir+'Surface_Trace_Gas/CH4/CH4_mbl_30S-30N.txt',data_start=72)
ch4_time = reform(dat.field1[3,*])
ch4_mbl = reform(dat.field1[4,*])

data = read_ascii(dir+'Surface_Trace_Gas/GMD/cfc12/combined/HATS_global_F12.txt',data_start=68)
f12_yr = reform(data.field01[0,*])
f12_mon = reform(data.field01[1,*])
f12_date = f12_yr + (f12_mon-0.5)/12.
f12_global = reform(data.field01[6,*])

restore,dir+'Balloon/Engel_mean_ages_opt.sav'
bal4 = bal
dates4 = dates
nf4 = n_elements(dates4)

;  Read in the full OMS profiles.
restore,dir+'Balloon/OMS_common_gc.sav'

restore,dir+'Aircraft/Missions/idlsave_files/N2O_CH4_norm_curve_post2000.sav'


read1 = keyword_set(build_profiles)

if read1 then begin

  bal = hash()

  restore,dir+'Balloon/Engel_mean_ages_opt.sav'

  ;  Read the Ray et al updated Engel file.
  dat = read_ascii(dir+'Balloon/Engel_balloon_all5.txt',data_start=3)
  dates = reform(dat.field01[0,*])
  yrs = reform(dat.field01[1,*])
  alt = reform(dat.field01[2,*])
  p = reform(dat.field01[3,*])
  theta = reform(dat.field01[4,*])
  lat = reform(dat.field01[5,*])
  lon = reform(dat.field01[6,*])
  elat = reform(dat.field01[7,*])
  n2o = reform(dat.field01[8,*])
  f11 = reform(dat.field01[9,*])
  f12 = reform(dat.field01[10,*])
  f113 = reform(dat.field01[11,*])
  ch4 = reform(dat.field01[12,*])
  co2 = reform(dat.field01[13,*])
  sf6 = reform(dat.field01[15,*])
  dum = ''
  openr,1,dir+'Balloon/Engel_balloon_all5.txt'
  for j = 0, 2 do readf,1,dum
  dat_str = strarr(n_elements(dates))
  for j = 0, n_elements(dates)-1 do begin
    readf,1,dum
    dat_str[j] = strmid(dum,1,8)
  endfor
  close,1
  di = uniq(dat_str)
  dates = dat_str[di]
  dates5 = dates
  yr_frac = yrs[di]

  i_jp = [5,indgen(5)+7,13,14,16,19,21,23,24,26,30,33,34,35,39,40]

  nd = n_elements(dates)
  for i = 0, nd-1 do begin
    ii = where(dat_str eq dates[i],ni)
    bal['lon_'+dates[i]] = lon[ii]
    bal['lat_'+dates[i]] = lat[ii]
    bal['alt_'+dates[i]] = alt[ii]
    bal['pres_'+dates[i]] = p[ii]
    bal['sf6_'+dates[i]] = sf6[ii]
    bal['n2o_'+dates[i]] = n2o[ii]
    bal['co2_'+dates[i]] = co2[ii]
    bal['ch4_'+dates[i]] = ch4[ii]
    bal['f12_'+dates[i]] = f12[ii]
    bal['f11_'+dates[i]] = f11[ii]

    ;  **********   Adjust CO2 from Japanese flask data as in Sugawara et al., 2025   *****************
    chk = where(i_jp eq i,nchk)
    if nchk gt 0 then co2[ii] -= 0.2

    ti = interpol(indgen(n_elements(n2o_date)),n2o_date,yr_frac[i])
    n2o_surface = interpolate(n2o_global,ti)
    n2o_norm = n2o[ii]/n2o_surface
    bal['n2o_norm_'+dates[i]] = n2o_norm
    ti = interpol(indgen(n_elements(f12_date)),f12_date,yr_frac[i])
    f12_surface = interpolate(f12_global,ti)
    f12_norm = f12[ii]/f12_surface
    bal['f12_norm_'+dates[i]] = f12_norm
    ti = interpol(indgen(n_elements(ch4_time)),ch4_time,yr_frac[i])
    ch4_surface = interpolate(ch4_mbl,ti)

    co2_err = replicate(0.05,ni)

    calculate_opt_co2_mean_age,yr_frac[i],lat[ii],co2[ii],co2_err,ch4[ii],n2o_norm,co2_adj,co2_age,co2_age_range,ratio,lat_source,ch4_entry
    bal['co2_adj_'+dates[i]] = co2_adj
    bal['co2_age_'+dates[i]] = co2_age
    diff1 = co2_age - reform(co2_age_range[0,*])
    diff2 = reform(co2_age_range[1,*]) - co2_age
    chk_d = where(diff1 gt 2,nchkd)
    if nchkd gt 0 then co2_age_range[0,chk_d] = co2_age[chk_d] - diff2[chk_d]
    bal['co2_age_range_'+dates[i]] = co2_age_range
    bal['co2_age_ratio_'+dates[i]] = ratio
    bal['ch4_norm_'+dates[i]] = ch4[ii] / ch4_entry

    calculate_opt_sf6_age,yr_frac[i],lat[ii],sf6[ii],n2o_norm,sf6_age,sf6_age_range,ratio,lat_source,sf6_adj,ch4_entry
    bal['sf6_adj_'+dates[i]] = sf6_adj
    bal['sf6_age_'+dates[i]] = sf6_age
    bal['sf6_age_range_'+dates[i]] = sf6_age_range
    bal['sf6_age_ratio_'+dates[i]] = ratio
    gd = where(finite(co2_age),ngd)
    if ngd eq 0 then bal['ch4_norm_'+dates[i]] = ch4[ii] / ch4_entry

    chk = where(finite(n2o_norm),nchk)
    if nchk eq 0 then begin
      chk = where(finite(n2o_norm_ch4_norm_avg[*,0]),nchk)
      ii = interpol(indgen(nchk),n2o_norm_ch4_norm_avg[chk,0],ch4[ii]/ch4_surface)
      n2o_norm = interpolate(n2o_norm_grid[chk],ii)
    endif

    calculate_opt_age,yr_frac[i],lat[ii],sf6[ii],co2[ii],ch4[ii],n2o_norm,co2_adj,sf6_adj,age,age_range,ratio,lat_source,ch4_entry,sf6_scale,co2_scale,sf6_bias_opt
    bal['sf6_adj_opt_'+dates[i]] = sf6_adj
    bal['co2_adj_opt_'+dates[i]] = co2_adj
    bal['opt_age_'+dates[i]] = age
    bal['opt_age_range_'+dates[i]] = age_range
    bal['opt_ratio_'+dates[i]] = ratio
    bal['sf6_scale_'+dates[i]] = sf6_scale
    bal['co2_scale_'+dates[i]] = co2_scale
    gd = where(finite(co2_age),ngd)
    if ngd eq 0 then bal['ch4_norm_opt_'+dates[i]] = ch4[ii] / ch4_entry
  endfor

  save,dates,yr_frac,bal,filename=dir+'Balloon/Engel_profiles_v2025.sav'

endif else begin
  bal = 0  ; avoid IDL 9.1 crash when restoring a hash over an existing hash
  restore,dir+'Balloon/Engel_profiles_v2025.sav'
endelse


;  Add in full OMS data.
read2 = keyword_set(build_profiles)

if read2 then begin

  restore,dir+'Balloon/Balloon_mean_ages_sweep_sf6_co2_early.sav'
  restore,dir+'Balloon/Balloon_mean_ages_co2_early.sav'
  restore,dir+'Balloon/Balloon_mean_ages_sf6_early.sav'
  restore,dir+'Balloon/OMS_common_gc_merge.sav'

  chk_c = where(finite(co2) and finite(sf6),nchk_c)
  ico2 = where(finite(co2),nco2)
  isf6 = where(finite(sf6),nsf6)

  mi = [0,2,5,8,9,10]
  dates = [dates[0:17],oms_dates[mi[1]],dates[18:-1]]
  yr_frac = [yr_frac[0:17],flight_year_frac[mi[1]],yr_frac[18:-1]]

  for i = 0, n_elements(mi)-1 do begin
    chk = where(mission eq mi[i]+1,nchk)
    bal['lon_'+oms_dates[mi[i]]] = lon[chk]
    bal['lat_'+oms_dates[mi[i]]] = lat[chk]
    bal['alt_'+oms_dates[mi[i]]] = alt[chk]
    bal['pres_'+oms_dates[mi[i]]] = pres[chk]
    bal['elat_'+oms_dates[mi[i]]] = elat_m[chk]
    bal['sf6_'+oms_dates[mi[i]]] = sf6[chk]
    bal['co2_'+oms_dates[mi[i]]] = co2_all[chk]
    bal['n2o_'+oms_dates[mi[i]]] = n2o[chk]
    bal['n2o_norm_'+oms_dates[mi[i]]] = n2o_norm[chk]
    bal['ch4_'+oms_dates[mi[i]]] = ch4[chk]
    bal['ch4_norm_'+oms_dates[mi[i]]] = ch4_all_norm[chk]
    bal['f12_norm_'+oms_dates[mi[i]]] = f12_norm[chk]

    age_both = replicate(!values.f_nan,nchk) & sf6_age = age_both & sf6_age_range = replicate(!values.f_nan,2,nchk) & co2_age = age_both & co2_age_range = sf6_age_range
    ai1 = where(dat.mission eq mi[i]+1,nai1)
    ai = where(finite(co2_all[chk]) and finite(sf6[chk]),nai)
    if nai gt 0 then age_both[ai] = dat.age_opt_all[ai1]
    bal['both_age_'+oms_dates[mi[i]]] = age_both
    ai1 = where(dat_co2.mission eq mi[i]+1,nai1)
    ai = where(finite(co2[chk]),nai)
    if nai gt 0 then begin
      co2_age[ai] = dat_co2.age_opt_all[ai1]
      co2_age_range[*,ai] = dat_co2.age_range_all[*,ai1]
    endif
    bal['co2_age_'+oms_dates[mi[i]]] = co2_age
    diff1 = co2_age - reform(co2_age_range[0,*])
    diff2 = reform(co2_age_range[1,*]) - co2_age
    chk_d = where(diff1 gt 2,nchkd)
    if nchkd gt 0 then co2_age_range[0,chk_d] = co2_age[chk_d] - diff2[chk_d]
    bal['co2_age_range_'+oms_dates[mi[i]]] = co2_age_range
    ai1 = where(dat_sf6.mission eq mi[i]+1,nai1)
    ai = where(finite(sf6[chk]),nai)
    if nai gt 0 then begin
      sf6_age[ai] = dat_sf6.age_opt_all[ai1]
      sf6_age_range[*,ai] = dat_sf6.age_range_all[*,ai1]
    endif
    bal['sf6_age_'+oms_dates[mi[i]]] = sf6_age
    bal['sf6_age_range_'+oms_dates[mi[i]]] = sf6_age_range
    gd = where(finite(sf6_age),ngd)
    gd2 = where(finite(sf6_age_range),ngd2)
  endfor

  save,dates,yr_frac,bal,filename=dir+'Balloon/Engel_profiles_v2025.sav'

endif else begin
  bal = 0  ; avoid IDL 9.1 crash when restoring a hash over an existing hash
  restore,dir+'Balloon/Engel_profiles_v2025.sav'
endelse


;  Add the six AirCore profiles in 2015-17.
read3 = keyword_set(build_profiles)

if read3 then begin

  year_frac = fltarr(nf4)

  for i = nf4-6, nf4-1 do begin
    year_frac[i] = float(strmid(dates4[i],0,4)) + (float(strmid(dates4[i],4,2))-1.)/12. + (float(strmid(dates4[i],6,2))-1.)/365.25
    bal['lon_'+dates4[i]] = bal4['lon_'+dates4[i]]
    bal['lat_'+dates4[i]] = bal4['lat_'+dates4[i]]
    bal['alt_'+dates4[i]] = bal4['alt_'+dates4[i]]
    bal['pres_'+dates4[i]] = bal4['pres_'+dates4[i]]
    bal['sf6_'+dates4[i]] = bal4['sf6_'+dates4[i]]
    bal['n2o_'+dates4[i]] = bal4['n2o_'+dates4[i]]
    bal['co2_'+dates4[i]] = bal4['co2_'+dates4[i]]
    bal['ch4_'+dates4[i]] = bal4['ch4_'+dates4[i]]
    np = n_elements(bal['alt_'+dates4[i]])
    bal['f12_norm_'+dates4[i]] = replicate(!values.f_nan,np)
    bal['sf6_age_'+dates4[i]] = replicate(!values.f_nan,np)


    chk = where(bal['alt_'+dates4[i]] lt 0,nchk)
    if nchk gt 0 then bal['alt_'+dates4[i],chk] = !values.f_nan
    chk = where(bal['ch4_'+dates4[i]] lt 0,nchk)
    if nchk gt 0 then bal['ch4_'+dates4[i],chk] = !values.f_nan

    ti = interpol(indgen(n_elements(n2o_date)),n2o_date,year_frac[i])
    n2o_surface = interpolate(n2o_global,ti)
    bal['n2o_norm_'+dates4[i]] = bal['n2o_'+dates4[i]] / n2o_surface

    co2_err = replicate(0.05,np)

    calculate_opt_co2_mean_age,year_frac[i],bal['lat_'+dates4[i]],bal['co2_'+dates4[i]],co2_err,bal['ch4_'+dates4[i]],bal['n2o_norm_'+dates4[i]],co2_adj,co2_age,co2_age_range,ratio,lat_source, $
      ch4_entry
    bal['co2_adj_'+dates4[i]] = co2_adj
    bal['co2_age_'+dates4[i]] = co2_age
    diff1 = co2_age - reform(co2_age_range[0,*])
    diff2 = reform(co2_age_range[1,*]) - co2_age
    chk_d = where(diff1 gt 2,nchkd)
    if nchkd gt 0 then co2_age_range[0,chk_d] = co2_age[chk_d] - diff2[chk_d]
    bal['co2_age_range_'+dates4[i]] = co2_age_range
    bal['co2_age_ratio_'+dates4[i]] = ratio
    bal['ch4_norm_'+dates4[i]] = bal['ch4_'+dates4[i]] / ch4_entry
  endfor

  dates = [dates[0:40],dates4[-6:-1],dates[-1]]
  yr_frac = [yr_frac[0:40],year_frac[-6:-1],yr_frac[-1]]

  save,dates,yr_frac,bal,filename=dir+'Balloon/Engel_profiles_v2025.sav'

endif else begin
  bal = 0  ; avoid IDL 9.1 crash when restoring a hash over an existing hash
  restore,dir+'Balloon/Engel_profiles_v2025.sav'
endelse


;  Add in NOAA AirCores.
read4 = keyword_set(build_profiles)

if read4 then begin

  restore,dir+'Balloon/Balloon_mean_ages_sweep_sf6_co2_late_v2025.sav'
  restore,dir+'Balloon/Balloon_mean_ages_co2_late_v2025.sav'
  restore,dir+'Balloon/Balloon_mean_ages_sf6_late.sav'
  restore,dir+'Balloon/Aircore_common_merge.sav'

  chk_c = where(finite(co2) and finite(sf6),nchk_c)
  ico2 = where(finite(co2),nco2)
  isf6 = where(finite(sf6),nsf6)

  nf = n_elements(flight_dates)

  for i = 0, nf-1 do begin
    chk = where(mission eq i+1,nchk)
    bal['lon_'+vdt_all[i]] = lon[chk]
    bal['lat_'+vdt_all[i]] = lat[chk]
    bal['alt_'+vdt_all[i]] = alt[chk]
    bal['pres_'+vdt_all[i]] = pres[chk]
    bal['elat_'+vdt_all[i]] = elat_m[chk]
    bal['sf6_'+vdt_all[i]] = sf6[chk]
    bal['co2_'+vdt_all[i]] = co2_all[chk]
    bal['n2o_'+vdt_all[i]] = n2o[chk]
    bal['n2o_norm_'+vdt_all[i]] = n2o_norm[chk]
    bal['n2o_ch4_norm_'+vdt_all[i]] = n2o_ch4_norm[chk]
    bal['ch4_'+vdt_all[i]] = ch4[chk]
    bal['ch4_norm_'+vdt_all[i]] = ch4_all_norm[chk]
    bal['f12_norm_'+vdt_all[i]] = f12_norm[chk]

    age_both = replicate(!values.f_nan,nchk) & sf6_age = age_both & sf6_age_range = replicate(!values.f_nan,2,nchk) & co2_age = age_both & co2_age_range = sf6_age_range
    ai1 = where(dat.mission eq i+1,nai1)
    ai = where(finite(co2_all[chk]) and finite(sf6[chk]),nai)
    if nai gt 0 then age_both[ai] = dat.age_opt_all[ai1]
    bal['both_age_'+vdt_all[i]] = age_both
    ai1 = where(dat_co2.mission eq i+1,nai1)
    ai = where(finite(co2[chk]),nai)
    if nai gt 0 then begin
      co2_age[ai] = dat_co2.age_opt_all[ai1]
      co2_age_range[*,ai] = dat_co2.age_range_all[*,ai1]
      diff1 = co2_age - reform(co2_age_range[0,*])
      diff2 = reform(co2_age_range[1,*]) - co2_age
      chk_d = where(diff1 gt 2,nchkd)
      if nchkd gt 0 then co2_age_range[0,chk_d] = co2_age[chk_d] - diff2[chk_d]
    endif
    bal['co2_age_'+vdt_all[i]] = co2_age
    bal['co2_age_range_'+vdt_all[i]] = co2_age_range
    ai1 = where(dat_sf6.mission eq i+1,nai1)
    ai = where(finite(sf6[chk]),nai)
    if nai gt 0 then begin
      sf6_age[ai] = dat_sf6.age_opt_all[ai1]
      sf6_age_range[*,ai] = dat_sf6.age_range_all[*,ai1]
    endif
    bal['sf6_age_'+vdt_all[i]] = sf6_age
    bal['sf6_age_range_'+vdt_all[i]] = sf6_age_range

    ;  Replace the older data from the first two AirCore flights with the updated v2025 version.
    if i le 1 then begin
      bal['co2_'+flight_dates[i+38]] = bal['co2_'+vdt_all[i]]
      bal['ch4_'+flight_dates[i+38]] = bal['ch4_'+vdt_all[i]]
      bal['ch4_norm_'+flight_dates[i+38]] = bal['ch4_norm_'+vdt_all[i]]
      bal['co2_age_'+flight_dates[i+38]] = bal['co2_age_'+vdt_all[i]]
      bal['co2_age_range_'+flight_dates[i+38]] = bal['co2_age_range_'+vdt_all[i]]
    endif

  endfor

  dates = [dates[0:39],vdt_all[2:24],dates[40:42],vdt_all[25:28],dates[43:44],vdt_all[29:45],dates[45:46],vdt_all[46:115],dates[47],vdt_all[116:-1]]
  yr_frac = [yr_frac[0:39],flight_year_frac[2:24],yr_frac[40:42],flight_year_frac[25:28],yr_frac[43:44],flight_year_frac[29:45],yr_frac[45:46],flight_year_frac[46:115],yr_frac[47], $
    flight_year_frac[116:-1]]

  ;  Quality control
  chk = where(bal['co2_age_'+dates[37]] gt 4.5)
  bal['co2_age_'+dates[37],chk] = !values.f_nan
  bal['co2_age_'+dates[54],0] = !values.f_nan
  save,dates,yr_frac,vdt_all,bal,filename=dir+'Balloon/Engel_profiles_v2025.sav'

endif else begin
  bal = 0  ; avoid IDL 9.1 crash when restoring a hash over an existing hash
  restore,dir+'Balloon/Engel_profiles_v2025.sav'
endelse


;  Add in WAS profiles from 2004-5.
read5 = keyword_set(build_profiles)

if read5 then begin

  restore,dir+'Balloon/WAS_mean_ages.sav'

  was_dates = ['20040929','20051001']
  was_yr_frac = float(strmid(was_dates,0,4)) + (float(strmid(was_dates,4,2))-1.)/12. + (float(strmid(was_dates,6,2))-1.)/365.

  dates = [dates[0:32],was_dates[0],dates[33:-1]]
  yr_frac = [yr_frac[0:32],was_yr_frac[0],yr_frac[33:-1]]

  for i = 0, 0 do begin
    np = n_elements(was_data[was_dates[i]+'_lon'])
    bal['lon_'+was_dates[i]] = was_data[was_dates[i]+'_lon']
    bal['lat_'+was_dates[i]] = was_data[was_dates[i]+'_lat']
    bal['alt_'+was_dates[i]] = was_data[was_dates[i]+'_alt']
    bal['pres_'+was_dates[i]] = was_data[was_dates[i]+'_pres']
    bal['elat_'+was_dates[i]] = was_data[was_dates[i]+'_lat']
    if i eq 0 then bal['sf6_'+was_dates[i]] = was_data[was_dates[i]+'_sf6'] else bal['sf6_'+was_dates[i]] = was_data[was_dates[i]+'_sf6_avg']
    if i eq 0 then bal['co2_'+was_dates[i]] = replicate(!values.f_nan,np) else bal['co2_'+was_dates[i]] = was_data[was_dates[i]+'_co2_eng']
    bal['n2o_'+was_dates[i]] = was_data[was_dates[i]+'_n2o']
    bal['n2o_norm_'+was_dates[i]] = was_data[was_dates[i]+'_n2o_norm']
    bal['ch4_'+was_dates[i]] = was_data[was_dates[i]+'_ch4']
    bal['ch4_norm_'+was_dates[i]] = was_data[was_dates[i]+'_ch4_norm']
    if i eq 0 then begin
      bal['co2_age_'+was_dates[i]] = replicate(!values.f_nan,np)
      bal['co2_age_range_'+was_dates[i]] = replicate(!values.f_nan,np)
      bal['sf6_age_'+was_dates[i]] = was_data[was_dates[i]+'_sf6_age_opt_all']
      bal['sf6_age_range_'+was_dates[i]] = was_data[was_dates[i]+'_sf6_age_range_all']
    endif else begin
      bal['co2_age_'+was_dates[i]] = was_data[was_dates[i]+'_co2_age_opt_all']
      bal['co2_age_range_'+was_dates[i]] = was_data[was_dates[i]+'_co2_age_range_all']
      bal['sf6_age_'+was_dates[i]] = was_data[was_dates[i]+'_sf6_avg_age_opt_all']
      bal['sf6_age_range_'+was_dates[i]] = was_data[was_dates[i]+'_sf6_avg_age_range_all']
    endelse
  endfor
  save,dates,yr_frac,bal,filename=dir+'Balloon/Engel_profiles_v2025.sav'

endif else begin
  bal = 0  ; avoid IDL 9.1 crash when restoring a hash over an existing hash
  restore,dir+'Balloon/Engel_profiles_v2025.sav'
endelse


n_70s = 3
i_70s = indgen(n_70s)
n_80s = 7
i_80s = indgen(n_80s)+3
n_90s = 11
i_90s = [indgen(5)+10,16,17,indgen(4)+20]
n_00s = 9
i_00s = [indgen(5)+24,31,indgen(4)+33]
n_10s = 2
i_10s = [37,64]
n_20s = 1
i_20s = [162]


i_jp = [5,indgen(5)+7,13,14,16,20,22,24,25,27,31,35,36,37,64,162]
i_flask = [i_70s,i_80s,i_90s,i_00s,i_10s,i_20s]
i_was = [33,34]
i_oms = [15,18,19,29,30,32]
i_aircore = [indgen(26)+38,indgen(97)+65,indgen(84)+163]
n_flask = n_elements(i_flask)
n_oms = n_elements(i_oms)
n_aircore = n_elements(i_aircore)
np = n_elements(yr_frac)


;  Quality control based on age-N2O or age-CH4 relationships or altitude profiles compared to OMS ranges.
for t = 0, np-1 do begin
  chk = where(i_flask eq t,nchk)
  if nchk eq 1 then bal['samp_type_'+dates[t]] = 0
  chk = where(i_oms eq t,nchk)
  if nchk eq 1 then bal['samp_type_'+dates[t]] = 1
  chk = where(i_aircore eq t,nchk)
  if nchk eq 1 then bal['samp_type_'+dates[t]] = 2
  chk = where(i_was eq t,nchk)
  if nchk eq 1 then bal['samp_type_'+dates[t]] = 0
  bal['n2o_norm_q_'+dates[t]] = bal['n2o_norm_'+dates[t]]
  bal['ch4_norm_q_'+dates[t]] = bal['ch4_norm_'+dates[t]]
  if ~bal.haskey('f12_norm_'+dates[t]) then bal['f12_norm_'+dates[t]] = replicate(!values.f_nan,n_elements(bal['ch4_norm_'+dates[t]]))
  if bal.haskey('f12_norm_'+dates[t]) then bal['f12_norm_q_'+dates[t]] = bal['f12_norm_'+dates[t]]
  bal['co2_age_q_'+dates[t]] = bal['co2_age_'+dates[t]]
  bal['sf6_age_q_'+dates[t]] = bal['sf6_age_'+dates[t]]

  if t eq 0 or t eq 1 then bal['f12_norm_q_'+dates[t]] = replicate(!values.f_nan,n_elements(bal['ch4_norm_'+dates[t]]))
  if t eq 2 then bal['n2o_norm_q_'+dates[t],3:-1] = bal['n2o_norm_'+dates[t],3:-1] - 0.05
  if t eq 2 then bal['ch4_norm_q_'+dates[t],3:-2] = bal['ch4_norm_'+dates[t],3:-2] + 0.05
  if t eq 3 then bal['ch4_norm_q_'+dates[t]] = replicate(!values.f_nan,n_elements(bal['ch4_norm_'+dates[t]]))
  if t eq 4 then begin
    bal['f12_norm_q_'+dates[t],8] = !values.f_nan
    bal['f12_norm_q_'+dates[t],9] = !values.f_nan
    bal['f12_norm_q_'+dates[t],11] = !values.f_nan
    bal['f12_norm_q_'+dates[t],12] = !values.f_nan
  endif

  ;  Subjective removal of highest point on several flights due to negative vertical gradient in age similar to 19810920 flight.
  if t eq 1 then bal['sf6_age_q_'+dates[t],1] = !values.f_nan
  if t eq 3 then begin
    bal['sf6_age_q_'+dates[t],1] = !values.f_nan
  endif
  if t eq 4 then bal['sf6_age_q_'+dates[t],0] = !values.f_nan
  if t eq 6 then bal['sf6_age_q_'+dates[t],0] = !values.f_nan
  if t eq 16 then begin
    bal['sf6_age_q_'+dates[t],0] = !values.f_nan
    bal['sf6_age_q_'+dates[t],2:3] = !values.f_nan
  endif
  if t eq 26 then begin
    bal['co2_age_q_'+dates[t],9] = !values.f_nan
  endif
  if t eq i_oms[3] or t eq i_oms[4] or t eq i_oms[5] then bal['sf6_age_q_'+dates[t],*] = !values.f_nan
  if t eq i_oms[5] then bal['both_age_'+dates[t],*] = !values.f_nan
endfor


nz2 = 30
dz2 = 1.0
alt_grid2 = findgen(nz2)*dz2+6.0
nz3 = 7
dz3 = 5.0
alt_grid3 = findgen(nz3)*dz3+12.5
nn = 20
dn = 0.05
norm_grid = findgen(nn)*dn+0.025
da1 = 0.5
age_grid1 = findgen(nn)*da1
gg = replicate(!values.f_nan,np,nz2,5) & co2_age_alt_grid = gg & sf6_age_alt_grid = gg & hh = replicate(!values.f_nan,np,nn,5) & co2_age_ch4_bin = hh & sf6_age_ch4_bin = hh
sf6_age_n2o_bin = hh & co2_age_n2o_bin = hh & age_n2o_bin = hh & age_ch4_bin = hh & co2_age_q_alt_grid = gg & sf6_age_q_alt_grid = gg & co2_age_q_ch4_bin = hh & sf6_age_q_ch4_bin = hh
co2_age_q_n2o_bin = hh & sf6_age_q_n2o_bin = hh & age_combo_n2o_bin = hh & age_combo_n2o_bin_nofl = hh & age_combo_n2o_bin_fl = hh & age_q_combo_n2o_bin_fl = hh & age_combo_ch4_bin = hh
age_combo_ch4_bin_nofl = hh & age_q_combo_ch4_bin_fl = hh & age_combo_ch4_bin_fl = hh & age_alt_grid = gg & age_combo_alt_grid = gg & age_combo_alt_grid_nofl = gg & age_combo_alt_grid_fl = gg
n2o_norm_alt_grid = gg & ch4_norm_alt_grid = gg & n2o_norm_alt_grid_nofl = gg & n2o_norm_alt_grid_fl = gg & ch4_norm_alt_grid_nofl = gg & ch4_norm_alt_grid_fl = gg & samp_type_b = intarr(np)
vv = replicate(!values.f_nan,np) & ss = replicate(!values.f_nan,np,nz3,5) & max_alt_b = vv & min_n_norm_b = vv & min_c_norm_b = vv & co2_age_q_alt_grid3 = ss & sf6_age_q_alt_grid3 = ss & age_alt_grid3 = ss
age_combo_alt_grid_nofl3 = ss & n2o_ch4_norm_alt_grid = gg & n2o_from_ch4 = intarr(np) & co2_age_alt_grid_fl = gg & sf6_age_alt_grid_fl = gg & f12_norm_alt_grid = gg & f12_norm_alt_grid_fl = gg
co2_age_f12_bin = hh & sf6_age_f12_bin = hh & co2_age_q_f12_bin = hh & sf6_age_q_f12_bin = hh & age_combo_f12_bin = hh & age_combo_f12_bin_nofl = hh & age_f12_bin = hh & age_q_combo_f12_bin_fl = hh
age_combo_f12_bin_fl = hh & n2o_ch4_norm_alt_grid_nofl = gg & age_combo_n2o_bin_nofl_noch4n2o = hh & age_sf6_n2o_bin_nofl = hh & n2o_on_age_q_bin = hh
age_co2_b_fl = 0. & age_sf6_b_fl = 0. & alt_b_fl = 0. & n2o_b_fl = 0. & ch4_b_fl = 0. & age_co2_b_nofl = 0. & age_sf6_b_nofl = 0. & alt_b_nofl = 0. & n2o_b_nofl = 0. & ch4_b_nofl = 0.
yrfrac_b_fl = 0. & yrfrac_b_nofl = 0. & f12_b_fl = 0.
for t = 0, np-1 do begin
  fi = where(i_flask eq t,nfi)
  if nfi gt 0 then fl = 1 else fl = 0
  if dates[t] eq '20040929' then fl = 0
  alts = bal['alt_'+dates[t]]
  max_alt_b[t] = max(alts)
  npts = n_elements(alts)
  yrfrac_b = replicate(yr_frac[t],npts)
  c_age = bal['co2_age_'+dates[t]]
  s_age = bal['sf6_age_'+dates[t]]
  c_age_q = bal['co2_age_q_'+dates[t]]
  s_age_q = bal['sf6_age_q_'+dates[t]]
  c_age_err = 0.5 * (bal['co2_age_range_'+dates[t],1,*] - bal['co2_age_range_'+dates[t],0,*])
  if bal.haskey('sf6_age_range_'+dates[t]) then s_age_err = 0.5 * (bal['sf6_age_range_'+dates[t],1,*] - bal['sf6_age_range_'+dates[t],0,*]) else s_age_err = replicate(!values.f_nan,npts)
  c_age_weights = 1. / c_age_err^2
  s_age_weights = 1. / s_age_err^2
  samp_type_b[t] = bal['samp_type_'+dates[t]]
  if bal.haskey('both_age_'+dates[t]) then begin
    b_age = bal['both_age_'+dates[t]]
  endif else begin
    if bal.haskey('opt_age_'+dates[t]) then b_age = bal['opt_age_'+dates[t]] else b_age = replicate(!values.f_nan,npts)
  endelse
  f_norm = bal['f12_norm_q_'+dates[t]]
  c_norm = bal['ch4_norm_q_'+dates[t]]
  if bal.haskey('n2o_ch4_norm_'+dates[t]) then n2o_ch4_norm = bal('n2o_ch4_norm_'+dates[t]) else n2o_ch4_norm = replicate(!values.f_nan,npts)
  gd = where(~finite(n2o_ch4_norm),ngd)
  if ngd eq npts then n2o_from_ch4[t] = 0 else n2o_from_ch4[t] = 1
  n_norm = bal['n2o_norm_'+dates[t]]
  n2o_b = bal['n2o_'+dates[t]]
  gd = where(n_norm gt 0,ngd)
  if ngd gt 0 then min_n_norm_b[t] = min(n_norm[gd])
  gd = where(c_norm gt 0)
  min_c_norm_b[t] = min(c_norm[gd])
  if fl then begin
    alt_b_fl = [alt_b_fl,alts]
    yrfrac_b_fl = [yrfrac_b_fl,yrfrac_b]
    age_co2_b_fl = [age_co2_b_fl,c_age]
    age_sf6_b_fl = [age_sf6_b_fl,s_age]
    n2o_b_fl = [n2o_b_fl,n_norm]
    ch4_b_fl = [ch4_b_fl,c_norm]
    f12_b_fl = [f12_b_fl,f_norm]
  endif else begin
    alt_b_nofl = [alt_b_nofl,alts]
    yrfrac_b_nofl = [yrfrac_b_nofl,yrfrac_b]
    age_co2_b_nofl = [age_co2_b_nofl,c_age]
    age_sf6_b_nofl = [age_sf6_b_nofl,s_age_q]
    n2o_b_nofl = [n2o_b_nofl,n_norm]
    ch4_b_nofl = [ch4_b_nofl,c_norm]
  endelse
  for z = 0, nz2-1 do begin
    zi = where(alts ge alt_grid2[z]-dz2/2. and alts lt alt_grid2[z]+dz2/2.,nzi)
    if nzi gt 0 then begin
      gd = where(finite(c_age[zi]),ngd)
      if ngd gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[zi[gd]]))
        wmean = total(c_age_weights[zi[gd]] * c_age[zi[gd]]) / total(c_age_weights[zi[gd]])
        stats = moment(c_age[zi],sd=sd,/nan)
        co2_age_alt_grid[t,z,*] = [wmean,err,stats[0],sd,ngd]
      endif
      gd = where(finite(s_age[zi]),ngd)
      if ngd gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[zi[gd]]))
        wmean = total(s_age_weights[zi[gd]] * s_age[zi[gd]]) / total(s_age_weights[zi[gd]])
        stats = moment(s_age[zi],sd=sd,/nan)
        sf6_age_alt_grid[t,z,*] = [wmean,err,stats[0],sd,ngd]
      endif
      gd = where(finite(c_age_q[zi]),ngd)
      if ngd gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[zi[gd]]))
        wmean = total(c_age_weights[zi[gd]] * c_age_q[zi[gd]]) / total(c_age_weights[zi[gd]])
        stats = moment(c_age_q[zi],sd=sd,/nan)
        co2_age_q_alt_grid[t,z,*] = [wmean,err,stats[0],sd,ngd]
      endif
      gd = where(finite(s_age_q[zi]),ngd)
      if ngd gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[zi[gd]]))
        wmean = total(s_age_weights[zi[gd]] * s_age_q[zi[gd]]) / total(s_age_weights[zi[gd]])
        stats = moment(s_age_q[zi],sd=sd,/nan)
        sf6_age_q_alt_grid[t,z,*] = [wmean,err,stats[0],sd,ngd]
      endif
      gd = where(finite(b_age[zi]),ngd)
      if ngd gt 0 then begin
        stats = moment(b_age[zi],sd=sd,/nan)
        age_alt_grid[t,z,*] = [stats[0],sd,stats[0],sd,ngd]
      endif
      if finite(age_alt_grid[t,z,0]) and ~finite(age_alt_grid[t,z,1]) or age_alt_grid[t,z,1] lt 0.25 then age_alt_grid[t,z,1] = 0.25
      chk = where(n_norm[zi] gt 0 and finite(n_norm[zi]),nchk)
      if nchk gt 0 then begin
        stats = moment(n_norm[zi[chk]],sd=sd,/nan)
        n2o_norm_alt_grid[t,z,0:1] = [stats[0],sd]
      endif
      chk = where(n2o_ch4_norm[zi] gt 0 and finite(n2o_ch4_norm[zi]),nchk)
      if nchk gt 0 then begin
        stats = moment(n2o_ch4_norm[zi[chk]],sd=sd,/nan)
        n2o_ch4_norm_alt_grid[t,z,0:1] = [stats[0],sd]
      endif
      gd = where(c_norm[zi] gt 0 and c_norm[zi] lt 1.05,ngd)
      if ngd gt 0 then begin
        stats = moment(c_norm[zi[gd]],sd=sd,/nan)
        ch4_norm_alt_grid[t,z,0:1] = [stats[0],sd]
      endif
      stats = moment(f_norm[zi],sd=sd,/nan)
      f12_norm_alt_grid[t,z,0:1] = [stats[0],sd]
    endif
    if finite(sf6_age_q_alt_grid[t,z,0]) then age_combo_alt_grid[t,z,*] = sf6_age_q_alt_grid[t,z,*]
    if finite(co2_age_q_alt_grid[t,z,0]) then age_combo_alt_grid[t,z,*] = co2_age_q_alt_grid[t,z,*]
    if finite(age_alt_grid[t,z,0]) then age_combo_alt_grid[t,z,*] = age_alt_grid[t,z,*]
    if fl eq 0 then begin
      if finite(sf6_age_q_alt_grid[t,z,0]) then age_combo_alt_grid_nofl[t,z,*] = sf6_age_q_alt_grid[t,z,*]
      if finite(co2_age_q_alt_grid[t,z,0]) then age_combo_alt_grid_nofl[t,z,*] = co2_age_q_alt_grid[t,z,*]
      if finite(age_alt_grid[t,z,0]) then age_combo_alt_grid_nofl[t,z,*] = age_alt_grid[t,z,*]
      if n2o_norm_alt_grid[t,z,0] gt 0 then n2o_norm_alt_grid_nofl[t,z,*] = n2o_norm_alt_grid[t,z,*]
      if n2o_ch4_norm_alt_grid[t,z,0] gt 0 then n2o_ch4_norm_alt_grid_nofl[t,z,*] = n2o_ch4_norm_alt_grid[t,z,*]
      if ch4_norm_alt_grid[t,z,0] gt 0 and ch4_norm_alt_grid[t,z,0] lt 1.05 then ch4_norm_alt_grid_nofl[t,z,*] = ch4_norm_alt_grid[t,z,*]
    endif else begin
      if finite(sf6_age_q_alt_grid[t,z,0]) then age_combo_alt_grid_fl[t,z,*] = sf6_age_q_alt_grid[t,z,*]
      if finite(sf6_age_q_alt_grid[t,z,0]) then sf6_age_alt_grid_fl[t,z,*] = sf6_age_q_alt_grid[t,z,*]
      if finite(co2_age_q_alt_grid[t,z,0]) then age_combo_alt_grid_fl[t,z,*] = co2_age_q_alt_grid[t,z,*]
      if finite(co2_age_q_alt_grid[t,z,0]) then co2_age_alt_grid_fl[t,z,*] = co2_age_q_alt_grid[t,z,*]
      if finite(age_alt_grid[t,z,0]) then age_combo_alt_grid_fl[t,z,*] = age_alt_grid[t,z,*]
      if n2o_norm_alt_grid[t,z,0] gt 0 then n2o_norm_alt_grid_fl[t,z,*] = n2o_norm_alt_grid[t,z,*]
      if ch4_norm_alt_grid[t,z,0] gt 0 and ch4_norm_alt_grid[t,z,0] lt 1.05 then ch4_norm_alt_grid_fl[t,z,*] = ch4_norm_alt_grid[t,z,*]
      if f12_norm_alt_grid[t,z,0] gt 0 then f12_norm_alt_grid_fl[t,z,*] = f12_norm_alt_grid[t,z,*]
    endelse
  endfor
  for z = 0, nz3-1 do begin
    zi = where(alts ge alt_grid3[z]-dz3/2. and alts lt alt_grid3[z]+dz3/2.,nzi)
    if nzi gt 0 then begin
      gd = where(finite(c_age_q[zi]),ngd)
      if ngd gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[zi[gd]]))
        wmean = total(c_age_weights[zi[gd]] * c_age_q[zi[gd]]) / total(c_age_weights[zi[gd]])
        stats = moment(c_age_q[zi],sd=sd,/nan)
        co2_age_q_alt_grid3[t,z,*] = [wmean,err,stats[0],sd,ngd]
      endif
      gd = where(finite(s_age_q[zi]),ngd)
      if ngd gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[zi[gd]]))
        wmean = total(s_age_weights[zi[gd]] * s_age_q[zi[gd]]) / total(s_age_weights[zi[gd]])
        stats = moment(s_age_q[zi],sd=sd,/nan)
        sf6_age_q_alt_grid3[t,z,*] = [wmean,err,stats[0],sd,ngd]
      endif
      gd = where(finite(b_age[zi]),ngd)
      if ngd gt 0 then begin
        stats = moment(b_age[zi],sd=sd,/nan)
        age_alt_grid3[t,z,*] = [stats[0],sd,stats[0],sd,ngd]
      endif
      if finite(age_alt_grid3[t,z,0]) and ~finite(age_alt_grid3[t,z,1]) or age_alt_grid3[t,z,1] lt 0.25 then age_alt_grid3[t,z,1] = 0.25
    endif
    if fl eq 0 then begin
      if finite(sf6_age_q_alt_grid3[t,z,0]) then age_combo_alt_grid_nofl3[t,z,*] = sf6_age_q_alt_grid3[t,z,*]
      if finite(co2_age_q_alt_grid3[t,z,0]) then age_combo_alt_grid_nofl3[t,z,*] = co2_age_q_alt_grid3[t,z,*]
      if finite(age_alt_grid3[t,z,0]) then age_combo_alt_grid_nofl3[t,z,*] = age_alt_grid3[t,z,*]
    endif
  endfor
  for i = 0, nn-1 do begin
    gd = where(c_norm ge norm_grid[i] - dn/2. and c_norm lt norm_grid[i] + dn/2.,ngd)
    if ngd ge 1 then begin
      gd2 = where(finite(c_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[gd[gd2]]))
        wmean = total(c_age_weights[gd[gd2]] * c_age[gd[gd2]]) / total(c_age_weights[gd[gd2]])
        stats = moment(c_age[gd],sdev=sdev,/nan)
        co2_age_ch4_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(s_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[gd[gd2]]))
        wmean = total(s_age_weights[gd[gd2]] * s_age[gd[gd2]]) / total(s_age_weights[gd[gd2]])
        stats = moment(s_age[gd],sdev=sdev,/nan)
        sf6_age_ch4_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(c_age_q[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[gd[gd2]]))
        wmean = total(c_age_weights[gd[gd2]] * c_age_q[gd[gd2]]) / total(c_age_weights[gd[gd2]])
        stats = moment(c_age_q[gd],sdev=sdev,/nan)
        co2_age_q_ch4_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(s_age_q[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[gd[gd2]]))
        wmean = total(s_age_weights[gd[gd2]] * s_age_q[gd[gd2]]) / total(s_age_weights[gd[gd2]])
        stats = moment(s_age_q[gd],sdev=sdev,/nan)
        sf6_age_q_ch4_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(b_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = min([co2_age_q_ch4_bin[t,i,1],sf6_age_q_ch4_bin[t,i,1]])
        stats = moment(b_age[gd],sdev=sdev,/nan)
        age_ch4_bin[t,i,*] = [stats[0],err,stats[0],sdev,ngd2]
      endif
    endif
    gd = where(n_norm ge norm_grid[i] - dn/2. and n_norm lt norm_grid[i] + dn/2.,ngd)
    if ngd ge 1 then begin
      gd2 = where(finite(c_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[gd[gd2]]))
        wmean = total(c_age_weights[gd[gd2]] * c_age[gd[gd2]]) / total(c_age_weights[gd[gd2]])
        stats = moment(c_age[gd],sdev=sdev,/nan)
        co2_age_n2o_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(s_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[gd[gd2]]))
        wmean = total(s_age_weights[gd[gd2]] * s_age[gd[gd2]]) / total(s_age_weights[gd[gd2]])
        stats = moment(s_age[gd],sdev=sdev,/nan)
        sf6_age_n2o_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(c_age_q[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[gd[gd2]]))
        wmean = total(c_age_weights[gd[gd2]] * c_age_q[gd[gd2]]) / total(c_age_weights[gd[gd2]])
        stats = moment(c_age_q[gd],sdev=sdev,/nan)
        co2_age_q_n2o_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(s_age_q[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[gd[gd2]]))
        wmean = total(s_age_weights[gd[gd2]] * s_age_q[gd[gd2]]) / total(s_age_weights[gd[gd2]])
        stats = moment(s_age_q[gd],sdev=sdev,/nan)
        sf6_age_q_n2o_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(b_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = min([co2_age_q_n2o_bin[t,i,1],sf6_age_q_n2o_bin[t,i,1]])
        stats = moment(b_age[gd],sdev=sdev,/nan)
        age_n2o_bin[t,i,*] = [stats[0],err,stats[0],sdev,ngd2]
      endif
    endif
    gd = where(f_norm ge norm_grid[i] - dn/2. and f_norm lt norm_grid[i] + dn/2.,ngd)
    if ngd ge 1 then begin
      gd2 = where(finite(c_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[gd[gd2]]))
        wmean = total(c_age_weights[gd[gd2]] * c_age[gd[gd2]]) / total(c_age_weights[gd[gd2]])
        stats = moment(c_age[gd],sdev=sdev,/nan)
        co2_age_f12_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(s_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[gd[gd2]]))
        wmean = total(s_age_weights[gd[gd2]] * s_age[gd[gd2]]) / total(s_age_weights[gd[gd2]])
        stats = moment(s_age[gd],sdev=sdev,/nan)
        sf6_age_f12_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(c_age_q[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(c_age_weights[gd[gd2]]))
        wmean = total(c_age_weights[gd[gd2]] * c_age_q[gd[gd2]]) / total(c_age_weights[gd[gd2]])
        stats = moment(c_age_q[gd],sdev=sdev,/nan)
        co2_age_q_f12_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(s_age_q[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = 1. / sqrt(total(s_age_weights[gd[gd2]]))
        wmean = total(s_age_weights[gd[gd2]] * s_age_q[gd[gd2]]) / total(s_age_weights[gd[gd2]])
        stats = moment(s_age_q[gd],sdev=sdev,/nan)
        sf6_age_q_f12_bin[t,i,*] = [wmean,err,stats[0],sd,ngd2]
      endif
      gd2 = where(finite(b_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        err = min([co2_age_q_f12_bin[t,i,1],sf6_age_q_f12_bin[t,i,1]])
        stats = moment(b_age[gd],sdev=sdev,/nan)
        age_f12_bin[t,i,*] = [stats[0],err,stats[0],sdev,ngd2]
      endif
    endif
    if finite(sf6_age_q_ch4_bin[t,i,0]) then age_combo_ch4_bin[t,i,*] = sf6_age_q_ch4_bin[t,i,*]
    if finite(co2_age_q_ch4_bin[t,i,0]) then age_combo_ch4_bin[t,i,*] = co2_age_q_ch4_bin[t,i,*]
    if finite(age_ch4_bin[t,i,0]) then age_combo_ch4_bin[t,i,*] = age_ch4_bin[t,i,*]
    if finite(sf6_age_q_n2o_bin[t,i,0]) then age_combo_n2o_bin[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
    if finite(co2_age_q_n2o_bin[t,i,0]) then age_combo_n2o_bin[t,i,*] = co2_age_q_n2o_bin[t,i,*]
    if finite(age_n2o_bin[t,i,0]) then age_combo_n2o_bin[t,i,*] = age_n2o_bin[t,i,*]
    if finite(sf6_age_q_f12_bin[t,i,0]) then age_combo_f12_bin[t,i,*] = sf6_age_q_f12_bin[t,i,*]
    if finite(co2_age_q_f12_bin[t,i,0]) then age_combo_f12_bin[t,i,*] = co2_age_q_f12_bin[t,i,*]
    if finite(age_f12_bin[t,i,0]) then age_combo_f12_bin[t,i,*] = age_f12_bin[t,i,*]
    if fl eq 0 then begin
      if finite(sf6_age_q_ch4_bin[t,i,0]) then age_combo_ch4_bin_nofl[t,i,*] = sf6_age_q_ch4_bin[t,i,*]
      if finite(co2_age_q_ch4_bin[t,i,0]) then age_combo_ch4_bin_nofl[t,i,*] = co2_age_q_ch4_bin[t,i,*]
      if finite(age_ch4_bin[t,i,0]) then age_combo_ch4_bin_nofl[t,i,*] = age_ch4_bin[t,i,*]
      if ~finite(age_combo_ch4_bin_nofl[t,i,1]) then age_combo_ch4_bin_nofl[t,i,1] = 0.1
      if finite(sf6_age_q_n2o_bin[t,i,0]) then age_combo_n2o_bin_nofl[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
      if finite(co2_age_q_n2o_bin[t,i,0]) then age_combo_n2o_bin_nofl[t,i,*] = co2_age_q_n2o_bin[t,i,*]
      if finite(age_n2o_bin[t,i,0]) then age_combo_n2o_bin_nofl[t,i,*] = age_n2o_bin[t,i,*]
      if finite(sf6_age_q_n2o_bin[t,i,0]) then age_sf6_n2o_bin_nofl[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
      if finite(sf6_age_q_n2o_bin[t,i,0]) and n2o_from_ch4[t] eq 0 then age_combo_n2o_bin_nofl_noch4n2o[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
      if finite(co2_age_q_n2o_bin[t,i,0]) and n2o_from_ch4[t] eq 0 then age_combo_n2o_bin_nofl_noch4n2o[t,i,*] = co2_age_q_n2o_bin[t,i,*]
      if finite(age_n2o_bin[t,i,0]) and n2o_from_ch4[t] eq 0 then age_combo_n2o_bin_nofl_noch4n2o[t,i,*] = age_n2o_bin[t,i,*]
      if finite(sf6_age_q_f12_bin[t,i,0]) then age_combo_f12_bin_nofl[t,i,*] = sf6_age_q_f12_bin[t,i,*]
      if finite(co2_age_q_f12_bin[t,i,0]) then age_combo_f12_bin_nofl[t,i,*] = co2_age_q_f12_bin[t,i,*]
      if finite(age_f12_bin[t,i,0]) then age_combo_f12_bin_nofl[t,i,*] = age_f12_bin[t,i,*]
    endif else begin
      if finite(sf6_age_q_ch4_bin[t,i,0]) then age_q_combo_ch4_bin_fl[t,i,*] = sf6_age_q_ch4_bin[t,i,*]
      if finite(co2_age_q_ch4_bin[t,i,0]) then age_q_combo_ch4_bin_fl[t,i,*] = co2_age_q_ch4_bin[t,i,*]
      if finite(age_ch4_bin[t,i,0]) then age_q_combo_ch4_bin_fl[t,i,*] = age_ch4_bin[t,i,*]
      if finite(sf6_age_ch4_bin[t,i,0]) then age_combo_ch4_bin_fl[t,i,*] = sf6_age_ch4_bin[t,i,*]
      if finite(co2_age_ch4_bin[t,i,0]) then age_combo_ch4_bin_fl[t,i,*] = co2_age_ch4_bin[t,i,*]
      if finite(age_ch4_bin[t,i,0]) then age_combo_ch4_bin_fl[t,i,*] = age_ch4_bin[t,i,*]
      if finite(sf6_age_q_n2o_bin[t,i,0]) then age_q_combo_n2o_bin_fl[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
      if finite(co2_age_q_n2o_bin[t,i,0]) then age_q_combo_n2o_bin_fl[t,i,*] = co2_age_q_n2o_bin[t,i,*]
      if finite(age_n2o_bin[t,i,0]) then age_q_combo_n2o_bin_fl[t,i,*] = age_n2o_bin[t,i,*]
      if finite(sf6_age_n2o_bin[t,i,0]) then age_combo_n2o_bin_fl[t,i,*] = sf6_age_n2o_bin[t,i,*]
      if finite(co2_age_n2o_bin[t,i,0]) then age_combo_n2o_bin_fl[t,i,*] = co2_age_n2o_bin[t,i,*]
      if finite(age_n2o_bin[t,i,0]) then age_combo_n2o_bin_fl[t,i,*] = age_n2o_bin[t,i,*]
      if finite(sf6_age_q_f12_bin[t,i,0]) then age_q_combo_f12_bin_fl[t,i,*] = sf6_age_q_f12_bin[t,i,*]
      if finite(co2_age_q_f12_bin[t,i,0]) then age_q_combo_f12_bin_fl[t,i,*] = co2_age_q_f12_bin[t,i,*]
      if finite(age_f12_bin[t,i,0]) then age_q_combo_f12_bin_fl[t,i,*] = age_f12_bin[t,i,*]
      if finite(sf6_age_f12_bin[t,i,0]) then age_combo_f12_bin_fl[t,i,*] = sf6_age_f12_bin[t,i,*]
      if finite(co2_age_f12_bin[t,i,0]) then age_combo_f12_bin_fl[t,i,*] = co2_age_f12_bin[t,i,*]
      if finite(age_f12_bin[t,i,0]) then age_combo_f12_bin_fl[t,i,*] = age_f12_bin[t,i,*]
    endelse
  endfor
  for i = 0, nn-1 do begin
    gd = where((c_age_q ge age_grid1[i] - da1/2. and c_age_q lt age_grid1[i] + da1/2.) or (s_age_q ge age_grid1[i] - da1/2. and s_age_q lt age_grid1[i] + da1/2.),ngd)
    if ngd ge 1 then begin
      gd2 = where(finite(n2o_b[gd]),ngd2)
      if ngd2 gt 0 then begin
        stats = moment(n2o_b[gd],sdev=sdev,/nan)
        n2o_on_age_q_bin[t,i,0:2] = [stats[0],sdev,ngd2]
      endif
    endif
  endfor
endfor

age_co2_b_fl = age_co2_b_fl[1:-1] & age_sf6_b_fl = age_sf6_b_fl[1:-1] & alt_b_fl = alt_b_fl[1:-1] & n2o_b_fl = n2o_b_fl[1:-1] & ch4_b_fl = ch4_b_fl[1:-1] & f12_b_fl = f12_b_fl[1:-1]
age_co2_b_nofl = age_co2_b_nofl[1:-1] & age_sf6_b_nofl = age_sf6_b_nofl[1:-1] & alt_b_nofl = alt_b_nofl[1:-1] & n2o_b_nofl = n2o_b_nofl[1:-1] & ch4_b_nofl = ch4_b_nofl[1:-1]


;  Make flask average profiles from 1970s-2000s.
aa = replicate(!values.f_nan,nz2,5) & n2o_norm_alt_grid_fl_avg = aa & ch4_norm_alt_grid_fl_avg = aa & f12_norm_alt_grid_fl_avg = aa & co2_age_alt_grid_fl_avg = aa & sf6_age_alt_grid_fl_avg = aa
bb = replicate(!values.f_nan,nn,5) & co2_age_q_n2o_bin_fl_avg = bb & co2_age_q_ch4_bin_fl_avg = bb & age_q_combo_n2o_bin_fl_avg = bb & age_q_combo_ch4_bin_fl_avg = bb
co2_age_q_f12_bin_fl_avg = bb & age_q_combo_f12_bin_fl_avg = bb
ti = where(yr_frac lt 2010)
for z = 0, nz2-1 do begin
  stats = moment(n2o_norm_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  n2o_norm_alt_grid_fl_avg[z,0:1] = [stats[0],sdev]
  stats = moment(ch4_norm_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  ch4_norm_alt_grid_fl_avg[z,0:1] = [stats[0],sdev]
  stats = moment(f12_norm_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  f12_norm_alt_grid_fl_avg[z,0:1] = [stats[0],sdev]
  gd = where(finite(co2_age_alt_grid_fl[ti,z,0]),ngd)
  err = 1. / sqrt(total(co2_age_alt_grid_fl[ti[gd],z,1]^2))
  wmean = total(1. / co2_age_alt_grid_fl[ti[gd],z,1]^2 * co2_age_alt_grid_fl[ti[gd],z,0]) / total(1. / co2_age_alt_grid_fl[ti[gd],z,1]^2)
  npts = total(co2_age_alt_grid_fl[ti[gd],z,4])
  stats = moment(co2_age_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  co2_age_alt_grid_fl_avg[z,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(sf6_age_alt_grid_fl[ti,z,0]),ngd)
  err = 1. / sqrt(total(sf6_age_alt_grid_fl[ti[gd],z,1]^2))
  wmean = total(1. / sf6_age_alt_grid_fl[ti[gd],z,1]^2 * sf6_age_alt_grid_fl[ti[gd],z,0]) / total(1. / sf6_age_alt_grid_fl[ti[gd],z,1]^2)
  npts = total(sf6_age_alt_grid_fl[ti[gd],z,4])
  stats = moment(sf6_age_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  sf6_age_alt_grid_fl_avg[z,*] = [wmean,err,stats[0],sdev,npts]
endfor
for i = 0, nn-1 do begin
  gd = where(finite(co2_age_q_n2o_bin[ti,i,0]),ngd)
  err = 1. / sqrt(total(co2_age_q_n2o_bin[ti[gd],i,1]^2))
  wmean = total(1. / co2_age_q_n2o_bin[ti[gd],i,1]^2 * co2_age_q_n2o_bin[ti[gd],i,0]) / total(1. / co2_age_q_n2o_bin[ti[gd],i,1]^2)
  npts = total(co2_age_q_n2o_bin[ti[gd],i,4])
  stats = moment(co2_age_q_n2o_bin[ti,i,0],sdev=sdev,/nan)
  co2_age_q_n2o_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(age_q_combo_n2o_bin_fl[ti,i,0]),ngd)
  err = 1. / sqrt(total(age_q_combo_n2o_bin_fl[ti[gd],i,1]^2))
  wmean = total(1. / age_q_combo_n2o_bin_fl[ti[gd],i,1]^2 * age_q_combo_n2o_bin_fl[ti[gd],i,0]) / total(1. / age_q_combo_n2o_bin_fl[ti[gd],i,1]^2)
  npts = total(age_q_combo_n2o_bin_fl[ti[gd],i,4])
  stats = moment(age_q_combo_n2o_bin_fl[ti,i,0],sdev=sdev,/nan)
  age_q_combo_n2o_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(co2_age_q_ch4_bin[ti,i,0]),ngd)
  err = 1. / sqrt(total(co2_age_q_ch4_bin[ti[gd],i,1]^2))
  wmean = total(1. / co2_age_q_ch4_bin[ti[gd],i,1]^2 * co2_age_q_ch4_bin[ti[gd],i,0]) / total(1. / co2_age_q_ch4_bin[ti[gd],i,1]^2)
  npts = total(co2_age_q_ch4_bin[ti[gd],i,4])
  stats = moment(co2_age_q_ch4_bin[ti,i,0],sdev=sdev,/nan)
  co2_age_q_ch4_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(age_q_combo_ch4_bin_fl[ti,i,0]),ngd)
  err = 1. / sqrt(total(age_q_combo_ch4_bin_fl[ti[gd],i,1]^2))
  wmean = total(1. / age_q_combo_ch4_bin_fl[ti[gd],i,1]^2 * age_q_combo_ch4_bin_fl[ti[gd],i,0]) / total(1. / age_q_combo_ch4_bin_fl[ti[gd],i,1]^2)
  npts = total(age_q_combo_ch4_bin_fl[ti[gd],i,4])
  stats = moment(age_q_combo_ch4_bin_fl[ti,i,0],sdev=sdev,/nan)
  age_q_combo_ch4_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(co2_age_q_f12_bin[ti,i,0]),ngd)
  err = 1. / sqrt(total(co2_age_q_f12_bin[ti[gd],i,1]^2))
  wmean = total(1. / co2_age_q_f12_bin[ti[gd],i,1]^2 * co2_age_q_f12_bin[ti[gd],i,0]) / total(1. / co2_age_q_f12_bin[ti[gd],i,1]^2)
  npts = total(co2_age_q_f12_bin[ti[gd],i,4])
  stats = moment(co2_age_q_f12_bin[ti,i,0],sdev=sdev,/nan)
  co2_age_q_f12_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(age_q_combo_f12_bin_fl[ti,i,0]),ngd)
  err = 1. / sqrt(total(age_q_combo_f12_bin_fl[ti[gd],i,1]^2))
  wmean = total(1. / age_q_combo_f12_bin_fl[ti[gd],i,1]^2 * age_q_combo_f12_bin_fl[ti[gd],i,0]) / total(1. / age_q_combo_f12_bin_fl[ti[gd],i,1]^2)
  npts = total(age_q_combo_f12_bin_fl[ti[gd],i,4])
  stats = moment(age_q_combo_f12_bin_fl[ti,i,0],sdev=sdev,/nan)
  age_q_combo_f12_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
endfor

nsamp_flask = fltarr(nz2) & nsamp_insitu = nsamp_flask
for z = 0, nz2-1 do begin
  nsamp_flask[z] = total(age_combo_alt_grid_fl[*,z,4],/nan)
endfor


;  Read in the AirCore Ch4-N2O correlation
restore,dir+'Balloon/Aircore_n2o_ch4_v2024.sav'
chk = where(~finite(ch4_norm_n2o_norm_ac_sm),nchk)
if nchk gt 0 then ch4_norm_n2o_norm_ac_sm[chk] = 0.2

nnorm = n_elements(n2o_norm_grid)
ni = interpol(indgen(nnorm),n2o_norm_grid,norm_grid)
ch4_norm_n2o_norm_ac_sm_coarse = interpolate(ch4_norm_n2o_norm_ac_sm,ni)
ch4_norm_n2o_norm_ac_sm_coarse[0:10] = [replicate(0.01,6),0.04,0.07,0.13,0.21,0.29]

ni = interpol(indgen(nn),ch4_norm_n2o_norm_ac_sm_coarse,norm_grid)
n2o_norm_ch4_norm_ac_sm_coarse = interpolate(norm_grid,ni)


nt = 2644
dt = 1./48.
years = findgen(nt)*dt+1975.
kk = replicate(!values.f_nan,nt,nn,5) & co2_age_n2o_bin_yrs = kk & sf6_age_n2o_bin_yrs = kk & age_combo_n2o_bin_yrs = kk & co2_age_ch4_bin_yrs = kk & sf6_age_ch4_bin_yrs = kk
age_combo_n2o_bin_nofl_yrs = kk & age_combo_n2o_bin_fl_yrs = kk & age_q_combo_n2o_bin_fl_yrs = kk & age_combo_ch4_bin_yrs = kk & age_combo_ch4_bin_nofl_yrs = kk
mm = replicate(!values.f_nan,nt,nz2,5) & age_combo_ch4_bin_fl_yrs = kk & age_q_combo_ch4_bin_fl_yrs = kk & age_combo_alt_grid_yrs = mm & age_combo_alt_grid_nofl_yrs = mm
age_combo_alt_grid_fl_yrs = mm & n2o_norm_alt_grid_fl_yrs = mm & ch4_norm_alt_grid_fl_yrs = mm & n2o_norm_alt_grid_nofl_yrs = mm & ch4_norm_alt_grid_nofl_yrs = mm
rr = replicate(!values.f_nan,nt) & ee = replicate(!values.f_nan,nt,nz3,5) & samp_type_b_yrs = rr & max_alt_b_yrs = rr & min_n_norm_b_yrs = rr & min_c_norm_b_yrs = rr & n2o_from_ch4_yrs = rr
age_combo_alt_grid_nofl3_yrs = ee & n2o_ch4_norm_alt_grid_nofl_yrs = mm & age_combo_n2o_bin_nofl_noch4n2o_yrs = kk & age_sf6_n2o_bin_nofl_yrs = kk & n2o_norm_alt_grid_yrs = mm
n2o_on_age_bin_nofl_yrs = kk
for t = 0, nt-1 do begin
  gd = where(yr_frac ge years[t]-dt/2. and yr_frac lt years[t]+dt/2.,ni1)
  if ni1 ge 1 then begin
    samp_type_b_yrs[t] = samp_type_b[gd[0]]
    max_alt_b_yrs[t] = max(max_alt_b[gd])
    min_n_norm_b_yrs[t] = max(min_n_norm_b[gd])
    min_c_norm_b_yrs[t] = max(min_c_norm_b[gd])
    n2o_from_ch4_yrs[t] = n2o_from_ch4[gd[0]]
    for i = 0, nn-1 do begin
      gd2 = where(finite(co2_age_ch4_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./co2_age_ch4_bin[gd[gd2],i,1]^2))
      wmean = total(1. / co2_age_ch4_bin[gd[gd2],i,1]^2 * co2_age_ch4_bin[gd[gd2],i,0]) / total(1. / co2_age_ch4_bin[gd[gd2],i,1]^2)
      npts = total(sf6_age_ch4_bin[gd[gd2],i,4])
      stats = moment(co2_age_ch4_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(co2_age_ch4_bin[gd,i,1],sdev=sdev,/nan)
      co2_age_ch4_bin_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(sf6_age_ch4_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./sf6_age_ch4_bin[gd[gd2],i,1]^2))
      wmean = total(1. / sf6_age_ch4_bin[gd[gd2],i,1]^2 * sf6_age_ch4_bin[gd[gd2],i,0]) / total(1. / sf6_age_ch4_bin[gd[gd2],i,1]^2)
      npts = total(sf6_age_ch4_bin[gd[gd2],i,4])
      stats = moment(sf6_age_ch4_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(sf6_age_ch4_bin[gd,i,1],sdev=sdev,/nan)
      sf6_age_ch4_bin_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(co2_age_q_n2o_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./co2_age_q_n2o_bin[gd[gd2],i,1]^2))
      wmean = total(1. / co2_age_q_n2o_bin[gd[gd2],i,1]^2 * co2_age_q_n2o_bin[gd[gd2],i,0]) / total(1. / co2_age_q_n2o_bin[gd[gd2],i,1]^2)
      npts = total(co2_age_q_n2o_bin[gd[gd2],i,4])
      stats = moment(co2_age_q_n2o_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(co2_age_q_n2o_bin[gd,i,1],sdev=sdev,/nan)
      co2_age_n2o_bin_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(sf6_age_q_n2o_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./sf6_age_q_n2o_bin[gd[gd2],i,1]^2))
      wmean = total(1. / sf6_age_q_n2o_bin[gd[gd2],i,1]^2 * sf6_age_q_n2o_bin[gd[gd2],i,0]) / total(1. / sf6_age_q_n2o_bin[gd[gd2],i,1]^2)
      npts = total(sf6_age_q_n2o_bin[gd[gd2],i,4])
      stats = moment(sf6_age_q_n2o_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(sf6_age_q_n2o_bin[gd,i,1],sdev=sdev,/nan)
      sf6_age_n2o_bin_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_ch4_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_ch4_bin[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_ch4_bin[gd[gd2],i,1]^2 * age_combo_ch4_bin[gd[gd2],i,0]) / total(1. / age_combo_ch4_bin[gd[gd2],i,1]^2)
      npts = total(age_combo_ch4_bin[gd[gd2],i,4])
      stats = moment(age_combo_ch4_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_ch4_bin[gd,i,1],sdev=sdev,/nan)
      age_combo_ch4_bin_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_ch4_bin_nofl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_ch4_bin_nofl[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_ch4_bin_nofl[gd[gd2],i,1]^2 * age_combo_ch4_bin_nofl[gd[gd2],i,0]) / total(1. / age_combo_ch4_bin_nofl[gd[gd2],i,1]^2)
      npts = total(age_combo_ch4_bin_nofl[gd[gd2],i,4])
      stats = moment(age_combo_ch4_bin_nofl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_ch4_bin_nofl[gd,i,1],sdev=sdev,/nan)
      age_combo_ch4_bin_nofl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_ch4_bin_fl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_ch4_bin_fl[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_ch4_bin_fl[gd[gd2],i,1]^2 * age_combo_ch4_bin_fl[gd[gd2],i,0]) / total(1. / age_combo_ch4_bin_fl[gd[gd2],i,1]^2)
      npts = total(age_combo_ch4_bin_fl[gd[gd2],i,4])
      stats = moment(age_combo_ch4_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_ch4_bin_fl[gd,i,1],sdev=sdev,/nan)
      age_combo_ch4_bin_fl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_q_combo_ch4_bin_fl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_q_combo_ch4_bin_fl[gd[gd2],i,1]^2))
      wmean = total(1. / age_q_combo_ch4_bin_fl[gd[gd2],i,1]^2 * age_q_combo_ch4_bin_fl[gd[gd2],i,0]) / total(1. / age_q_combo_ch4_bin_fl[gd[gd2],i,1]^2)
      npts = total(age_q_combo_ch4_bin_fl[gd[gd2],i,4])
      stats = moment(age_q_combo_ch4_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_q_combo_ch4_bin_fl[gd,i,1],sdev=sdev,/nan)
      age_q_combo_ch4_bin_fl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_n2o_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_n2o_bin[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_n2o_bin[gd[gd2],i,1]^2 * age_combo_n2o_bin[gd[gd2],i,0]) / total(1. / age_combo_n2o_bin[gd[gd2],i,1]^2)
      npts = total(age_combo_n2o_bin[gd[gd2],i,4])
      stats = moment(age_combo_n2o_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin[gd,i,1],sdev=sdev,/nan)
      age_combo_n2o_bin_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_n2o_bin_nofl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_n2o_bin_nofl[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_n2o_bin_nofl[gd[gd2],i,1]^2 * age_combo_n2o_bin_nofl[gd[gd2],i,0]) / total(1. / age_combo_n2o_bin_nofl[gd[gd2],i,1]^2)
      npts = total(age_combo_n2o_bin_nofl[gd[gd2],i,4])
      stats = moment(age_combo_n2o_bin_nofl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin_nofl[gd,i,1],sdev=sdev,/nan)
      age_combo_n2o_bin_nofl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_sf6_n2o_bin_nofl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_sf6_n2o_bin_nofl[gd[gd2],i,1]^2))
      wmean = total(1. / age_sf6_n2o_bin_nofl[gd[gd2],i,1]^2 * age_sf6_n2o_bin_nofl[gd[gd2],i,0]) / total(1. / age_sf6_n2o_bin_nofl[gd[gd2],i,1]^2)
      npts = total(age_sf6_n2o_bin_nofl[gd[gd2],i,4])
      stats = moment(age_sf6_n2o_bin_nofl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_sf6_n2o_bin_nofl[gd,i,1],sdev=sdev,/nan)
      age_sf6_n2o_bin_nofl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_n2o_bin_nofl_noch4n2o[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_n2o_bin_nofl_noch4n2o[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_n2o_bin_nofl_noch4n2o[gd[gd2],i,1]^2 * age_combo_n2o_bin_nofl_noch4n2o[gd[gd2],i,0]) / total(1. / age_combo_n2o_bin_nofl_noch4n2o[gd[gd2],i,1]^2)
      npts = total(age_combo_n2o_bin_nofl_noch4n2o[gd[gd2],i,4])
      stats = moment(age_combo_n2o_bin_nofl_noch4n2o[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin_nofl_noch4n2o[gd,i,1],sdev=sdev,/nan)
      age_combo_n2o_bin_nofl_noch4n2o_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_n2o_bin_fl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_n2o_bin_fl[gd[gd2],i,1]^2))
      wmean = total(1. / age_combo_n2o_bin_fl[gd[gd2],i,1]^2 * age_combo_n2o_bin_fl[gd[gd2],i,0]) / total(1. / age_combo_n2o_bin_fl[gd[gd2],i,1]^2)
      npts = total(age_combo_n2o_bin_fl[gd[gd2],i,4])
      stats = moment(age_combo_n2o_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin_fl[gd,i,1],sdev=sdev,/nan)
      age_combo_n2o_bin_fl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_q_combo_n2o_bin_fl[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./age_q_combo_n2o_bin_fl[gd[gd2],i,1]^2))
      wmean = total(1. / age_q_combo_n2o_bin_fl[gd[gd2],i,1]^2 * age_q_combo_n2o_bin_fl[gd[gd2],i,0]) / total(1. / age_q_combo_n2o_bin_fl[gd[gd2],i,1]^2)
      npts = total(age_q_combo_n2o_bin_fl[gd[gd2],i,4])
      stats = moment(age_q_combo_n2o_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_q_combo_n2o_bin_fl[gd,i,1],sdev=sdev,/nan)
      age_q_combo_n2o_bin_fl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]

      gd2 = where(finite(n2o_on_age_q_bin[gd,i,0]),ngd2)
      err = 1. / sqrt(total(1./n2o_on_age_q_bin[gd[gd2],i,1]^2))
      wmean = total(1. / n2o_on_age_q_bin[gd[gd2],i,1]^2 * n2o_on_age_q_bin[gd[gd2],i,0]) / total(1. / n2o_on_age_q_bin[gd[gd2],i,1]^2)
      npts = total(n2o_on_age_q_bin[gd[gd2],i,2])
      stats = moment(n2o_on_age_q_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(n2o_on_age_q_bin[gd,i,1],sdev=sdev,/nan)
      n2o_on_age_bin_nofl_yrs[t,i,*] = [wmean,err,stats[0],stats2[0],npts]
    endfor
    for z = 0, nz2-1 do begin
      gd2 = where(finite(age_combo_alt_grid[gd,z,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_alt_grid[gd[gd2],z,1]^2))
      wmean = total(1. / age_combo_alt_grid[gd[gd2],z,1]^2 * age_combo_alt_grid[gd[gd2],z,0]) / total(1. / age_combo_alt_grid[gd[gd2],z,1]^2)
      npts = total(age_combo_alt_grid[gd[gd2],z,4])
      stats = moment(age_combo_alt_grid[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid[gd,z,1],sdev=sdev,/nan)
      age_combo_alt_grid_yrs[t,z,*] = [wmean,err,stats[0],stats2[0],npts]
      gd2 = where(finite(age_combo_alt_grid_nofl[gd,z,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_alt_grid_nofl[gd[gd2],z,1]^2))
      wmean = total(1. / age_combo_alt_grid_nofl[gd[gd2],z,1]^2 * age_combo_alt_grid_nofl[gd[gd2],z,0]) / total(1. / age_combo_alt_grid_nofl[gd[gd2],z,1]^2)
      npts = total(age_combo_alt_grid_nofl[gd[gd2],z,4])
      stats = moment(age_combo_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      age_combo_alt_grid_nofl_yrs[t,z,*] = [wmean,err,stats[0],stats2[0],npts]
      if finite(age_combo_alt_grid_nofl_yrs[t,z,0]) and ~finite(age_combo_alt_grid_nofl_yrs[t,z,3]) or age_combo_alt_grid_nofl_yrs[t,z,3] lt 0.25 then age_combo_alt_grid_nofl_yrs[t,z,3] = 0.25
      gd2 = where(finite(age_combo_alt_grid_fl[gd,z,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_alt_grid_fl[gd[gd2],z,1]^2))
      wmean = total(1. / age_combo_alt_grid_fl[gd[gd2],z,1]^2 * age_combo_alt_grid_fl[gd[gd2],z,0]) / total(1. / age_combo_alt_grid_fl[gd[gd2],z,1]^2)
      npts = total(age_combo_alt_grid_fl[gd[gd2],z,4])
      stats = moment(age_combo_alt_grid_fl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid_fl[gd,z,1],sdev=sdev,/nan)
      age_combo_alt_grid_fl_yrs[t,z,*] = [wmean,err,stats[0],stats2[0],npts]
      if finite(age_combo_alt_grid_fl_yrs[t,z,0]) and ~finite(age_combo_alt_grid_fl_yrs[t,z,3]) or age_combo_alt_grid_fl_yrs[t,z,3] lt 0.25 then age_combo_alt_grid_fl_yrs[t,z,3] = 0.25
      stats = moment(n2o_norm_alt_grid[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_norm_alt_grid[gd,z,1],sdev=sdev,/nan)
      n2o_norm_alt_grid_yrs[t,z,0:1] = [stats[0],stats2[0]]
      stats = moment(n2o_norm_alt_grid_fl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_norm_alt_grid_fl[gd,z,1],sdev=sdev,/nan)
      n2o_norm_alt_grid_fl_yrs[t,z,0:1] = [stats[0],stats2[0]]
      stats = moment(n2o_norm_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_norm_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      n2o_norm_alt_grid_nofl_yrs[t,z,0:1] = [stats[0],stats2[0]]
      stats = moment(n2o_ch4_norm_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_ch4_norm_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      n2o_ch4_norm_alt_grid_nofl_yrs[t,z,0:1] = [stats[0],stats2[0]]
      stats = moment(ch4_norm_alt_grid_fl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(ch4_norm_alt_grid_fl[gd,z,1],sdev=sdev,/nan)
      ch4_norm_alt_grid_fl_yrs[t,z,0:1] = [stats[0],stats2[0]]
      stats = moment(ch4_norm_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(ch4_norm_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      ch4_norm_alt_grid_nofl_yrs[t,z,0:1] = [stats[0],stats2[0]]
    endfor
    for z = 0, nz3-1 do begin
      gd2 = where(finite(age_combo_alt_grid_nofl3[gd,z,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_alt_grid_nofl3[gd[gd2],z,1]^2))
      wmean = total(1. / age_combo_alt_grid_nofl3[gd[gd2],z,1]^2 * age_combo_alt_grid_nofl3[gd[gd2],z,0]) / total(1. / age_combo_alt_grid_nofl3[gd[gd2],z,1]^2)
      npts = total(age_combo_alt_grid_nofl3[gd[gd2],z,4])
      stats = moment(age_combo_alt_grid_nofl3[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid_nofl3[gd,z,1],sdev=sdev,/nan)
      age_combo_alt_grid_nofl3_yrs[t,z,*] = [wmean,err,stats[0],stats2[0],npts]
      if finite(age_combo_alt_grid_nofl3_yrs[t,z,0]) and ~finite(age_combo_alt_grid_nofl3_yrs[t,z,3]) or age_combo_alt_grid_nofl3_yrs[t,z,3] lt 0.25 then $
        age_combo_alt_grid_nofl3_yrs[t,z,3] = 0.25
    endfor
  endif
endfor

save,samp_type_b_yrs,filename=dir+'Balloon/Sample_type_time_series.sav'


restore,dir+'Balloon/OMS_tracer_profile_avgs.sav'
restore,dir+'Balloon/Mean_age_relationships.sav'
age_grid_hist = age_grid
restore,dir+'Balloon/Balloon_Tracers_vs_mean_age.sav'
restore,dir+'Balloon/Balloon_mean_ages_co2_late_v2025_bins.sav'
bins_co2_b2_v2025 = bins_co2
restore,dir+'Balloon/Balloon_mean_ages_sf6_late_bins.sav'
bins_sf6_b2 = bins_sf6
restore,dir+'Balloon/Balloon_mean_ages_optimum_sf6_co2_early_bins.sav'
bins_b0 = bins
restore,dir+'Balloon/Balloon_mean_ages_sf6_early.sav'
dat_sf6_b0 = dat_sf6
restore,dir+'Balloon/Balloon_mean_ages_sf6_early_bins.sav'
bins_sf6_b0 = bins_sf6
restore,dir+'Balloon/Balloon_mean_ages_co2_early_bins.sav'
bins_co2_b0 = bins_co2
restore,dir+'Aircraft/Missions/idlsave_files//Mean_age_relationships_aircraft.sav'
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_optimum_sf6_co2_early_bins.sav'
bins_a0 = bins
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_co2_early_bins.sav'
bins_co2_a0 = bins_co2
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_early_bins.sav'
bins_sf6_a0 = bins_sf6
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


restore,dir+'Models/WACCM/FWSD_early_aircraft_balloon_subsample_mean_age_N2O.sav'


;  Rebin the aircraft data into coarser balloon bins.
nz = 60
dz = 0.5
alt_grid_a = findgen(nz)*dz+6.
nt_n2o = n_elements(n2o_years)
jj = replicate(!values.f_nan,nn,5,nt_n2o) & bins_co2_a0_n2o_norm_age_all_nh_years = jj & bins_sf6_a2_n2o_norm_age_all_nh_years = jj & bins_co2_a2_n2o_norm_age_all_nh_years = jj
bins_a0_n2o_norm_age_all_nh_years = jj & bins_a2_n2o_norm_age_all_nh_years = jj & bins_a_combo_n2o_norm_age_all_nh_years = jj & bins_a0_ch4_norm_age_all_nh_years = jj
bins_co2_a0_ch4_norm_age_all_nh_years = jj & bins_a2_ch4_norm_age_all_nh_years = jj & bins_sf6_a2_ch4_norm_age_all_nh_years = jj & bins_co2_a2_ch4_norm_age_all_nh_years = jj
uu = replicate(!values.f_nan,nz2,5,nt_n2o) & bins_a_combo_ch4_norm_age_all_nh_years = jj & bins_a0_co2_age_alt_nh_years = uu & bins_a2_co2_age_alt_nh_years = uu
bins_a_combo_age_alt_nh_years = uu & bins_a0_co2_n2o_alt_nh_years = uu & bins_a2_co2_n2o_alt_nh_years = uu & bins_a_combo_n2o_alt_nh_years = uu & bins_a2_sf6_age_alt_nh_years = uu
bins_a2_sf6_n2o_alt_nh_years = uu & bins_a2_sf6_ch4_alt_nh_years = uu & bins_a0_co2_ch4_alt_nh_years = uu & bins_a2_co2_ch4_alt_nh_years = uu & bins_a_combo_ch4_alt_nh_years = uu
ww = replicate(!values.f_nan,nt_n2o) & xx = replicate(!values.f_nan,nz3,5,nt_n2o) & max_alts_a = ww & min_n_norm_a = ww & min_c_norm_a = ww & bins_a0_co2_age_alt3_nh_years = xx
bins_a2_co2_age_alt3_nh_years = xx & bins_a2_sf6_age_alt3_nh_years = xx & bins_a_combo_age_alt3_nh_years = xx & bins_co2_a1_ch4_norm_age_all_nh_years = jj
bins_sf6_a1_ch4_norm_age_all_nh_years = jj & bins_co2_a1_n2o_norm_age_all_nh_years = jj & bins_sf6_a1_n2o_norm_age_all_nh_years = jj & bins_a_sf6_n2o_norm_age_all_nh_years = jj
for t = 0, nt_n2o-1 do begin

  for i = 0, nn-1 do begin

    ii = where(n2o_norm_grid ge norm_grid[i]-dn/2. and n2o_norm_grid le norm_grid[i]+dn/2.)

    npts = total(bins_co2_a0.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_a0.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_a0_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],sdev,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a0.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_co2_a0.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_co2_a0.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_co2_a0_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a1.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_co2_a1.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_co2_a1.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_co2_a1_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_sf6_a1.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_sf6_a1.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_sf6_a1.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_sf6_a1_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    npts = total(bins_sf6_a2.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_a2.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_a2_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],sdev,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_sf6_a2.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_sf6_a2.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_sf6_a2.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_sf6_a2_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a2.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_co2_a2.ch4_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_co2_a2.ch4_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_co2_a2_ch4_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    if finite(bins_co2_a0_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_co2_a0_ch4_norm_age_all_nh_years[i,*,t]
    if finite(bins_a0_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_a0_ch4_norm_age_all_nh_years[i,*,t]
    if finite(bins_co2_a1_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_co2_a1_ch4_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a1_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_sf6_a1_ch4_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a1_ch4_norm_age_all_nh_years[i,0,t]) and finite(bins_co2_a1_ch4_norm_age_all_nh_years[i,0,t]) then begin
      for j = 0, 3 do bins_a_combo_ch4_norm_age_all_nh_years[i,j,t] = mean([bins_sf6_a1_ch4_norm_age_all_nh_years[i,j,t],bins_co2_a1_ch4_norm_age_all_nh_years[i,j,t]])
      bins_a_combo_ch4_norm_age_all_nh_years[i,4,t] = total([bins_sf6_a1_ch4_norm_age_all_nh_years[i,4,t],bins_co2_a1_ch4_norm_age_all_nh_years[i,4,t]])
    endif
    if finite(bins_co2_a2_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_co2_a2_ch4_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a2_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_sf6_a2_ch4_norm_age_all_nh_years[i,*,t]
    if finite(bins_a2_ch4_norm_age_all_nh_years[i,0,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,*,t] = bins_a2_ch4_norm_age_all_nh_years[i,*,t]
    if ~finite(bins_a_combo_ch4_norm_age_all_nh_years[i,2,t]) then bins_a_combo_ch4_norm_age_all_nh_years[i,2,t] = 0.1

    ii = where(n2o_norm_grid ge norm_grid[i]-dn/2. and n2o_norm_grid le norm_grid[i]+dn/2.)

    npts = total(bins_co2_a0.n2o_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_a0.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_a0_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],sdev,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a0.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_co2_a0.n2o_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_co2_a0.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_co2_a0_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    if finite(stats[0]) and ~finite(sdev) then bins_co2_a0_n2o_norm_age_all_nh_years[i,3,t] = 0.25
    err = 1. / sqrt(total(1./bins_co2_a1.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_co2_a1.n2o_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_co2_a1.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_co2_a1_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    ;  Leave out GLOPAC mission for now.
    if years[t] lt 2010 or years[t] gt 2011 then begin
      err = 1. / sqrt(total(1./bins_sf6_a1.n2o_norm_age_uncert_nh_years[ii,t]^2))
      npts = total(bins_sf6_a1.n2o_norm_age_npts_nh_years[ii,t])
      stats = moment(bins_sf6_a1.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
      bins_sf6_a1_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    endif
    npts = total(bins_sf6_a2.n2o_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_a2.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_a2_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],sdev,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_sf6_a2.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_sf6_a2.n2o_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_sf6_a2.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_sf6_a2_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a2.n2o_norm_age_uncert_nh_years[ii,t]^2))
    npts = total(bins_co2_a2.n2o_norm_age_npts_nh_years[ii,t])
    stats = moment(bins_co2_a2.n2o_norm_age_all_nh_years[ii,0,t],sdev=sdev,/nan)
    bins_co2_a2_n2o_norm_age_all_nh_years[i,*,t] = [stats[0],err,stats[0],sdev,npts]
    if finite(bins_a0_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_a0_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_co2_a0_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_co2_a0_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_co2_a1_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_co2_a1_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a1_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_sf6_a1_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a1_n2o_norm_age_all_nh_years[i,0,t]) and finite(bins_co2_a1_n2o_norm_age_all_nh_years[i,0,t]) then begin
      for j = 0, 3 do bins_a_combo_n2o_norm_age_all_nh_years[i,j,t] = mean([bins_sf6_a1_n2o_norm_age_all_nh_years[i,j,t],bins_co2_a1_n2o_norm_age_all_nh_years[i,j,t]])
      bins_a_combo_n2o_norm_age_all_nh_years[i,4,t] = total([bins_sf6_a1_n2o_norm_age_all_nh_years[i,4,t],bins_co2_a1_n2o_norm_age_all_nh_years[i,4,t]])
    endif
    if finite(bins_co2_a2_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_co2_a2_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a2_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_sf6_a2_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_a2_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_combo_n2o_norm_age_all_nh_years[i,*,t] = bins_a2_n2o_norm_age_all_nh_years[i,*,t]

    if finite(bins_a0_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_sf6_n2o_norm_age_all_nh_years[i,*,t] = bins_a0_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a1_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_sf6_n2o_norm_age_all_nh_years[i,*,t] = bins_sf6_a1_n2o_norm_age_all_nh_years[i,*,t]
    if finite(bins_sf6_a2_n2o_norm_age_all_nh_years[i,0,t]) then bins_a_sf6_n2o_norm_age_all_nh_years[i,*,t] = bins_sf6_a2_n2o_norm_age_all_nh_years[i,*,t]
  endfor

  chk = where(finite(bins_co2_a0.n2o_norm_age_all_nh_years[*,0,t]),nchk)
  if nchk gt 0 then min_n_norm_a[t] = n2o_norm_grid[chk[0]]
  chk = where(finite(bins_co2_a2.n2o_norm_age_all_nh_years[*,0,t]),nchk)
  if nchk gt 0 then min_n_norm_a[t] = n2o_norm_grid[chk[0]]
  chk = where(finite(bins_sf6_a2.n2o_norm_age_all_nh_years[*,0,t]),nchk)
  if nchk gt 0 then min_n_norm_a[t] = n2o_norm_grid[chk[0]]
  chk = where(finite(bins_co2_a0.ch4_norm_age_all_nh_years[*,0,t]),nchk)
  if nchk gt 0 then min_c_norm_a[t] = n2o_norm_grid[chk[0]]
  chk = where(finite(bins_co2_a2.ch4_norm_age_all_nh_years[*,0,t]),nchk)
  if nchk gt 0 then min_c_norm_a[t] = n2o_norm_grid[chk[0]]
  chk = where(finite(bins_sf6_a2.ch4_norm_age_all_nh_years[*,0,t]),nchk)
  if nchk gt 0 then min_c_norm_a[t] = n2o_norm_grid[chk[0]]

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
    if finite(bins_a0_co2_age_alt_nh_years[z,0,t]) and n2o_years[t] lt 2000 then bins_a_combo_age_alt_nh_years[z,*,t] = bins_a0_co2_age_alt_nh_years[z,*,t]
    if finite(bins_a2_co2_age_alt_nh_years[z,0,t]) then bins_a_combo_age_alt_nh_years[z,*,t] = bins_a2_co2_age_alt_nh_years[z,*,t]
    if finite(bins_a2_sf6_age_alt_nh_years[z,0,t]) then bins_a_combo_age_alt_nh_years[z,*,t] = bins_a2_sf6_age_alt_nh_years[z,*,t]
    stats = moment(bins_co2_a0.n2o_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a0_co2_n2o_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    stats = moment(bins_co2_a2.n2o_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_co2_n2o_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    stats = moment(bins_sf6_a2.n2o_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_sf6_n2o_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    if finite(bins_a0_co2_n2o_alt_nh_years[z,0,t]) and n2o_years[t] lt 2000 then bins_a_combo_n2o_alt_nh_years[z,*,t] = bins_a0_co2_n2o_alt_nh_years[z,*,t]
    if finite(bins_a2_co2_n2o_alt_nh_years[z,0,t]) then bins_a_combo_n2o_alt_nh_years[z,*,t] = bins_a2_co2_n2o_alt_nh_years[z,*,t]
    if finite(bins_a2_sf6_n2o_alt_nh_years[z,0,t]) then bins_a_combo_n2o_alt_nh_years[z,*,t] = bins_a2_sf6_n2o_alt_nh_years[z,*,t]
    stats = moment(bins_co2_a0.ch4_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a0_co2_ch4_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    stats = moment(bins_co2_a2.ch4_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_co2_ch4_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    stats = moment(bins_sf6_a2.ch4_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_sf6_ch4_alt_nh_years[z,0:1,t] = [stats[0],sdev]
    if finite(bins_a0_co2_ch4_alt_nh_years[z,0,t]) and n2o_years[t] lt 2000 then bins_a_combo_ch4_alt_nh_years[z,*,t] = bins_a0_co2_ch4_alt_nh_years[z,*,t]
    if finite(bins_a2_co2_ch4_alt_nh_years[z,0,t]) then bins_a_combo_ch4_alt_nh_years[z,*,t] = bins_a2_co2_ch4_alt_nh_years[z,*,t]
    if finite(bins_a2_sf6_ch4_alt_nh_years[z,0,t]) then bins_a_combo_ch4_alt_nh_years[z,*,t] = bins_a2_sf6_ch4_alt_nh_years[z,*,t]
  endfor

  for z = 0, nz3-1 do begin
    zi = where(alt_grid_a ge alt_grid3[z]-dz3/2. and alt_grid_a le alt_grid3[z]+dz3/2.)
    err = 1. / sqrt(total(1./bins_co2_a0.age_uncert_alt_nh_years[zi,t]^2))
    npts = total(bins_co2_a0.age_npts_alt_nh_years[zi,t])
    stats = moment(bins_co2_a0.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a0_co2_age_alt3_nh_years[z,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_co2_a2.age_uncert_alt_nh_years[zi,t]^2))
    npts = total(bins_co2_a2.age_npts_alt_nh_years[zi,t])
    stats = moment(bins_co2_a2.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_co2_age_alt3_nh_years[z,*,t] = [stats[0],err,stats[0],sdev,npts]
    err = 1. / sqrt(total(1./bins_sf6_a2.age_uncert_alt_nh_years[zi,t]^2))
    npts = total(bins_sf6_a2.age_npts_alt_nh_years[zi,t])
    stats = moment(bins_sf6_a2.age_alt_nh_years[zi,0,t],sdev=sdev,/nan)
    bins_a2_sf6_age_alt3_nh_years[z,*,t] = [stats[0],err,stats[0],sdev,npts]
    if finite(bins_a0_co2_age_alt3_nh_years[z,0,t]) and n2o_years[t] lt 2000 then bins_a_combo_age_alt3_nh_years[z,*,t] = bins_a0_co2_age_alt3_nh_years[z,*,t]
    if finite(bins_a2_co2_age_alt3_nh_years[z,0,t]) then bins_a_combo_age_alt3_nh_years[z,*,t] = bins_a2_co2_age_alt3_nh_years[z,*,t]
    if finite(bins_a2_sf6_age_alt3_nh_years[z,0,t]) then bins_a_combo_age_alt3_nh_years[z,*,t] = bins_a2_sf6_age_alt3_nh_years[z,*,t]
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


;  Combined aircraft/balloon time series.

bins_all_combo_n2o_norm_age_all_nh_years = kk & bins_all_combo_ch4_norm_age_all_nh_years = kk & bins_fl_combo_ch4_norm_age_all_nh_years = kk & bins_fl_combo_n2o_norm_age_all_nh_years = kk
bins_all_combo_age_alt_nh_years = mm & bins_all_combo_n2o_norm_alt_nh_years = mm & bins_all_combo_ch4_norm_alt_nh_years = mm & bins_all_combo_fl_age_alt_nh_years = mm
bins_all_combo_fl_n2o_alt_nh_years = mm & bins_all_combo_fl_ch4_alt_nh_years = mm & bins_all_combo_nofl_ch4_norm_age_all_nh_years = kk & bins_all_combo_nofl_n2o_norm_age_all_nh_years = kk
bins_all_combo_age_alt3_nh_years = ee & samp_type_a_yrs = replicate(!values.f_nan,nn,nt) & samp_type_a_yrs_alt = replicate(!values.f_nan,nz2,nt)
bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years = kk & bins_sf6_nofl_n2o_norm_age_all_nh_years = kk & bins_all_combo_nofl_n2o_norm_alt_nh_years = mm
;for t = 0, nt-nt1-1 do begin
for t = 0, nt-1 do begin
  for i = 0, nn-1 do begin
    for j = 0, 3 do begin
      bins_all_combo_nofl_ch4_norm_age_all_nh_years[t,i,j] = mean([age_combo_ch4_bin_nofl_yrs[t,i,j],bins_a_combo_ch4_norm_age_all_nh_years[i,j,t]],/nan)
      if finite(age_combo_ch4_bin_fl_yrs[t,i,j]) then bins_all_combo_ch4_norm_age_all_nh_years[t,i,j] = mean([age_combo_ch4_bin_fl_yrs[t,i,j], $
        bins_a_combo_ch4_norm_age_all_nh_years[i,j,t]],/nan)
      if finite(age_combo_ch4_bin_nofl_yrs[t,i,j]) then bins_all_combo_ch4_norm_age_all_nh_years[t,i,j] = mean([age_combo_ch4_bin_nofl_yrs[t,i,j], $
        bins_a_combo_ch4_norm_age_all_nh_years[i,j,t]],/nan)
      bins_all_combo_nofl_n2o_norm_age_all_nh_years[t,i,j] = mean([age_combo_n2o_bin_nofl_yrs[t,i,j],bins_a_combo_n2o_norm_age_all_nh_years[i,j,t]],/nan)
      if j eq 3 and bins_all_combo_nofl_n2o_norm_age_all_nh_years[t,i,j] lt 0.15 then bins_all_combo_nofl_n2o_norm_age_all_nh_years[t,i,j] = 0.15
      if j eq 3 and ~finite(bins_all_combo_nofl_n2o_norm_age_all_nh_years[t,i,j]) then bins_all_combo_nofl_n2o_norm_age_all_nh_years[t,i,j] = 0.15
      bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years[t,i,j] = mean([age_combo_n2o_bin_nofl_noch4n2o_yrs[t,i,j],bins_a_combo_n2o_norm_age_all_nh_years[i,j,t]],/nan)
      if j eq 1 and ~finite(bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years[t,i,j]) then bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years[t,i,j] = 0.1

      bins_sf6_nofl_n2o_norm_age_all_nh_years[t,i,j] = mean([age_sf6_n2o_bin_nofl_yrs[t,i,j],bins_a_sf6_n2o_norm_age_all_nh_years[i,j,t]],/nan)

      if finite(bins_a_combo_n2o_norm_age_all_nh_years[i,j,t]) then samp_type_a_yrs[i,t] = 1
      if finite(age_combo_n2o_bin_fl_yrs[t,i,j]) then bins_all_combo_n2o_norm_age_all_nh_years[t,i,j] = mean([age_combo_n2o_bin_fl_yrs[t,i,j], $
        bins_a_combo_n2o_norm_age_all_nh_years[i,j,t]],/nan)
      if finite(age_combo_n2o_bin_nofl_yrs[t,i,j]) then bins_all_combo_n2o_norm_age_all_nh_years[t,i,j] = mean([age_combo_n2o_bin_nofl_yrs[t,i,j], $
        bins_a_combo_n2o_norm_age_all_nh_years[i,j,t]],/nan)
      if j eq 3 and bins_all_combo_n2o_norm_age_all_nh_years[t,i,j] lt 0.15 then bins_all_combo_n2o_norm_age_all_nh_years[t,i,j] = 0.15
      if j eq 3 and ~finite(bins_all_combo_n2o_norm_age_all_nh_years[t,i,j]) then bins_all_combo_n2o_norm_age_all_nh_years[t,i,j] = 0.15
    endfor
    for j = 0, 3 do bins_fl_combo_ch4_norm_age_all_nh_years[t,i,j] = age_q_combo_ch4_bin_fl_yrs[t,i,j]
    for j = 0, 3 do bins_fl_combo_n2o_norm_age_all_nh_years[t,i,j] = age_q_combo_n2o_bin_fl_yrs[t,i,j]
  endfor
  for z = 0, nz2-1 do begin
    for j = 0, 3 do begin
      bins_all_combo_age_alt_nh_years[t,z,j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      if j eq 0 and finite(bins_a_combo_age_alt_nh_years[z,j,t]) then samp_type_a_yrs_alt[z,t] = 1
      bins_all_combo_fl_age_alt_nh_years[t,z,j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      if finite(age_combo_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z,j] = mean([age_combo_alt_grid_fl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      if finite(age_combo_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z,j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      bins_all_combo_n2o_norm_alt_nh_years[t,z,j] = mean([n2o_norm_alt_grid_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
      bins_all_combo_nofl_n2o_norm_alt_nh_years[t,z,j] = mean([n2o_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
      bins_all_combo_fl_n2o_alt_nh_years[t,z,j] = mean([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
      if finite(n2o_norm_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_n2o_alt_nh_years[t,z,j] = mean([n2o_norm_alt_grid_fl_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
      if finite(n2o_norm_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_n2o_alt_nh_years[t,z,j] = mean([n2o_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
      bins_all_combo_ch4_norm_alt_nh_years[t,z,j] = mean([ch4_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
      bins_all_combo_fl_ch4_alt_nh_years[t,z,j] = mean([ch4_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
      if finite(ch4_norm_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_ch4_alt_nh_years[t,z,j] = mean([ch4_norm_alt_grid_fl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
      if finite(ch4_norm_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_ch4_alt_nh_years[t,z,j] = mean([ch4_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
    endfor
    j = 4
    bins_all_combo_age_alt_nh_years[t,z,j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_fl_age_alt_nh_years[t,z,j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    if finite(age_combo_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z,j] = total([age_combo_alt_grid_fl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    if finite(age_combo_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_age_alt_nh_years[t,z,j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_n2o_norm_alt_nh_years[t,z,j] = total([n2o_norm_alt_grid_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_nofl_n2o_norm_alt_nh_years[t,z,j] = total([n2o_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_fl_n2o_alt_nh_years[t,z,j] = total([age_combo_alt_grid_nofl_yrs[t,z,j],bins_a_combo_age_alt_nh_years[z,j,t]],/nan)
    if finite(n2o_norm_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_n2o_alt_nh_years[t,z,j] = total([n2o_norm_alt_grid_fl_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
    if finite(n2o_norm_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_n2o_alt_nh_years[t,z,j] = total([n2o_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_n2o_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_ch4_norm_alt_nh_years[t,z,j] = total([ch4_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
    bins_all_combo_fl_ch4_alt_nh_years[t,z,j] = total([ch4_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
    if finite(ch4_norm_alt_grid_fl_yrs[t,z,j]) then bins_all_combo_fl_ch4_alt_nh_years[t,z,j] = total([ch4_norm_alt_grid_fl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
    if finite(ch4_norm_alt_grid_nofl_yrs[t,z,j]) then bins_all_combo_fl_ch4_alt_nh_years[t,z,j] = total([ch4_norm_alt_grid_nofl_yrs[t,z,j],bins_a_combo_ch4_alt_nh_years[z,j,t]],/nan)
  endfor

  for z = 0, nz3-1 do for j = 0, 1 do bins_all_combo_age_alt3_nh_years[t,z,j] = mean([age_combo_alt_grid_nofl3_yrs[t,z,j],bins_a_combo_age_alt3_nh_years[z,j,t]],/nan)
endfor


restore,dir+'Balloon/Mean_age_and_N2O_trend_alt_profile_elat_adj.sav'
restore,dir+'Balloon/Age_N2O_gridded_N2O_elat_adj_time_series.sav'


; ------------------------  Trends  ---------------------------------------------------------------------------------------

zi = where(alt_grid2 ge 24)
bins_all_combo_fl_age_midstrat_nh_years = replicate(!values.f_nan,nt,5) & bins_all_combo_fl_age_maxalt_nh_years = rr
bins_all_combo_nofl_age_midstrat_nh_years = bins_all_combo_fl_age_midstrat_nh_years
for t = 0, nt-1 do begin
  gd = where(finite(bins_all_combo_fl_age_alt_nh_years[t,zi,0]),ngd)
  if ngd gt 0 then begin
    err_fl = 1. / sqrt(total(1./bins_all_combo_fl_age_alt_nh_years[t,zi[gd],1]^2))
    wmean = total(1. / bins_all_combo_fl_age_alt_nh_years[t,zi[gd],1]^2 * bins_all_combo_fl_age_alt_nh_years[t,zi[gd],0]) / $
      total(1. / bins_all_combo_fl_age_alt_nh_years[t,zi[gd],1]^2)
    bins_all_combo_fl_age_midstrat_nh_years[t,0:1] = [wmean,err_fl]
  endif
  gd = where(finite(bins_all_combo_age_alt_nh_years[t,zi,1]),ngd)
  if ngd gt 0 then begin
    err_nofl = 1. / sqrt(total(1./bins_all_combo_age_alt_nh_years[t,zi[gd],1]^2))
    wmean = total(1. / bins_all_combo_age_alt_nh_years[t,zi[gd],1]^2 * bins_all_combo_age_alt_nh_years[t,zi[gd],0]) / $
      total(1. / bins_all_combo_age_alt_nh_years[t,zi[gd],1]^2)
    bins_all_combo_nofl_age_midstrat_nh_years[t,0:1] = [wmean,err_nofl]
  endif
  chk = where(finite(bins_all_combo_fl_age_alt_nh_years[t,*,0]),nchk)
  if nchk gt 0 then bins_all_combo_fl_age_maxalt_nh_years[t] = alt_grid2[chk[-1]]
  for j = 2, 3 do begin
    bins_all_combo_fl_age_midstrat_nh_years[t,j] = mean(bins_all_combo_fl_age_alt_nh_years[t,zi,j],/nan)
    bins_all_combo_nofl_age_midstrat_nh_years[t,j] = mean(bins_all_combo_age_alt_nh_years[t,zi,j],/nan)
  endfor
  bins_all_combo_fl_age_midstrat_nh_years[t,4] = total(bins_all_combo_fl_age_alt_nh_years[t,zi,4],/nan)
  bins_all_combo_nofl_age_midstrat_nh_years[t,4] = total(bins_all_combo_age_alt_nh_years[t,zi,4],/nan)
endfor
gd = where(finite(bins_all_combo_fl_age_midstrat_nh_years[*,0]))
age_fl_midstrat_linfit = linfit(years[gd],bins_all_combo_fl_age_midstrat_nh_years[gd,0],sigma=sigma)
age_fl_midstrat_trends = [age_fl_midstrat_linfit[1],sigma[1]]
stats = moment(bins_all_combo_fl_age_midstrat_nh_years[gd,0],sdev=sdev)
age_fl_midstrat_trends_per = 1e3*age_fl_midstrat_trends/stats[0]
gd = where(finite(bins_all_combo_nofl_age_midstrat_nh_years[*,0]))
age_nofl_midstrat_linfit = linfit(years[gd],bins_all_combo_nofl_age_midstrat_nh_years[gd,0],sigma=sigma)
age_nofl_midstrat_trends = [age_nofl_midstrat_linfit[1],sigma[1]]
stats = moment(bins_all_combo_nofl_age_midstrat_nh_years[gd,0],sdev=sdev)
age_nofl_midstrat_trends_per = 1e3*age_nofl_midstrat_trends/stats[0]

vv = replicate(!values.f_nan,nn,2) & n2o_age_trends = vv & n2o_age_trends_sigma = vv & ch4_age_trends = vv & ch4_age_trends_sigma = vv & ch4_fl_age_trends = vv & ch4_fl_age_trends_sigma = vv
ww = replicate(!values.f_nan,nz2,2) & n2o_fl_age_trends = vv & n2o_fl_age_trends_sigma = vv & age_alt_trends = ww & age_alt_trends_sigma = ww & n2o_alt_trends = ww & n2o_alt_trends_sigma = ww
ch4_alt_trends = ww & ch4_alt_trends_sigma = ww & age_fl_alt_trends = ww & age_fl_alt_trends_sigma = ww & n2o_fl_alt_trends = ww & n2o_fl_alt_trends_sigma = ww & age_alt_means = ww
ch4_fl_alt_trends = ww & ch4_fl_alt_trends_sigma = ww & age_fl_alt_means = ww & n2o_nofl_age_trends = vv & n2o_nofl_age_trends_sigma = vv & ch4_nofl_age_trends = vv & ch4_nofl_age_trends_sigma = vv
age_ac_alt_trends = ww & age_ac_alt_trends_sigma = ww & age_ac_alt_means = ww & age_n2o_means = vv & aa = replicate(!values.f_nan,nz3,2) & age_alt3_trends = aa & age_alt3_means = aa
age_alt3_trends_sigma = aa & n2o_nofl_age_trends_count = fltarr(nn) & age_nofl_alt_trends_count = fltarr(nz2) & n2o_nofl_alt_trends_count = fltarr(nz2) & n2o_wfl_alt_trends_count = fltarr(nz2)
ch4_nofl_alt_trends_count = fltarr(nz2) & ch4_wfl_alt_trends_count = fltarr(nz2) & ch4_nofl_age_trends_count = fltarr(nz2) & age_fl_alt_trends_count = fltarr(nz2)
n2o_nofl_noch4n2o_age_trends_count = fltarr(nz2) & n2o_nofl_noch4n2o_age_trends = vv & n2o_nofl_noch4n2o_age_trends_sigma = vv & n2o_nofl_sf6_age_trends = vv & n2o_nofl_sf6_age_trends_sigma = vv
age_sf6_n2o_means = vv & age_fl_alt_trends_recent = ww & age_fl_alt_trends_sigma_recent = ww & age_fl_alt_means_recent = ww & age_on_ch4_trends = vv
age_on_ch4_trends_sigma = vv & age_ch4_means = vv & age_on_ch4_trends_count = fltarr(nn) & n2o_nofl_age_trends_ac = vv & n2o_nofl_age_trends_ac_sigma = vv
for i = 0, nn-1 do begin
  gd = where(finite(bins_all_combo_n2o_norm_age_all_nh_years[*,i,0]),ngd)
  if ngd ge 15 then begin
    n2o_age_trends[i,*] = linfit(years[gd],bins_all_combo_n2o_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    n2o_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_n2o_adj_tseries[2,i,*]),ngd)
  n2o_nofl_age_trends_count[i] = ngd
  if ngd ge 20 then begin
    n2o_nofl_age_trends[i,*] = linfit(years[gd],age_on_n2o_adj_tseries[2,i,gd],sigma=sigma)
    n2o_nofl_age_trends_sigma[i,*] = sigma
    stats = moment(age_on_n2o_adj_tseries[2,i,gd],sdev=sdev)
    age_n2o_means[i,*] = [stats[0],sdev]
  endif
  gd = where(finite(age_on_n2o_adj_tseries[2,i,*]) and age_on_n2o_adj_tseries[5,i,*] eq 3,ngd)
  if ngd ge 10 then begin
    n2o_nofl_age_trends_ac[i,*] = linfit(years[gd],age_on_n2o_adj_tseries[2,i,gd],sigma=sigma)
    n2o_nofl_age_trends_ac_sigma[i,*] = sigma
  endif
  gd = where(finite(bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years[*,i,0]),ngd)
  n2o_nofl_noch4n2o_age_trends_count[i] = ngd
  if ngd ge 25 then begin
    n2o_nofl_noch4n2o_age_trends[i,*] = linfit(years[gd],bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    n2o_nofl_noch4n2o_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(bins_sf6_nofl_n2o_norm_age_all_nh_years[*,i,0]),ngd)
  if ngd ge 10 then begin
    n2o_nofl_sf6_age_trends[i,*] = linfit(years[gd],bins_sf6_nofl_n2o_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    n2o_nofl_sf6_age_trends_sigma[i,*] = sigma
    stats = moment(bins_sf6_nofl_n2o_norm_age_all_nh_years[gd,i,0],sdev=sdev)
    age_sf6_n2o_means[i,*] = [stats[0],sdev]
  endif
  gd = where(finite(bins_fl_combo_n2o_norm_age_all_nh_years[*,i,0]),ngd)
  if ngd ge 15 then begin
    n2o_fl_age_trends[i,*] = linfit(years[gd],bins_fl_combo_n2o_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    n2o_fl_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(bins_all_combo_ch4_norm_age_all_nh_years[*,i,0]),ngd)
  if ngd ge 15 then begin
    ch4_age_trends[i,*] = linfit(years[gd],bins_all_combo_ch4_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    ch4_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(bins_all_combo_nofl_ch4_norm_age_all_nh_years[*,i,0]),ngd)
  ch4_nofl_age_trends_count[i] = ngd
  if ngd ge 15 then begin
    ch4_nofl_age_trends[i,*] = linfit(years[gd],bins_all_combo_nofl_ch4_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    ch4_nofl_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(bins_fl_combo_ch4_norm_age_all_nh_years[*,i,0]),ngd)
  if ngd ge 15 then begin
    ch4_fl_age_trends[i,*] = linfit(years[gd],bins_fl_combo_ch4_norm_age_all_nh_years[gd,i,0],sigma=sigma)
    ch4_fl_age_trends_sigma[i,*] = sigma
  endif
  gd = where(finite(age_on_ch4_adj_tseries[2,i,*]),ngd)
  age_on_ch4_trends_count[i] = ngd
  if ngd ge 20 then begin
    age_on_ch4_trends[i,*] = linfit(years[gd],age_on_ch4_adj_tseries[2,i,gd],sigma=sigma)
    age_on_ch4_trends_sigma[i,*] = sigma
    stats = moment(age_on_ch4_adj_tseries[2,i,gd],sdev=sdev)
    age_ch4_means[i,*] = [stats[0],sdev]
  endif
endfor

age_nh_coarse_wflask = age_nh_coarse

for z = 0, nz2-1 do begin
  gd = where(finite(bins_all_combo_age_alt_nh_years[*,z,0]),ngd)
  age_nofl_alt_trends_count[z] = ngd
  if ngd ge 15 then begin
    chk = where(~finite(bins_all_combo_age_alt_nh_years[gd,z,1]),nchk)
    if nchk gt 0 then bins_all_combo_age_alt_nh_years[gd[chk],z,1] = 0.25
    age_alt_trends[z,*] = linfit(years[gd],bins_all_combo_age_alt_nh_years[gd,z,0],sigma=sigma)
    age_alt_trends_sigma[z,*] = sigma
    stats = moment(bins_all_combo_age_alt_nh_years[gd,z,0],sdev=sdev)
    age_alt_means[z,*] = [stats[0],sdev]
  endif
  gd = where(years gt 2010 and finite(bins_all_combo_age_alt_nh_years[*,z,0]),ngd)
  if ngd ge 8 then begin
    age_ac_alt_trends[z,*] = linfit(years[gd],bins_all_combo_age_alt_nh_years[gd,z,0],sigma=sigma)
    age_ac_alt_trends_sigma[z,*] = sigma
    stats = moment(bins_all_combo_age_alt_nh_years[gd,z,0],sdev=sdev)
    age_ac_alt_means[z,*] = [stats[0],sdev]
  endif

  ;  Trends including flask data.  Add in the elat adjusted non-flask data to the time series.
  si = where(samp_type_b_yrs eq 0,nsi)
  for j = 0, 1 do age_nh_coarse_wflask[j,z,si] = bins_all_combo_fl_age_alt_nh_years[si,z,j]
  for j = 2, 3 do age_nh_coarse_wflask[j,z,si] = bins_all_combo_fl_age_alt_nh_years[si,z,j-2]
  age_nh_coarse_wflask[4,z,si] = bins_all_combo_fl_age_alt_nh_years[si,z,4]
  chk = where(samp_type_b_yrs eq 0 and age_nh_coarse[5,z,*] gt 0,nchk)
  if nchk eq 0 then age_nh_coarse_wflask[5,z,si] = 0
  gd = where(finite(age_nh_coarse_wflask[2,z,*]),ngd)
  age_fl_alt_trends_count[z] = ngd
  if ngd ge 15 then begin
    age_fl_alt_trends[z,*] = linfit(years[gd],age_nh_coarse_wflask[2,z,gd],sigma=sigma)
    age_fl_alt_trends_sigma[z,*] = sigma
    stats = moment(age_nh_coarse_wflask[2,z,gd],sdev=sdev)
    age_fl_alt_means[z,*] = [stats[0],sdev]
  endif
  gd = where(finite(age_nh_coarse_wflask[2,z,*]) and years ge 1993,ngd)
  if ngd ge 15 then begin
    age_fl_alt_trends_recent[z,*] = linfit(years[gd],age_nh_coarse_wflask[2,z,gd],sigma=sigma)
    age_fl_alt_trends_sigma_recent[z,*] = sigma
    stats = moment(age_nh_coarse_wflask[2,z,gd],sdev=sdev)
    age_fl_alt_means_recent[z,*] = [stats[0],sdev]
  endif
  gd = where(finite(bins_all_combo_n2o_norm_alt_nh_years[*,z,0]) and years ge 1990,ngd)
  n2o_nofl_alt_trends_count[z] = ngd
  if ngd ge 15 then begin
    n2o_alt_trends[z,*] = linfit(years[gd],bins_all_combo_n2o_norm_alt_nh_years[gd,z,0],sigma=sigma)
    n2o_alt_trends_sigma[z,*] = sigma
  endif
  gd = where(finite(bins_all_combo_fl_n2o_alt_nh_years[*,z,0]),ngd)
  n2o_wfl_alt_trends_count[z] = ngd
  if ngd ge 15 then begin
    n2o_fl_alt_trends[z,*] = linfit(years[gd],bins_all_combo_fl_n2o_alt_nh_years[gd,z,0],sigma=sigma)
    n2o_fl_alt_trends_sigma[z,*] = sigma
  endif
  gd = where(finite(bins_all_combo_ch4_norm_alt_nh_years[*,z,0]),ngd)
  ch4_nofl_alt_trends_count[z] = ngd
  if ngd ge 10 then begin
    ch4_alt_trends[z,*] = linfit(years[gd],bins_all_combo_ch4_norm_alt_nh_years[gd,z,0],sigma=sigma)
    ch4_alt_trends_sigma[z,*] = sigma
  endif
  gd = where(finite(bins_all_combo_fl_ch4_alt_nh_years[*,z,0]),ngd)
  ch4_wfl_alt_trends_count[z] = ngd
  if ngd ge 10 then begin
    ch4_fl_alt_trends[z,*] = linfit(years[gd],bins_all_combo_fl_ch4_alt_nh_years[gd,z,0],sigma=sigma)
    ch4_fl_alt_trends_sigma[z,*] = sigma
  endif
endfor
for z = 0, nz3-1 do begin
  gd = where(finite(bins_all_combo_age_alt3_nh_years[*,z,0]),ngd)
  if ngd ge 8 then begin
    chk = where(~finite(bins_all_combo_age_alt3_nh_years[gd,z,1]),nchk)
    if nchk gt 0 then bins_all_combo_age_alt3_nh_years[gd[chk],z,1] = 0.25
    age_alt3_trends[z,*] = linfit(years[gd],bins_all_combo_age_alt3_nh_years[gd,z,0],sigma=sigma)
    age_alt3_trends_sigma[z,*] = sigma
    stats = moment(bins_all_combo_age_alt3_nh_years[gd,z,0],sdev=sdev)
    age_alt3_means[z,*] = [stats[0],sdev]
  endif
endfor


;  Find seasonal cycles.
seas_frac = findgen(12)/12.+1./24.
mean_age_ac_seas = replicate(!values.f_nan,12,nn) & bins_all_combo_n2o_norm_age_all_nh_years_noseas = bins_all_combo_n2o_norm_age_all_nh_years
bins_all_combo_nofl_n2o_norm_age_all_nh_years_noseas = bins_all_combo_nofl_n2o_norm_age_all_nh_years
bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years_noseas = bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years
si = where(samp_type_b_yrs eq 2,nsi)
si1 = where(samp_type_b_yrs eq 1,nsi)
seas = years[si]-fix(years[si])
seas1 = years[si1]-fix(years[si1])
ti = interpol(findgen(12),seas_frac,seas)
ti1 = interpol(findgen(12),seas_frac,seas1)
for i = nn-1, 0, -1 do begin
  seas_anom = bins_all_combo_n2o_norm_age_all_nh_years[si,i,0] - (n2o_nofl_age_trends[i,0] + n2o_nofl_age_trends[i,1]*years[si])
  for t = 0, 11 do begin
    chk = where(seas ge seas_frac[t]-1./24. and seas le seas_frac[t]+1./24.,nchk)
    if nchk gt 0 then begin
      mean_age_ac_seas[t,i] = mean(seas_anom[chk],/nan)
    endif
  endfor
  chk = where(finite(mean_age_ac_seas[*,i]),nchk)
  if nchk eq 12 then mean_age_ac_seas[*,i] = lowpass_cfc(mean_age_ac_seas[*,i], BOX=3, EDGE_PFCAST=1)
  if i eq 8 then mean_age_ac_seas[*,i] = 1.8 * mean_age_ac_seas[*,i+2]
  if i eq 7 or i eq 6 then mean_age_ac_seas[*,i] = 1.4 * mean_age_ac_seas[*,i+1]

  tot = total(mean_age_ac_seas[*,i])
  mean_age_ac_seas[*,i] -= tot/12.

  ;  Remove seasonal cycle from N2O time series.
  bins_all_combo_n2o_norm_age_all_nh_years_noseas[si,i,0] -= interpolate(mean_age_ac_seas[*,i],ti)
  bins_all_combo_nofl_n2o_norm_age_all_nh_years_noseas[si,i,0] -= interpolate(mean_age_ac_seas[*,i],ti)
  bins_all_combo_nofl_noch4n2o_n2o_norm_age_all_nh_years_noseas[si,i,0] -= interpolate(mean_age_ac_seas[*,i],ti)
  bins_all_combo_n2o_norm_age_all_nh_years_noseas[si1,i,0] -= interpolate(mean_age_ac_seas[*,i],ti1)

  sia = where(samp_type_a_yrs[i,*] eq 1)
  seas_a = years[sia]-fix(years[sia])
  tia = interpol(findgen(12),seas_frac,seas_a)
  bins_all_combo_nofl_n2o_norm_age_all_nh_years_noseas[sia,i,0] -= interpolate(mean_age_ac_seas[*,i],tia)
endfor

;  Interpolate the seasonal cycle to the fine N2O grid.
ni = interpol(indgen(nn),norm_grid,n2o_norm_grid)
mean_age_ac_seas_hires = replicate(!values.f_nan,12,nnorm)
for t = 0, 11 do begin
  mean_age_ac_seas_hires[t,*] = interpolate(mean_age_ac_seas[t,*],ni)
endfor


;  Trends w/out seasonal cycle.
n2o_nofl_age_noseas_trends = vv & n2o_nofl_age_noseas_trends_sigma = vv
for i = 0, nn-1 do begin
  gd = where(finite(bins_all_combo_nofl_n2o_norm_age_all_nh_years_noseas[*,i,0]),ngd)
  n2o_nofl_age_trends_count[i] = ngd
  if ngd ge 20 then begin
    n2o_nofl_age_noseas_trends[i,*] = linfit(years[gd],bins_all_combo_nofl_n2o_norm_age_all_nh_years_noseas[gd,i,0],sigma=sigma)
    n2o_nofl_age_noseas_trends_sigma[i,*] = sigma
  endif
endfor


;  Read in the full aircraft and balloon data to make histograms.
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sweep_sf6_co2_early.sav'
dat_a0 = dat
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sweep_sf6_co2_late.sav'
dat_a2 = dat
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_early.sav'
dat_co2_a0 = dat_co2
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_co2_late.sav'
dat_co2_a2 = dat_co2
restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_sf6_early.sav'
dat_sf6_a0 = dat_sf6
restore,dir+'Aircraft/Airborne_Save_Files/aircraft_mean_ages_sf6_late.sav'
dat_sf6_a2 = dat_sf6
restore,dir+'Aircraft/Airborne_Save_Files/Aircraft_mean_ages_early_bins_comb.sav'
restore,dir+'Balloon/Balloon_mean_ages_co2_late_v2025.sav'
dat_co2_b2 = dat_co2

elat_a = [dat_co2_a0.elat,dat_co2_a2.elat]
alt_a = [dat_co2_a0.alt,dat_co2_a2.alt]
yrfrac_a = [dat_co2_a0.yrfrac,dat_co2_a2.yrfrac]
age_a = [dat_co2_a0.age_opt_all,dat_co2_a2.age_opt_all]
n2o_norm_a = [dat_co2_a0.n2o_norm,dat_co2_a2.n2o_norm]
ch4_norm_a = [dat_co2_a0.ch4_norm,dat_co2_a2.ch4_norm]

nna = 11
dna = 0.1
norm_grid_a = findgen(nna)*dna+0.05
na = 32
da = 0.25
age_grid = findgen(na)*da
qq = replicate(!values.f_nan,na,nna) & H_age_n2o_a_norm = qq & H_age_ch4_a_norm = qq & H_age_co2_n2o_bfl_norm = qq & H_age_co2_ch4_bfl_norm = qq & H_age_sf6_n2o_bfl_norm = qq
H_age_sf6_ch4_bfl_norm = qq & H_age_n2o_bfl_norm = qq & H_age_ch4_bfl_norm = qq & H_age_co2_n2o_bnofl_norm = qq & H_age_co2_ch4_bnofl_norm = qq & H_age_sf6_n2o_bnofl_norm = qq
H_age_sf6_ch4_bnofl_norm = qq & H_age_n2o_bnofl_norm = qq & H_age_ch4_bnofl_norm = qq
for i = 0, nna-1 do begin
  gd = where(elat_a ge 25 and yrfrac_a le 2010 and n2o_norm_a ge norm_grid_a[i] - dn/2. and n2o_norm_a lt norm_grid_a[i] + dn/2.,ngd)
  if ngd gt 1 then begin
    H = histogram(age_a[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_n2o_a_norm[*,i] = H / total(H)
  endif
  gd = where(elat_a ge 25 and yrfrac_a le 2010 and ch4_norm_a ge norm_grid_a[i] - dn/2. and ch4_norm_a lt norm_grid_a[i] + dn/2.,ngd)
  if ngd gt 1 then begin
    H = histogram(age_a[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_ch4_a_norm[*,i] = H / total(H)
  endif
  gd = where(yrfrac_b_fl le 2010 and n2o_b_fl ge norm_grid_a[i] - dn/2. and n2o_b_fl lt norm_grid_a[i] + dn/2.,ngd)
  if ngd ge 1 then begin
    H = histogram(age_co2_b_fl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_co2_n2o_bfl_norm[*,i] = H / total(H)
    H = histogram(age_sf6_b_fl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_sf6_n2o_bfl_norm[*,i] = H / total(H)
  endif
  H_age_n2o_bfl_norm[*,i] = 0.5 * (H_age_co2_n2o_bfl_norm[*,i] + H_age_sf6_n2o_bfl_norm[*,i])
  gd = where(yrfrac_b_nofl le 2010 and n2o_b_nofl ge norm_grid_a[i] - dn/2. and n2o_b_nofl lt norm_grid_a[i] + dn/2.,ngd)
  if ngd ge 1 then begin
    H = histogram(age_co2_b_nofl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_co2_n2o_bnofl_norm[*,i] = H / total(H)
    H = histogram(age_sf6_b_nofl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_sf6_n2o_bnofl_norm[*,i] = H / total(H)
  endif
  H_age_n2o_bnofl_norm[*,i] = 0.5 * (H_age_co2_n2o_bnofl_norm[*,i] + H_age_sf6_n2o_bnofl_norm[*,i])
  gd = where(yrfrac_b_fl le 2010 and ch4_b_fl ge norm_grid_a[i] - dn/2. and ch4_b_fl lt norm_grid_a[i] + dn/2.,ngd)
  if ngd ge 1 then begin
    H = histogram(age_co2_b_fl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_co2_ch4_bfl_norm[*,i] = H / total(H)
    H = histogram(age_sf6_b_fl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_sf6_ch4_bfl_norm[*,i] = H / total(H)
  endif
  H_age_ch4_bfl_norm[*,i] = 0.5 * (H_age_co2_ch4_bfl_norm[*,i] + H_age_sf6_ch4_bfl_norm[*,i])
  gd = where(yrfrac_b_nofl le 2010 and ch4_b_nofl ge norm_grid_a[i] - dn/2. and ch4_b_nofl lt norm_grid_a[i] + dn/2.,ngd)
  if ngd ge 1 then begin
    H = histogram(age_co2_b_nofl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_co2_ch4_bnofl_norm[*,i] = H / total(H)
    H = histogram(age_sf6_b_nofl[gd],binsize=da,max=age_grid[-1],min=0)
    H_age_sf6_ch4_bnofl_norm[*,i] = H / total(H)
  endif
  H_age_ch4_bnofl_norm[*,i] = 0.5 * (H_age_co2_ch4_bnofl_norm[*,i] + H_age_sf6_ch4_bnofl_norm[*,i])
endfor



;  The seasonal 1990s grids are built by make_balloon_seas_grid.pro (run before this program by run_all.pro).
restore,dir+'Balloon/balloon_mean_age_seas_grid_early.sav'
restore,dir+'Balloon/OMS_grid_elat_adj.sav'
restore,dir+'Balloon/Aircore_grid_elat_adj.sav'
nyg = 36
dy = 5.
nz = 60
dz = 0.5
lat_grid_a = findgen(nyg)*dy-85.
nza = 14
dza = 1.0
alt_grid_a0 = findgen(nz)*dz+6.
alt_grid_a02 = findgen(nza)*dza+12.
yi = where(lat_grid_a ge 30 and lat_grid_a le 60)


restore,dir+'Aircraft/Missions/idlsave_files/aircraft_mean_ages_early.sav'


;  Changes in seasonal averages.
www = replicate(!values.f_nan,nyg,nz,4) & yyy = replicate(!values.f_nan,nyg,nz) & age_diffs_seas = www & n2o_diffs_seas = www & age_diffs = yyy & n2o_diffs = yyy & nseas_n2o = yyy & nseas_age = yyy
for s = 0, 3 do for z = 0, nz-1 do begin
  ii = where(finite(n2o_norm_grid_seas_both[0,*,z,s]) and finite(n2o_norm_combo_grid_seas_late[0,*,z,s]),nii)
  n2o_diffs_seas[ii,z,s] = n2o_norm_combo_grid_seas_late[0,ii,z,s]-n2o_norm_grid_seas_both[0,ii,z,s]
  ii = where(finite(age_combo_grid_seas_sm_early[0,*,z,s]) and finite(age_combo_grid_seas_late[0,*,z,s]),nii)
  age_diffs_seas[ii,z,s] = age_combo_grid_seas_late[0,ii,z,s]-age_combo_grid_seas_sm_early[0,ii,z,s]
endfor
for z = 0, nz-1 do for y = 0, nyg-1 do begin
  chk = where(finite(n2o_diffs_seas[y,z,*]),nchk)
  nseas_n2o[y,z] = nchk
  n2o_diffs[y,z] = mean(n2o_diffs_seas[y,z,*],/nan)
  chk = where(finite(age_diffs_seas[y,z,*]),nchk)
  nseas_age[y,z] = nchk
  age_diffs[y,z] = mean(age_diffs_seas[y,z,*],/nan)
endfor


;  Read in mean age vs N2O from elat adj program.
restore,dir+'Balloon/Mean_age_vs_N2O_1990s.sav'


;  Calculate 2020s mean age vs N2O from 1990s avg and trend.
nic = interpol(indgen(nnorm),n2o_norm_grid,norm_grid)
bins_a0_n2o_norm_age_nh_coarse = interpolate(bins_a0_n2o_norm_age_nh,nic)
n2o_nofl_noch4n2o_age_trends[-1,1] = -8e-3
bins_a2_n2o_norm_age_nh_trend = bins_a0_n2o_norm_age_nh_coarse + 25.*n2o_nofl_age_trends[*,1]

mean_age_on_n2o_90s = mean_age_all_on_n2o_0[*,2]
mean_age_on_n2o_90s[0:-4] = smooth(mean_age_on_n2o_90s[0:-4],5,/edge_truncate)
mean_age_on_n2o_90s[0:-2] = smooth(mean_age_on_n2o_90s[0:-2],3,/edge_truncate)
ni = interpol(indgen(n_elements(norm_grid)),norm_grid,norm_grid2)
n2o_age_trend_interp = interpolate(n2o_nofl_age_trends[*,1],ni)
n2o_age_trend_interp[-3:-1] = [-0.012,-0.01,-0.005]
mean_age_on_n2o_20s_from_trend = mean_age_on_n2o_90s + 25.*n2o_age_trend_interp
mean_age_on_n2o_20s_from_trend[19:-3] = smooth(mean_age_on_n2o_20s_from_trend[19:-3],3,/edge_truncate)
mean_age_on_n2o_20s_from_trend_coarse = interpolate(mean_age_on_n2o_20s_from_trend,nic)

;  Calculate mean age trend based on N2O trend and constant mean age vs. N2O relationship.
aaa = replicate(!values.f_nan,nn) & bbb = replicate(!values.f_nan,2,nz2) & mean_age_from_n2o_trend = aaa & n2o_norm_alt_profile_avg = bbb & age_alt_profile_avg = bbb
n2o_norm_alt_profile_avg_new = bbb & age_alt_profile_avg_new = bbb
ti = where(years lt 2000)
ti2 = where(years ge 2015)
for z = 0, nz2-1 do begin
  gd = where(finite(bins_all_combo_n2o_norm_alt_nh_years[ti,z,0]),ngd)
  if ngd gt 1 then for j = 0, 1 do n2o_norm_alt_profile_avg[j,z] = mean(bins_all_combo_n2o_norm_alt_nh_years[ti[gd],z,j])
  gd = where(finite(bins_all_combo_n2o_norm_alt_nh_years[ti2,z,0]),ngd)
  if ngd gt 1 then for j = 0, 1 do n2o_norm_alt_profile_avg_new[j,z] = mean(bins_all_combo_n2o_norm_alt_nh_years[ti2[gd],z,j])
  gd = where(finite(bins_all_combo_age_alt_nh_years[ti,z,0]),ngd)
  if ngd gt 1 then for j = 0, 1 do age_alt_profile_avg[j,z] = mean(bins_all_combo_age_alt_nh_years[ti[gd],z,j])
  gd = where(finite(bins_all_combo_age_alt_nh_years[ti2,z,0]),ngd)
  if ngd gt 1 then for j = 0, 1 do age_alt_profile_avg_new[j,z] = mean(bins_all_combo_age_alt_nh_years[ti2[gd],z,j])
endfor
gd = where(finite(n2o_norm_alt_profile_avg))
n2o_norm_alt_profile_avg[gd] = smooth(n2o_norm_alt_profile_avg[gd],5,/edge_truncate)
gd = where(finite(n2o_norm_alt_profile_avg_new))
n2o_norm_alt_profile_avg_new[gd] = smooth(n2o_norm_alt_profile_avg_new[gd],3,/edge_truncate)
gd = where(finite(age_alt_profile_avg[0,*]))
age_alt_profile_avg[0,gd] = smooth(age_alt_profile_avg[0,gd],5,/edge_truncate)
gd = where(finite(age_alt_profile_avg_new[0,*]))
age_alt_profile_avg_new[0,gd] = smooth(age_alt_profile_avg_new[0,gd],3,/edge_truncate)

age_alt_profile_from_trend = 25.*age_alt_trends[*,1] + reform(age_alt_profile_avg[0,*])

zi = interpol(indgen(nz2),n2o_norm_alt_profile_avg[0,*],norm_grid[5:-1])
norm_grid_n2o_alts = interpolate(alt_grid2,zi)

;  Add ACE N2O trends onto top of in situ N2O trend profiles.
n2o_alt_trends_combo = [n2o_alt_trends[0:15,1],-5e-4,0,3e-4,6e-4]
n2o_alt_trends_combo = n2o_nh_lat_adj_trends_coarse[*,1]
n2o_alt_trends_combo[22] = 0.
n2o_norm_alt_profile_from_trend = 25.*n2o_alt_trends_combo + reform(n2o_norm_alt_profile_avg[0,*])
n2o_trends_interp = interpolate(n2o_alt_trends_combo,zi)

;  25 year mean age change from N2O trend based on constant mean age vs. N2O relationship.
ni = interpol(indgen(nn-5),norm_grid[5:-1],norm_grid[5:-1]+25.*n2o_trends_interp)
mean_age_from_n2o = interpolate(bins_a0_n2o_norm_age_nh_coarse[5:-1],ni)
mean_age_from_n2o_trend = 0.4 * (mean_age_from_n2o - bins_a0_n2o_norm_age_nh_coarse[5:-1])
mean_age_from_n2o_trend_per = 1e2 * 0.4 * (mean_age_from_n2o - bins_a0_n2o_norm_age_nh_coarse[5:-1]) / bins_a0_n2o_norm_age_nh_coarse[5:-1]


;  Categorize the balloon sampling type.
bal['samp_type_'+dates[i_flask]] = 'flask'
bal['samp_type_'+dates[i_jp]] = 'flask_jp'
bal['samp_type_'+dates[i_oms]] = 'in_situ'
bal['samp_type_'+dates[i_aircore]] = 'aircore'

restore,dir+'Balloon/Mean_age_and_N2O_trend_alt_profile_elat_adj.sav'


;  Read in the gridded ACE data.
restore,dir+'/Satellite/ACE/ACE_gridding_5p3.sav'
yim = where(lat_grid gt 30 and lat_grid lt 60)
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


restore,dir+'swoosh/n2o_merge.sav'
yic = where(slat gt 30 and slat lt 70)
nzc = n_elements(level)
ntc = n_elements(swoosh_year)
alt_s = -7.0*alog(level/1e3)
swoosh_n2o_nh = replicate(!values.f_nan,nzc,ntc) & swoosh_n2o_full_nh = swoosh_n2o_nh & swoosh_n2o_yrsm_nh = swoosh_n2o_nh & swoosh_n2o_nh_trends = replicate(!values.f_nan,2,nzc)
swoosh_n2o_nh_trends_chi = replicate(!values.f_nan,nzc)
for t = 0, ntc-1 do for zz = 0, nzc-1 do begin
  swoosh_n2o_nh[zz,t] = mean(combn2oq_norm_sm[yic,zz,t],/nan)
  swoosh_n2o_yrsm_nh[zz,t] = mean(combn2oq_norm_yrsm[yic,zz,t],/nan)
  swoosh_n2o_full_nh[zz,t] = mean(combn2oq_norm[yic,zz,t],/nan)
endfor
seas_cyc = findgen(12)*1./12.+1./24.
swoosh_n2o_nh_seas = replicate(!values.f_nan,nzc,12)
for t = 0, 11 do begin
  ti = where(swoosh_year-fix(swoosh_year) ge seas_cyc[t]-1./24. and swoosh_year-fix(swoosh_year) lt seas_cyc[t]+1./24.,nti)
  for zz = 0, nzc-1 do swoosh_n2o_nh_seas[zz,t] = mean(swoosh_n2o_full_nh[zz,ti] - swoosh_n2o_yrsm_nh[zz,ti],/nan)
endfor
for zz = 0, nzc-1 do swoosh_n2o_nh_seas[zz,*] = smooth(swoosh_n2o_nh_seas[zz,*],3)
for zz = 0, nzc-1 do begin
  chk = where(finite(swoosh_n2o_nh[zz,*]))
  swoosh_n2o_nh_trends[*,zz] = linfit(swoosh_year[chk],swoosh_n2o_nh[zz,chk],sigma=sigma)
  swoosh_n2o_nh_trends_chi[zz] = sigma[1]
endfor


n2o_nh_alt_trends_combo = ace_n2o_nh_trends
n2o_nh_alt_trends_combo[0,*] = ace_n2o_nh_trends[1,*]
n2o_nh_alt_trends_combo[1,*] = ace_n2o_nh_trends_chi
zia = interpol(indgen(nz2),alt_grid2,ace_alts)
n2o_nh_alt_trends_combo[0,10:26] = interpolate(n2o_alt_trends[0:16,1],zia[10:26])
n2o_nh_alt_trends_combo[0,26] = 0.5 * (n2o_nh_alt_trends_combo[0,26] + ace_n2o_nh_trends[1,26])
n2o_nh_alt_trends_combo[0,27:28] += [5e-4,2e-4]

norm_grid_obs = norm_grid


restore,dir+'Models/CCMI/CCMI-1_refc1_means.sav'
restore,dir+'Models/CCMI/CCMI-2022_refd1_means.sav'
restore,dir+'Models/CCMI/CCMI-2022_refd1_trends.sav'
restore,dir+'Models/CCMI/CCMI-2022_refd2_trends.sav'
ny_c = n_elements(lat_c)
nz_c = n_elements(pres_c)
nnm = n_elements(n2o_norm_grid_model)
alt_c = -7.*alog(pres_c/1e3)
nt_c = n_elements(years_c)
nm = n_elements(models)
ny_c1 = n_elements(lat_c1)
nz_c1 = n_elements(pres_c1)
pres_c1 /= 1e2
alt_c1 = -7.*alog(pres_c1/1e3)
nm1 = n_elements(models_c1)
nm1_n2o = n_elements(models_n2o_c1)
zz = fltarr(nz_c,2) & aoa_trends_post2000_per_minmax = zz & aoa_trends_pre2000_per_minmax = zz & n2o_norm_trends_post2000_minmax = zz & n2o_norm_trends_pre2000_minmax = zz
aoa_trends_post2000_minmax = zz & aoa_trends_pre2000_minmax = zz & ch4_norm_trends_pre2000_minmax = zz & ch4_norm_trends_post2000_minmax = zz & aoa_trends_minmax = zz
zzz = fltarr(nz_c) & wstar_trends_pre2000_minmax = zz & aoa_trends_d1_minmax = zz & aoa_trends_d1_avg = zzz
w_trends_d1_minmax = zz & aoa_trends_d2_minmax = zz & aoa_trends_d2_avg = zzz & n2o_norm_trends_d1_minmax = zz & n2o_norm_trends_d1_avg = zzz & w_trends_d1_avg = zzz
n2o_norm_trends_d2_minmax = zz & n2o_norm_trends_d2_avg = zzz & w_trends_d2_minmax = zz & w_trends_d2_avg = zzz & w_trends_d1_d2sub_minmax = zz
w_trends_d1_d2sub_avg = zzz & ch4_norm_trends_d1_minmax = zz & ch4_norm_trends_d1_avg = zzz & ch4_norm_trends_d2_minmax = zz & ch4_norm_trends_d2_avg = zzz
w_turn_trends_per_d2_minmax = zz & w_turn_trends_per_d2_avg = zzz & w_turn_trends_per_d1_minmax = zz & w_turn_trends_per_d1_avg = zzz
mi_gd = [indgen(22),indgen(3)+23]
wi_d2sub = [0,1,2,7,8]
for z = 0, nz_c-1 do begin
  wstar_trends_pre2000_minmax[z,0] = min(wstar_trends_pre2000[0,1,z,mi_d1_common],/nan)
  wstar_trends_pre2000_minmax[z,1] = max(wstar_trends_pre2000[0,1,z,mi_d1_common],/nan)
  w_trends_d1_minmax[z,0] = min(w_trends_d1[0,1,z,mi_d1_common,10:14],/nan)
  w_trends_d1_minmax[z,1] = max(w_trends_d1[0,1,z,mi_d1_common,10:14],/nan)
  w_trends_d1_avg[z] = mean(w_trends_avg_d1[0,1,z,10:14])
  w_trends_d1_d2sub_minmax[z,0] = min(w_trends_d1[0,1,z,mi_d1_common[wi_d2sub],10:14],/nan)
  w_trends_d1_d2sub_minmax[z,1] = max(w_trends_d1[0,1,z,mi_d1_common[wi_d2sub],10:14],/nan)
  w_trends_d1_d2sub_avg[z] = mean(w_trends_avg_d1_d2sub[0,1,z,10:14])
  w_trends_d2_minmax[z,0] = min(w_trends_d2[0,1,z,mi_d2_common,10:15],/nan)
  w_trends_d2_minmax[z,1] = max(w_trends_d2[0,1,z,mi_d2_common,10:15],/nan)
  w_trends_d2_avg[z] = mean(w_trends_avg_d2[0,1,z,10:15])
  aoa_trends_pre2000_per_minmax[z,0] = min(aoa_trends_pre2000_per[0,2,z,mi_gd],/nan)
  aoa_trends_pre2000_per_minmax[z,1] = max(aoa_trends_pre2000_per[0,2,z,mi_gd],/nan)
  aoa_trends_post2000_per_minmax[z,0] = min(aoa_trends_post2000_per[0,2,z,mi_gd],/nan)
  aoa_trends_post2000_per_minmax[z,1] = max(aoa_trends_post2000_per[0,2,z,mi_gd],/nan)
  aoa_trends_minmax[z,0] = min(aoa_trends[0,2,z,*],/nan)
  aoa_trends_minmax[z,1] = max(aoa_trends[0,2,z,*],/nan)
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
  aoa_trends_pre2000_minmax[z,0] = min(aoa_trends_pre2000[0,2,z,mi_d1_common],/nan)
  aoa_trends_pre2000_minmax[z,1] = max(aoa_trends_pre2000[0,2,z,mi_d1_common],/nan)
  aoa_trends_post2000_minmax[z,0] = min(aoa_trends_post2000[0,2,z,mi_gd],/nan)
  aoa_trends_post2000_minmax[z,1] = max(aoa_trends_post2000[0,2,z,mi_gd],/nan)
  n2o_norm_trends_pre2000_minmax[z,0] = min(n2o_norm_trends_pre2000[0,2,z,mi_d1_common],/nan)
  n2o_norm_trends_pre2000_minmax[z,1] = max(n2o_norm_trends_pre2000[0,2,z,mi_d1_common],/nan)
  n2o_norm_trends_post2000_minmax[z,0] = min(n2o_norm_trends_post2000[0,2,z,*],/nan)
  n2o_norm_trends_post2000_minmax[z,1] = max(n2o_norm_trends_post2000[0,2,z,*],/nan)
  n2o_norm_trends_d1_minmax[z,0] = min(n2o_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  n2o_norm_trends_d1_minmax[z,1] = max(n2o_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  n2o_norm_trends_d1_avg[z] = mean(n2o_norm_trends_avg_d1[1,2,z,10:14])
  n2o_norm_trends_d2_minmax[z,0] = min(n2o_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  n2o_norm_trends_d2_minmax[z,1] = max(n2o_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  n2o_norm_trends_d2_avg[z] = mean(n2o_norm_trends_avg_d2[1,2,z,10:15])
  ch4_norm_trends_pre2000_minmax[z,0] = min(ch4_norm_trends_pre2000[0,2,z,mi_d1_common],/nan)
  ch4_norm_trends_pre2000_minmax[z,1] = max(ch4_norm_trends_pre2000[0,2,z,mi_d1_common],/nan)
  ch4_norm_trends_post2000_minmax[z,0] = min(ch4_norm_trends_post2000[0,2,z,*],/nan)
  ch4_norm_trends_post2000_minmax[z,1] = max(ch4_norm_trends_post2000[0,2,z,*],/nan)
  ch4_norm_trends_d1_minmax[z,0] = min(ch4_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  ch4_norm_trends_d1_minmax[z,1] = max(ch4_norm_trends_d1[1,2,z,mi_d1_common,10:14],/nan)
  ch4_norm_trends_d1_avg[z] = mean(ch4_norm_trends_avg_d1[1,2,z,10:14])
  ch4_norm_trends_d2_minmax[z,0] = min(ch4_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  ch4_norm_trends_d2_minmax[z,1] = max(ch4_norm_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  ch4_norm_trends_d2_avg[z] = mean(ch4_norm_trends_avg_d2[1,2,z,10:15])
endfor

yy = fltarr(nnm,2) & aoa_on_n2o_trends_pre2000_per_minmax = yy & aoa_on_n2o_trends_post2000_per_minmax = yy & aoa_on_n2o_trends_pre2000_minmax = yy & aoa_on_n2o_trends_post2000_minmax = yy
aoa_on_ch4_trends_pre2000_per_minmax = yy & aoa_on_ch4_trends_post2000_per_minmax = yy & aoa_on_ch4_trends_pre2000_minmax = yy & aoa_on_ch4_trends_post2000_minmax = yy
aoa_on_n2o_trends_minmax = yy & aoa_on_n2o_trends_d1_avg = fltarr(nnm) & aoa_on_n2o_trends_d2_avg = fltarr(nnm) & aoa_on_n2o_trends_d2_minmax = yy
aoa_on_ch4_trends_d1_minmax = yy & aoa_on_ch4_trends_d1_avg = fltarr(nnm) & aoa_on_ch4_trends_d2_minmax = yy & aoa_on_ch4_trends_d2_avg = fltarr(nnm)
ww = fltarr(nnm,2,3) & aoa_on_n2o_trends_d1_minmax = ww
for z = 0, nnm-1 do begin
  aoa_on_n2o_trends_minmax[z,0] = min(aoa_on_n2o_trends[0,2,z,*],/nan)
  aoa_on_n2o_trends_minmax[z,1] = max(aoa_on_n2o_trends[0,2,z,*],/nan)
  aoa_on_n2o_trends_pre2000_minmax[z,0] = min(aoa_on_n2o_trends_pre2000[0,2,z,mi_d1_common],/nan)
  aoa_on_n2o_trends_pre2000_minmax[z,1] = max(aoa_on_n2o_trends_pre2000[0,2,z,mi_d1_common],/nan)
  aoa_on_n2o_trends_post2000_minmax[z,0] = min(aoa_on_n2o_trends_d1[0,2,z,mi_d1_common,10:14],/nan)
  aoa_on_n2o_trends_post2000_minmax[z,1] = max(aoa_on_n2o_trends_d1[0,2,z,mi_d1_common,10:14],/nan)
  for j = 0, 2 do begin
    aoa_on_n2o_trends_d1_minmax[z,0,j] = min(aoa_on_n2o_trends_d1[j,2,z,mi_d1_common,10:14],/nan)
    aoa_on_n2o_trends_d1_minmax[z,1,j] = max(aoa_on_n2o_trends_d1[j,2,z,mi_d1_common,10:14],/nan)
  endfor
  aoa_on_n2o_trends_d1_avg[z] = mean(aoa_on_n2o_trends_avg_d1[1,2,z,10:14])
  aoa_on_n2o_trends_d2_minmax[z,0] = min(aoa_on_n2o_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_on_n2o_trends_d2_minmax[z,1] = max(aoa_on_n2o_trends_d2[1,2,z,mi_d2_common,10:15],/nan)
  aoa_on_n2o_trends_d2_avg[z] = mean(aoa_on_n2o_trends_avg_d2[1,2,z,10:15])
  aoa_on_n2o_trends_pre2000_per_minmax[z,0] = min(aoa_on_n2o_trends_pre2000_per[0,2,z,*],/nan)
  aoa_on_n2o_trends_pre2000_per_minmax[z,1] = max(aoa_on_n2o_trends_pre2000_per[0,2,z,*],/nan)
  aoa_on_n2o_trends_post2000_per_minmax[z,0] = min(aoa_on_n2o_trends_post2000_per[0,2,z,*],/nan)
  aoa_on_n2o_trends_post2000_per_minmax[z,1] = max(aoa_on_n2o_trends_post2000_per[0,2,z,*],/nan)
  gd = where(abs(aoa_on_ch4_trends_pre2000[0,2,z,*]) lt 0.04)
  aoa_on_ch4_trends_pre2000_minmax[z,0] = min(aoa_on_ch4_trends_pre2000[0,2,z,mi_d1_common],/nan)
  aoa_on_ch4_trends_pre2000_minmax[z,1] = max(aoa_on_ch4_trends_pre2000[0,2,z,mi_d1_common],/nan)
  gd = where(abs(aoa_on_ch4_trends_post2000[0,2,z,*]) lt 0.01)
  aoa_on_ch4_trends_post2000_minmax[z,0] = min(aoa_on_ch4_trends_post2000[0,2,z,gd],/nan)
  aoa_on_ch4_trends_post2000_minmax[z,1] = max(aoa_on_ch4_trends_post2000[0,2,z,gd],/nan)
  aoa_on_ch4_trends_d1_minmax[z,0] = min(aoa_on_ch4_trends_d1[0,2,z,mi_d1_common,10:14],/nan)
  aoa_on_ch4_trends_d1_minmax[z,1] = max(aoa_on_ch4_trends_d1[0,2,z,mi_d1_common,10:14],/nan)
  aoa_on_ch4_trends_d1_avg[z] = mean(aoa_on_ch4_trends_avg_d1[0,2,z,10:14])
  aoa_on_ch4_trends_d2_minmax[z,0] = min(aoa_on_ch4_trends_d2[0,2,z,mi_d2_common,10:15],/nan)
  aoa_on_ch4_trends_d2_minmax[z,1] = max(aoa_on_ch4_trends_d2[0,2,z,mi_d2_common,10:15],/nan)
  aoa_on_ch4_trends_d2_avg[z] = mean(aoa_on_ch4_trends_avg_d2[0,2,z,10:15])
  aoa_on_ch4_trends_pre2000_per_minmax[z,0] = min(aoa_on_ch4_trends_pre2000_per[0,2,z,mi_d1_common],/nan)
  aoa_on_ch4_trends_pre2000_per_minmax[z,1] = max(aoa_on_ch4_trends_pre2000_per[0,2,z,mi_d1_common],/nan)
  aoa_on_ch4_trends_post2000_per_minmax[z,0] = min(aoa_on_ch4_trends_post2000_per[0,2,z,*],/nan)
  aoa_on_ch4_trends_post2000_per_minmax[z,1] = max(aoa_on_ch4_trends_post2000_per[0,2,z,*],/nan)
endfor


n2o_norm_aoa_d1_minmax = replicate(!values.f_nan,2,3,nnm,nt_c)
for t = 0, nt_c-1 do for z = 0, nnm-1 do for y = 0, 2 do begin
  n2o_norm_aoa_d1_minmax[0,y,z,t] = min(n2o_norm_aoa_latavgs_c[y,z,t,mi_d1_common],/nan)
  n2o_norm_aoa_d1_minmax[1,y,z,t] = max(n2o_norm_aoa_latavgs_c[y,z,t,mi_d1_common],/nan)
endfor


;  Make CCMI time averages.
ccc = replicate(!values.f_nan,ny_c,nz_c,nnm) & fff = replicate(!values.f_nan,ny_c1,nz_c1,nm1) & n2o_norm_c_1990s = ccc & age_c_1990s = ccc & age_c1_1990s = fff
ggg = replicate(!values.f_nan,ny_c1,nz_c1,nm1_n2o) & n2o_norm_c1_1990s = ggg
ti1 = where(years_c ge 1990 and years_c lt 2000)
for m = 0, nnm-1 do for z = 0, nz_c-1 do for y = 0, ny_c-1 do begin
  n2o_norm_c_1990s[y,z,m] = mean(n2o_norm_d1[y,z,ti1,m])
  age_c_1990s[y,z,m] = mean(age_d1[y,z,ti1,m])
endfor
ti1 = where(years_c1 ge 1990 and years_c1 lt 2000)
for m = 0, nm1-1 do for z = 0, nz_c1-1 do for y = 0, ny_c1-1 do begin
  age_c1_1990s[y,z,m] = mean(age_c1[y,z,ti1,m])
endfor
ti1 = where(years_n2o_c1 ge 1990 and years_n2o_c1 lt 2000)
for m = 0, nm1_n2o-1 do for z = 0, nz_c1-1 do for y = 0, ny_c1-1 do begin
  n2o_norm_c1_1990s[y,z,m] = mean(n2o_norm_c1[y,z,ti1,m])
endfor

; Interpolate to obs vertical grid.
ddd = replicate(!values.f_nan,ny_c,nz,nnm) & n2o_norm_c_1990s_zobs = ddd & age_c_1990s_zobs = ddd & eee = replicate(!values.f_nan,ny_c1,nz,nm1) & age_c1_1990s_zobs = eee
hhh = replicate(!values.f_nan,ny_c1,nz,nm1_n2o) & n2o_norm_c1_1990s_zobs = hhh
zic = interpol(indgen(nz_c),alt_c,alt_grid_a0)
for m = 0, nnm-1 do for y = 0, ny_c-1 do begin
  n2o_norm_c_1990s_zobs[y,*,m] = interpolate(n2o_norm_c_1990s[y,*,m],zic)
  age_c_1990s_zobs[y,*,m] = interpolate(age_c_1990s[y,*,m],zic)
endfor
zic1 = interpol(indgen(nz_c1),alt_c1,alt_grid_a0)
for m = 0, nm1-1 do for y = 0, ny_c1-1 do age_c1_1990s_zobs[y,*,m] = interpolate(age_c1_1990s[y,*,m],zic1)
for m = 0, nm1_n2o-1 do for y = 0, ny_c1-1 do n2o_norm_c1_1990s_zobs[y,*,m] = interpolate(n2o_norm_c1_1990s[y,*,m],zic1)


save,norm_grid_obs,bins_a0_n2o_norm_age_nh_coarse,bins_a2_n2o_norm_age_nh_trend,filename=dir+'Aircraft/Airborne_Save_Files/N2O_vs_mean_age_coarse.sav'

;  Write out the coarse mean age vs. N2O average relationships into an ascii file.
;openw,1,dir+'Aircraft/Airborne_Save_Files/N2O_vs_mean_age_coarse.txt'
;printf,1,'Average NH midlatitude relationships based on aircraft and balloon data'
;printf,1,'Normalized N2O    Mean Age 1990s    Mean Age 2020s'
;for i = 0, nn-1 do printf,1,norm_grid[i],bins_a0_n2o_norm_age_nh_coarse[i],bins_a2_n2o_norm_age_nh_trend[i]
;close,1


nseas_old = fltarr(nyg)
for y = 0, nyg-1 do begin
  chk = where(finite(age_co2_old_grid_seas_sm[0,y,z,*]),nchk)
  nseas_old[y] = nchk
endfor


;  Read in TLP model output

restore,dir+'Models/TLP/ideal_sweep_w_LowStrat_UpStrat_Budget.sav'
tlp_lu = tlp & run_key_lu = run_key

nz_tlp = n_elements(alt_tlp)
nnt = 21
dn2ot = 0.05
n2o_norm_grid_tlp = findgen(nnt)*dn2ot
ntlp = n_elements(run_key_lu)-1

alt_sw = -7.*alog(level/1e3)

age_nh_diff_alt_profile = 30.*age_nh_lat_adj_trends_coarse[*,1]
n2o_nh_diff_alt_profile = 30.*n2o_nh_lat_adj_trends_per_coarse[*,1]
n2o_nh_diff_alt_profile_sigma = 30.*n2o_nh_lat_adj_trends_per_coarse_sigma[*,1]

n2o_nh_alt_trends_per_combo = n2o_nh_alt_trends_combo
for i = 0, 1 do n2o_nh_alt_trends_per_combo[i,*] = 1e2*reform(n2o_nh_alt_trends_combo[i,*])/ace_n2o_nh_avg

for z = 23, 24 do begin
  n2o_nh_diff_alt_profile[z] = 30 * 0.5 * total(n2o_nh_alt_trends_per_combo[0,z+5:z+6])
  n2o_nh_diff_alt_profile_sigma[z] = 30 * 0.5 * total(n2o_nh_alt_trends_per_combo[1,z+5:z+6])
endfor

tlp_diffs = hash()

ni = interpol(indgen(nnt),n2o_norm_grid_tlp,n2o_norm_grid)
mm = fltarr(ntlp+1) & n2o_mean_age_tot = mm 

for ii = 1, ntlp do begin
  tlp_diffs[run_key_lu[ii]+'_mean_age'] = tlp_lu[run_key_lu[ii]+'_mean_age'] - tlp_lu[run_key_lu[0]+'_mean_age']
  tlp_diffs[run_key_lu[ii]+'_n2o'] = tlp_lu[run_key_lu[ii]+'_n2o'] - tlp_lu[run_key_lu[0]+'_n2o']
  tlp_diffs[run_key_lu[ii]+'_n2o_per'] = 1e2*tlp_diffs[run_key_lu[ii]+'_n2o'] / tlp_lu[run_key_lu[0]+'_n2o']

  n2o_mean_age_diffs = replicate(!values.f_nan,3,n_elements(n2o_norm_grid))
  for y = 0, 2 do n2o_mean_age_diffs[y,*] = interpolate(tlp_lu[run_key_lu[ii]+'_mean_age_n2o',y,*]-tlp_lu[run_key_lu[0]+'_mean_age_n2o',y,*],ni)
  tlp_diffs[run_key_lu[ii]+'_n2o_mean_age'] = n2o_mean_age_diffs
endfor


aa = replicate(!values.f_nan,ntlp+1) & diffs_age_nh = aa & rmse_age_nh = aa & diffs_n2o_nh = aa & rmse_n2o_nh = aa & diffs_up_n2o_nh = aa
rmse_up_n2o_nh = aa & diffs_mid_age_nh = aa & rmse_mid_age_nh = aa & mid_w = aa

;  Mean age profile differences from 16-17km.
nage = 2.
age_nh_i = indgen(nage)
;  Mean age profile differences from 22-27km.
nage_mid = 6.
age_mid_nh_i = indgen(nage_mid)+6
;  N2O profile differences from 20-27km.
nn2o = 10.
n2o_nh_i = indgen(nn2o)+5
;  N2O profile differences from 26-29km.
nu_n2o = 8.
n2o_up_nh_i = indgen(nu_n2o)+7

for i = 1, ntlp do begin
  mid_w[i] = float(strmid(run_key_lu[i],3,4,/reverse_offset))
  mean_age_tlp_nh_diffs = reform(tlp_diffs[run_key_lu[i]+'_mean_age',2,*])
  diffs_age_nh[i] = total(abs(mean_age_tlp_nh_diffs[age_nh_i] - age_nh_diff_alt_profile[age_nh_i+10]))
  rmse_age_nh[i] = sqrt(total(mean_age_tlp_nh_diffs[age_nh_i] - age_nh_diff_alt_profile[age_nh_i+10])^2/nage)
  diffs_mid_age_nh[i] = total(abs(mean_age_tlp_nh_diffs[age_mid_nh_i] - age_nh_diff_alt_profile[age_mid_nh_i+10]))
  rmse_mid_age_nh[i] = sqrt(total(mean_age_tlp_nh_diffs[age_mid_nh_i] - age_nh_diff_alt_profile[age_mid_nh_i+10])^2/nage)

  n2o_tlp_nh_diffs = reform(tlp_diffs[run_key_lu[i]+'_n2o_per',2,*])
  diffs_n2o_nh[i] = total(abs(n2o_tlp_nh_diffs[n2o_nh_i] - n2o_nh_diff_alt_profile[n2o_nh_i+10]))
  rmse_n2o_nh[i] = sqrt(total(n2o_tlp_nh_diffs[n2o_nh_i] - n2o_nh_diff_alt_profile[n2o_nh_i+10])^2/nn2o)
  diffs_up_n2o_nh[i] = total(abs(n2o_tlp_nh_diffs[n2o_up_nh_i] - n2o_nh_diff_alt_profile[n2o_up_nh_i+10]))
  rmse_up_n2o_nh[i] = sqrt(total(n2o_tlp_nh_diffs[n2o_up_nh_i] - n2o_nh_diff_alt_profile[n2o_up_nh_i+10])^2/nu_n2o)
endfor

rmse_age_thresh = 0.16
rmse_n2o_thresh = 1
rmse_up_n2o_thresh = 15
best_age_nh_i = where(rmse_age_nh le rmse_age_thresh,nb_ra)
best_age_mid_nh_i = where(rmse_mid_age_nh le rmse_age_thresh,nb_rma)
best_n2o_nh_i = where(rmse_n2o_nh le rmse_n2o_thresh,nb_rn)

;  Best N2O within the age NH set of solutions.
best_age_n2o_nh_i = where(rmse_n2o_nh[best_age_nh_i] le rmse_up_n2o_thresh and rmse_mid_age_nh[best_age_nh_i] le 2.5*rmse_age_thresh,nbu_rn)



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

yi = interpol(findgen(ny_tp),lat_tp,lat_grid_a)
tp_alt_grid = interpolate(tpause_alt,yi)



restore,dir+'Aircraft/Airborne_Save_Files/Andrews_N2O_age.sav'
restore,dir+'Models/TLP/ideal_sweep_w_LowStrat.sav'
tlp_low = tlp & run_key_low = run_key
restore,dir+'Models/TLP/ideal_sweep_w_UpStrat_budget.sav'


mission_names = ['aaoe','aase','aase2','ashoe','spade','strat','polaris','solve','accent','crystalf','pre_ave','ave_0506','cr_ave','tc4','glopac','attrex', $
  'seac4rs','ATom1','ATom2','ATom3','ATom4','dcotts1','dcotts2','sabre']


colors = ['sky blue','Blue','dark blue','Purple','Medium Purple','Magenta','Lime Green','Green','gold','Orange','orange red','Red','Brown','Grey','Dark Grey','dodger blue','cornflower', $
  'Purple','Medium Purple','Magenta']

LOADCT, 39, RGB_TABLE = rgb
LOADCT, 70, RGB_TABLE = rgb_70

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb']
seas = ['DJF','MAM','JJA','SON']

ia1 = 101
ia2 = n_aircore-1
z = 13
ii = z*5+0
tti = 15

restore,dir+'Balloon/NSamp_OMS.sav'
tot_samp = fltarr(nz2) & tot_samp_aircraft = tot_samp & tot_samp_aircore = tot_samp & tot_samp_flask = tot_samp
for z = 0, nz2-1 do begin
  si = where(age_nh_coarse[5,z,*] eq 1)
  tot_samp_aircraft[z] = total(age_nh_coarse[4,z,si],/nan)
  si = where(age_nh_coarse[5,z,*] eq 3)
  tot_samp_aircore[z] = total(age_nh_coarse[4,z,si],/nan)
  tot_samp[z] = tot_samp_aircraft[z] + tot_samp_aircore[z] + nsamp_flask[z] + nsamp_oms[z]
endfor


;  Read in aircraft and surface data for Supplement Fig. 5
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


n_levels = 21
ll = findgen(n_levels)*0.015
lls = findgen(n_levels/2+1)*0.03
tickname = strmid(strcompress(string(lls),/r),0,4)

y = 26
z = 12
z2 = 12

tickname = strmid(strcompress(string(findgen(7)*0.1+0.4),/r),0,3)
tickname1 = strmid(strcompress(string(findgen(12)*0.5),/r),0,3)
tickname4 = strmid(strcompress(string(findgen(11)*0.02-0.1),/r),0,5)
tickname4[5] = '0'
tickname5 = strmid(strcompress(string(findgen(9)*0.2-0.8),/r),0,4)

m = 8
s = 3


;  Ray et al 2026 paper plots.

;  Figure 1
  
;  p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[0,7.3],ytitle='Normalized N$_2$O',xtitle='Mean Age (years)',title='Mean Age vs. N$_2$O',font_size=11,position=[0.1,0.1,0.97,0.9], $
;    dimensions=[700,500])
;  c = contour(mean_age_all_on_n2o_hist,age_grid_h,norm_grid2,c_value=[0.03,0.03],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['n2o_norm_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75, $
;    linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_'+dates[i_flask[i]]],bal['n2o_norm_'+dates[i_flask[i]]],color='lime green',symbol='tu',sym_size=1, $
;    linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['n2o_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o', $
;    sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_'+dates[i_flask[i]]],bal['n2o_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='tu', $
;    sym_size=1,linestyle=6,/overplot)
;  p = errorplot(age_q_combo_n2o_bin_fl_avg[*,2],norm_grid,age_q_combo_n2o_bin_fl_avg[*,3],replicate(0,nn),symbol='s',/sym_filled,thick=2,color='dark green', $
;    linestyle=2,errorbar_capsize=0,/overplot)
;  for i = 26, 26 do p = plot(tlp_low[run_key_low[i]+'_mean_age',2,0:-20],tlp_low[run_key_low[i]+'_n2o',2,0:-20],thick=3,color='orange',linestyle=1,/overplot)
;  for i = 47, 47 do p = plot(tlp[run_key[i]+'_mean_age',2,0:-16],tlp[run_key[i]+'_n2o',2,0:-16],thick=3,color='orange',linestyle=2,/overplot)
;  p = errorplot(mean_age_all_on_n2o_a0[*,2],norm_grid2,mean_age_all_on_n2o_a0[*,3],replicate(0,n_elements(norm_grid2)),symbol='s',color='dodger blue',sym_size=1, $
;    /sym_filled,linestyle=6,errorbar_capsize=0,/overplot)
;  p = errorplot(mean_age_all_on_n2o_b0[*,2],norm_grid2,mean_age_all_on_n2o_b0[*,3],replicate(0,n_elements(norm_grid2)),symbol='s',color='purple',sym_size=1,/sym_filled, $
;    linestyle=6,errorbar_capsize=0,/overplot)
;  s = symbol(0.69,0.28,'s',sym_size=1,/sym_filled,sym_color='purple',/norm)
;  t = text(0.71,0.27,'In situ balloon 1990s',color='purple',font_size=11,/norm)
;  s = symbol(0.69,0.24,'s',sym_size=1,/sym_filled,sym_color='dodger blue',/norm)
;  t = text(0.71,0.23,'In situ aircraft 1990s',color='dodger blue',font_size=11,/norm)
;  p = plot([4.7,5],[0.9,0.9],thick=3,linestyle=1,color='orange',/overplot)
;  t = text(0.71,0.19,'Idealized model',color='orange',font_size=11,/norm)
;  p = plot([4.7,5],[0.95,0.95],thick=3,linestyle=2,color='orange',/overplot)
;  t = text(0.71,0.15,'w+20%',color='orange',font_size=11,/norm)
;  t = text(0.16,0.82,'Flask balloon 1970s-2000s',color='dark green',font_size=11,/norm)
;  s = symbol(0.17,0.79,'s',sym_size=1,/sym_filled,sym_color='dark green',/norm)
;  t = text(0.18,0.78,' = Weighted mean $\pm$ uncert. on mean',color='dark green',font_size=11,/norm)
;  s = symbol(0.17,0.75,'td',sym_size=1.5,/sym_filled,sym_color='lime green',/norm)
;  t = text(0.18,0.74,' = CO$_2$',color='lime green',font_size=11,/norm)
;  s = symbol(0.17,0.71,'o',sym_size=1,/sym_filled,sym_color='lime green',/norm)
;  t = text(0.18,0.7,' = SF$_6$',color='lime green',font_size=11,/norm)
;;  p.save,dir+'Plots/Age Trend Paper/Figure_1.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_1.pdf'

;  Figure 2

;  save,alt_c,aoa_trends_d2_minmax,aoa_trends_d1_minmax,aoa_trends_d1_avg,aoa_trends_d2_avg,age_nh_lat_adj_trends_coarse,age_nh_lat_adj_trends_coarse_sigma,nz2, $
;    alt_grid2,filename=dir+'Plots/Age Trend Paper/Fig_2a.sav'
;  
;  p = plot(indgen(2),/nodata,xrange=[-0.3,0.07],yrange=[13.5,30],ytitle='Altitude (km)',xtitle='Years/Decade',title='Mean Age',font_size=10,margin=[0.1,0.1,0.04,0.08], $
;    dimensions=[450,500])
;  p = plot([0,0],[13,35],linestyle=2,/overplot)
;;  poly = polygon([10*aoa_trends_d2_minmax[7:30,0],reverse(10*aoa_trends_d2_minmax[7:30,1])],[alt_c[7:30],reverse(alt_c[7:30])],/data,/fill_background, $
;;    fill_color='violet',fill_transparency=70,color='pink')
;;  poly = polygon([10*aoa_trends_d2_minmax[7:30,0],reverse(10*aoa_trends_d2_minmax[7:30,1])],[alt_c[7:30],reverse(alt_c[7:30])],/data,/fill_background, $
;;    fill_color='light pink',pattern_orientation=45,pattern_spacing=10,color='pink')
;;  poly = polygon([10*aoa_trends_d1_minmax[*,0],reverse(10*aoa_trends_d1_minmax[*,1])],[alt_c,reverse(alt_c)],/data,/fill_background,fill_color='lime green', $
;;    fill_transparency=65,color='light green')
;;  for ii = 0, n_elements(mi_d2_common)-1 do for i = 10, 15 do p = plot(10*aoa_trends_d2[1,2,*,mi_d2_common[ii],i],alt_c,color='pink',thick=1,/overplot)
;  p = plot(10.*aoa_trends_d1_minmax[*,0],alt_c,color='lime green',thick=1,/overplot)
;  p = plot(10.*aoa_trends_d1_minmax[*,1],alt_c,color='lime green',thick=1,/overplot)
;  p = plot(10.*aoa_trends_d2_minmax[*,0],alt_c,color='violet',thick=1,/overplot)
;  p = plot(10.*aoa_trends_d2_minmax[*,1],alt_c,color='violet',thick=1,/overplot)
;  p = plot(10.*aoa_trends_d1_avg,alt_c,color='dark green',thick=3,/overplot)
;  p = plot(10.*aoa_trends_d2_avg,alt_c,color='medium orchid',thick=3,/overplot)
;;  p = plot(10*aoa_trends_d1_minmax[*,0],alt_c,color='lime green',thick=3,/overplot)
;;  p = plot(10*aoa_trends_d1_minmax[*,1],alt_c,color='lime green',thick=3,/overplot)
;;  p = plot(10*aoa_trends_d2_minmax[*,0],alt_c,color='pink',thick=3,/overplot)
;;  p = plot(10*aoa_trends_d2_minmax[*,1],alt_c,color='pink',thick=3,/overplot)
;  p = errorplot(10.*age_nh_lat_adj_trends_coarse[*,1],alt_grid2,10.*age_nh_lat_adj_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1.5, $
;    /sym_filled,color='blue',linestyle=6,errorbar_capsize=0,/overplot)
;  p = errorplot([-0.14,-0.1],[15.5,18.5],[0.06,0.05],replicate(!values.f_nan,2),symbol='s',sym_size=1.25,color='black',linestyle=6,errorbar_capsize=0,/overplot)
;  s = symbol(0.17,0.51,'o',sym_size=1.5,sym_color='blue',/sym_filled,/norm)
;  t = text(0.2,0.5,'In situ',color='blue',font_size=10,/norm)
;  s = symbol(0.17,0.47,'s',sym_size=1.25,sym_color='black',/norm)
;  t = text(0.2,0.46,'ACE-FTS',color='black',font_size=10,/norm)
;;;  p.save,dir+'Plots/Age Trend Paper/Figure_2a.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_2a.pdf',/vector,resolution=300
;  p.save,dir+'Plots/Age Trend Paper/Figure_2a.svg';,bitmap=1
  
;save,alt_c,n2o_norm_trends_d2_minmax,n2o_norm_trends_d1_minmax,n2o_norm_trends_d1_avg,n2o_norm_trends_d2_avg,n2o_nh_lat_adj_trends_coarse, $
;  n2o_nh_lat_adj_trends_coarse_sigma,nz2,alt_grid2,ace_n2o_nh_trends,ace_n2o_nh_trends_chi,ace_alts,nza,filename=dir+'Plots/Age Trend Paper/Fig_2b.sav'
;
;  p = plot(indgen(2),/nodata,yrange=[13.5,30],xrange=[-0.04,0.019],ytitle='Altitude (km)',xtitle='Normalized Tracer Fraction/Decade',title='N$_2$O', $
;    font_size=10,margin=[0.1,0.1,0.04,0.08],dimensions=[450,500])
;  p = plot([0,0],[15,36],linestyle=2,/overplot)
;  poly = polygon([10*n2o_norm_trends_d2_minmax[7:30,0],reverse(10*n2o_norm_trends_d2_minmax[7:30,1])],[alt_c[7:30],reverse(alt_c[7:30])],/data,/fill_background, $
;    fill_color='violet',fill_transparency=70,color='pink')
;  poly = polygon([10*n2o_norm_trends_d1_minmax[*,0],reverse(10*n2o_norm_trends_d1_minmax[*,1])],[alt_c,reverse(alt_c)],/data,/fill_background,fill_color='lime green', $
;    fill_transparency=65,color='light green')
;;  p = errorplot(10.*ace_n2o_nh_trends[1,*],ace_alts,10.*ace_n2o_nh_trends_chi,replicate(!values.f_nan,nza),symbol='s',sym_size=1.25,color='black',errorbar_capsize=0,linestyle=6,/overplot)
;;  p = errorplot(10.*n2o_nh_lat_adj_trends_coarse[*,1],alt_grid2,10.*n2o_nh_lat_adj_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1.5, $
;;    /sym_filled,color='blue',linestyle=6,errorbar_capsize=0,/overplot)
;;  p = plot(10.*n2o_norm_trends_d1_minmax[*,0],alt_c,color='lime green',thick=1,/overplot)
;;  p = plot(10.*n2o_norm_trends_d1_minmax[*,1],alt_c,color='lime green',thick=1,/overplot)
;;  p = plot(10.*n2o_norm_trends_d2_minmax[*,0],alt_c,color='violet',thick=1,/overplot)
;;  p = plot(10.*n2o_norm_trends_d2_minmax[*,1],alt_c,color='violet',thick=1,/overplot)
;;  p = plot(10.*n2o_norm_trends_d1_avg,alt_c,color='dark green',thick=3,/overplot)
;;  p = plot(10.*n2o_norm_trends_d2_avg,alt_c,color='medium orchid',thick=3,/overplot)
;  t = text(0.175,0.31,'CCMI-2022',font_size=10,/norm)
;  p = plot([-0.035,-0.031],[17.1,17.1],thick=2,color='dark green',/overplot)
;  poly = polygon([-0.035,-0.031,-0.031,-0.035],[17.3,17.3,16.9,16.9],/data,/fill_background,fill_color='lime green',fill_transparency=70,color='light green')
;  t = text(0.25,0.27,'refD1',color='dark green',font_size=10,/norm)
;  p = plot([-0.035,-0.031],[16.3,16.3],thick=2,color='medium orchid',/overplot)
;  poly = polygon([-0.035,-0.031,-0.031,-0.035],[16.5,16.5,16.1,16.1],/data,/fill_background,fill_color='violet',fill_transparency=70,color='pink')
;  t = text(0.25,0.23,'refD2',color='medium orchid',font_size=10,/norm)
;;;  p.save,dir+'Plots/Age Trend Paper/Figure_2b.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_2b.pdf',/vector,resolution=300
;
;  p = plot(indgen(2),/nodata,xrange=[0,85],xtickinterval=20,xminor=1,yrange=[9,26],position=[0.09,0.1,0.87,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)', $
;    font_size=14,xticklen=0.04,yticklen=0.03,dimensions=[700,450],title='N$_2$O Changes 1990s to 2015-25')
;  for z = 0, nz-1 do begin
;    ii = where(finite(n2o_diffs[0:-3,z]) and nseas_n2o[0:-3,z] ge 3 and tp_alt_grid[0:-3]-0.5 le alt_grid[z],nii)
;    if nii gt 0 then begin
;      vcols = (-(n2o_diffs[ii,z]) + 0.1)*1275.
;      test = where(vcols ge 250.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 250.)] = 250.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;      p = plot(lat_grid_a[ii]-1.6,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;      p = plot(lat_grid_a[ii]+1.6,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1.7,/overplot)
;    endif
;    ii = where(finite(n2o_diffs[0:-3,z]) and nseas_n2o[0:-3,z] ge 1 and nseas_n2o[0:-3,z] le 2 and tp_alt_grid[0:-3]-0.5 le alt_grid[z],nii)
;    if nii gt 0 then begin
;      vcols = (-(n2o_diffs[ii,z]) + 0.1)*1275.
;      test = where(vcols ge 250.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 250.)] = 250.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1,/overplot)
;      p = plot(lat_grid_a[ii]-1,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1,/overplot)
;      p = plot(lat_grid_a[ii]+1,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=70,vert_colors=vcols,sym_size=1,/overplot)
;    endif
;  endfor
;  p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;  cb = COLORBAR(POSITION=[0.88,0.2,0.9,0.75],tickname=tickname4,tickvalues=findgen(11)*25.,RGB_TABLE=reverse(rgb_70),BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='dN$_2$O')
;  p.save,dir+'Plots/Age Trend Paper/Figure_2d.png'

  z = 12

;  Figure 3
    
;  p = plot(indgen(2),/nodata,xrange=[1990,2028],yrange=[3.4,5.2],xtickinterval=5,ytitle='Mean Age (years)',xtitle='Years',title='Mean Age in N$_2$O$_{norm}$ = '+ $
;    strmid(strcompress(string(norm_grid[z]-dn/2.),/r),0,4)+'-'+strmid(strcompress(string(norm_grid[z]+dn/2.),/r),0,4),font_size=12, $
;    margin=[0.1,0.1,0.04,0.08],dimensions=[700,400])
;  si = where(reform(age_on_n2o_adj_tseries[5,z,*]) eq 3 and n2o_from_ch4_yrs,nsi)
;  p = errorplot(years[si],age_on_n2o_adj_tseries[2,z,si],replicate(0,nsi),age_on_n2o_adj_tseries[3,z,si],color='magenta',symbol='d',sym_size=2,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  si = where(reform(age_on_n2o_adj_tseries[5,z,*]) eq 3 and n2o_from_ch4_yrs eq 0,nsi)
;  p = errorplot(years[si],age_on_n2o_adj_tseries[2,z,si],replicate(0,nsi),age_on_n2o_adj_tseries[3,z,si],color='magenta',symbol='d',/sym_filled,sym_size=2,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  si = where(reform(age_on_n2o_adj_tseries[5,z,*]) eq 1,nsi)
;  p = errorplot(years[si],age_on_n2o_adj_tseries[2,z,si],replicate(0,nsi),age_on_n2o_adj_tseries[3,z,si],color='blue',symbol='o',/sym_filled,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  si = where(reform(age_on_n2o_adj_tseries[5,z,*]) eq 2,nsi)
;  p = errorplot(years[si],age_on_n2o_adj_tseries[2,z,si],replicate(0,nsi),age_on_n2o_adj_tseries[3,z,si],color='purple',symbol='tu',sym_size=2,/sym_filled,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  p = plot([1992,2026],n2o_nofl_age_trends[z,0] + n2o_nofl_age_trends[z,1]*[1992,2026],color='blue',linestyle=2,/overplot)
;  p = plot([1990,2019],aoa_on_n2o_trends_avg_d1[0,2,z,12]+aoa_on_n2o_trends_avg_d1[0,2,z,12]*[1990,2019],color='lime green',thick=2,linestyle=2,/overplot)
;  s = symbol(0.17,0.24,'tu',sym_size=2,sym_color='purple',/sym_filled,/norm)
;  t = text(0.2,0.23,'In situ balloon',color='purple',font_size=11,/norm)
;  s = symbol(0.17,0.2,'o',sym_size=1.5,sym_color='blue',/sym_filled,/norm)
;  t = text(0.2,0.19,'In situ aircraft',color='blue',font_size=11,/norm)
;  s = symbol(0.65,0.83,'d',sym_size=2,sym_color='magenta',/sym_filled,/norm)
;  t = text(0.67,0.82,'AirCore N$_2$O measured',color='magenta',font_size=11,/norm)
;  s = symbol(0.65,0.78,'d',sym_size=2,sym_color='magenta',/norm)
;  t = text(0.67,0.77,'AirCore N$_2$O derived',color='magenta',font_size=11,/norm)
;;;  p.save,dir+'Plots/Age Trend Paper/Figure_3a.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_3a.pdf',/vector,resolution=300

;save,norm_grid,n2o_norm_grid_model,n2o_nofl_age_trends,n2o_nofl_age_trends_sigma,aoa_on_n2o_trends_d2_minmax,aoa_on_n2o_trends_d1_minmax, $
;  aoa_on_n2o_trends_d1_avg,nz2,aoa_on_n2o_trends_d2_avg,filename=dir+'Plots/Age Trend Paper/Fig_3b.sav'
;
;  p = plot(indgen(2),/nodata,xrange=[-0.31,0.11],yrange=[1.02,0.2],ytitle='Normalized N$_2$O',xtitle='Years/Decade',title='NH Mean Age Trends',font_size=10,axis_style=1, $
;    margin=[0.12,0.09,0.12,0.07],dimensions=[450,500])
;  p = plot([0,0],[1,0.1],linestyle=2,/overplot)
;  p = errorplot(10.*n2o_nofl_age_trends[7:-2,1],norm_grid[7:-2],10.*n2o_nofl_age_trends_sigma[7:-2,1],replicate(!values.f_nan,nn-8),symbol='o',sym_size=1.5, $
;    /sym_filled,color='blue',errorbar_capsize=0,/overplot)
;  poly = polygon([10*aoa_on_n2o_trends_d2_minmax[*,0],reverse(10*aoa_on_n2o_trends_d2_minmax[*,1])],[n2o_norm_grid_model,reverse(n2o_norm_grid_model)],/data, $
;    /fill_background,fill_color='violet',fill_transparency=70,color='pink')
;  poly = polygon([10*aoa_on_n2o_trends_d1_minmax[*,0,1],reverse(10*aoa_on_n2o_trends_d1_minmax[*,1,1])],[n2o_norm_grid_model,reverse(n2o_norm_grid_model)], $
;    /data,/fill_background,fill_color='lime green',fill_transparency=70,color='light green')
;;  poly = polygon([10*aoa_on_n2o_trends_post2000_minmax[*,0],reverse(10*aoa_on_n2o_trends_post2000_minmax[*,1])],[n2o_norm_grid_model,reverse(n2o_norm_grid_model)], $
;;    /data,/fill_background,fill_color='lime green',fill_transparency=70,color='light green')
;;  p = plot(10.*aoa_on_n2o_trends_d1_minmax[*,0,1],n2o_norm_grid_model,color='lime green',thick=1,/overplot)
;;  p = plot(10.*aoa_on_n2o_trends_d1_minmax[*,1,1],n2o_norm_grid_model,color='lime green',thick=1,/overplot)
;;  p = plot(10.*aoa_on_n2o_trends_d2_minmax[*,0],n2o_norm_grid_model,color='violet',thick=1,/overplot)
;;  p = plot(10.*aoa_on_n2o_trends_d2_minmax[*,1],n2o_norm_grid_model,color='violet',thick=1,/overplot)
;  p = plot(10.*aoa_on_n2o_trends_d1_avg,n2o_norm_grid_model,color='dark green',thick=3,/overplot)
;  p = plot(10.*aoa_on_n2o_trends_d2_avg[0:-2],n2o_norm_grid_model[0:-2],color='medium orchid',thick=3,/overplot)
;  xaxis = axis('X',location='top',tickfont_size=0)
;  yaxis = axis('Y',location='right',coord_transform=[0,1],minor=0,tickfont_size=10,tickname=['35','27','22','19','12'],title='Approximate Altitude (km)')
;  s = symbol(0.17,0.87,'o',sym_size=1.25,sym_color='blue',/sym_filled,/norm)
;  t = text(0.2,0.86,'In situ',color='blue',font_size=10,/norm)
;  t = text(0.42,0.865,'CCMI-2022',font_size=10,alignment=0.5,/norm)
;  p = plot([-0.19,-0.17],[0.29,0.29],thick=3,color='dark green',/overplot)
;  poly = polygon([-0.19,-0.17,-0.17,-0.19],[0.28,0.28,0.3,0.3],/data,/fill_background,fill_color='lime green',fill_transparency=70,color='light green')
;  t = text(0.44,0.83,'refD1',color='dark green',font_size=10,alignment=0.5,/norm)
;  p = plot([-0.19,-0.17],[0.33,0.33],thick=3,color='medium orchid',/overplot)
;  poly = polygon([-0.19,-0.17,-0.17,-0.19],[0.32,0.32,0.34,0.34],/data,/fill_background,fill_color='violet',fill_transparency=70,color='pink')
;  t = text(0.44,0.79,'refD2',color='medium orchid',font_size=10,alignment=0.5,/norm)
;;;  p.save,dir+'Plots/Age Trend Paper/Figure_3b.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_3b.pdf',/vector,resolution=300
;  p.save,dir+'Plots/Age Trend Paper/Figure_3b.svg'

;  p = plot(indgen(2),/nodata,yrange=[1.03,0.2],xrange=[0,6],ytitle='Normalized N$_2$O',xtitle='Mean Age (years)',title='NH Extratropical Mean Age vs. N$_2$O', $
;    font_size=12,margin=[0.1,0.1,0.04,0.08],dimensions=[700,400])
;  p = plot(mean_age_on_n2o_90s,norm_grid2,color='blue',thick=3,/overplot)
;  p = plot(mean_age_on_n2o_20s_from_trend,norm_grid2,color='red',thick=3,/overplot)
;  for z = 11, 11 do a = arrow([age_alt_profile_avg[0,z],age_alt_profile_from_trend[z]],0.01+[n2o_norm_alt_profile_avg[0,z],n2o_norm_alt_profile_from_trend[z]],/data, $
;    head_size=0.75,color='blue violet',/current)
;  for z = 13, 13 do a = arrow([age_alt_profile_avg[0,z],age_alt_profile_from_trend[z]],0.02+[n2o_norm_alt_profile_avg[0,z],n2o_norm_alt_profile_from_trend[z]],/data, $
;    head_size=0.75,color='blue violet',/current)
;  for z = 14, 14 do a = arrow([age_alt_profile_avg[0,z],age_alt_profile_from_trend[z]],0.01+[n2o_norm_alt_profile_avg[0,z],n2o_norm_alt_profile_from_trend[z]],/data, $
;    head_size=0.75,color='blue violet',/current)
;  for z = 15, 15 do a = arrow([age_alt_profile_avg[0,z],age_alt_profile_from_trend[z]],0.01+[n2o_norm_alt_profile_avg[0,z],n2o_norm_alt_profile_from_trend[z]],/data, $
;    head_size=0.75,color='blue violet',/current)
;  for z = 17, 17 do a = arrow([age_alt_profile_avg[0,z],age_alt_profile_from_trend[z]],-0.01+[n2o_norm_alt_profile_avg[0,z],n2o_norm_alt_profile_from_trend[z]],/data, $
;    head_size=0.75,color='blue violet',/current)
;  for z = 19, 19 do a = arrow([age_alt_profile_avg[0,z],age_alt_profile_from_trend[z]],[n2o_norm_alt_profile_avg[0,z],n2o_norm_alt_profile_from_trend[z]],/data, $
;    head_size=0.75,color='blue violet',/current)
;  t = text(0.18,0.8,'1990s',color='blue',font_size=11,/norm)
;  t = text(0.18,0.76,'2020s',color='red',font_size=11,/norm)
;;;  p.save,dir+'Plots/Age Trend Paper/Figure_3c.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_3c.pdf',/vector,resolution=300


;  Figure 4

  wdiffs = fltarr(nbu_rn,nz_tlp)
  for i = 0, nbu_rn-1 do wdiffs[i,*] = 1e2*(tlp_lu[run_key_lu[best_age_nh_i[best_age_n2o_nh_i[i]]]+'_w',1,*]-tlp_lu[run_key_lu[0]+'_w',1,*])/tlp_lu[run_key_lu[0]+'_w',1,*]
  diffs_minmax = fltarr(nz_tlp,2)
  for z = 0, nz_tlp-1 do diffs_minmax[z,*] = [min(wdiffs[1:-1,z]),max(wdiffs[1:-1,z])]
  diffs = 0.5 * (rmse_age_nh[best_age_nh_i[best_age_n2o_nh_i]] + rmse_n2o_nh[best_age_nh_i[best_age_n2o_nh_i]])
  weights = 1. / diffs[1:-1]^2
  wmean = fltarr(nz_tlp)
  for z = 0, nz_tlp-1 do wmean[z] = total(weights * wdiffs[1:-1,z]) / total(weights)

;  save,diffs_minmax,alt_tlp,w_turn_trends_per_d2_minmax,alt_c,w_turn_trends_per_d1_minmax,w_turn_trends_per_d2_avg,w_turn_trends_per_d1_avg, $
;     wmean,filename=dir+'Plots/Age Trend Paper/Fig_4.sav'
;  
;  p = plot(indgen(2),/nodata,yrange=[16,40],xrange=[-3,12],xtitle='%/Decade',ytitle='Altitude',title='Tropical Upwelling Trends',font_size=11, $
;    dimensions=[350,500],margin=[0.15,0.07,0.05,0.07])
;  poly = polygon([smooth(diffs_minmax[*,0],3,/edge_truncate),smooth(reverse(diffs_minmax[*,1]),3,/edge_truncate)]/3.,[alt_tlp,reverse(alt_tlp)],/data, $
;    /fill_background,fill_color='sky blue',fill_transparency=80,color='sky blue')
;  poly = polygon([10*w_turn_trends_per_d2_minmax[7:30,0],reverse(10*w_turn_trends_per_d2_minmax[7:30,1])],[alt_c[7:30],reverse(alt_c[7:30])],/data,/fill_background, $
;    fill_color='violet',fill_transparency=80,color='pink')
;  poly = polygon([10*w_turn_trends_per_d1_minmax[10:30,0],reverse(10*w_turn_trends_per_d1_minmax[10:30,1])],[alt_c[10:30],reverse(alt_c[10:30])],/data, $
;    /fill_background,fill_color='lime green',fill_transparency=80,color='light green')
;;  p = plot(10.*w_turn_trends_per_d1_minmax[*,0],alt_c,color='lime green',thick=1,/overplot)
;;  p = plot(10.*w_turn_trends_per_d1_minmax[*,1],alt_c,color='lime green',thick=1,/overplot)
;;  p = plot(10.*w_turn_trends_per_d2_minmax[*,0],alt_c,color='violet',thick=1,/overplot)
;;  p = plot(10.*w_turn_trends_per_d2_minmax[*,1],alt_c,color='violet',thick=1,/overplot)
;  p = plot(10.*w_turn_trends_per_d2_avg,alt_c,color='medium orchid',thick=3,/overplot)
;  p = plot(10.*w_turn_trends_per_d1_avg,alt_c,color='dark green',thick=3,/overplot)
;  p = plot(smooth(diffs_minmax[*,0]/3.,3,/edge_truncate),alt_tlp,thick=1,color='sky blue',/overplot)
;  p = plot(smooth(diffs_minmax[*,1]/3.,3,/edge_truncate),alt_tlp,thick=1,color='sky blue',/overplot)
;  p = plot(smooth(wmean/3.,3,/edge_truncate),alt_tlp,thick=3,color='blue',/overplot)
;  p = plot([0,0],[16,50],linestyle=2,thick=2,/overplot)
;  t = text(0.76,0.33,'CCMI-2022',font_size=10,alignment=0.5,/norm)
;  p = plot([6.5,7.5],[22.5,22.5],thick=3,color='dark green',/overplot)
;  poly = polygon([6.5,7.5,7.5,6.5],[22.3,22.3,22.7,22.7],/data,/fill_background,fill_color='lime green',fill_transparency=80,color='light green')
;  t = text(0.79,0.29,'refD1',color='dark green',font_size=10,alignment=0.5,/norm)
;  p = plot([6.5,7.5],[21.6,21.6],thick=3,color='medium orchid',/overplot)
;  poly = polygon([6.5,7.5,7.5,6.5],[21.4,21.4,21.8,21.8],/data,/fill_background,fill_color='violet',fill_transparency=80,color='pink')
;  t = text(0.79,0.26,'refD2',color='medium orchid',font_size=10,alignment=0.5,/norm)
;  p = plot([6.3,7.3],[19.4,19.4],thick=3,color='blue',/overplot)
;  poly = polygon([6.3,7.3,7.3,6.3],[19.2,19.2,19.6,19.6],/data,/fill_background,fill_color='sky blue',fill_transparency=80,color='sky blue')
;  t = text(0.81,0.18,'TLP model',color='blue',font_size=10,alignment=0.5,/norm)
;  ;p.save,dir_base+'Plots/TLP_trajecs/TLP_CCMI_tropical_w_trend_profiles.png'
;  p.save,dir+'Plots/Age Trend Paper/Figure_4.pdf',/vector,resolution=300


;  Extended Data Figure 1

;  p = plot(indgen(2),/nodata,yrange=[1.03,0.1],xrange=[0,7.5],ytitle='Normalized CH$_4$',xtitle='Mean Age (years)',title='Mean Age vs. CH$_4$',font_size=11,position=[0.1,0.1,0.97,0.9], $
;    dimensions=[700,500])
;  c = contour(mean_age_all_on_ch4_hist,age_grid_hist,norm_grid2,c_value=[0.03,0.03],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['ch4_norm_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_'+dates[i_flask[i]]],bal['ch4_norm_'+dates[i_flask[i]]],color='lime green',symbol='tu',sym_size=1.25,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['ch4_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_'+dates[i_flask[i]]],bal['ch4_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='tu',sym_size=1.25,linestyle=6,/overplot)
;  for i = 26, 26 do p = plot(tlp_low[run_key_low[i]+'_mean_age',2,0:-18],tlp_low[run_key_low[i]+'_ch4',2,0:-18],thick=3,color='orange',linestyle=1,/overplot)
;  for i = 47, 47 do p = plot(tlp[run_key[i]+'_mean_age',2,0:-14],tlp[run_key[i]+'_ch4',2,0:-14],thick=3,color='orange',linestyle=2,/overplot)
;  p = errorplot(age_q_combo_ch4_bin_fl_avg[*,2],norm_grid,age_q_combo_ch4_bin_fl_avg[*,3],replicate(0,nn),symbol='s',/sym_filled,thick=3,color='green',linestyle=2,errorbar_capsize=0,/overplot)
;  p = errorplot(mean_age_all_on_ch4_a0[*,2],norm_grid2,mean_age_all_on_ch4_a0[*,3],replicate(0,n_elements(norm_grid2)),symbol='s',color='dodger blue',sym_size=1,/sym_filled, $
;    linestyle=6,errorbar_capsize=0,/overplot)
;  p = errorplot(mean_age_all_on_ch4_b0[*,2],norm_grid2,mean_age_all_on_ch4_b0[*,3],replicate(0,n_elements(norm_grid2)),symbol='s',color='purple',sym_size=1,/sym_filled, $
;    linestyle=6,errorbar_capsize=0,/overplot)
;  t = text(0.16,0.82,'Flask balloon 1970s-2000s',color='dark green',font_size=11,/norm)
;  s = symbol(0.17,0.79,'s',sym_size=1,/sym_filled,sym_color='dark green',/norm)
;  t = text(0.18,0.78,' = Avg',color='dark green',font_size=11,/norm)
;  s = symbol(0.17,0.75,'td',sym_size=1.5,/sym_filled,sym_color='lime green',/norm)
;  t = text(0.18,0.74,' = CO$_2$',color='lime green',font_size=11,/norm)
;  s = symbol(0.17,0.71,'o',sym_size=1,/sym_filled,sym_color='lime green',/norm)
;  t = text(0.18,0.7,' = SF$_6$',color='lime green',font_size=11,/norm)
;;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_1a.png'
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_1a.jpg'
;
;  p = plot(indgen(2),/nodata,yrange=[1.03,-0.03],xrange=[0,7.5],ytitle='Normalized CFC-12',xtitle='Mean Age (years)',title='Mean Age vs. CFC-12',font_size=11, $
;    position=[0.1,0.1,0.97,0.9],dimensions=[700,500])
;  c = contour(mean_age_all_on_f12_hist,age_grid_hist,norm_grid2,c_value=[0.03,0.03],c_color=['plum','plum'],c_label_show=0,transparency=50,/fill,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['f12_norm_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_'+dates[i_flask[i]]],bal['f12_norm_'+dates[i_flask[i]]],color='lime green',symbol='tu',sym_size=1,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['f12_norm_q_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,/sym_filled,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_'+dates[i_flask[i]]],bal['f12_norm_q_'+dates[i_flask[i]]],color='lime green',symbol='tu',sym_size=1,/sym_filled,linestyle=6,/overplot)
;  p = errorplot(age_q_combo_f12_bin_fl_avg[*,2],norm_grid,age_q_combo_f12_bin_fl_avg[*,3],replicate(0,nn),symbol='s',/sym_filled,thick=3,color='green',linestyle=2,errorbar_capsize=0,/overplot)
;  for i = 26, 26 do p = plot(tlp_low[run_key_low[i]+'_mean_age',2,0:-20],tlp_low[run_key_low[i]+'_f12',2,0:-20],thick=3,color='orange',linestyle=1,/overplot)
;  for i = 47, 47 do p = plot(tlp[run_key[i]+'_mean_age',2,0:-16],tlp[run_key[i]+'_f12',2,0:-16],thick=3,color='orange',linestyle=2,/overplot)
;  p = errorplot(mean_age_all_on_f12_a0[*,2],norm_grid2,mean_age_all_on_f12_a0[*,3],replicate(0,n_elements(norm_grid2)),symbol='s',color='dodger blue',sym_size=1,/sym_filled,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  p = errorplot(mean_age_all_on_f12_b0[*,2],norm_grid2,mean_age_all_on_f12_b0[*,3],replicate(0,n_elements(norm_grid2)),symbol='s',color='purple',sym_size=1,/sym_filled,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  s = symbol(0.69,0.28,'s',sym_size=1,/sym_filled,sym_color='purple',/norm)
;  t = text(0.71,0.27,'In situ balloon 1990s',color='purple',font_size=11,/norm)
;  s = symbol(0.69,0.24,'s',sym_size=1,/sym_filled,sym_color='dodger blue',/norm)
;  t = text(0.71,0.23,'In situ aircraft 1990s',color='dodger blue',font_size=11,/norm)
;  p = plot([4.8,5.1],[0.9,0.9],thick=3,linestyle=1,color='orange',/overplot)
;  t = text(0.71,0.19,'Idealized model',color='orange',font_size=11,/norm)
;  p = plot([4.8,5.1],[0.95,0.95],thick=3,linestyle=2,color='orange',/overplot)
;  t = text(0.71,0.15,'w+20%',color='orange',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_1b.png'

;  Extended Data Figure 2

;  p = plot(indgen(2),/nodata,yrange=[14,36],xrange=[1.03,0],ytitle='Altitude (km)',xtitle='Normalized N$_2$O',font_size=11,title='Normalized N$_2$O Profiles',margin=[0.09,0.09,0.04,0.08], $
;    dimensions=[550,500])
;  p = plot(alt_n2o_norm_avg_tr_percent_oms[*,3],alt_grid,thick=3,color='red',/overplot)
;  p = plot(alt_n2o_norm_avg_percent_oms[*,6],alt_grid,thick=1,color='sky blue',/overplot)
;  p = plot(alt_n2o_norm_avg_percent_oms[*,3],alt_grid,thick=3,color='blue',/overplot)
;  p = plot(alt_n2o_norm_avg_percent_oms[*,0],alt_grid,thick=1,color='sky blue',/overplot)
;  p = plot(alt_n2o_norm_vx_oms[*,3],alt_grid,thick=3,color='purple',/overplot)
;  p = plot(alt_n2o_norm_vx_oms[*,6],alt_grid,thick=1,color='orchid',/overplot)
;  p = plot(alt_n2o_norm_vx_oms[*,0],alt_grid,thick=1,color='orchid',/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['n2o_norm_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['n2o_norm_q_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  p = errorplot(n2o_norm_alt_grid_fl_avg[*,0],alt_grid2,n2o_norm_alt_grid_fl_avg[*,1],replicate(0,nz2),symbol='s',/sym_filled,thick=3,color='green',linestyle=2,errorbar_capsize=0,/overplot)
;  t = text(0.13,0.85,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.13,0.81,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.25,0.81,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.35,0.81,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.13,0.76,'Symbols = Flask balloon',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_2a.png'
;
;  p = plot(indgen(2),/nodata,yrange=[14,36],xrange=[1.05,0.1],ytitle='Altitude (km)',xtitle='Normalized CH$_4$',font_size=11,title='Normalized CH$_4$ Profiles',margin=[0.09,0.09,0.04,0.08], $
;    dimensions=[550,500])
;  p = plot(alt_ch4_norm_avg_tr_percent_oms[*,3],alt_grid,thick=3,color='red',/overplot)
;  p = plot(alt_ch4_norm_avg_percent_oms[*,3],alt_grid,thick=3,color='blue',/overplot)
;  p = plot(alt_ch4_norm_avg_percent_oms[*,6],alt_grid,thick=1,color='sky blue',/overplot)
;  p = plot(alt_ch4_norm_avg_percent_oms[*,0],alt_grid,thick=1,color='sky blue',/overplot)
;  p = plot(alt_ch4_norm_vx_oms[*,3],alt_grid,thick=3,color='purple',/overplot)
;  p = plot(alt_ch4_norm_vx_oms[*,6],alt_grid,thick=1,color='orchid',/overplot)
;  p = plot(alt_ch4_norm_vx_oms[*,0],alt_grid,thick=1,color='orchid',/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['ch4_norm_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['ch4_norm_q_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  p = errorplot(ch4_norm_alt_grid_fl_avg[*,0],alt_grid2,ch4_norm_alt_grid_fl_avg[*,1],replicate(0,nz2),symbol='s',/sym_filled,thick=3,color='green',linestyle=2,errorbar_capsize=0,/overplot)
;  t = text(0.13,0.85,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.13,0.81,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.25,0.81,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.35,0.81,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.55,0.17,'Symbols = Flask balloon',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_2b.png'
;
;  p = plot(indgen(2),/nodata,yrange=[14,36],xrange=[1.03,0],ytitle='Altitude (km)',xtitle='Normalized CFC-12',font_size=11,title='Normalized CFC-12 Profiles',margin=[0.09,0.09,0.04,0.08], $
;    dimensions=[550,500])
;  p = plot(alt_f12_norm_vx_oms[*,3],alt_grid,thick=3,color='purple',/overplot)
;  p = plot(alt_f12_norm_vx_oms[*,0],alt_grid,thick=1,color='orchid',/overplot)
;  p = plot(alt_f12_norm_vx_oms[*,6],alt_grid,thick=1,color='orchid',/overplot)
;  p = plot(alt_f12_norm_sf6_tr_percent_oms[*,3],alt_grid,thick=3,color='red',/overplot)
;  p = plot(alt_f12_norm_sf6_percent_oms[*,3],alt_grid,thick=3,color='blue',/overplot)
;  p = plot(alt_f12_norm_sf6_percent_oms[*,6],alt_grid,thick=1,color='sky blue',/overplot)
;  p = plot(alt_f12_norm_sf6_percent_oms[*,0],alt_grid,thick=1,color='sky blue',/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['f12_norm_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['f12_norm_q_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  p = errorplot(f12_norm_alt_grid_fl_avg[*,0],alt_grid2,f12_norm_alt_grid_fl_avg[*,1],replicate(0,nz2),symbol='s',/sym_filled,thick=3,color='green',linestyle=2,errorbar_capsize=0,/overplot)
;  t = text(0.13,0.85,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.13,0.81,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.25,0.81,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.35,0.81,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.13,0.76,'Symbols = Flask balloon',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_2c.png'
;
;  p = plot(indgen(2),/nodata,yrange=[14,36],xrange=[0,7.5],xtitle='Mean Age (years)',ytitle='Altitude (km)',font_size=11,title='Mean Age Profiles',margin=[0.09,0.09,0.04,0.08],dimensions=[550,500])
;  p = plot(alt_age_avg_tr_percent_oms[*,3],alt_grid,thick=3,color='red',/overplot)
;  p = plot(alt_age_vx_oms[*,0],alt_grid,thick=1,color='orchid',/overplot)
;  p = plot(alt_age_vx_oms[*,3],alt_grid,thick=2,color='purple',/overplot)
;  p = plot(alt_age_vx_oms[*,6],alt_grid,thick=1,color='orchid',/overplot)
;  p = plot(alt_age_avg_percent_oms[*,0],alt_grid,thick=1,color='sky blue',/overplot)
;  p = plot(alt_age_avg_percent_oms[*,3],alt_grid,thick=3,color='blue',/overplot)
;  p = plot(alt_age_avg_percent_oms[*,6],alt_grid,thick=1,color='sky blue',/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['sf6_age_q_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['co2_age_q_'+dates[i_flask[i]]],bal['alt_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='tu',sym_size=1,linestyle=6,/overplot)
;  p = errorplot(co2_age_alt_grid_fl_avg[*,0],alt_grid2,co2_age_alt_grid_fl_avg[*,1],replicate(0,nz2),symbol='tu',/sym_filled,sym_size=1.5,thick=3,color='green',linestyle=2,errorbar_capsize=0,/overplot)
;  p = errorplot(sf6_age_alt_grid_fl_avg[*,0],alt_grid2,sf6_age_alt_grid_fl_avg[*,1],replicate(0,nz2),symbol='o',/sym_filled,thick=3,color='green',linestyle=1,errorbar_capsize=0,/overplot)
;  t = text(0.13,0.85,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.13,0.81,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.25,0.81,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.35,0.81,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.65,0.17,'Flask balloon ',color='lime green',font_size=11,/norm)
;  t = text(0.56,0.13,'Circles=SF$_6$, Triangles=CO$_2$',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_2d.png'

;  Extended Data Figure 3

;  p = plot(indgen(2),/nodata,yrange=[1.1,0],xrange=[1.1,0],ytitle='Normalized CH$_4$',xtitle='Normalized N$_2$O',title='CH$_4$ vs. N$_2$O',font_size=11,position=[0.1,0.1,0.97,0.9], $
;    dimensions=[600,500])
;  p = plot([1.1,0],[1.1,0],linestyle=2,/overplot)
;  p = plot(norm_grid,n2o_norm_ch4_norm_ac_sm_coarse,symbol='o',color='green',/sym_filled,/overplot)
;  p = plot(alt_n2o_norm_avg_tr_percent_oms[*,3],alt_ch4_norm_avg_tr_percent_oms[*,3],thick=3,color='red',/overplot)
;  p = plot(alt_n2o_norm_vx_oms[*,3],alt_ch4_norm_vx_oms[*,3],thick=3,color='purple',/overplot)
;  p = plot(alt_n2o_norm_avg_percent_oms[*,3],alt_ch4_norm_avg_percent_oms[*,3],thick=3,color='blue',/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['n2o_norm_'+dates[i_flask[i]]],bal['ch4_norm_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 0, n_flask-1 do p = plot(bal['n2o_norm_q_'+dates[i_flask[i]]],bal['ch4_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  t = text(0.15,0.82,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.15,0.78,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.27,0.78,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.37,0.78,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.15,0.73,'Symbols = Flask balloon',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_3a.png'
;
;  p = plot(indgen(2),/nodata,yrange=[1.01,0.05],xrange=[1.04,-0.03],xtitle='Normalized CFC-12',ytitle='Normalized N$_2$O',title='N$_2$O vs. CFC-12',font_size=11,position=[0.1,0.1,0.97,0.9], $
;    dimensions=[600,500])
;  p = plot([0,1],[0,1],linestyle=2,/overplot)
;  p = plot(alt_f12_norm_sf6_tr_percent_oms[*,3],alt_n2o_norm_avg_tr_percent_oms[*,3],thick=3,color='red',/overplot)
;  p = plot(alt_f12_norm_vx_oms[*,3],alt_n2o_norm_vx_oms[*,3],thick=3,color='purple',/overplot)
;  p = plot(alt_f12_norm_sf6_percent_oms[*,3],alt_n2o_norm_avg_percent_oms[*,3],thick=3,color='blue',/overplot)
;  for i = 2, n_flask-1 do p = plot(bal['f12_norm_'+dates[i_flask[i]]],bal['n2o_norm_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 2, n_flask-1 do p = plot(bal['f12_norm_q_'+dates[i_flask[i]]],bal['n2o_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  t = text(0.15,0.82,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.15,0.78,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.27,0.78,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.37,0.78,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.15,0.73,'Symbols = Flask balloon',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_3b.png'
;
;  p = plot(indgen(2),/nodata,yrange=[1.1,0],xrange=[1.1,0],ytitle='Normalized CH$_4$',xtitle='Normalized CFC-12',title='CH$_4$ vs. CFC-12',font_size=11,position=[0.1,0.1,0.97,0.9], $
;    dimensions=[600,500])
;  p = plot([1.1,0],[1.1,0],linestyle=2,/overplot)
;  p = plot(alt_f12_norm_sf6_tr_percent_oms[*,3],alt_ch4_norm_avg_tr_percent_oms[*,3],thick=3,color='red',/overplot)
;  p = plot(alt_f12_norm_vx_oms[*,3],alt_ch4_norm_vx_oms[*,3],thick=3,color='purple',/overplot)
;  p = plot(alt_f12_norm_sf6_percent_oms[*,3],alt_ch4_norm_avg_percent_oms[*,3],thick=3,color='blue',/overplot)
;  for i = 2, n_flask-1 do p = plot(bal['f12_norm_'+dates[i_flask[i]]],bal['ch4_norm_'+dates[i_flask[i]]],color='lime green',symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  for i = 2, n_flask-1 do p = plot(bal['f12_norm_q_'+dates[i_flask[i]]],bal['ch4_norm_q_'+dates[i_flask[i]]],color='lime green',/sym_filled,symbol='o',sym_size=0.75,linestyle=6,/overplot)
;  t = text(0.15,0.82,'Solid lines = In situ balloon',font_size=11,/norm)
;  t = text(0.15,0.78,'Tropics',color='red',font_size=11,/norm)
;  t = text(0.27,0.78,'Midlat',color='blue',font_size=11,/norm)
;  t = text(0.37,0.78,'Vortex',color='purple',font_size=11,/norm)
;  t = text(0.15,0.73,'Symbols = Flask balloon',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_3c.png'
  
;  Extended Data Figure 4

;  p = plot(indgen(2),/nodata,xrange=[1970,2030],yrange=[14.5,36],ytitle='Altitude (km)',xtitle='Years',title='Sampling Altitude',font_size=11,margin=[0.08,0.1,0.04,0.08],dimensions=[700,400])
;  p = plot([1970,2030],[24,24],linestyle=2,/overplot)
;  for t = 0, nt-1 do if finite(max_alt_b_yrs[t]) and samp_type_b_yrs[t] eq 0 then p = plot(replicate(years[t],2),[15,max_alt_b_yrs[t]],color='lime green',thick=2,/overplot)
;  for t = 0, nt-1 do if finite(max_alt_b_yrs[t]) and samp_type_b_yrs[t] eq 1 then p = plot(replicate(years[t],2),[15,max_alt_b_yrs[t]],color='purple',thick=2,/overplot)
;  for t = 0, nt-1 do if finite(max_alt_b_yrs[t]) and samp_type_b_yrs[t] eq 2 then p = plot(replicate(years[t],2),[15,max_alt_b_yrs[t]],color='magenta',thick=2,/overplot)
;  for t = 0, nt_n2o-1 do if finite(max_alts_a[t]) then p = plot(replicate(n2o_years[t],2),[15,max_alts_a[t]],color='blue',thick=2,/overplot)
;  t = text(0.93,0.83,'Flask',color='lime green',font_size=11,alignment=1,/norm)
;  t = text(0.93,0.79,'In situ balloon',color='purple',font_size=11,alignment=1,/norm)
;  t = text(0.93,0.75,'In situ aircraft',color='blue',font_size=11,alignment=1,/norm)
;  t = text(0.93,0.71,'AirCore',color='magenta',font_size=11,alignment=1,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_4a.png'
;
;  p = plot(indgen(2),/nodata,xrange=[1970,2030],yrange=[3,7],ytitle='Mean Age (years)',xtitle='Years',title='Mean Age in 24-35km Altitude Range',font_size=11, $
;    margin=[0.08,0.1,0.04,0.08],dimensions=[700,400])
;  ti = where(samp_type_b_yrs eq 0 or samp_type_b_yrs eq 3,nti)
;  p = errorplot(years[ti],bins_all_combo_fl_age_midstrat_nh_years[ti,0],replicate(0,nti),bins_all_combo_fl_age_midstrat_nh_years[ti,1],color='lime green',symbol='o',linestyle=6,/sym_filled, $
;    errorbar_capsize=0,/overplot)
;  si = where(samp_type_b_yrs eq 2,nsi)
;  p = errorplot(years[si],bins_all_combo_nofl_age_midstrat_nh_years[si,0],replicate(0,nsi),bins_all_combo_nofl_age_midstrat_nh_years[si,1],color='magenta',symbol='s',/sym_filled,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  si = where(samp_type_b_yrs eq 1,nsi)
;  p = errorplot(years[si],bins_all_combo_nofl_age_midstrat_nh_years[si,0],replicate(0,nsi),bins_all_combo_nofl_age_midstrat_nh_years[si,1],color='purple',symbol='s',/sym_filled,linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  p = plot([1995,2026],age_nofl_midstrat_linfit[0] + age_nofl_midstrat_linfit[1]*[1995,2026],thick=2,linestyle=2,color='purple',/overplot)
;  p = plot([1975,2026],age_fl_midstrat_linfit[0] + age_fl_midstrat_linfit[1]*[1975,2026],thick=2,linestyle=2,color='green',/overplot)
;  t = text(0.13,0.82,'In situ balloon',color='purple',font_size=11,/norm)
;  t = text(0.13,0.78,'AirCore',color='magenta',font_size=11,/norm)
;  t = text(0.13,0.74,'Flask',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_4b.png'

;  Extended Data Figure 5

;  p = plot(indgen(2),/nodata,yrange=[13.5,35],xrange=[-0.35,0.4],ytitle='Altitude (km)',xtitle='Years/Decade',title='NH Midlatitude Mean Age Trends',font_size=11, $
;    margin=[0.13,0.1,0.08,0.08],dimensions=[450,500])
;  p = plot([0,0],[14,35],linestyle=2,/overplot)
;  p = errorplot(10.*age_nh_trends_coarse[*,1],alt_grid2,10.*age_nh_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1.5, $
;    /sym_filled,color='sky blue',linestyle=6,errorbar_capsize=0,/overplot)
;  p = errorplot(10.*age_nh_lat_adj_trends_coarse[*,1],alt_grid2,10.*age_nh_lat_adj_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1.5, $
;    /sym_filled,color='blue',linestyle=6,errorbar_capsize=0,/overplot)
;  p = errorplot(10.*age_fl_alt_trends[*,1],alt_grid2,10.*age_fl_alt_trends_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1,color='lime green',linestyle=6, $
;    errorbar_capsize=0,/overplot)
;  p = errorplot([-0.07,0.09,0.2,0.22],findgen(4)*5+17.5,[0.037,0.075,0.075,0.082],replicate(!values.f_nan,4),symbol='D',sym_size=1.5,color='black',linestyle=6,errorbar_capsize=0,/overplot)
;  p = errorplot([0.083,0.083],[28.5,28.5],[0.156,0.156],replicate(!values.f_nan,2),symbol='s',sym_size=1.5,color='purple',linestyle=6,errorbar_capsize=0,/overplot)
;  t = text(0.17,0.86,'In situ w/o flasks lat adj',color='blue',font_size=11,/norm)
;  t = text(0.17,0.82,'In situ w/o flasks',color='sky blue',font_size=11,/norm)
;  t = text(0.17,0.78,'In situ w/flasks',color='lime green',font_size=11,/norm)
;  t = text(0.17,0.74,'Ray et al., 2014',color='black',font_size=11,/norm)
;  t = text(0.17,0.7,'Fritsch et al., 2020',color='purple',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_5a.png'
;
;  p = plot(indgen(2),/nodata,/xlog,xrange=[0.7,1e6],yrange=[13.5,35],ytitle='Altitude (km)',xtitle='Number of measurements',title='NH Midlatitude number of measurements',font_size=11, $
;    margin=[0.13,0.1,0.08,0.08],dimensions=[450,500])
;  p = plot(tot_samp,alt_grid2,color='red',thick=2,/overplot)
;  p = plot(tot_samp_aircraft,alt_grid2,symbol='s',/sym_filled,color='blue',linestyle=6,/overplot)
;  p = plot(nsamp_oms,alt_grid2,symbol='s',/sym_filled,color='purple',linestyle=6,/overplot)
;  p = plot(tot_samp_aircore,alt_grid2,symbol='s',/sym_filled,color='magenta',linestyle=6,/overplot)
;  p = plot(nsamp_flask,alt_grid2,symbol='s',/sym_filled,color='lime green',linestyle=6,/overplot)
;  t = text(0.67,0.82,'Aircraft in situ',color='blue',font_size=11,/norm)
;  t = text(0.67,0.78,'Balloon in situ',color='purple',font_size=11,/norm)
;  t = text(0.67,0.74,'AirCore',color='magenta',font_size=11,/norm)
;  t = text(0.67,0.7,'Flask',color='lime green',font_size=11,/norm)
;  t = text(0.67,0.66,'Total',color='red',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_5b.png'
;
;  Extended Data Figure 6
;
;  p = plot(indgen(2),/nodata,yrange=[13.5,31],xrange=[-0.028,0.018],ytitle='Altitude (km)',xtitle='Normalized Tracer Fraction/Decade',title='NH CH$_4$ Trend Profiles', $
;    font_size=11,margin=[0.12,0.09,0.04,0.08],dimensions=[450,500])
;  p = plot([0,0],[15,36],linestyle=2,/overplot)
;  poly = polygon([10*ch4_norm_trends_d2_minmax[7:30,0],reverse(10*ch4_norm_trends_d2_minmax[7:30,1])],[alt_c[7:30],reverse(alt_c[7:30])],/data,/fill_background, $
;    fill_color='violet',fill_transparency=70,color='pink')
;  poly = polygon([10*ch4_norm_trends_d1_minmax[*,0],reverse(10*ch4_norm_trends_d1_minmax[*,1])],[alt_c,reverse(alt_c)],/data,/fill_background,fill_color='lime green', $
;    fill_transparency=70,color='light green')
;  p = plot(10.*ch4_norm_trends_d1_avg,alt_c,color='dark green',thick=2,/overplot)
;  p = plot(10.*ch4_norm_trends_d2_avg,alt_c,color='medium orchid',thick=2,/overplot)
;  p = errorplot(10.*ace_ch4_nh_trends[1,*],ace_alts,10.*ace_ch4_nh_trends_chi,replicate(!values.f_nan,nza),symbol='s',sym_size=1.25,color='black',errorbar_capsize=0, $
;    linestyle=6,/overplot)
;  p = errorplot(10.*ch4_nh_trends_coarse[*,1],alt_grid2,10.*ch4_nh_trends_coarse_sigma[*,1],replicate(!values.f_nan,nz2),symbol='o',sym_size=1.5, $
;    /sym_filled,color='blue',linestyle=6,errorbar_capsize=0,/overplot)
;  s = symbol(0.19,0.86,'o',sym_size=1.25,sym_color='blue',/sym_filled,/norm)
;  t = text(0.22,0.85,'In situ',color='blue',font_size=10,/norm)
;  s = symbol(0.19,0.82,'s',sym_size=1.25,sym_color='black',/norm)
;  t = text(0.22,0.81,'ACE-FTS',color='black',font_size=10,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_6a.png'
;
;  p = plot(indgen(2),/nodata,yrange=[1.01,0.4],xrange=[-0.25,0.11],ytitle='Normalized CH$_4$',xtitle='Years/Decade',font_size=11,title='Mean Age Trends', $
;    axis_style=1,margin=[0.13,0.09,0.12,0.08],dimensions=[450,500])
;  p = plot([0,0],[1,0.4],linestyle=2,/overplot)
;  poly = polygon([10*aoa_on_ch4_trends_d2_minmax[*,0],reverse(10*aoa_on_ch4_trends_d2_minmax[*,1])],[n2o_norm_grid_model,reverse(n2o_norm_grid_model)],/data, $
;    /fill_background,fill_color='violet',fill_transparency=70,color='pink')
;  poly = polygon([10*aoa_on_ch4_trends_d1_minmax[*,0],reverse(10*aoa_on_ch4_trends_d1_minmax[*,1])],[n2o_norm_grid_model,reverse(n2o_norm_grid_model)],/data, $
;    /fill_background,fill_color='lime green',fill_transparency=70,color='light green')
;  p = errorplot(10.*age_on_ch4_trends[5:-2,1],norm_grid[5:-2],10.*age_on_ch4_trends_sigma[5:-2,1],replicate(!values.f_nan,nn-6),symbol='o',sym_size=1.5, $
;    /sym_filled,color='blue',errorbar_capsize=0,/overplot)
;  p = plot(10.*aoa_on_ch4_trends_d1_avg,n2o_norm_grid_model,color='dark green',thick=2,/overplot)
;  p = plot(10.*aoa_on_ch4_trends_d2_avg,n2o_norm_grid_model,color='medium orchid',thick=2,/overplot)
;  xaxis = axis('X',location='top',tickfont_size=0)
;  yaxis = axis('Y',location='right',coord_transform=[0,1],minor=0,tickfont_size=11,tickname=['36','31','27','23','21','19','15'],title='Approximate Altitude (km)')
;  t = text(0.31,0.86,'CCMI-2022',font_size=10,alignment=0.5,/norm)
;  p = plot([-0.19,-0.17],[0.465,0.465],thick=2,color='dark green',/overplot)
;  poly = polygon([-0.19,-0.17,-0.17,-0.19],[0.455,0.455,0.475,0.475],/data,/fill_background,fill_color='lime green',fill_transparency=70,color='light green')
;  t = text(0.34,0.82,'refD1',color='dark green',font_size=10,alignment=0.5,/norm)
;  p = plot([-0.19,-0.17],[0.495,0.495],thick=2,color='medium orchid',/overplot)
;  poly = polygon([-0.19,-0.17,-0.17,-0.19],[0.485,0.485,0.505,0.505],/data,/fill_background,fill_color='violet',fill_transparency=70,color='pink')
;  t = text(0.34,0.78,'refD2',color='medium orchid',font_size=10,alignment=0.5,/norm)
;  p.save,dir+'Plots/Age Trend Paper/ED_Figure_6b.png'

;  Supplement Figure 2

;  p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;    yticklen=0.03,dimensions=[700,450],title='Mean Ages 1990s')
;  for z = 0, nz-1 do begin
;    ii = where(finite(age_combo_grid_early[0,0:-3,z]) and age_combo_grid_early[2,0:-3,z] ge 3 and tp_alt_grid[0:-3]-1.5 le alt_grid[z],nii)
;    if nii gt 0 then begin
;      vcols = age_combo_grid_early[0,ii,z]*46.
;      test = where(vcols ge 240.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;      p = plot(lat_grid_a[ii]-0.6,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;      p = plot(lat_grid_a[ii]+0.6,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;    endif
;    ii = where(finite(age_combo_grid_early[0,*,z]) and age_combo_grid_early[2,*,z] ge 1 and age_combo_grid_early[2,*,z] le 2 and tp_alt_grid-1.5 le alt_grid[z],nii)
;    if nii gt 0 then begin
;      vcols = age_combo_grid_early[0,ii,z]*46.
;      test = where(vcols ge 240.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]-0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]+0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;    endif
;    ii = where(finite(age_combo_grid_early[0,-2:-1,z]) and age_combo_grid_early[2,-2:-1,z] ge 1,nii) + nyg-2
;    if nii gt 0 then begin
;      vcols = age_combo_grid_early[0,ii,z]*46.
;      test = where(vcols ge 240.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]-0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]+0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;    endif
;    if alt_grid[z] le 21 then begin
;      ii = where(finite(age_combo_grid_early[0,-2:-1,z]) and age_combo_grid_early[2,-2:-1,z] ge 3,nii) + nyg-2
;      if nii gt 0 then begin
;        vcols = age_combo_grid_early[0,ii,z]*46.
;        test = where(vcols ge 240.,ntest)
;        if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;        test = where(vcols lt 0.,ntest)
;        if ntest gt 0 then vcols[test] = 0.
;        p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;        p = plot(lat_grid_a[ii]-0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;        p = plot(lat_grid_a[ii]+0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;      endif
;    endif
;  endfor
;  p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;  cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname1,tickvalues=findgen(12)*23.,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='years')
;  p.save,dir+'Plots/Age Trend Paper/Supp_Figure_2a.png'
;
;  p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xminor=1,yrange=[9,25],position=[0.1,0.1,0.89,0.9],xtitle='Equivalent Latitude',ytitle='Altitude (km)',font_size=11,xticklen=0.04, $
;    yticklen=0.03,dimensions=[700,450],title='Normalized N$_2$O 1990s')
;  for z = 0, nz-1 do begin
;    ii = where(finite(n2o_norm_grid_early[0,0:-3,z]) and n2o_norm_grid_early[2,0:-3,z] ge 3 and tp_alt_grid[0:-3]-1.5 le alt_grid[z],nii)
;    if nii gt 0 then begin
;      vcols = (n2o_norm_grid_early[0,ii,z]-0.4)*416.
;      test = where(vcols ge 240.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;      p = plot(lat_grid_a[ii]-0.6,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;      p = plot(lat_grid_a[ii]+0.6,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;    endif
;    ii = where(finite(n2o_norm_grid_early[0,0:-2,z]) and n2o_norm_grid_early[2,0:-2,z] ge 1 and n2o_norm_grid_early[2,0:-2,z] le 2 and tp_alt_grid[0:-2]-1.5 le alt_grid[z],nii)
;    if nii gt 0 then begin
;      vcols = (n2o_norm_grid_early[0,ii,z]-0.4)*416.
;      test = where(vcols ge 240.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]-0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]+0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;    endif
;    ii = where(finite(n2o_norm_grid_early[0,-2,z]) and n2o_norm_grid_early[2,-2,z] ge 1,nii) + nyg-2
;    if nii gt 0 then begin
;      vcols = (n2o_norm_grid_early[0,ii,z]-0.4)*416.
;      test = where(vcols ge 240.,ntest)
;      if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;      test = where(vcols lt 0.,ntest)
;      if ntest gt 0 then vcols[test] = 0.
;      p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]-0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;      p = plot(lat_grid_a[ii]+0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=0.8,/overplot)
;    endif
;    if alt_grid[z] le 21 then begin
;      ii = where(finite(n2o_norm_grid_early[0,-2,z]) and n2o_norm_grid_early[2,-2,z] ge 3,nii) + nyg-2
;      if nii gt 0 then begin
;        vcols = (n2o_norm_grid_early[0,ii,z]-0.4)*416.
;        test = where(vcols ge 240.,ntest)
;        if ntest gt 0 then vcols[where(vcols ge 240.)] = 240.
;        test = where(vcols lt 0.,ntest)
;        if ntest gt 0 then vcols[test] = 0.
;        p = plot(lat_grid_a[ii],replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;        p = plot(lat_grid_a[ii]-0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;        p = plot(lat_grid_a[ii]+0.5,replicate(alt_grid[z],nii),symbol='s',/sym_filled,linestyle=6,rgb_table=39,vert_colors=vcols,sym_size=1.6,/overplot)
;      endif
;    endif
;  endfor
;  p = plot(lat_tp,tpause_alt,thick=2,color='magenta',/overplot)
;  cb = COLORBAR(POSITION=[0.905,0.2,0.925,0.75],tickname=tickname,tickvalues=findgen(7)*41.6,RGB_TABLE=rgb,BORDER=1, ORIENTATION=1,/textpos,/tickdir,font_size=10,title='Normalized N$_2$O')
;  p.save,dir+'Plots/Age Trend Paper/Supp_Figure_2b.png'
;
;  z = 28
;  
;  p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xtitle='Latitude',yrange=[0,6],ytitle='Mean Age (years)',xticklen=0.02,yticklen=0.02,title='Mean age at 20 km  1990s', $
;    margin=[0.08,0.1,0.04,0.08],font_size=11,dimensions=[700,400])
;  for m = 0, nnm-1 do p = plot(lat_c,age_c_1990s_zobs[*,z,m],color='lime green',thick=1,/overplot)
;  gd = where(nseas_old ge 3)
;  p = plot(lat_grid_a[0:-2],age_combo_grid_early[0,0:-2,z],color='blue',thick=3,/overplot)
;  t = text(0.12,0.28,'In situ',color='blue',font_size=11,/norm)
;  t = text(0.12,0.23,'CCMI-2022',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/Supp_Figure_2c.png'
;
;  p = plot(indgen(2),/nodata,xrange=[-90,90],xtickinterval=20,xtitle='Latitude',yrange=[0.4,1],ytitle='Normalized N$_2$O',xticklen=0.02,yticklen=0.02,title='N$_2$O at 20 km  1990s', $
;    margin=[0.09,0.1,0.03,0.08],font_size=11,dimensions=[700,400])
;  for m = 3, 14 do p = plot(lat_c,n2o_norm_c_1990s_zobs[*,z,m],color='lime green',thick=1,/overplot)
;  for m = 16, nnm-1 do p = plot(lat_c,n2o_norm_c_1990s_zobs[*,z,m],color='lime green',thick=1,/overplot)
;  p = plot(lat_grid_a[0:-2],n2o_norm_grid_early[0,0:-2,z],color='blue',thick=3,/overplot)
;  t = text(0.12,0.82,'In situ',color='blue',font_size=11,/norm)
;  t = text(0.12,0.77,'CCMI-2022',color='lime green',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/Supp_Figure_2d.png'

;  Supplement Figure 4

;  p = plot(indgen(2),/nodata,yrange=[1.03,0],xrange=[0,7],ytitle='Normalized N$_2$O',xtitle='Mean Age (years)',title='Mean Age vs. N$_2$O',font_size=11)
;  p = plot(bins_a0_n2o_norm_age_nh_coarse,norm_grid,symbol='s',color='blue',sym_size=1,/sym_filled,thick=2,/overplot)
;  p = plot(bins_a2_n2o_norm_age_nh_trend,norm_grid,symbol='s',color='red',sym_size=1,/sym_filled,thick=2,/overplot)
;  p = plot(n2o_age_a,n2o_norm_grid,thick=2,color='lime green',/overplot)
;  p = plot(n2o_age_engel,n2o_norm_grid,thick=2,color='magenta',/overplot)
;  t = text(0.18,0.8,'1990s this study',color='blue',font_size=11,/norm)
;  t = text(0.18,0.76,'2020s this study',color='red',font_size=11,/norm)
;  t = text(0.18,0.72,'1990s (Andrews et al., 2001)',color='lime green',font_size=11,/norm)
;  t = text(0.18,0.68,'1990s (Engel et al., 2002)',color='magenta',font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/Supp_Figure_4.png'


;  Supplement Figure 5

;i = 29  ; Example date index
;
;p = plot(indgen(2),/nodata,xrange=[1988.5,1995.5],xtickinterval=1,xminor=3,yrange=[350.2,360],ytitle='CO$_2$ (ppmv)',xtitle='Year',font_size=11, $
;  title='Surface and ASHOE  '+flight_track[i].date+' ER-2 CO$_2$',dimensions=[700,400])
;p = plot(years_sm+2./12.,co2_sm,thick=2,color='lime green',/overplot)
;p = plot(years_sm+2./12.,co2_sm_smooth,thick=2,color='green',/overplot)
;t = text(0.16,0.8,'MLO-SMO Avg w/2 mon delay',color='lime green',font_size=11,/norm)
;t = text(0.16,0.76,'MLO-SMO Avg smooth',color='green',font_size=11,/norm)
;t = text(0.16,0.72,'MLO-SMO Avg const growth rate',color='black',font_size=11,/norm)
;t = text(0.16,0.68,'Adj CO$_2$ and lag age',color='red',font_size=11,/norm)
;t = text(0.16,0.64,'Archive Adj CO$_2$ and age',color='blue',font_size=11,/norm)
;t = text(0.16,0.6,'Flight data',color='purple',font_size=11,/norm)
;p = plot(co2_steady_growth_time,co2_steady_growth,thick=2,color='black',/overplot)
;p = plot(flight_track[i].date_float-*flight_track[i].co2_age,*flight_track[i].co2 - (1.725-*flight_track[i].ch4/1e3),symbol='o',color='blue',/overplot)
;p = plot(flight_track[i].date_float-*flight_track[i].co2_lag_age,*flight_track[i].co2_adj,symbol='o',color='red',linestyle=6,/overplot)
;p = plot(replicate(flight_track[i].date_float,n_elements(*flight_track[i].co2_adj)),*flight_track[i].co2_adj,symbol='o',sym_size=0.75,linestyle=6,color='purple', $
;  /overplot)
;p = plot(replicate(flight_track[i].date_float,n_elements(*flight_track[i].co2_adj)),*flight_track[i].co2 - (1.725-*flight_track[i].ch4/1e3),symbol='o', $
;  sym_size=0.75,linestyle=6,color='purple',/overplot)
;p.save,dir+'Plots/Age Trend Paper/Supp_Figure_5.png'

;  Supplement Figure 6

;  diffs_avg = replicate(!values.f_nan,nnorm)
;  for n = 0, nnorm-1 do begin
;    diff_tot = 0. & num = 0
;    for m = 3, 7 do begin
;      if finite(bins_co2_a0.n2o_norm_co2_age_m[n,0,m]) and finite(bins_co2_a0.n2o_norm_age_all_m[n,0,m]) then begin
;        diff_tot += (bins_co2_a0.n2o_norm_co2_age_m[n,0,m]-bins_co2_a0.n2o_norm_age_all_m[n,0,m])
;        num += 1
;      endif
;    endfor
;    if num gt 0 then diffs_avg[n] = diff_tot / num
;  endfor
;  gd = where(finite(diffs_avg))
;  diffs_avg[gd] = smooth(diffs_avg[gd],7,/edge_truncate)
;
;  p = plot(indgen(2),/nodata,yrange=[1,0.1],xrange=[-1.4,0.7],ytitle='Normalized N$_2$O',xtitle='Mean Age Difference(years)',title='Mean Age Difference vs. N$_2$O', $
;    font_size=11,dimensions=[450,550],margin=[0.13,0.08,0.03,0.08])
;  p = plot([0,0],[1,0],linestyle=2,/overplot)
;  for m = 3, 7 do p = errorplot(bins_co2_a0.n2o_norm_co2_age_m[*,0,m]-bins_co2_a0.n2o_norm_age_all_m[*,0,m],n2o_norm_grid,bins_co2_a0.n2o_norm_age_all_m[*,1,m], $
;    replicate(0,nnorm),symbol='o',sym_size=0.75,/sym_filled,linestyle=6,color=colors[m],errorbar_capsize=0,/overplot)
;  p = plot(diffs_avg,n2o_norm_grid,thick=5,color='blue',/overplot)
;  t = text(0.7,0.8,'Average',color='blue',font_size=11,/norm)
;  t = text(0.7,0.76,'1992-3 ('+mission_names[4]+')',color=colors[4],font_size=11,/norm)
;  t = text(0.7,0.72,'1994   ('+mission_names[3]+')',color=colors[3],font_size=11,/norm)
;  t = text(0.7,0.68,'1995-6 ('+mission_names[5]+')',color=colors[5],font_size=11,/norm)
;  t = text(0.7,0.64,'1997   ('+mission_names[6]+')',color=colors[6],font_size=11,/norm)
;  t = text(0.7,0.6,'1999-',color=colors[7],font_size=11,/norm)
;  t = text(0.7,0.565,'2000   ('+mission_names[7]+')',color=colors[7],font_size=11,/norm)
;  p.save,dir+'Plots/Age Trend Paper/Supp_Figure_6.png'
  

;  Write the data plotted in each main-text and Extended Data figure to NetCDF
;  (read by the plotting routines in code/figures/).
@export_figure_data_main
@export_figure_data_supp_ats

end
