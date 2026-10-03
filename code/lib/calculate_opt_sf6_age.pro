pro calculate_opt_sf6_age,year_frac,lat,sf6,n2o_norm,age,age_range,ratio,lat_source,sf6_adj_all,ch4_entry

dir = mean_age_data_dir()

restore,dir+'Surface_Trace_Gas/Tracer_entry_time_series.sav'
ny = n_elements(lat_grid_entry)
nr = n_elements(ratios)
nr_trop = n_elements(ratios_trop)
n_ages_trop = n_elements(grn_mean_ages_trop)
nte = n_elements(time_entry)

d_age = grn_mean_ages_trop[1] - grn_mean_ages_trop[0]
d_age_fine = 0.05
n_ages_trop_fine = 2*n_ages_trop
grn_mean_ages_trop_fine = findgen(n_ages_trop_fine)*d_age_fine
ai = interpol(indgen(n_ages_trop),grn_mean_ages_trop,grn_mean_ages_trop_fine)

d_ratios_fine = 0.025
nr_fine = 2.5/d_ratios_fine-1
ratios_fine = findgen(nr_fine)*d_ratios_fine+0.1
rii = interpol(indgen(nr_trop),ratios_trop,ratios_fine)

tii = interpol(indgen(nte),time_entry,year_frac)

;  Interpolate the surface convolutions in time and to a finer grid.
aa = replicate(!values.f_nan,n_ages_trop,ny,nr_trop) & cc = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_trop) & ee = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine)
sf6_et = aa & ch4_e = aa & sf6_e_fine = cc & ch4_e_fine = cc & lat_e_fine = cc & frac_extratr_e_fine = cc
conv = {sf6_e_fine2:ee, ch4_e_fine2:ee, lat_e_fine2:ee, frac_extratr_e_fine2:ee}
for r = 0, nr_trop-1 do for y = 0, ny-1 do for a = 0, n_ages_trop-1 do begin
  sf6_et[a,y,r] = interpolate(sf6_entry_trop[a,y,r,*],tii)
  ch4_e[a,y,r] = interpolate(ch4_entry_trop[a,y,r,*],tii)
endfor
for r = 0, nr_trop-1 do for y = 0, ny-1 do begin
  sf6_e_fine[*,y,r] = interpolate(sf6_et[*,y,r],ai)
  ch4_e_fine[*,y,r] = interpolate(ch4_e[*,y,r],ai)
  lat_e_fine[*,y,r] = interpolate(lat_entry_trop[*,y,r],ai)
  frac_extratr_e_fine[*,y,r] = interpolate(frac_extratropics_trop[*,r],ai)
endfor
for y = 0, ny-1 do for a = 0, n_ages_trop_fine-1 do begin
  conv.sf6_e_fine2[a,y,*] = interpolate(sf6_e_fine[a,y,*],rii)
  conv.ch4_e_fine2[a,y,*] = interpolate(ch4_e_fine[a,y,*],rii)
  conv.lat_e_fine2[a,y,*] = interpolate(lat_e_fine[a,y,*],rii)
  conv.frac_extratr_e_fine2[a,y,*] = interpolate(frac_extratr_e_fine[a,y,*],rii)
endfor


;  Read in Garny SF6 age correction curve.
data = read_ascii(dir+'Models/EMAC/SF6_age_bias_correction_Garny.txt',data_start=5)
year_b = reform(data.field1[0,*])
ft_b = reform(data.field1[1,*])

sf6_bias_corr = fltarr(n_ages_trop_fine,n_elements(year_b)) & sf6_bias_corr_trop = fltarr(n_ages_trop,n_elements(year_b))
for t = 0, n_elements(year_b)-1 do begin
  sf6_bias_corr[*,t] = grn_mean_ages_trop_fine * (ft_b[t] / (5.20991983e+02*exp(-3.80974431e-01*grn_mean_ages_trop_fine)))
  sf6_bias_corr_trop[*,t] = grn_mean_ages_trop * (ft_b[t] / (5.20991983e+02*exp(-3.80974431e-01*grn_mean_ages_trop)))
endfor

sf6_e = reform(sf6_entry_trop[n_ages_trop/2,ny/2,nr/2,*])
sf6_yrsm = lowpass_cfc(sf6_e, BOX=12, EDGE_PFCAST=1)
sf6_gr = 12*sg_smooth(sf6_yrsm, NLEFT=12/2, NRIGHT=12/2-1, DERIV=1, DELTA=1.0, EDGE_PFCAST=1)

tii = interpol(indgen(nte),time_entry,year_frac)
flight_sf6_gr = interpolate(sf6_gr,tii)


aa = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_trop) & bb = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine)
dd = replicate(!values.f_nan,2,n_ages_trop_fine,ny) & ee = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine)

ni = n_elements(sf6)
  
cc = replicate(!values.f_nan,ni) & age = cc & ratio = cc & lat_source = cc & lat_source_c = cc & sf6_adj_all = cc & age_orig = cc & ratio_orig = cc & ch4_entry = cc
age_range = replicate(!values.f_nan,2,ni)

bii = where(year_b eq fix(year_frac))
sf6_bias_yr = sf6_bias_corr[*,bii]

for i = 0, ni-1 do begin

  if finite(sf6[i]) and sf6[i] gt 0 then begin

    range_ratio_trop_sf6 = dd & sf6_diffs = ee & ages_full = ee & source_lats_full = ee & source_lats_c_full = ee & ratios_full = ee & frac_extratr_full = ee & sf6_adj = ee
    ch4_entry_full = ee
  
    yi1 = 0
    yi2 = ny-1
    if lat[i] ge 20 then yi1 = 4
    if lat[i] le -20 then yi2 = ny-5
  
    ri1 = 0
  
    for y = yi1, yi2 do begin
      for a = 0, n_ages_trop_fine-1 do begin
        sf6_bias = sf6_bias_yr[a] * flight_sf6_gr
        sf6_adj_b = sf6[i] + sf6_bias
        sf6_min = reform(sf6_adj_b-0.03*sf6[i])
        sf6_max = reform(sf6_adj_b+0.03*sf6[i])
        sf6_min = sf6_min[0]
        sf6_max = sf6_max[0]
        chk_sf6 = where(conv.sf6_e_fine2[a,y,*] ge sf6_min and conv.sf6_e_fine2[a,y,*] le sf6_max,nchk_sf6)
        if nchk_sf6 gt 0 then begin
          range_ratio_trop_sf6[*,a,y] = [ratios_fine[chk_sf6[0]+ri1],ratios_fine[chk_sf6[-1]+ri1]]
          sf6_diffs[a,y,chk_sf6+ri1] = replicate(sf6_adj_b,nchk_sf6) - reform(conv.sf6_e_fine2[a,y,chk_sf6+ri1])
          sf6_diffs[a,y,chk_sf6+ri1] = abs(sf6_diffs[a,y,chk_sf6+ri1])
          ratios_full[a,y,chk_sf6+ri1] = ratios_fine[chk_sf6+ri1]
          ages_full[a,y,chk_sf6+ri1] = replicate(grn_mean_ages_trop_fine[a],nchk_sf6)
          source_lats_full[a,y,chk_sf6+ri1] = replicate(lat_grid_entry[y],nchk_sf6)
          source_lats_c_full[a,y,chk_sf6+ri1] = conv.lat_e_fine2[a,y,chk_sf6+ri1]
          frac_extratr_full[a,y,chk_sf6+ri1] = conv.frac_extratr_e_fine2[a,y,chk_sf6+ri1]
          sf6_adj[a,y,chk_sf6+ri1] = replicate(sf6_adj_b,nchk_sf6)
        endif
      endfor
    endfor
  
    zz = where(finite(ages_full),nzz)
    diffs_rev = max(sf6_diffs[zz])-sf6_diffs[zz]
    tot = total(diffs_rev)
    diffs_norm = diffs_rev / tot
    age_range[*,i] = [min(ages_full[zz]),max(ages_full[zz])]
    age[i] = total(ages_full[zz] * diffs_norm)
    ratio[i] = total(ratios_full[zz] * diffs_norm)
    lat_source[i] = total(source_lats_full[zz] * diffs_norm)
    lat_source_c[i] = total(source_lats_c_full[zz] * diffs_norm)
    sf6_adj_all[i] = total(sf6_adj[zz] * diffs_norm)
  endif
endfor

;  Repeat calculation for original SF6 mixing ratios without loss adjustment.
for i = 0, ni-1 do begin

  if finite(sf6[i]) and sf6[i] gt 0 then begin

    range_ratio_trop_sf6 = dd & sf6_diffs = ee & ages_full = ee & source_lats_full = ee & source_lats_c_full = ee & ratios_full = ee & frac_extratr_full = ee & ch4_entry_full = ee
  
    yi1 = 0
    yi2 = ny-1
    if lat[i] ge 20 then yi1 = 4
    if lat[i] le -20 then yi2 = ny-5
    for y = yi1, yi2 do begin
      for a = 0, n_ages_trop_fine-1 do begin
        sf6_min = reform(sf6[i]-0.03*sf6[i])
        sf6_max = reform(sf6[i]+0.03*sf6[i])
        sf6_min = sf6_min[0]
        sf6_max = sf6_max[0]
        chk_sf6 = where(conv.sf6_e_fine2[a,y,*] ge sf6_min and conv.sf6_e_fine2[a,y,*] le sf6_max,nchk_sf6)
        if nchk_sf6 gt 0 then begin
          range_ratio_trop_sf6[*,a,y] = [ratios_fine[chk_sf6[0]+ri1],ratios_fine[chk_sf6[-1]+ri1]]
          sf6_diffs[a,y,chk_sf6+ri1] = replicate(sf6[i],nchk_sf6) - reform(conv.sf6_e_fine2[a,y,chk_sf6+ri1])
          sf6_diffs[a,y,chk_sf6+ri1] = abs(sf6_diffs[a,y,chk_sf6+ri1])
          ratios_full[a,y,chk_sf6+ri1] = ratios_fine[chk_sf6+ri1]
          ages_full[a,y,chk_sf6+ri1] = replicate(grn_mean_ages_trop_fine[a],nchk_sf6)
          source_lats_full[a,y,chk_sf6+ri1] = replicate(lat_grid_entry[y],nchk_sf6)
          source_lats_c_full[a,y,chk_sf6+ri1] = conv.lat_e_fine2[a,y,chk_sf6+ri1]
          frac_extratr_full[a,y,chk_sf6+ri1] = conv.frac_extratr_e_fine2[a,y,chk_sf6+ri1]
          ch4_entry_full[a,y,chk_sf6+ri1] = conv.ch4_e_fine2[a,y,chk_sf6+ri1]
        endif
      endfor
    endfor
    zz = where(finite(ages_full),nzz)
    diffs_rev = max(sf6_diffs[zz])-sf6_diffs[zz]
    tot = total(diffs_rev)
    diffs_norm = diffs_rev / tot
    age_orig[i] = total(ages_full[zz] * diffs_norm)
    ratio_orig[i] = total(ratios_full[zz] * diffs_norm)
    ch4_entry[i] = total(ch4_entry_full[zz] * diffs_norm)
  endif
endfor

end