pro calculate_opt_co2_mean_age,year_frac,lat,co2,co2_err,ch4,n2o_norm,co2_adj,age,age_range,ratio,lat_source,ch4_entry

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
co2_e = aa & ch4_e = aa & co2_e_fine = cc & ch4_e_fine = cc & lat_e_fine = cc & frac_extratr_e_fine = cc
conv = {co2_e_fine2:ee, ch4_e_fine2:ee, lat_e_fine2:ee, frac_extratr_e_fine2:ee}
for r = 0, nr_trop-1 do for y = 0, ny-1 do for a = 0, n_ages_trop-1 do begin
  co2_e[a,y,r] = interpolate(co2_entry_trop[a,y,r,*],tii)
  ch4_e[a,y,r] = interpolate(ch4_entry_trop[a,y,r,*],tii)
endfor
for r = 0, nr_trop-1 do for y = 0, ny-1 do begin
  co2_e_fine[*,y,r] = interpolate(co2_e[*,y,r],ai)
  ch4_e_fine[*,y,r] = interpolate(ch4_e[*,y,r],ai)
  lat_e_fine[*,y,r] = interpolate(lat_entry_trop[*,y,r],ai)
  frac_extratr_e_fine[*,y,r] = interpolate(frac_extratropics_trop[*,r],ai)
endfor
for y = 0, ny-1 do for a = 0, n_ages_trop_fine-1 do begin
  conv.co2_e_fine2[a,y,*] = interpolate(co2_e_fine[a,y,*],rii)
  conv.ch4_e_fine2[a,y,*] = interpolate(ch4_e_fine[a,y,*],rii)
  conv.lat_e_fine2[a,y,*] = interpolate(lat_e_fine[a,y,*],rii)
  conv.frac_extratr_e_fine2[a,y,*] = interpolate(frac_extratr_e_fine[a,y,*],rii)
endfor

;co2_err = 0.05  ;  In ppm
ch4_err = 0.05  ;  In percent based on Weinheimer et al., 1998
ch4_conversion = 0.95   ;  Conversion fraction of CH4 to CO2
ch4_conv_uncert = 0.05  ;  Addition to adjusted CO2 uncertainty due to CH4 conversion uncertainty.

ni = n_elements(co2)

dd = replicate(!values.f_nan,2,n_ages_trop_fine,ny) & ee = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine) & oo = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine)
cc = replicate(!values.f_nan,ni) & age = cc & ratio = cc & lat_source = cc & lat_source_c = cc & co2_adj = cc & ch4_entry = cc & age_range = replicate(!values.f_nan,2,ni)

for i = 0, ni-1 do begin

  if finite(co2[i]) and co2[i] gt 0 then begin
  
    range_ratio_trop_co2 = dd & co2_diffs = ee & ch4_loss_full = ee & ages_full = ee & source_lats_full = ee & source_lats_c_full = ee & ratios_full = ee & co2_err_full = oo & co2_adj_full = oo
    ch4_entry_full = ee
  
    yi1 = 0
    yi2 = ny-1
    if lat[i] ge 20 then yi1 = 4
    if lat[i] le -20 then yi2 = ny-5
  
    ri1 = 0
  
    for y = yi1, yi2 do begin
      for a = 0, n_ages_trop_fine-1 do begin
        if finite(ch4[i]) then begin
          co2_adj_sc = co2[i] - 1e-3 * ch4_conversion * abs(mean(conv.ch4_e_fine2[a,y,*]) - ch4[i])
          co2_adj_full[a,y,*] = co2[i] - 1e-3 * ch4_conversion * abs(conv.ch4_e_fine2[a,y,*] - ch4[i])
          co2_err_sc = co2_err[i] + 1e-3 * ch4_conv_uncert * abs(mean(conv.ch4_e_fine2[a,y,*]) - (ch4[i]-ch4_err*ch4[i]))
          co2_err_full[a,y,*] = co2_err[i] + 1e-3 * ch4_conv_uncert * abs(conv.ch4_e_fine2[a,y,*] - (ch4[i]-ch4_err*ch4[i]))
        endif else begin
          co2_adj_sc = co2[i] - 1.1 * (1. - n2o_norm[i])
          co2_adj_full[a,y,*] = replicate(co2_adj_sc,nr_fine)
          co2_err_sc = co2_err[i]
          co2_err_full[a,y,*] = replicate(co2_err_sc,nr_fine)
        endelse
        chk_co2 = where(conv.co2_e_fine2[a,y,*] ge co2_adj_sc-co2_err_sc and conv.co2_e_fine2[a,y,*] le co2_adj_sc+co2_err_sc,nchk_co2)
        if nchk_co2 gt 0 then begin
          range_ratio_trop_co2[*,a,y] = [ratios_fine[chk_co2[0]+ri1],ratios_fine[chk_co2[-1]+ri1]]
          co2_diffs[a,y,chk_co2+ri1] = abs(co2_adj_sc - conv.co2_e_fine2[a,y,chk_co2+ri1])
          if finite(ch4[i]) then begin
            ch4_loss_full[a,y,chk_co2+ri1] = ch4_conversion * abs(conv.ch4_e_fine2[a,y,chk_co2+ri1] - ch4[i])   ;  Assume 95% conversion of CH4 to CO2.
            ch4_entry_full[a,y,chk_co2+ri1] = conv.ch4_e_fine2[a,y,chk_co2+ri1]
          endif else begin
            ch4_loss_full[a,y,chk_co2+ri1] = 1.1e3 * (1. - n2o_norm[i])
            ch4_entry_full[a,y,chk_co2+ri1] = !values.f_nan
          endelse
          ratios_full[a,y,chk_co2+ri1] = ratios_fine[chk_co2+ri1]
          ages_full[a,y,chk_co2+ri1] = replicate(grn_mean_ages_trop_fine[a],nchk_co2)
          source_lats_full[a,y,chk_co2+ri1] = replicate(lat_grid_entry[y],nchk_co2)
          source_lats_c_full[a,y,chk_co2+ri1] = conv.lat_e_fine2[a,y,chk_co2+ri1]
        endif
      endfor
    endfor
  
    zz = where(finite(ages_full),nzz)
    diffs_rev = max(co2_diffs[zz])-co2_diffs[zz]
    tot = total(diffs_rev)
    diffs_norm = diffs_rev / tot
    age_range[*,i] = [min(ages_full[zz]),max(ages_full[zz])]
    age[i] = total(ages_full[zz] * diffs_norm)
    ratio[i] = total(ratios_full[zz] * diffs_norm)
    lat_source[i] = total(source_lats_full[zz] * diffs_norm)
    lat_source_c[i] = total(source_lats_c_full[zz] * diffs_norm)
    co2_adj[i] = total(co2_adj_full[zz] * diffs_norm)
    ch4_entry[i] = total(ch4_entry_full[zz] * diffs_norm)
  endif

endfor

end