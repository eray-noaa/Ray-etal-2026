pro calculate_opt_age,year_frac,lat,sf6,co2,ch4,n2o_norm,co2_adj,sf6_opt,age,age_range,ratio,lat_source,ch4_entry,sf6_scale,co2_scale,sf6_bias_opt,lag_age=lag_age

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
sf6_et = aa & co2_e = aa & ch4_e = aa & sf6_e_fine = cc & co2_e_fine = cc & ch4_e_fine = cc & lat_e_fine = cc & frac_extratr_e_fine = cc
conv = {sf6_e_fine2:ee, co2_e_fine2:ee, ch4_e_fine2:ee, lat_e_fine2:ee, frac_extratr_e_fine2:ee}
for r = 0, nr_trop-1 do for y = 0, ny-1 do for a = 0, n_ages_trop-1 do begin
  sf6_et[a,y,r] = interpolate(sf6_entry_trop[a,y,r,*],tii)
  co2_e[a,y,r] = interpolate(co2_entry_trop[a,y,r,*],tii)
  ch4_e[a,y,r] = interpolate(ch4_entry_trop[a,y,r,*],tii)
endfor
for r = 0, nr_trop-1 do for y = 0, ny-1 do begin
  sf6_e_fine[*,y,r] = interpolate(sf6_et[*,y,r],ai)
  co2_e_fine[*,y,r] = interpolate(co2_e[*,y,r],ai)
  ch4_e_fine[*,y,r] = interpolate(ch4_e[*,y,r],ai)
  lat_e_fine[*,y,r] = interpolate(lat_entry_trop[*,y,r],ai)
  frac_extratr_e_fine[*,y,r] = interpolate(frac_extratropics_trop[*,r],ai)
endfor
for y = 0, ny-1 do for a = 0, n_ages_trop_fine-1 do begin
  conv.sf6_e_fine2[a,y,*] = interpolate(sf6_e_fine[a,y,*],rii)
  conv.co2_e_fine2[a,y,*] = interpolate(co2_e_fine[a,y,*],rii)
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

nn = 1e2
dn2o = 1e-2
n2o_norm_grid = findgen(nn)*dn2o+3e-2

sf6_scale_inc = 0.03
sf6_diffs_both_scale = 1.0
co2_scale_inc = 0.02

ns = 41
nsc = 15
sf6_scale_all = dindgen(ns)*sf6_scale_inc-ns/2*sf6_scale_inc
co2_scale_all = dindgen(nsc)*co2_scale_inc-nsc/2*co2_scale_inc
sf6_scale_full = fltarr(ns,nsc) & co2_scale_full = sf6_scale_full
for s = 0, nsc-1 do sf6_scale_full[*,s] = sf6_scale_all
for s = 0, ns-1 do co2_scale_full[s,*] = co2_scale_all


co2_err = 0.05  ;  In ppm
ch4_err = 0.05  ;  In percent based on Weinheimer et al., 1998
ch4_conversion = 0.95   ;  Conversion fraction of CH4 to CO2
ch4_conv_uncert = 0.05  ;  Addition to adjusted CO2 uncertainty due to CH4 conversion uncertainty.

ni = n_elements(co2)

aa = replicate(!values.f_nan,n_ages_trop_fine,ny)
dd = replicate(!values.f_nan,2,n_ages_trop_fine,ny) & ee = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine) & oo = replicate(!values.f_nan,n_ages_trop_fine,ny,nr_fine)
cc = replicate(!values.f_nan,ni) & gg = replicate(!values.f_nan,ns,nsc) & age = cc & ratio = cc & lat_source = cc & lat_source_c = cc & co2_adj = cc & ch4_entry = cc & sf6_scale = cc
ff = replicate(!values.f_nan,2,ny) & pp = replicate(!values.f_nan,n_ages_trop_fine,ny) & qq = replicate(!values.f_nan,2,n_ages_trop_fine,ny) & co2_scale = cc & sf6_bias_opt = cc & sf6_opt = cc
age_range = replicate(!values.f_nan,2,ni)

for i = 0, ni-1 do begin

  if finite(sf6[i]) and finite(co2[i]) and co2[i] gt 0 and finite(n2o_norm[i]) then begin

    print,i,lat[i],sf6[i],co2[i],n2o_norm[i],ch4[i]
    
    range_ratio_trop_co2 = dd & co2_diffs = ee & ch4_loss_full = ee & ages_full = ee & source_lats_full = ee & source_lats_c_full = ee & ratios_full = ee & co2_err_full = oo & co2_adj_full = oo
    ch4_entry_full = ee & diffs_both_min_tot = gg & opt_ratio_diffs_min = gg & opt_lat_source_diffs_min = gg & opt_age_diffs_min = gg & opt_age_diffs_all = gg & opt_ratio_diffs_all = gg
    opt_lat_source_diffs_all = gg & opt_lat_source_c_diffs_all = gg & opt_ch4_loss_diffs_all = gg & opt_sf6_all = gg & opt_sf6_bias_all = gg & opt_co2_adj_all = gg & diffs_both_min_tot_num = gg

    bii = where(year_b eq fix(year_frac))
    sf6_bias_yr = sf6_bias_corr[*,bii]

    yi1 = 0
    yi2 = ny-1
    if lat[i] ge 20 then yi1 = 4
    if lat[i] le -20 then yi2 = ny-5

    if n2o_norm[i] lt 0.5 then begin
      sb1 = 0
      sbb1 = 0
      sb2 = ns-1
      sbb2 = nsc-1
    endif else begin
      n2o_scale = n2o_norm_grid[-1] - 2. * (n2o_norm_grid[-1] - n2o_norm[i])
      sb1 = fix(n2o_scale*(ns/4.))
      sbb1 = fix(n2o_scale*(nsc/3+1.))
      sb2 = ns-1 - sb1
      sbb2 = nsc-1 - sbb1
    endelse

    cutit = 0

    for b = 0, ns-1 do begin

      sf6_adj_s = sf6[i] + sf6_scale_all[b]

      for bb = 3, nsc-4 do begin

        co2_diffs = oo & sf6_diffs = oo & range_ratio_trop_co2 = qq & range_ratio_trop_sf6 = qq & diffs_both_all = oo & diffs_both_min = pp & range_ratio_both = qq & opt_ratio_both = aa
        ages_full = oo & ages_full_sf6 = oo & source_lats_full = oo & ratios_full = oo & opt_lat_source_both = dd & diffs_both_min_lat = dd & opt_ratio_min_diff_lat = dd & opt_age_min_diff_lat = dd
        opt_age_sf6_min_diff_lat = dd & age_range_both = ff & ch4_loss = oo & ch4_loss_min = pp & opt_ch4_loss_min_diff_lat = dd & ch4_loss_full = oo & ages_full_co2 = oo
        co2_adj_full = oo & co2_err_full = oo & sf6_full = oo & sf6_e2 = oo & sf6_both = oo & co2_adj_both = oo & source_lats_c_full = oo & sf6_bias_all = oo & sf6_bias_both = oo

        for y = yi1, yi2 do begin

          for a = 0, n_ages_trop_fine-1 do begin

            sf6_bias = sf6_bias_yr[a] * flight_sf6_gr
            sf6_adj_b = sf6_adj_s + sf6_bias

            ri1 = where(ratios_fine ge grn_mean_ages_trop_fine[a]/8.)
            ri1 = ri1[0]

            if finite(ch4[i]) then begin
              co2_adj_sc = co2[i] + co2_scale_all[bb] - 1e-3 * ch4_conversion * abs(mean(conv.ch4_e_fine2[a,y,ri1:-1]) - ch4[i])
              co2_adj_full[a,y,*] = co2[i] + co2_scale_all[bb] - 1e-3 * ch4_conversion * abs(conv.ch4_e_fine2[a,y,*] - ch4[i])
              co2_err_sc = co2_err + 1e-3 * ch4_conv_uncert * abs(mean(conv.ch4_e_fine2[a,y,ri1:-1]) - (ch4[i]-ch4_err*ch4[i]))
              co2_err_full[a,y,*] = co2_err + 1e-3 * ch4_conv_uncert * abs(conv.ch4_e_fine2[a,y,*] - (ch4[i]-ch4_err*ch4[i]))
            endif else begin
              co2_adj_sc = co2[i] + co2_scale_all[bb] - 1.1 * (1. - n2o_norm[i])
              co2_adj_full[a,y,*] = replicate(co2_adj_sc,nr_fine)
              co2_err_sc = co2_err
              co2_err_full[a,y,*] = replicate(co2_err_sc,nr_fine)
            endelse
            chk_co2 = where(conv.co2_e_fine2[a,y,ri1:-1] ge co2_adj_sc-co2_err_sc and conv.co2_e_fine2[a,y,ri1:-1] le co2_adj_sc+co2_err_sc,nchk_co2)
            if nchk_co2 gt 0 then begin
              range_ratio_trop_co2[*,a,y] = [ratios_fine[chk_co2[0]+ri1],ratios_fine[chk_co2[-1]+ri1]]
              co2_diffs[a,y,chk_co2+ri1] = abs(co2_adj_sc - conv.co2_e_fine2[a,y,chk_co2+ri1])
              if finite(ch4[i]) then begin
                ch4_loss[a,y,chk_co2+ri1] = ch4_conversion * abs(conv.ch4_e_fine2[a,y,chk_co2+ri1] - ch4[i])  ;  Assume 95% conversion of CH4 to CO2.
                ch4_entry_full[a,y,chk_co2+ri1] = conv.ch4_e_fine2[a,y,chk_co2+ri1]
              endif else begin
                ch4_loss[a,y,chk_co2+ri1] = 1.1e3 * (1. - n2o_norm[i])
                ch4_entry_full[a,y,chk_co2+ri1] = !values.f_nan
              endelse
              ages_full_co2[a,y,chk_co2+ri1] = replicate(grn_mean_ages_trop_fine[a],nchk_co2)
            endif
            sf6_min = reform(sf6_adj_b-0.03*sf6[i])
            sf6_max = reform(sf6_adj_b+0.03*sf6[i])
            sf6_min = sf6_min[0]
            sf6_max = sf6_max[0]
            chk_sf6 = where(conv.sf6_e_fine2[a,y,ri1:-1] ge sf6_min and conv.sf6_e_fine2[a,y,ri1:-1] le sf6_max,nchk_sf6)
            if nchk_sf6 gt 0 then begin
              range_ratio_trop_sf6[*,a,y] = [ratios_fine[chk_sf6[0]+ri1],ratios_fine[chk_sf6[-1]+ri1]]
              sf6_diffs[a,y,chk_sf6+ri1] = replicate(sf6_adj_b,nchk_sf6) - reform(conv.sf6_e_fine2[a,y,chk_sf6+ri1])
              sf6_diffs[a,y,chk_sf6+ri1] = abs(sf6_diffs[a,y,chk_sf6+ri1])
              ages_full_sf6[a,y,chk_sf6+ri1] = replicate(grn_mean_ages_trop_fine[a],nchk_sf6)
              sf6_full[a,y,chk_sf6+ri1] = replicate(sf6_adj_b,nchk_sf6)
              sf6_e2[a,y,chk_sf6+ri1] = reform(conv.sf6_e_fine2[a,y,chk_sf6+ri1])
              sf6_bias_all[a,y,chk_sf6+ri1] = replicate(sf6_bias,nchk_sf6)
            endif
          endfor

          for a = 0, n_ages_trop_fine-1 do begin
            ri1 = where(ratios_fine ge grn_mean_ages_trop_fine[a]/8.)
            ri1 = ri1[0]
            if finite(range_ratio_trop_co2[0,a,y]) then begin
              if finite(range_ratio_trop_sf6[0,a,y]) then begin
                si = where(finite(co2_diffs[a,y,*]) and finite(sf6_diffs[a,y,*]),nsi)
                if nsi gt 0 then begin
                  diffs_both_all[a,y,si] = 0.5 * (co2_diffs[a,y,si]/co2_err_full[a,y,si] + sf6_diffs[a,y,si]/0.03)
                  si_min = where(diffs_both_all[a,y,si] eq min(diffs_both_all[a,y,si]))
                  diffs_both_min[a,y] = diffs_both_all[a,y,si[si_min[0]]]
                  range_ratio_both[*,a,y] = [ratios_fine[si[0]],ratios_fine[si[-1]]]
                  opt_ratio_both[a,y] = median(ratios_fine[si])
                  ratios_full[a,y,si] = ratios_fine[si]
                  ages_full[a,y,si] = replicate(grn_mean_ages_trop_fine[a],nsi)
                  source_lats_full[a,y,si] = replicate(lat_grid_entry[y],nsi)
                  source_lats_c_full[a,y,si] = conv.lat_e_fine2[a,y,si]
                  ch4_loss_full[a,y,si] = ch4_loss[a,y,si]
                  ch4_loss_min[a,y] = ch4_loss[a,y,si[si_min[0]]]
                  sf6_both[a,y,si] = sf6_full[a,y,si]
                  sf6_bias_both[a,y,si] = sf6_bias_all[a,y,si]
                  co2_adj_both[a,y,si] = co2_adj_full[a,y,si]
                endif
              endif
            endif
          endfor

          chk = where(finite(opt_ratio_both[*,y]),nchk)
          if nchk gt 0 then begin
            age_range_both[*,y] = [grn_mean_ages_trop_fine[chk[0]],grn_mean_ages_trop_fine[chk[-1]]]
            opt_lat_source_both[y] = 1
            dmin = min(diffs_both_min[chk,y],mi)
            diffs_both_min_lat[y] = dmin
            opt_age_min_diff_lat[y] = grn_mean_ages_trop_fine[chk[mi]]
            opt_age_sf6_min_diff_lat[y] = grn_mean_ages_trop_fine[chk[mi]]
            opt_ch4_loss_min_diff_lat[y] = ch4_loss_min[chk[mi],y]
          endif

        endfor

        chk = where(finite(opt_lat_source_both),nchk)
        if nchk gt 0 then begin
          dmin = min(diffs_both_min_lat[chk],mi)
          diffs_both_min_tot[b,bb] = dmin
          opt_lat_source_diffs_min[b,bb] = lat_grid_entry[chk[mi]]
          opt_age_diffs_min[b,bb] = opt_age_min_diff_lat[chk[mi]]

          ;  Check on size of diffs_both_min_tot for no CO2 scaling to see if the SF6 scaling can be stopped.
          if b ge ns/2+3 then begin
            chk = where(diffs_both_min_tot[b-6:b,bb] le 0.5,nchk)
            if nchk eq 7 then cutit = 1
          endif

          chk = where(finite(diffs_both_all) and diffs_both_all le 0.15,nchk)  ;  Number of low difference scaling.
          if nchk gt 0 then diffs_both_min_tot_num[b,bb] = nchk

          ;  Calculate the optimal ages and ratios based on scaled diffs.
          zz = where(finite(ages_full),nzz)
          diffs_rev = max(diffs_both_all[zz])-diffs_both_all[zz]
          tot = total(diffs_rev)
          diffs_norm = diffs_rev / tot
          opt_age_diffs_all[b,bb] = total(ages_full[zz] * diffs_norm)
          opt_ratio_diffs_all[b,bb] = total(ratios_full[zz] * diffs_norm)
          opt_lat_source_diffs_all[b,bb] = total(source_lats_full[zz] * diffs_norm)
          opt_lat_source_c_diffs_all[b,bb] = total(source_lats_c_full[zz] * diffs_norm)
          opt_ch4_loss_diffs_all[b,bb] = total(ch4_loss_full[zz] * diffs_norm)
          opt_sf6_all[b,bb] = total(sf6_both[zz] * diffs_norm)
          opt_sf6_bias_all[b,bb] = total(sf6_bias_both[zz] * diffs_norm)
          opt_co2_adj_all[b,bb] = total(co2_adj_both[zz] * diffs_norm)
        endif

      endfor

      if cutit then b = ns-1

    endfor

    diffs = diffs_both_min_tot
    counts = diffs_both_min_tot_num
    ages = opt_age_diffs_all

    ;  Set the conditions for the scalings that are included in the optimization.
    if n2o_norm[i] lt 0.8 then zz = where(finite(ages) and counts ge 5 and diffs le 0.1,nzz)
    if n2o_norm[i] ge 0.8 then begin
      if diffs[ns/2,nsc/2] le 0.15 or diffs[ns/2-1,nsc/2] le 0.15 or diffs[ns/2+1,nsc/2] le 0.15 then begin
        diffs[0:16,*] = !values.f_nan
        diffs[-17:-1,*] = !values.f_nan
        zz = where(finite(ages) and counts ge 5 and diffs le 0.15,nzz)
      endif else zz = where(finite(ages) and counts ge 5 and diffs le 0.15,nzz)
    endif
    print,nzz

    diffs_rev = (max(diffs[zz])-diffs[zz])^2
    tot = total(diffs_rev)
    diffs_norm = diffs_rev / tot
    age_range[*,i] = [min(ages[zz]),max(ages[zz])]
    age[i] = total(ages[zz] * diffs_norm)
    ratio[i] = total(opt_ratio_diffs_all[zz] * diffs_norm)
    lat_source[i] = total(opt_lat_source_diffs_all[zz] * diffs_norm)
    lat_source_c[i] = total(opt_lat_source_c_diffs_all[zz] * diffs_norm)
    sf6_scale[i] = total(sf6_scale_full[zz] * diffs_norm)
    co2_scale[i] = total(co2_scale_full[zz] * diffs_norm)
    ch4_loss[i] = total(opt_ch4_loss_diffs_all[zz] * diffs_norm)
    sf6_opt[i] = total(opt_sf6_all[zz] * diffs_norm)
    sf6_bias_opt[i] = total(opt_sf6_bias_all[zz] * diffs_norm)
    co2_adj[i] = total(opt_co2_adj_all[zz] * diffs_norm)
    ch4_entry[i] = total(ch4_entry_full[zz] * diffs_norm)

  endif

endfor

end