;+
;  Step 4: balloon profiles.
;
;  Quality-controls the cryo-flask, OMS in situ, WAS and AirCore balloon profiles,
;  bins them by altitude and normalized N2O, CH4 and CFC-12, averages the 1970s-2000s
;  flask profiles, and makes the balloon time series on the common time grid.
;
;  Reads:  Balloon/Engel_profiles_v2025.sav (all balloon profiles with mean ages),
;          Balloon/OMS_tracer_profile_avgs.sav, Mean_age_relationships.sav, Aircore_n2o_ch4_v2024.sav
;  Writes: <output>/intermediate/balloon_profiles.nc, flask_profiles.nc;
;          <output>/figure_data/ed_fig2_3.nc
;
;  Generated from the original analysis programs (code/original/) by keeping only the code
;  that contributes to the paper figures; results are identical to the originals.
;-
pro step4_balloon_profiles

dir = mean_age_data_dir()

restore,dir+'Balloon/Engel_profiles_v2025.sav'

n_70s = 3
i_70s = indgen(n_70s)
n_80s = 7
i_80s = indgen(n_80s)+3
i_90s = [indgen(5)+10,16,17,indgen(4)+20]
i_00s = [indgen(5)+24,31,indgen(4)+33]
i_10s = [37,64]
i_20s = [162]
i_flask = [i_70s,i_80s,i_90s,i_00s,i_10s,i_20s]
i_was = [33,34]
i_oms = [15,18,19,29,30,32]
i_aircore = [indgen(26)+38,indgen(97)+65,indgen(84)+163]
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
  if ~bal.haskey('f12_norm_'+dates[t]) then bal['f12_norm_'+dates[t]] = replicate(!values.f_nan, $
    n_elements(bal['ch4_norm_'+dates[t]]))
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
gg = replicate(!values.f_nan,np,nz2,5)
hh = replicate(!values.f_nan,np,nn,5)
co2_age_ch4_bin = hh
sf6_age_ch4_bin = hh
sf6_age_n2o_bin = hh
co2_age_n2o_bin = hh
age_n2o_bin = hh
age_ch4_bin = hh
co2_age_q_alt_grid = gg
sf6_age_q_alt_grid = gg
co2_age_q_ch4_bin = hh
sf6_age_q_ch4_bin = hh
co2_age_q_n2o_bin = hh
sf6_age_q_n2o_bin = hh
age_combo_n2o_bin = hh
age_combo_n2o_bin_nofl = hh
age_combo_n2o_bin_fl = hh
age_q_combo_n2o_bin_fl = hh
age_combo_ch4_bin = hh
age_combo_ch4_bin_nofl = hh
age_q_combo_ch4_bin_fl = hh
age_combo_ch4_bin_fl = hh
age_alt_grid = gg
age_combo_alt_grid = gg
age_combo_alt_grid_nofl = gg
age_combo_alt_grid_fl = gg
n2o_norm_alt_grid = gg
ch4_norm_alt_grid = gg
n2o_norm_alt_grid_nofl = gg
n2o_norm_alt_grid_fl = gg
ch4_norm_alt_grid_nofl = gg
ch4_norm_alt_grid_fl = gg
samp_type_b = intarr(np)
vv = replicate(!values.f_nan,np)
ss = replicate(!values.f_nan,np,nz3,5)
max_alt_b = vv
co2_age_q_alt_grid3 = ss
sf6_age_q_alt_grid3 = ss
age_alt_grid3 = ss
age_combo_alt_grid_nofl3 = ss
n2o_ch4_norm_alt_grid = gg
n2o_from_ch4 = intarr(np)
co2_age_alt_grid_fl = gg
sf6_age_alt_grid_fl = gg
f12_norm_alt_grid = gg
f12_norm_alt_grid_fl = gg
co2_age_q_f12_bin = hh
sf6_age_q_f12_bin = hh
age_f12_bin = hh
age_q_combo_f12_bin_fl = hh
n2o_ch4_norm_alt_grid_nofl = gg
age_combo_n2o_bin_nofl_noch4n2o = hh
age_sf6_n2o_bin_nofl = hh
n2o_on_age_q_bin = hh
for t = 0, np-1 do begin
  fi = where(i_flask eq t,nfi)
  if nfi gt 0 then begin
    fl = 1
  endif else begin
    fl = 0
  endelse
  if dates[t] eq '20040929' then fl = 0
  alts = bal['alt_'+dates[t]]
  max_alt_b[t] = max(alts)
  npts = n_elements(alts)
  c_age = bal['co2_age_'+dates[t]]
  s_age = bal['sf6_age_'+dates[t]]
  c_age_q = bal['co2_age_q_'+dates[t]]
  s_age_q = bal['sf6_age_q_'+dates[t]]
  c_age_err = 0.5 * (bal['co2_age_range_'+dates[t],1,*] - bal['co2_age_range_'+dates[t],0,*])
  if bal.haskey('sf6_age_range_'+dates[t]) then begin
    s_age_err = 0.5 * (bal['sf6_age_range_'+dates[t],1,*] - bal['sf6_age_range_'+dates[t],0,*])
  endif else begin
    s_age_err = replicate(!values.f_nan,npts)
  endelse
  c_age_weights = 1. / c_age_err^2
  s_age_weights = 1. / s_age_err^2
  samp_type_b[t] = bal['samp_type_'+dates[t]]
  if bal.haskey('both_age_'+dates[t]) then begin
    b_age = bal['both_age_'+dates[t]]
  endif else begin
    if bal.haskey('opt_age_'+dates[t]) then begin
      b_age = bal['opt_age_'+dates[t]]
    endif else begin
      b_age = replicate(!values.f_nan,npts)
    endelse
  endelse
  f_norm = bal['f12_norm_q_'+dates[t]]
  c_norm = bal['ch4_norm_q_'+dates[t]]
  if bal.haskey('n2o_ch4_norm_'+dates[t]) then begin
    n2o_ch4_norm = bal('n2o_ch4_norm_'+dates[t])
  endif else begin
    n2o_ch4_norm = replicate(!values.f_nan,npts)
  endelse
  gd = where(~finite(n2o_ch4_norm),ngd)
  if ngd eq npts then begin
    n2o_from_ch4[t] = 0
  endif else begin
    n2o_from_ch4[t] = 1
  endelse
  n_norm = bal['n2o_norm_'+dates[t]]
  n2o_b = bal['n2o_'+dates[t]]
  for z = 0, nz2-1 do begin
    zi = where(alts ge alt_grid2[z]-dz2/2. and alts lt alt_grid2[z]+dz2/2.,nzi)
    if nzi gt 0 then begin
      gd = where(finite(c_age[zi]),ngd)
      if ngd gt 0 then begin
        stats = moment(c_age[zi],sd=sd,/nan)
      endif
      gd = where(finite(s_age[zi]),ngd)
      if ngd gt 0 then begin
        stats = moment(s_age[zi],sd=sd,/nan)
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
      if finite(age_alt_grid[t,z,0]) and ~finite(age_alt_grid[t,z,1]) or age_alt_grid[t,z,1] lt 0.25 then age_alt_grid[t,z, $
        1] = 0.25
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
      if ch4_norm_alt_grid[t,z,0] gt 0 and ch4_norm_alt_grid[t,z,0] lt 1.05 then ch4_norm_alt_grid_nofl[t,z, $
        *] = ch4_norm_alt_grid[t,z,*]
    endif else begin
      if finite(sf6_age_q_alt_grid[t,z,0]) then age_combo_alt_grid_fl[t,z,*] = sf6_age_q_alt_grid[t,z,*]
      if finite(sf6_age_q_alt_grid[t,z,0]) then sf6_age_alt_grid_fl[t,z,*] = sf6_age_q_alt_grid[t,z,*]
      if finite(co2_age_q_alt_grid[t,z,0]) then age_combo_alt_grid_fl[t,z,*] = co2_age_q_alt_grid[t,z,*]
      if finite(co2_age_q_alt_grid[t,z,0]) then co2_age_alt_grid_fl[t,z,*] = co2_age_q_alt_grid[t,z,*]
      if finite(age_alt_grid[t,z,0]) then age_combo_alt_grid_fl[t,z,*] = age_alt_grid[t,z,*]
      if n2o_norm_alt_grid[t,z,0] gt 0 then n2o_norm_alt_grid_fl[t,z,*] = n2o_norm_alt_grid[t,z,*]
      if ch4_norm_alt_grid[t,z,0] gt 0 and ch4_norm_alt_grid[t,z,0] lt 1.05 then ch4_norm_alt_grid_fl[t,z, $
        *] = ch4_norm_alt_grid[t,z,*]
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
      if finite(age_alt_grid3[t,z,0]) and ~finite(age_alt_grid3[t,z,1]) or age_alt_grid3[t,z,1] lt 0.25 then age_alt_grid3[t,z, $
        1] = 0.25
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
        stats = moment(c_age[gd],sdev=sdev,/nan)
      endif
      gd2 = where(finite(s_age[gd]),ngd2)
      if ngd2 gt 0 then begin
        stats = moment(s_age[gd],sdev=sdev,/nan)
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
    if fl eq 0 then begin
      if finite(sf6_age_q_ch4_bin[t,i,0]) then age_combo_ch4_bin_nofl[t,i,*] = sf6_age_q_ch4_bin[t,i,*]
      if finite(co2_age_q_ch4_bin[t,i,0]) then age_combo_ch4_bin_nofl[t,i,*] = co2_age_q_ch4_bin[t,i,*]
      if finite(age_ch4_bin[t,i,0]) then age_combo_ch4_bin_nofl[t,i,*] = age_ch4_bin[t,i,*]
      if ~finite(age_combo_ch4_bin_nofl[t,i,1]) then age_combo_ch4_bin_nofl[t,i,1] = 0.1
      if finite(sf6_age_q_n2o_bin[t,i,0]) then age_combo_n2o_bin_nofl[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
      if finite(co2_age_q_n2o_bin[t,i,0]) then age_combo_n2o_bin_nofl[t,i,*] = co2_age_q_n2o_bin[t,i,*]
      if finite(age_n2o_bin[t,i,0]) then age_combo_n2o_bin_nofl[t,i,*] = age_n2o_bin[t,i,*]
      if finite(sf6_age_q_n2o_bin[t,i,0]) then age_sf6_n2o_bin_nofl[t,i,*] = sf6_age_q_n2o_bin[t,i,*]
      if finite(sf6_age_q_n2o_bin[t,i,0]) and n2o_from_ch4[t] eq 0 then age_combo_n2o_bin_nofl_noch4n2o[t,i, $
        *] = sf6_age_q_n2o_bin[t,i,*]
      if finite(co2_age_q_n2o_bin[t,i,0]) and n2o_from_ch4[t] eq 0 then age_combo_n2o_bin_nofl_noch4n2o[t,i, $
        *] = co2_age_q_n2o_bin[t,i,*]
      if finite(age_n2o_bin[t,i,0]) and n2o_from_ch4[t] eq 0 then age_combo_n2o_bin_nofl_noch4n2o[t,i,*] = age_n2o_bin[t,i,*]
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

;  Make flask average profiles from 1970s-2000s.
aa = replicate(!values.f_nan,nz2,5)
n2o_norm_alt_grid_fl_avg = aa
ch4_norm_alt_grid_fl_avg = aa
f12_norm_alt_grid_fl_avg = aa
co2_age_alt_grid_fl_avg = aa
sf6_age_alt_grid_fl_avg = aa
bb = replicate(!values.f_nan,nn,5)
age_q_combo_n2o_bin_fl_avg = bb
age_q_combo_ch4_bin_fl_avg = bb
age_q_combo_f12_bin_fl_avg = bb
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
  wmean = total(1. / co2_age_alt_grid_fl[ti[gd],z,1]^2 * co2_age_alt_grid_fl[ti[gd],z, $
    0]) / total(1. / co2_age_alt_grid_fl[ti[gd],z,1]^2)
  npts = total(co2_age_alt_grid_fl[ti[gd],z,4])
  stats = moment(co2_age_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  co2_age_alt_grid_fl_avg[z,*] = [wmean,err,stats[0],sdev,npts]
  gd = where(finite(sf6_age_alt_grid_fl[ti,z,0]),ngd)
  err = 1. / sqrt(total(sf6_age_alt_grid_fl[ti[gd],z,1]^2))
  wmean = total(1. / sf6_age_alt_grid_fl[ti[gd],z,1]^2 * sf6_age_alt_grid_fl[ti[gd],z, $
    0]) / total(1. / sf6_age_alt_grid_fl[ti[gd],z,1]^2)
  npts = total(sf6_age_alt_grid_fl[ti[gd],z,4])
  stats = moment(sf6_age_alt_grid_fl[ti,z,0],sdev=sdev,/nan)
  sf6_age_alt_grid_fl_avg[z,*] = [wmean,err,stats[0],sdev,npts]
endfor
for i = 0, nn-1 do begin
  stats = moment(co2_age_q_n2o_bin[ti,i,0],sdev=sdev,/nan)
  gd = where(finite(age_q_combo_n2o_bin_fl[ti,i,0]),ngd)
  err = 1. / sqrt(total(age_q_combo_n2o_bin_fl[ti[gd],i,1]^2))
  wmean = total(1. / age_q_combo_n2o_bin_fl[ti[gd],i,1]^2 * age_q_combo_n2o_bin_fl[ti[gd],i, $
    0]) / total(1. / age_q_combo_n2o_bin_fl[ti[gd],i,1]^2)
  npts = total(age_q_combo_n2o_bin_fl[ti[gd],i,4])
  stats = moment(age_q_combo_n2o_bin_fl[ti,i,0],sdev=sdev,/nan)
  age_q_combo_n2o_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  stats = moment(co2_age_q_ch4_bin[ti,i,0],sdev=sdev,/nan)
  gd = where(finite(age_q_combo_ch4_bin_fl[ti,i,0]),ngd)
  err = 1. / sqrt(total(age_q_combo_ch4_bin_fl[ti[gd],i,1]^2))
  wmean = total(1. / age_q_combo_ch4_bin_fl[ti[gd],i,1]^2 * age_q_combo_ch4_bin_fl[ti[gd],i, $
    0]) / total(1. / age_q_combo_ch4_bin_fl[ti[gd],i,1]^2)
  npts = total(age_q_combo_ch4_bin_fl[ti[gd],i,4])
  stats = moment(age_q_combo_ch4_bin_fl[ti,i,0],sdev=sdev,/nan)
  age_q_combo_ch4_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
  stats = moment(co2_age_q_f12_bin[ti,i,0],sdev=sdev,/nan)
  gd = where(finite(age_q_combo_f12_bin_fl[ti,i,0]),ngd)
  err = 1. / sqrt(total(age_q_combo_f12_bin_fl[ti[gd],i,1]^2))
  wmean = total(1. / age_q_combo_f12_bin_fl[ti[gd],i,1]^2 * age_q_combo_f12_bin_fl[ti[gd],i, $
    0]) / total(1. / age_q_combo_f12_bin_fl[ti[gd],i,1]^2)
  npts = total(age_q_combo_f12_bin_fl[ti[gd],i,4])
  stats = moment(age_q_combo_f12_bin_fl[ti,i,0],sdev=sdev,/nan)
  age_q_combo_f12_bin_fl_avg[i,*] = [wmean,err,stats[0],sdev,npts]
endfor

nsamp_flask = fltarr(nz2)
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
mm = replicate(!values.f_nan,nt,nz2,5)
age_combo_alt_grid_nofl_yrs = mm
age_combo_alt_grid_fl_yrs = mm
rr = replicate(!values.f_nan,nt)
samp_type_b_yrs = rr
max_alt_b_yrs = rr
n2o_from_ch4_yrs = rr
n2o_norm_alt_grid_yrs = mm
for t = 0, nt-1 do begin
  gd = where(yr_frac ge years[t]-dt/2. and yr_frac lt years[t]+dt/2.,ni1)
  if ni1 ge 1 then begin
    samp_type_b_yrs[t] = samp_type_b[gd[0]]
    max_alt_b_yrs[t] = max(max_alt_b[gd])
    n2o_from_ch4_yrs[t] = n2o_from_ch4[gd[0]]
    for i = 0, nn-1 do begin
      stats = moment(co2_age_ch4_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(co2_age_ch4_bin[gd,i,1],sdev=sdev,/nan)
      stats = moment(sf6_age_ch4_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(sf6_age_ch4_bin[gd,i,1],sdev=sdev,/nan)
      stats = moment(co2_age_q_n2o_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(co2_age_q_n2o_bin[gd,i,1],sdev=sdev,/nan)
      stats = moment(sf6_age_q_n2o_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(sf6_age_q_n2o_bin[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_ch4_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_ch4_bin[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_ch4_bin_nofl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_ch4_bin_nofl[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_ch4_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_ch4_bin_fl[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_q_combo_ch4_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_q_combo_ch4_bin_fl[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_n2o_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_n2o_bin_nofl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin_nofl[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_sf6_n2o_bin_nofl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_sf6_n2o_bin_nofl[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_n2o_bin_nofl_noch4n2o[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin_nofl_noch4n2o[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_combo_n2o_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_n2o_bin_fl[gd,i,1],sdev=sdev,/nan)
      stats = moment(age_q_combo_n2o_bin_fl[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(age_q_combo_n2o_bin_fl[gd,i,1],sdev=sdev,/nan)
      stats = moment(n2o_on_age_q_bin[gd,i,0],sdev=sdev,/nan)
      stats2 = moment(n2o_on_age_q_bin[gd,i,1],sdev=sdev,/nan)
    endfor
    for z = 0, nz2-1 do begin
      stats = moment(age_combo_alt_grid[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid[gd,z,1],sdev=sdev,/nan)
      gd2 = where(finite(age_combo_alt_grid_nofl[gd,z,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_alt_grid_nofl[gd[gd2],z,1]^2))
      wmean = total(1. / age_combo_alt_grid_nofl[gd[gd2],z,1]^2 * age_combo_alt_grid_nofl[gd[gd2],z, $
        0]) / total(1. / age_combo_alt_grid_nofl[gd[gd2],z,1]^2)
      npts = total(age_combo_alt_grid_nofl[gd[gd2],z,4])
      stats = moment(age_combo_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      age_combo_alt_grid_nofl_yrs[t,z,*] = [wmean,err,stats[0],stats2[0],npts]
      if finite(age_combo_alt_grid_nofl_yrs[t,z,0]) and ~finite(age_combo_alt_grid_nofl_yrs[t,z, $
        3]) or age_combo_alt_grid_nofl_yrs[t,z,3] lt 0.25 then age_combo_alt_grid_nofl_yrs[t,z,3] = 0.25
      gd2 = where(finite(age_combo_alt_grid_fl[gd,z,0]),ngd2)
      err = 1. / sqrt(total(1./age_combo_alt_grid_fl[gd[gd2],z,1]^2))
      wmean = total(1. / age_combo_alt_grid_fl[gd[gd2],z,1]^2 * age_combo_alt_grid_fl[gd[gd2],z, $
        0]) / total(1. / age_combo_alt_grid_fl[gd[gd2],z,1]^2)
      npts = total(age_combo_alt_grid_fl[gd[gd2],z,4])
      stats = moment(age_combo_alt_grid_fl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid_fl[gd,z,1],sdev=sdev,/nan)
      age_combo_alt_grid_fl_yrs[t,z,*] = [wmean,err,stats[0],stats2[0],npts]
      if finite(age_combo_alt_grid_fl_yrs[t,z,0]) and ~finite(age_combo_alt_grid_fl_yrs[t,z,3]) or age_combo_alt_grid_fl_yrs[t, $
        z,3] lt 0.25 then age_combo_alt_grid_fl_yrs[t,z,3] = 0.25
      stats = moment(n2o_norm_alt_grid[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_norm_alt_grid[gd,z,1],sdev=sdev,/nan)
      n2o_norm_alt_grid_yrs[t,z,0:1] = [stats[0],stats2[0]]
      stats = moment(n2o_norm_alt_grid_fl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_norm_alt_grid_fl[gd,z,1],sdev=sdev,/nan)
      stats = moment(n2o_norm_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_norm_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      stats = moment(n2o_ch4_norm_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(n2o_ch4_norm_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
      stats = moment(ch4_norm_alt_grid_fl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(ch4_norm_alt_grid_fl[gd,z,1],sdev=sdev,/nan)
      stats = moment(ch4_norm_alt_grid_nofl[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(ch4_norm_alt_grid_nofl[gd,z,1],sdev=sdev,/nan)
    endfor
    for z = 0, nz3-1 do begin
      stats = moment(age_combo_alt_grid_nofl3[gd,z,0],sdev=sdev,/nan)
      stats2 = moment(age_combo_alt_grid_nofl3[gd,z,1],sdev=sdev,/nan)
    endfor
  endif
endfor

restore,dir+'Balloon/OMS_tracer_profile_avgs.sav'
restore,dir+'Balloon/Mean_age_relationships.sav'

inter_dir__ = ma_output_dir() + 'intermediate' + path_sep()
fig_dir__ = ma_output_dir() + 'figure_data' + path_sep()
stat5__ = ['weighted_mean','uncertainty_of_weighted_mean','mean','standard_deviation','n_measurements']
flask_dates__ = dates[i_flask]
flask__ = ma_flatten_profiles(bal, flask_dates__, $
  ['alt','n2o_norm','n2o_norm_q','ch4_norm','ch4_norm_q','f12_norm','f12_norm_q','co2_age','co2_age_q','sf6_age','sf6_age_q'], $
  prefix='flask_')
flask__['flask_alt','units'] = 'km'
foreach k__, ['flask_co2_age','flask_co2_age_q','flask_sf6_age','flask_sf6_age_q'] do flask__[k__,'units'] = 'years'
flask__['flask_co2_age','long_name'] = 'cryo-flask CO2 mean age'
flask__['flask_sf6_age','long_name'] = 'cryo-flask SF6 mean age'
flask__['flask_co2_age_q','long_name'] = 'cryo-flask CO2 mean age after quality control (NaN = removed)'
flask__['flask_sf6_age_q','long_name'] = 'cryo-flask SF6 mean age after quality control (NaN = removed)'
foreach k__, ['n2o','ch4','f12'] do begin
  flask__['flask_'+k__+'_norm','long_name'] = 'cryo-flask normalized '+k__
  flask__['flask_'+k__+'_norm_q','long_name'] = 'cryo-flask normalized '+k__+' after quality control (NaN = removed)'
endforeach
ma_nc_write, inter_dir__+'flask_profiles.nc', flask__, $
  global=hash('title', 'Cryo-flask balloon profiles used in Fig. 1 and Extended Data Figs 1-3')
; ---------- Extended Data Figures 2 and 3: profiles and tracer correlations
v__ = orderedhash()
v__['alt_grid'] = ma_var(alt_grid, ['alt'], 'km', 'altitude grid of the OMS in situ balloon profiles')
pct__ = ', columns are percentiles; the figure uses column 3 (median) and columns 0 and 6 (range)'
foreach k__, ['n2o_norm','ch4_norm','f12_norm','age'] do begin
  case k__ of
    'n2o_norm': begin
      tr__ = alt_n2o_norm_avg_tr_percent_oms
      ml__ = alt_n2o_norm_avg_percent_oms
      vx__ = alt_n2o_norm_vx_oms
    end
    'ch4_norm': begin
      tr__ = alt_ch4_norm_avg_tr_percent_oms
      ml__ = alt_ch4_norm_avg_percent_oms
      vx__ = alt_ch4_norm_vx_oms
    end
    'f12_norm': begin
      tr__ = alt_f12_norm_sf6_tr_percent_oms
      ml__ = alt_f12_norm_sf6_percent_oms
      vx__ = alt_f12_norm_vx_oms
    end
    'age':      begin
      tr__ = alt_age_avg_tr_percent_oms
      ml__ = alt_age_avg_percent_oms
      vx__ = alt_age_vx_oms
    end
  endcase
  u__ = (k__ eq 'age') ? 'years' : ''
  v__['oms_tropics_'+k__] = ma_var(tr__, ['alt','pct'], u__, 'OMS in situ balloon tropical profile'+pct__)
  v__['oms_midlat_'+k__] = ma_var(ml__, ['alt','pct'], u__, 'OMS in situ balloon midlatitude profile'+pct__)
  v__['oms_vortex_'+k__] = ma_var(vx__, ['alt','pct'], u__, 'OMS in situ balloon vortex profile'+pct__)
endforeach
v__['alt_grid2'] = ma_var(alt_grid2, ['alt2'], 'km', '1 km altitude bins for flask averages')
v__['flask_avg_n2o_norm'] = ma_var(n2o_norm_alt_grid_fl_avg, ['alt2','stat5'], '', $
  'cryo-flask average normalized N2O (columns 0 mean, 1 sd)')
v__['flask_avg_ch4_norm'] = ma_var(ch4_norm_alt_grid_fl_avg, ['alt2','stat5'], '', $
  'cryo-flask average normalized CH4 (columns 0 mean, 1 sd)')
v__['flask_avg_f12_norm'] = ma_var(f12_norm_alt_grid_fl_avg, ['alt2','stat5'], '', $
  'cryo-flask average normalized CFC-12 (columns 0 mean, 1 sd)')
v__['flask_avg_co2_age'] = ma_var(co2_age_alt_grid_fl_avg, ['alt2','stat5'], 'years', 'cryo-flask CO2 mean age', $
  statistics='columns: '+strjoin(stat5__,', '))
v__['flask_avg_sf6_age'] = ma_var(sf6_age_alt_grid_fl_avg, ['alt2','stat5'], 'years', 'cryo-flask SF6 mean age', $
  statistics='columns: '+strjoin(stat5__,', '))
v__['norm_grid'] = ma_var(norm_grid, ['norm'], '', 'normalized N2O bin centres')
v__['aircore_ch4_norm_on_n2o'] = ma_var(n2o_norm_ch4_norm_ac_sm_coarse, ['norm'], '', $
  'AirCore CH4-N2O relationship (Extended Data Fig. 3a green symbols: x = norm_grid, y = this)')
v__ += flask__
ma_nc_write, fig_dir__+'ed_fig2_3.nc', v__, $
  global=hash('title', 'Extended Data Figures 2 and 3: tracer and mean age profiles and correlations')

;  Hand the results needed by later steps to them.
ma_save_nc, ma_output_dir() + 'intermediate' + path_sep() + 'balloon_profiles.nc', $
  ['age_combo_alt_grid_fl_yrs','age_combo_alt_grid_nofl_yrs','age_q_combo_ch4_bin_fl_avg','age_q_combo_f12_bin_fl_avg', $
  'age_q_combo_n2o_bin_fl_avg','alt_grid2','alt_grid3','dn','dz2','dz3','max_alt_b_yrs','mm','n2o_from_ch4_yrs', $
  'n2o_norm_alt_grid_yrs','nn','nnorm','norm_grid','nsamp_flask','nt','nz2','nz3','samp_type_b_yrs','sdev','years','z']
end

