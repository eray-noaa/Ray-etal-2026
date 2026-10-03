;+
; vars = MA_FLATTEN_PROFILES(bal, dates, keys [, PREFIX=prefix])
;
; Converts profile data stored in a HASH keyed by '<key>_<date>' (the 'bal'
; hash used by the balloon programs) into a contiguous ragged array suitable
; for NetCDF: all observations are stacked along an 'obs' dimension and
; PROFILE_INDEX gives the profile (0-based index into DATES) of each one.
; Missing keys are filled with NaN. Returns an ORDEREDHASH for MA_NC_WRITE.
;-
function ma_flatten_profiles, bal, dates, keys, PREFIX=prefix

  compile_opt idl2

  if n_elements(prefix) eq 0 then prefix = ''
  np = n_elements(dates)
  nobs = lonarr(np)
  for i = 0, np-1 do nobs[i] = n_elements(bal['alt_'+dates[i]])
  ntot = total(nobs, /integer)
  start = [0L, (total(nobs, /cumulative, /integer))[0:-2]]

  out = orderedhash()
  out[prefix+'profile_date'] = ma_var(string(dates), [prefix+'profile'], '', 'profile launch date (YYYYMMDD)')
  out[prefix+'profile_nobs'] = ma_var(nobs, [prefix+'profile'], '', 'number of observations in each profile')
  pidx = lonarr(ntot)
  for i = 0, np-1 do pidx[start[i]:start[i]+nobs[i]-1] = i
  out[prefix+'profile_index'] = ma_var(pidx, [prefix+'obs'], '', '0-based index of the profile each observation belongs to')
  foreach k, keys do begin
    a = replicate(!values.f_nan, ntot)
    for i = 0, np-1 do if bal.haskey(k+'_'+dates[i]) then begin
      v = float(reform(bal[k+'_'+dates[i]]))
      if n_elements(v) eq nobs[i] then a[start[i]:start[i]+nobs[i]-1] = v
    endif
    out[prefix+k] = ma_var(a, [prefix+'obs'])
  endforeach
  return, out

end
