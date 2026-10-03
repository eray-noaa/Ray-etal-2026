;+
; v = MA_VAR(data, dims [, units [, long_name]], _EXTRA=attrs)
;
; Convenience constructor for one entry of the VARS hash passed to MA_NC_WRITE.
; DIMS is a string array of dimension names in IDL order (omit or pass '' for
; a scalar). Extra keywords become attributes, e.g. DESCRIPTION='...'.
;-
function ma_var, data, dims, units, long_name, _EXTRA=attrs

  compile_opt idl2

  v = hash('data', data)
  if n_elements(dims) gt 0 && dims[0] ne '' then v['dims'] = dims
  if n_elements(units) gt 0 && units ne '' then v['units'] = units
  if n_elements(long_name) gt 0 then v['long_name'] = long_name
  if n_elements(attrs) gt 0 then begin
    tn = strlowcase(tag_names(attrs))
    for i = 0, n_elements(tn)-1 do v[tn[i]] = attrs.(i)
  endif
  return, v

end
