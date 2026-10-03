;+
; result = MA_NC_READ(file [, names])
;
; Reads variables from a NetCDF file written by MA_NC_WRITE into an
; ORDEREDHASH keyed by variable name. NAMES optionally selects variables.
; Attributes are not returned.
;-
function ma_nc_read, file, names

  compile_opt idl2

  if ~file_test(file) then message, 'Missing ' + file + ' (run the step that writes it first)'
  id = ncdf_open(file)
  info = ncdf_inquire(id)
  if n_elements(names) eq 0 then begin
    names = strarr(info.nvars)
    for i = 0, info.nvars-1 do names[i] = (ncdf_varinq(id, i)).name
  endif
  out = orderedhash()
  foreach n, names do begin
    ncdf_varget, id, ncdf_varid(id, n), v
    out[n] = v
  endforeach
  ncdf_close, id
  return, out

end
