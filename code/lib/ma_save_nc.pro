;+
; MA_SAVE_NC, file, names [, DESCRIPTION=description]
;
; Saves variables from the calling routine to a NetCDF file (the NetCDF
; equivalent of SAVE). NAMES is a string array (or comma-separated string)
; of variable names. Numeric and string arrays and scalars are supported;
; each array dimension is named <variable>_d<i>. Read the file back into a
; routine with MA_RESTORE_NC, or with any NetCDF reader.
;-
pro ma_save_nc, file, names, DESCRIPTION=description

  compile_opt idl2

  if n_elements(names) eq 1 then names = strtrim(strsplit(names, ',', /extract), 2)
  vars = orderedhash()
  foreach n, strlowcase(names) do begin
    v = scope_varfetch(n, level=-1)
    t = size(v, /type)
    if t eq 8 || t eq 10 || t eq 11 then message, n + ': structures, pointers and objects cannot be written with MA_SAVE_NC'
    nd = size(v, /n_dimensions)
    if nd eq 0 then vars[n] = hash('data', v) $
    else vars[n] = hash('data', v, 'dims', n + '_d' + strtrim(indgen(nd),2))
  endforeach
  g = hash('source', 'Intermediate file of the Ray et al. (2026) analysis code (read with ma_restore_nc)')
  if n_elements(description) gt 0 then g['title'] = description
  ma_nc_write, file, vars, global=g

end
