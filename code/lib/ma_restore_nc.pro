;+
; MA_RESTORE_NC, file [, names]
;
; Reads variables from a NetCDF file written by MA_SAVE_NC into the calling
; routine (the NetCDF equivalent of RESTORE). NAMES optionally selects variables.
;-
pro ma_restore_nc, file, names

  compile_opt idl2

  d = ma_nc_read(file, names)
  foreach v, d, k do (scope_varfetch(k, level=-1, /enter)) = v

end
