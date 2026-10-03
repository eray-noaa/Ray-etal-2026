;+
; vars = MA_NC_READ_FULL(file)
;
; Reads a NetCDF file into the ORDEREDHASH form used by MA_NC_WRITE: each
; entry is a HASH with 'data', 'dims' and the variable's attributes, so the
; variables can be written again (e.g. copied into another figure file).
;-
function ma_nc_read_full, file

  compile_opt idl2

  id = ncdf_open(file)
  info = ncdf_inquire(id)
  out = orderedhash()
  for i = 0, info.nvars-1 do begin
    vi = ncdf_varinq(id, i)
    ncdf_varget, id, i, v
    h = hash('data', v)
    if vi.ndims gt 0 then begin
      dn = strarr(vi.ndims)
      for j = 0, vi.ndims-1 do begin
        ncdf_diminq, id, vi.dim[j], name, len
        dn[j] = name
      endfor
      h['dims'] = dn
    endif
    for j = 0, vi.natts-1 do begin
      an = ncdf_attname(id, i, j)
      if an eq '_FillValue' then continue
      ncdf_attget, id, i, an, av
      if size(av, /type) eq 1 then av = string(av)
      h[an] = av
    endfor
    out[vi.name] = h
  endfor
  ncdf_close, id
  return, out

end
