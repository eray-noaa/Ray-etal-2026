;+
; MA_NC_WRITE, file, vars, GLOBAL=global
;
; Writes variables to a compressed NetCDF-4 file.
;
; VARS is an ORDEREDHASH keyed by variable name. Each value is a HASH with
;   'data'      the array (numeric or string)
;   'dims'      string array of dimension names, in IDL order (fastest varying
;               first); omit for a scalar
;   any other keys (e.g. 'units', 'long_name', 'description') are written as
;   variable attributes.
; Dimension lengths are taken from the data. A dimension name used by more
; than one variable must have the same length in each.
;
; GLOBAL is an optional HASH of global attributes. The paper DOI and the
; creation time are always added.
;-
pro ma_nc_write, file, vars, GLOBAL=global

  compile_opt idl2

  file_mkdir, file_dirname(file)
  id = ncdf_create(file, /clobber, /netcdf4_format)

  ; Define dimensions.
  dimids = hash()
  foreach v, vars, name do begin
    if ~v.haskey('dims') then continue
    d = v['dims']
    sz = size(v['data'], /dimensions)
    if n_elements(d) eq 1 && n_elements(v['data']) eq 1 then sz = [1L]
    if n_elements(sz) ne n_elements(d) then $
      message, name + ': data has ' + strtrim(n_elements(sz),2) + ' dimensions but ' + $
        strtrim(n_elements(d),2) + ' dimension names were given'
    for i = 0, n_elements(d)-1 do begin
      if dimids.haskey(d[i]) then begin
        ncdf_diminq, id, dimids[d[i]], dname, dlen
        if dlen ne sz[i] then message, name + ': dimension ' + d[i] + ' has length ' + $
          strtrim(sz[i],2) + ', already defined as ' + strtrim(dlen,2)
      endif else dimids[d[i]] = ncdf_dimdef(id, d[i], sz[i])
    endfor
  endforeach

  ; Define variables and attributes.
  varids = hash()
  foreach v, vars, name do begin
    t = size(v['data'], /type)
    dims = []
    if v.haskey('dims') then foreach d, v['dims'] do dims = [dims, dimids[d]]
    types = hash(1,'ubyte', 2,'short', 3,'long', 4,'float', 5,'double', 7,'string', $
      12,'ushort', 13,'ulong', 14,'int64', 15,'uint64')
    if ~types.haskey(t) then message, name + ': cannot write IDL type ' + strtrim(t,2) + ' to NetCDF'
    kw = create_struct(types[t], 1)
    if n_elements(dims) gt 0 && t ne 7 then kw = create_struct(kw, 'gzip', 4)
    vid = (n_elements(dims) gt 0) ? ncdf_vardef(id, name, dims, _extra=kw) : ncdf_vardef(id, name, _extra=kw)
    varids[name] = vid
    if t eq 4 || t eq 5 then ncdf_attput, id, vid, '_FillValue', (t eq 4) ? !values.f_nan : !values.d_nan
    foreach a, v, key do begin
      if key eq 'data' || key eq 'dims' then continue
      ncdf_attput, id, vid, key, a
    endforeach
  endforeach

  ; Global attributes.
  ncdf_attput, id, /global, 'references', 'Ray et al. (2026), Nature Geoscience, https://doi.org/10.1038/s41561-026-02011-3'
  ncdf_attput, id, /global, 'history', 'Created ' + systime() + ' by IDL ' + !version.release
  if n_elements(global) gt 0 then foreach a, global, key do ncdf_attput, id, /global, key, a
  ncdf_control, id, /endef

  foreach v, vars, name do ncdf_varput, id, varids[name], v['data']
  ncdf_close, id

end
