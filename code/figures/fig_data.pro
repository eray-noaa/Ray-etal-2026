;+
; d = FIG_DATA(name)
;
; Reads <output>/figure_data/<name>.nc (written by the analysis code) into a hash.
;-
function fig_data, name

  compile_opt idl2

  return, ma_nc_read(ma_output_dir() + 'figure_data' + path_sep() + name + '.nc')

end
