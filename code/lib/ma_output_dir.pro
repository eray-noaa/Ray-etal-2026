;+
; Returns the directory (with trailing separator) where the pipeline writes its
; NetCDF outputs. Set the environment variable MEAN_AGE_OUTPUT_DIR to change it.
; Defaults to the output/ folder at the root of this repository.
;-
function ma_output_dir

  compile_opt idl2

  dir = getenv('MEAN_AGE_OUTPUT_DIR')
  if dir eq '' then begin
    here = file_dirname(routine_filepath('ma_output_dir', /is_function))
    dir = file_dirname(file_dirname(here)) + path_sep() + 'output'
  endif
  if strmid(dir, strlen(dir)-1) ne path_sep() then dir = dir + path_sep()
  return, dir

end
