;+
; Returns the top-level data directory (with trailing separator) used by all
; routines in this repository. Set the environment variable MEAN_AGE_DATA_DIR
; to point somewhere else, e.g. setenv,'MEAN_AGE_DATA_DIR=/path/to/Data/'.
; Defaults to the data/ folder at the root of this repository.
;-
function mean_age_data_dir

  dir = getenv('MEAN_AGE_DATA_DIR')
  if dir eq '' then begin
    here = file_dirname(routine_filepath('mean_age_data_dir', /is_function))
    dir = file_dirname(file_dirname(here)) + path_sep() + 'data'
  endif
  if strmid(dir, strlen(dir)-1) ne path_sep() then dir = dir + path_sep()
  return, dir

end
