;+
; FIG_SAVE, p, name
;
; Saves graphic P as <output>/figures/<name>.pdf and .png and closes it.
;-
pro fig_save, p, name

  compile_opt idl2

  odir = ma_output_dir() + 'figures' + path_sep()
  file_mkdir, odir
  p.save, odir + name + '.png', resolution=150
  p.save, odir + name + '.pdf', /vector
  p.close
  print, 'Saved ', odir + name + '.pdf'

end
