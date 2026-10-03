;+
; FIG_FLASK_POINTS, d, xkey, ykey, _EXTRA=ex
;
; Overplots the cryo-flask balloon observations stored in a figure data hash:
; open symbols for all data and filled symbols for data passing quality
; control (the *_q variables). XKEY and YKEY name variables without the
; 'flask_' prefix, e.g. 'sf6_age' and 'n2o_norm'. SKIP_PROFILES drops the
; first N flask profiles (Extended Data Fig. 3b,c start at the third).
;-
pro fig_flask_points, d, xkey, ykey, xkey_q, ykey_q, SKIP_PROFILES=skip, _EXTRA=ex

  compile_opt idl2

  x = d['flask_'+xkey] & y = d['flask_'+ykey]
  keep = (n_elements(skip) gt 0) ? where(d['flask_profile_index'] ge skip) : lindgen(n_elements(x))
  p = plot(x[keep], y[keep], linestyle=6, color='lime green', /overplot, _EXTRA=ex)
  if n_elements(xkey_q) eq 0 then return
  xq = d['flask_'+xkey_q] & yq = d['flask_'+ykey_q]
  p = plot(xq[keep], yq[keep], linestyle=6, color='lime green', /sym_filled, /overplot, _EXTRA=ex)

end
