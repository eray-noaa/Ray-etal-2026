;+
; FIG_BAND, xrange2, y, _EXTRA=ex
;
; Draws a shaded band between XRANGE2[*,0] and XRANGE2[*,1] along Y
; (used for model min/max ranges). Points where either bound is NaN are dropped.
;-
pro fig_band, xr, y, _EXTRA=ex

  compile_opt idl2

  gd = where(finite(xr[*,0]) and finite(xr[*,1]), ngd)
  if ngd lt 2 then return
  poly = polygon([xr[gd,0], reverse(xr[gd,1])], [y[gd], reverse(y[gd])], /data, /fill_background, _EXTRA=ex)

end
