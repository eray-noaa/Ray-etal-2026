;+
; Supplementary Figure 8: mean age vs normalized N2O for 1990s and 2020s in situ
; data and two previously published 1990s relationships.
; Data: <output>/figure_data/supp_fig8.nc.
;-
pro supp_fig8

  compile_opt idl2

  d = fig_data('supp_fig8')
  p = plot(indgen(2), /nodata, /buffer, yrange=[1.03,0], xrange=[0,7], ytitle='Normalized N$_2$O', xtitle='Mean Age (years)', $
    title='Mean Age vs. N$_2$O', font_size=11)
  p = plot(d['age_1990s'], d['norm_grid'], symbol='s', color='blue', sym_size=1, /sym_filled, thick=2, /overplot)
  p = plot(d['age_2020s'], d['norm_grid'], symbol='s', color='red', sym_size=1, /sym_filled, thick=2, /overplot)
  p = plot(d['age_andrews2001'], d['n2o_norm_grid_fine'], thick=2, color='lime green', /overplot)
  p = plot(d['age_engel2002'], d['n2o_norm_grid_fine'], thick=2, color='magenta', /overplot)
  t = text(0.18, 0.8, '1990s this study', color='blue', font_size=11, /norm)
  t = text(0.18, 0.76, '2020s this study', color='red', font_size=11, /norm)
  t = text(0.18, 0.72, '1990s (Andrews et al., 2001)', color='lime green', font_size=11, /norm)
  t = text(0.18, 0.68, '1990s (Engel et al., 2002)', color='magenta', font_size=11, /norm)
  fig_save, p, 'Supplementary_Figure_8'

end
