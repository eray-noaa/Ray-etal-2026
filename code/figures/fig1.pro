;+
; Figure 1: mean age vs normalized N2O.
; Data: <output>/figure_data/fig1.nc (written by age_time_series.pro).
;-
pro fig1

  compile_opt idl2

  d = fig_data('fig1')

  p = plot(indgen(2), /nodata, /buffer, yrange=[1.03,0], xrange=[0,7.3], ytitle='Normalized N$_2$O', xtitle='Mean Age (years)', $
    title='Mean Age vs. N$_2$O', font_size=11, position=[0.1,0.1,0.97,0.9], dimensions=[700,500])
  c = contour(d['insitu_1990s_hist'], d['age_grid_hist'], d['norm_grid_hist'], c_value=[0.03,0.03], c_color=['plum','plum'], $
    c_label_show=0, transparency=50, /fill, /overplot)
  fig_flask_points, d, 'sf6_age', 'n2o_norm', 'sf6_age', 'n2o_norm_q', symbol='o', sym_size=0.75
  fig_flask_points, d, 'co2_age', 'n2o_norm', 'co2_age', 'n2o_norm_q', symbol='tu', sym_size=1
  fl = d['flask_avg_mean_age']
  p = errorplot(fl[*,2], d['norm_grid'], fl[*,3], replicate(0,n_elements(d['norm_grid'])), symbol='s', /sym_filled, thick=2, $
    color='dark green', linestyle=2, errorbar_capsize=0, /overplot)
  p = plot(d['tlp_base_mean_age_n2o'], d['tlp_base_n2o_norm'], thick=3, color='orange', linestyle=1, /overplot)
  p = plot(d['tlp_w20_mean_age_n2o'], d['tlp_w20_n2o_norm'], thick=3, color='orange', linestyle=2, /overplot)
  n2 = n_elements(d['norm_grid2'])
  a0 = d['aircraft_1990s_mean_age'] & b0 = d['balloon_1990s_mean_age']
  p = errorplot(a0[*,2], d['norm_grid2'], a0[*,3], replicate(0,n2), symbol='s', color='dodger blue', sym_size=1, /sym_filled, $
    linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(b0[*,2], d['norm_grid2'], b0[*,3], replicate(0,n2), symbol='s', color='purple', sym_size=1, /sym_filled, $
    linestyle=6, errorbar_capsize=0, /overplot)
  s = symbol(0.69, 0.28, 's', sym_size=1, /sym_filled, sym_color='purple', /norm)
  t = text(0.71, 0.27, 'In situ balloon 1990s', color='purple', font_size=11, /norm)
  s = symbol(0.69, 0.24, 's', sym_size=1, /sym_filled, sym_color='dodger blue', /norm)
  t = text(0.71, 0.23, 'In situ aircraft 1990s', color='dodger blue', font_size=11, /norm)
  p = plot([4.7,5], [0.9,0.9], thick=3, linestyle=1, color='orange', /overplot)
  t = text(0.71, 0.19, 'Idealized model', color='orange', font_size=11, /norm)
  p = plot([4.7,5], [0.95,0.95], thick=3, linestyle=2, color='orange', /overplot)
  t = text(0.71, 0.15, 'w+20%', color='orange', font_size=11, /norm)
  t = text(0.16, 0.82, 'Flask balloon 1970s-2000s', color='dark green', font_size=11, /norm)
  s = symbol(0.17, 0.79, 's', sym_size=1, /sym_filled, sym_color='dark green', /norm)
  t = text(0.18, 0.78, ' = Weighted mean $\pm$ uncert. on mean', color='dark green', font_size=11, /norm)
  s = symbol(0.17, 0.75, 'td', sym_size=1.5, /sym_filled, sym_color='lime green', /norm)
  t = text(0.18, 0.74, ' = CO$_2$', color='lime green', font_size=11, /norm)
  s = symbol(0.17, 0.71, 'o', sym_size=1, /sym_filled, sym_color='lime green', /norm)
  t = text(0.18, 0.7, ' = SF$_6$', color='lime green', font_size=11, /norm)
  fig_save, p, 'Figure_1'

end
