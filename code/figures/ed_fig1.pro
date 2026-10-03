;+
; Extended Data Figure 1: mean age vs normalized CH4 (a) and CFC-12 (b).
; Data: <output>/figure_data/ed_fig1a.nc and ed_fig1b.nc.
;-
pro ed_fig1

  compile_opt idl2

  foreach tr, ['ch4','f12'] do begin
    d = fig_data('ed_fig1' + ((tr eq 'ch4') ? 'a' : 'b'))
    lab = (tr eq 'ch4') ? 'CH$_4$' : 'CFC-12'
    yr = (tr eq 'ch4') ? [1.03,0.1] : [1.03,-0.03]
    p = plot(indgen(2), /nodata, /buffer, yrange=yr, xrange=[0,7.5], ytitle='Normalized '+lab, xtitle='Mean Age (years)', $
      title=((tr eq 'ch4') ? 'a' : 'b')+'  Mean Age vs. '+lab, font_size=11, position=[0.1,0.1,0.97,0.9], dimensions=[700,500])
    c = contour(d['insitu_1990s_hist'], d['age_grid_hist'], d['norm_grid_hist'], c_value=[0.03,0.03], c_color=['plum','plum'], $
      c_label_show=0, transparency=50, /fill, /overplot)
    fig_flask_points, d, 'sf6_age', tr+'_norm', 'sf6_age', tr+'_norm_q', symbol='o', sym_size=0.75
    fig_flask_points, d, 'co2_age', tr+'_norm', 'co2_age', tr+'_norm_q', symbol='tu', sym_size=(tr eq 'ch4') ? 1.25 : 1
    p = plot(d['tlp_base_mean_age_'+tr], d['tlp_base_'+tr+'_norm'], thick=3, color='orange', linestyle=1, /overplot)
    p = plot(d['tlp_w20_mean_age_'+tr], d['tlp_w20_'+tr+'_norm'], thick=3, color='orange', linestyle=2, /overplot)
    fl = d['flask_avg_mean_age']
    p = errorplot(fl[*,2], d['norm_grid'], fl[*,3], replicate(0,n_elements(d['norm_grid'])), symbol='s', /sym_filled, thick=3, $
      color='green', linestyle=2, errorbar_capsize=0, /overplot)
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
    t = text(0.16, 0.82, 'Flask balloon 1970s-2000s', color='dark green', font_size=11, /norm)
    s = symbol(0.17, 0.79, 's', sym_size=1, /sym_filled, sym_color='dark green', /norm)
    t = text(0.18, 0.78, ' = Avg', color='dark green', font_size=11, /norm)
    s = symbol(0.17, 0.75, 'td', sym_size=1.5, /sym_filled, sym_color='lime green', /norm)
    t = text(0.18, 0.74, ' = CO$_2$', color='lime green', font_size=11, /norm)
    s = symbol(0.17, 0.71, 'o', sym_size=1, /sym_filled, sym_color='lime green', /norm)
    t = text(0.18, 0.7, ' = SF$_6$', color='lime green', font_size=11, /norm)
    fig_save, p, 'Extended_Data_Figure_1' + ((tr eq 'ch4') ? 'a' : 'b')
  endforeach

end
