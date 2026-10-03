;+
; Figure 2: NH midlatitude mean age (a) and normalized N2O (b) trend profiles.
; Data: <output>/figure_data/fig2.nc. Trends are stored per year and plotted per decade.
;-
pro fig2

  compile_opt idl2

  d = fig_data('fig2')
  nz = n_elements(d['alt_grid2'])

  ; a: mean age
  p = plot(indgen(2), /nodata, /buffer, xrange=[-0.3,0.07], yrange=[13.5,30], ytitle='Altitude (km)', xtitle='Years/Decade', $
    title='a  Mean Age', font_size=10, margin=[0.1,0.1,0.04,0.08], dimensions=[450,500])
  p = plot([0,0], [13,35], linestyle=2, /overplot)
  fig_band, 10*d['ccmi_refd2_age_trend_range'], d['ccmi_alt'], fill_color='violet', fill_transparency=70, color='pink'
  fig_band, 10*d['ccmi_refd1_age_trend_range'], d['ccmi_alt'], fill_color='lime green', fill_transparency=65, color='light green'
  p = plot(10.*d['ccmi_refd1_age_trend_avg'], d['ccmi_alt'], color='dark green', thick=3, /overplot)
  p = plot(10.*d['ccmi_refd2_age_trend_avg'], d['ccmi_alt'], color='medium orchid', thick=3, /overplot)
  p = errorplot(10.*d['insitu_age_trend'], d['alt_grid2'], 10.*d['insitu_age_trend_sigma'], replicate(!values.f_nan,nz), symbol='o', $
    sym_size=1.5, /sym_filled, color='blue', linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(10.*d['ace_age_trend'], d['ace_age_trend_alt'], 10.*d['ace_age_trend_sigma'], replicate(!values.f_nan,2), symbol='s', $
    sym_size=1.25, color='black', linestyle=6, errorbar_capsize=0, /overplot)
  s = symbol(0.17, 0.51, 'o', sym_size=1.5, sym_color='blue', /sym_filled, /norm)
  t = text(0.2, 0.5, 'In situ', color='blue', font_size=10, /norm)
  s = symbol(0.17, 0.47, 's', sym_size=1.25, sym_color='black', /norm)
  t = text(0.2, 0.46, 'ACE-FTS', color='black', font_size=10, /norm)
  fig_save, p, 'Figure_2a'

  ; b: normalized N2O
  p = plot(indgen(2), /nodata, /buffer, yrange=[13.5,30], xrange=[-0.04,0.019], ytitle='Altitude (km)', $
    xtitle='Normalized Tracer Fraction/Decade', title='b  N$_2$O', font_size=10, margin=[0.1,0.1,0.04,0.08], dimensions=[450,500])
  p = plot([0,0], [15,36], linestyle=2, /overplot)
  fig_band, 10*d['ccmi_refd2_n2o_norm_trend_range'], d['ccmi_alt'], fill_color='violet', fill_transparency=70, color='pink'
  fig_band, 10*d['ccmi_refd1_n2o_norm_trend_range'], d['ccmi_alt'], fill_color='lime green', fill_transparency=65, color='light green'
  p = plot(10.*d['ccmi_refd1_n2o_norm_trend_avg'], d['ccmi_alt'], color='dark green', thick=3, /overplot)
  p = plot(10.*d['ccmi_refd2_n2o_norm_trend_avg'], d['ccmi_alt'], color='medium orchid', thick=3, /overplot)
  p = errorplot(10.*d['ace_n2o_trend'], d['ace_alt'], 10.*d['ace_n2o_trend_sigma'], replicate(!values.f_nan,n_elements(d['ace_alt'])), $
    symbol='s', sym_size=1.25, color='black', errorbar_capsize=0, linestyle=6, /overplot)
  p = errorplot(10.*d['insitu_n2o_norm_trend'], d['alt_grid2'], 10.*d['insitu_n2o_norm_trend_sigma'], replicate(!values.f_nan,nz), $
    symbol='o', sym_size=1.5, /sym_filled, color='blue', linestyle=6, errorbar_capsize=0, /overplot)
  t = text(0.175, 0.31, 'CCMI-2022', font_size=10, /norm)
  p = plot([-0.035,-0.031], [17.1,17.1], thick=2, color='dark green', /overplot)
  t = text(0.25, 0.27, 'refD1', color='dark green', font_size=10, /norm)
  p = plot([-0.035,-0.031], [16.3,16.3], thick=2, color='medium orchid', /overplot)
  t = text(0.25, 0.23, 'refD2', color='medium orchid', font_size=10, /norm)
  fig_save, p, 'Figure_2b'

end
