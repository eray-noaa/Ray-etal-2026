;+
; Extended Data Figure 5: (a) mean age trend profiles with and without flask
; data and from previous studies, (b) number of measurements by platform.
; Data: <output>/figure_data/ed_fig5.nc. Trends are plotted per decade.
;-
pro ed_fig5

  compile_opt idl2

  d = fig_data('ed_fig5')
  alt = d['alt_grid2'] & nz = n_elements(alt) & nanz = replicate(!values.f_nan, nz)

  p = plot(indgen(2), /nodata, /buffer, yrange=[13.5,35], xrange=[-0.35,0.4], ytitle='Altitude (km)', xtitle='Years/Decade', $
    title='a  NH Midlatitude Mean Age Trends', font_size=11, margin=[0.13,0.1,0.08,0.08], dimensions=[450,500])
  p = plot([0,0], [14,35], linestyle=2, /overplot)
  p = errorplot(10.*d['age_trend_no_flask'], alt, 10.*d['age_trend_no_flask_sigma'], nanz, symbol='o', sym_size=1.5, /sym_filled, $
    color='sky blue', linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(10.*d['age_trend_no_flask_lat_adj'], alt, 10.*d['age_trend_no_flask_lat_adj_sigma'], nanz, symbol='o', sym_size=1.5, $
    /sym_filled, color='blue', linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(10.*d['age_trend_with_flask'], alt, 10.*d['age_trend_with_flask_sigma'], nanz, symbol='o', sym_size=1, color='lime green', $
    linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(10.*d['ray2014_trend'], d['ray2014_alt'], 10.*d['ray2014_trend_sigma'], replicate(!values.f_nan,4), symbol='D', sym_size=1.5, $
    color='black', linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(10.*[d['fritsch2020_trend'],d['fritsch2020_trend']], [d['fritsch2020_alt'],d['fritsch2020_alt']], $
    10.*[d['fritsch2020_trend_sigma'],d['fritsch2020_trend_sigma']], replicate(!values.f_nan,2), symbol='s', sym_size=1.5, color='purple', $
    linestyle=6, errorbar_capsize=0, /overplot)
  t = text(0.17, 0.86, 'In situ w/o flasks lat adj', color='blue', font_size=11, /norm)
  t = text(0.17, 0.82, 'In situ w/o flasks', color='sky blue', font_size=11, /norm)
  t = text(0.17, 0.78, 'In situ w/flasks', color='lime green', font_size=11, /norm)
  t = text(0.17, 0.74, 'Ray et al., 2014', color='black', font_size=11, /norm)
  t = text(0.17, 0.7, 'Fritsch et al., 2020', color='purple', font_size=11, /norm)
  fig_save, p, 'Extended_Data_Figure_5a'

  p = plot(indgen(2), /nodata, /buffer, /xlog, xrange=[0.7,1e6], yrange=[13.5,35], ytitle='Altitude (km)', xtitle='Number of measurements', $
    title='b  NH Midlatitude number of measurements', font_size=11, margin=[0.13,0.1,0.08,0.08], dimensions=[450,500])
  p = plot(d['n_total'], alt, color='red', thick=2, /overplot)
  p = plot(d['n_aircraft'], alt, symbol='s', /sym_filled, color='blue', linestyle=6, /overplot)
  p = plot(d['n_insitu_balloon'], alt, symbol='s', /sym_filled, color='purple', linestyle=6, /overplot)
  p = plot(d['n_aircore'], alt, symbol='s', /sym_filled, color='magenta', linestyle=6, /overplot)
  p = plot(d['n_flask'], alt, symbol='s', /sym_filled, color='lime green', linestyle=6, /overplot)
  t = text(0.67, 0.82, 'Aircraft in situ', color='blue', font_size=11, /norm)
  t = text(0.67, 0.78, 'Balloon in situ', color='purple', font_size=11, /norm)
  t = text(0.67, 0.74, 'AirCore', color='magenta', font_size=11, /norm)
  t = text(0.67, 0.7, 'Flask', color='lime green', font_size=11, /norm)
  t = text(0.67, 0.66, 'Total', color='red', font_size=11, /norm)
  fig_save, p, 'Extended_Data_Figure_5b'

end
