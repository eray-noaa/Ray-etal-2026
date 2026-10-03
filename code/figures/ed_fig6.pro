;+
; Extended Data Figure 6: trend profiles of (a) normalized CH4 and (b) mean
; age in normalized CH4 bins. Data: <output>/figure_data/ed_fig6.nc.
;-
pro ed_fig6

  compile_opt idl2

  d = fig_data('ed_fig6')
  ac = d['ccmi_alt'] & alt = d['alt_grid2'] & nz = n_elements(alt)

  p = plot(indgen(2), /nodata, /buffer, yrange=[13.5,31], xrange=[-0.028,0.018], ytitle='Altitude (km)', $
    xtitle='Normalized Tracer Fraction/Decade', title='a  NH CH$_4$ Trend Profiles', font_size=11, margin=[0.12,0.09,0.04,0.08], dimensions=[450,500])
  p = plot([0,0], [15,36], linestyle=2, /overplot)
  r2 = d['ccmi_refd2_ch4_norm_trend_range']
  fig_band, 10*r2[7:30,*], ac[7:30], fill_color='violet', fill_transparency=70, color='pink'
  fig_band, 10*d['ccmi_refd1_ch4_norm_trend_range'], ac, fill_color='lime green', fill_transparency=70, color='light green'
  p = plot(10.*d['ccmi_refd1_ch4_norm_trend_avg'], ac, color='dark green', thick=2, /overplot)
  p = plot(10.*d['ccmi_refd2_ch4_norm_trend_avg'], ac, color='medium orchid', thick=2, /overplot)
  p = errorplot(10.*d['ace_ch4_trend'], d['ace_alt'], 10.*d['ace_ch4_trend_sigma'], replicate(!values.f_nan,n_elements(d['ace_alt'])), $
    symbol='s', sym_size=1.25, color='black', errorbar_capsize=0, linestyle=6, /overplot)
  p = errorplot(10.*d['insitu_ch4_norm_trend'], alt, 10.*d['insitu_ch4_norm_trend_sigma'], replicate(!values.f_nan,nz), symbol='o', $
    sym_size=1.5, /sym_filled, color='blue', linestyle=6, errorbar_capsize=0, /overplot)
  s = symbol(0.19, 0.86, 'o', sym_size=1.25, sym_color='blue', /sym_filled, /norm)
  t = text(0.22, 0.85, 'In situ', color='blue', font_size=10, /norm)
  s = symbol(0.19, 0.82, 's', sym_size=1.25, sym_color='black', /norm)
  t = text(0.22, 0.81, 'ACE-FTS', color='black', font_size=10, /norm)
  fig_save, p, 'Extended_Data_Figure_6a'

  ng = d['norm_grid'] & nn = n_elements(ng) & gm = d['ccmi_norm_grid']
  p = plot(indgen(2), /nodata, /buffer, yrange=[1.01,0.4], xrange=[-0.25,0.11], ytitle='Normalized CH$_4$', xtitle='Years/Decade', $
    font_size=11, title='b  Mean Age Trends', axis_style=1, margin=[0.13,0.09,0.12,0.08], dimensions=[450,500])
  p = plot([0,0], [1,0.4], linestyle=2, /overplot)
  fig_band, 10*d['ccmi_refd2_age_trend_on_ch4_range'], gm, fill_color='violet', fill_transparency=70, color='pink'
  fig_band, 10*d['ccmi_refd1_age_trend_on_ch4_range'], gm, fill_color='lime green', fill_transparency=70, color='light green'
  tr = d['insitu_age_trend_on_ch4'] & trs = d['insitu_age_trend_on_ch4_sigma']
  p = errorplot(10.*tr[5:-2], ng[5:-2], 10.*trs[5:-2], replicate(!values.f_nan,nn-6), symbol='o', sym_size=1.5, /sym_filled, color='blue', $
    errorbar_capsize=0, /overplot)
  p = plot(10.*d['ccmi_refd1_age_trend_on_ch4_avg'], gm, color='dark green', thick=2, /overplot)
  p = plot(10.*d['ccmi_refd2_age_trend_on_ch4_avg'], gm, color='medium orchid', thick=2, /overplot)
  xaxis = axis('X', location='top', tickfont_size=0)
  yaxis = axis('Y', location='right', coord_transform=[0,1], minor=0, tickfont_size=11, tickname=['36','31','27','23','21','19','15'], $
    title='Approximate Altitude (km)')
  t = text(0.31, 0.86, 'CCMI-2022', font_size=10, alignment=0.5, /norm)
  p = plot([-0.19,-0.17], [0.465,0.465], thick=2, color='dark green', /overplot)
  t = text(0.34, 0.82, 'refD1', color='dark green', font_size=10, alignment=0.5, /norm)
  p = plot([-0.19,-0.17], [0.495,0.495], thick=2, color='medium orchid', /overplot)
  t = text(0.34, 0.78, 'refD2', color='medium orchid', font_size=10, alignment=0.5, /norm)
  fig_save, p, 'Extended_Data_Figure_6b'

end
