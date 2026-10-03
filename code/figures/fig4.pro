;+
; Figure 4: tropical upwelling trends from CCMI-2022 and the TLP model.
; Data: <output>/figure_data/fig4.nc.
;-
pro fig4

  compile_opt idl2

  d = fig_data('fig4')
  alt = d['tlp_alt'] & rng = d['tlp_w_trend_range'] & ac = d['ccmi_alt']
  sm = [[smooth(rng[*,0],3,/edge_truncate)], [smooth(rng[*,1],3,/edge_truncate)]]

  p = plot(indgen(2), /nodata, /buffer, yrange=[16,40], xrange=[-3,12], xtitle='%/Decade', ytitle='Altitude (km)', $
    title='Tropical Upwelling Trends', font_size=11, dimensions=[350,500], margin=[0.15,0.07,0.05,0.07])
  fig_band, sm, alt, fill_color='sky blue', fill_transparency=80, color='sky blue'
  r2 = d['ccmi_refd2_w_trend_range'] & r1 = d['ccmi_refd1_w_trend_range']
  fig_band, r2[7:30,*], ac[7:30], fill_color='violet', fill_transparency=80, color='pink'
  fig_band, r1[10:30,*], ac[10:30], fill_color='lime green', fill_transparency=80, color='light green'
  p = plot(d['ccmi_refd2_w_trend_avg'], ac, color='medium orchid', thick=3, /overplot)
  p = plot(d['ccmi_refd1_w_trend_avg'], ac, color='dark green', thick=3, /overplot)
  p = plot(sm[*,0], alt, thick=1, color='sky blue', /overplot)
  p = plot(sm[*,1], alt, thick=1, color='sky blue', /overplot)
  p = plot(smooth(d['tlp_w_trend_wmean'],3,/edge_truncate), alt, thick=3, color='blue', /overplot)
  p = plot([0,0], [16,50], linestyle=2, thick=2, /overplot)
  t = text(0.76, 0.33, 'CCMI-2022', font_size=10, alignment=0.5, /norm)
  p = plot([6.5,7.5], [22.5,22.5], thick=3, color='dark green', /overplot)
  t = text(0.79, 0.29, 'refD1', color='dark green', font_size=10, alignment=0.5, /norm)
  p = plot([6.5,7.5], [21.6,21.6], thick=3, color='medium orchid', /overplot)
  t = text(0.79, 0.26, 'refD2', color='medium orchid', font_size=10, alignment=0.5, /norm)
  p = plot([6.3,7.3], [19.4,19.4], thick=3, color='blue', /overplot)
  t = text(0.81, 0.18, 'TLP model', color='blue', font_size=10, alignment=0.5, /norm)
  fig_save, p, 'Figure_4'

end
