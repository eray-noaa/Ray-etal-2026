;+
; Supplementary Figure 7: CCMI-2022 refD2 NH midlatitude mean age trend
; profiles for 1990s-2020s and 2010-2025.
; Data: <output>/figure_data/supp_fig7.nc. Trends are plotted per decade.
;-
pro supp_fig7

  compile_opt idl2

  d = fig_data('supp_fig7')
  ac = d['ccmi_alt'] & r1 = d['refd2_1990s_2020s_range'] & r2 = d['refd2_2010_2025_range']
  p = plot(indgen(2), /nodata, /buffer, xrange=[-0.4,0.1], yrange=[13.5,30], ytitle='Altitude (km)', xtitle='Years/Decade', $
    title='NH Midlatitude Mean Age Trend Profiles', font_size=10, margin=[0.1,0.1,0.04,0.08], dimensions=[450,500])
  p = plot([0,0], [13,35], linestyle=2, /overplot)
  fig_band, 10*r1[7:30,*], ac[7:30], fill_color='violet', fill_transparency=70, color='violet'
  p = plot(10.*d['refd2_1990s_2020s_avg'], ac, color='medium orchid', thick=3, /overplot)
  fig_band, 10*r2[7:30,*], ac[7:30], fill_color='lime green', fill_transparency=70, color='white'
  p = plot(10.*d['refd2_2010_2025_avg'], ac, color='green', thick=3, /overplot)
  t = text(0.2, 0.23, 'refD2 1990s-2020s', color='medium orchid', font_size=10, /norm)
  t = text(0.2, 0.19, 'refD2 2010-2020s', color='green', font_size=10, /norm)
  fig_save, p, 'Supplementary_Figure_7'

end
