;+
; Supplementary Figure 6: NH midlatitude mean age trends vs altitude (a) and
; normalized N2O (b) for 1993-2025 and the AirCore period 2011-2025.
; Data: <output>/figure_data/supp_fig6.nc. Trends are plotted per decade.
;-
pro supp_fig6

  compile_opt idl2

  d = fig_data('supp_fig6')
  alt = d['alt_grid2'] & nz = n_elements(alt)
  p = plot(indgen(2), /nodata, /buffer, xrange=[-0.4,0.2], yrange=[13.5,30], ytitle='Altitude (km)', xtitle='Years/Decade', $
    title='a  NH Midlatitude Mean Age Trend Profiles', font_size=10, margin=[0.1,0.1,0.04,0.08], dimensions=[450,500])
  p = plot([0,0], [13,35], linestyle=2, /overplot)
  p = errorplot(10.*d['trend_alt_1993_2025'], alt, 10.*d['trend_alt_1993_2025_sigma'], replicate(!values.f_nan,nz), symbol='o', sym_size=1.5, $
    /sym_filled, color='blue', linestyle=6, errorbar_capsize=0, /overplot)
  p = errorplot(10.*d['trend_alt_2011_2025'], alt, 10.*d['trend_alt_2011_2025_sigma'], replicate(!values.f_nan,nz), symbol='o', sym_size=1.5, $
    /sym_filled, color='sky blue', linestyle=6, errorbar_capsize=0, /overplot)
  t = text(0.16, 0.86, 'In situ (1993-2025)', color='blue', font_size=10, /norm)
  t = text(0.16, 0.82, 'In situ (2011-2025)', color='sky blue', font_size=10, /norm)
  fig_save, p, 'Supplementary_Figure_6a'

  ng = d['norm_grid'] & nn = n_elements(ng) & i = lindgen(nn-7) + 6
  p = plot(indgen(2), /nodata, /buffer, xrange=[-0.4,0.11], yrange=[1.02,0.2], ytitle='Normalized N$_2$O', xtitle='Years/Decade', $
    title='b  NH Mean Age Trends', font_size=11, margin=[0.13,0.09,0.12,0.07], dimensions=[450,500])
  p = plot([0,0], [1,0.1], linestyle=2, /overplot)
  a = d['trend_n2o_1993_2025'] & as = d['trend_n2o_1993_2025_sigma'] & b = d['trend_n2o_2011_2025'] & bs = d['trend_n2o_2011_2025_sigma']
  p = errorplot(10.*a[i], ng[i], 10.*as[i], replicate(!values.f_nan,n_elements(i)), symbol='o', sym_size=1.5, /sym_filled, color='blue', $
    errorbar_capsize=0, /overplot)
  p = errorplot(10.*b[i], ng[i], 10.*bs[i], replicate(!values.f_nan,n_elements(i)), symbol='o', sym_size=1.5, /sym_filled, color='sky blue', $
    errorbar_capsize=0, /overplot)
  t = text(0.16, 0.86, 'In situ (1993-2025)', color='blue', font_size=10, /norm)
  t = text(0.16, 0.82, 'In situ (2011-2025)', color='sky blue', font_size=10, /norm)
  fig_save, p, 'Supplementary_Figure_6b'

end
