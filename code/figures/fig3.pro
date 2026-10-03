;+
; Figure 3: (a) mean age time series in one normalized N2O bin, (b) mean age
; trends vs normalized N2O, (c) 1990s and 2020s mean age vs N2O relationships.
; Data: <output>/figure_data/fig3.nc.
;-
pro fig3

  compile_opt idl2

  d = fig_data('fig3')
  yrs = d['years'] & age = d['age_in_bin'] & sd = d['age_in_bin_sd'] & st = d['sample_type'] & fromch4 = d['n2o_from_ch4']
  bin = d['n2o_norm_bin']

  ; a
  p = plot(indgen(2), /nodata, /buffer, xrange=[1990,2028], yrange=[3.4,5.2], xtickinterval=5, ytitle='Mean Age (years)', xtitle='Years', $
    title='a  Mean Age in N$_2$O$_{norm}$ = '+string(bin[0],format='(f4.2)')+'-'+string(bin[1],format='(f4.2)'), font_size=12, $
    margin=[0.1,0.1,0.04,0.08], dimensions=[700,400])
  sets = list({k:3, ch4:1, color:'magenta', sym:'d', filled:0, size:2}, {k:3, ch4:0, color:'magenta', sym:'d', filled:1, size:2}, $
    {k:1, ch4:-1, color:'blue', sym:'o', filled:1, size:1}, {k:2, ch4:-1, color:'purple', sym:'tu', filled:1, size:2})
  foreach s, sets do begin
    si = (s.ch4 lt 0) ? where(st eq s.k, nsi) : where(st eq s.k and fromch4 eq s.ch4, nsi)
    if nsi gt 0 then p = errorplot(yrs[si], age[si], replicate(0,nsi), sd[si], color=s.color, symbol=s.sym, sym_filled=s.filled, $
      sym_size=s.size, linestyle=6, errorbar_capsize=0, /overplot)
  endforeach
  f = d['bin_trend_fit']
  p = plot([1992,2026], f[0] + f[1]*[1992,2026], color='blue', linestyle=2, /overplot)
  s = symbol(0.17, 0.24, 'tu', sym_size=2, sym_color='purple', /sym_filled, /norm)
  t = text(0.2, 0.23, 'In situ balloon', color='purple', font_size=11, /norm)
  s = symbol(0.17, 0.2, 'o', sym_size=1.5, sym_color='blue', /sym_filled, /norm)
  t = text(0.2, 0.19, 'In situ aircraft', color='blue', font_size=11, /norm)
  s = symbol(0.65, 0.83, 'd', sym_size=2, sym_color='magenta', /sym_filled, /norm)
  t = text(0.67, 0.82, 'AirCore N$_2$O measured', color='magenta', font_size=11, /norm)
  s = symbol(0.65, 0.78, 'd', sym_size=2, sym_color='magenta', /norm)
  t = text(0.67, 0.77, 'AirCore N$_2$O derived', color='magenta', font_size=11, /norm)
  fig_save, p, 'Figure_3a'

  ; b
  ng = d['norm_grid'] & nn = n_elements(ng) & gm = d['ccmi_n2o_norm_grid']
  p = plot(indgen(2), /nodata, /buffer, xrange=[-0.31,0.11], yrange=[1.02,0.2], ytitle='Normalized N$_2$O', xtitle='Years/Decade', $
    title='b  NH Mean Age Trends', font_size=10, axis_style=1, margin=[0.12,0.09,0.12,0.07], dimensions=[450,500])
  p = plot([0,0], [1,0.1], linestyle=2, /overplot)
  tr = d['insitu_age_trend_on_n2o'] & trs = d['insitu_age_trend_on_n2o_sigma']
  p = errorplot(10.*tr[7:-2], ng[7:-2], 10.*trs[7:-2], replicate(!values.f_nan,nn-8), symbol='o', sym_size=1.5, /sym_filled, color='blue', $
    errorbar_capsize=0, /overplot)
  fig_band, 10*d['ccmi_refd2_age_trend_on_n2o_range'], gm, fill_color='violet', fill_transparency=70, color='pink'
  fig_band, 10*d['ccmi_refd1_age_trend_on_n2o_range'], gm, fill_color='lime green', fill_transparency=70, color='light green'
  p = plot(10.*d['ccmi_refd1_age_trend_on_n2o_avg'], gm, color='dark green', thick=3, /overplot)
  d2 = d['ccmi_refd2_age_trend_on_n2o_avg']
  p = plot(10.*d2[0:-2], gm[0:-2], color='medium orchid', thick=3, /overplot)
  xaxis = axis('X', location='top', tickfont_size=0)
  yaxis = axis('Y', location='right', coord_transform=[0,1], minor=0, tickfont_size=10, tickname=['35','27','22','19','12'], $
    title='Approximate Altitude (km)')
  s = symbol(0.17, 0.87, 'o', sym_size=1.25, sym_color='blue', /sym_filled, /norm)
  t = text(0.2, 0.86, 'In situ', color='blue', font_size=10, /norm)
  t = text(0.42, 0.865, 'CCMI-2022', font_size=10, alignment=0.5, /norm)
  p = plot([-0.19,-0.17], [0.29,0.29], thick=3, color='dark green', /overplot)
  t = text(0.44, 0.83, 'refD1', color='dark green', font_size=10, alignment=0.5, /norm)
  p = plot([-0.19,-0.17], [0.33,0.33], thick=3, color='medium orchid', /overplot)
  t = text(0.44, 0.79, 'refD2', color='medium orchid', font_size=10, alignment=0.5, /norm)
  fig_save, p, 'Figure_3b'

  ; c
  p = plot(indgen(2), /nodata, /buffer, yrange=[1.03,0.2], xrange=[0,6], ytitle='Normalized N$_2$O', xtitle='Mean Age (years)', $
    title='c  NH Extratropical Mean Age vs. N$_2$O', font_size=12, margin=[0.1,0.1,0.04,0.08], dimensions=[700,400])
  p = plot(d['age_on_n2o_1990s'], d['norm_grid2'], color='blue', thick=3, /overplot)
  p = plot(d['age_on_n2o_2020s'], d['norm_grid2'], color='red', thick=3, /overplot)
  aa = d['arrow_age'] & an = d['arrow_n2o_norm']
  offs = [0.01, 0.02, 0.01, 0.01, -0.01, 0.]   ; vertical offsets used in the paper figure for legibility
  for i = 0, n_elements(offs)-1 do a = arrow(reform(aa[*,i]), offs[i] + reform(an[*,i]), /data, head_size=0.75, color='blue violet', /current)
  t = text(0.18, 0.8, '1990s', color='blue', font_size=11, /norm)
  t = text(0.18, 0.76, '2020s', color='red', font_size=11, /norm)
  fig_save, p, 'Figure_3c'

end
