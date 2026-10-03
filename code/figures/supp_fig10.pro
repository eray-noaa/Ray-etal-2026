;+
; Supplementary Figure 10: archived minus updated CO2 mean ages vs normalized
; N2O for 1990s ER-2 missions. Data: <output>/figure_data/supp_fig10.nc.
;-
pro supp_fig10

  compile_opt idl2

  d = fig_data('supp_fig10')
  ng = d['n2o_norm_grid'] & nn = n_elements(ng)
  diff = d['age_archived'] - d['age_updated'] & unc = d['age_updated_uncert']
  names = d['mission'] & cols = d['mission_color'] & nm = n_elements(names)

  ; Mission average of the differences, then a 7-point smooth.
  avg = replicate(!values.f_nan, nn)
  for n = 0, nn-1 do begin
    gd = where(finite(diff[n,*]), ngd)
    if ngd gt 0 then avg[n] = mean(diff[n,gd])
  endfor
  gd = where(finite(avg))
  avg[gd] = smooth(avg[gd], 7, /edge_truncate)

  p = plot(indgen(2), /nodata, /buffer, yrange=[1,0.1], xrange=[-1.4,0.7], ytitle='Normalized N$_2$O', xtitle='Mean Age Difference (years)', $
    title='Mean Age Difference vs. N$_2$O', font_size=11, dimensions=[450,550], margin=[0.13,0.08,0.03,0.08])
  p = plot([0,0], [1,0], linestyle=2, /overplot)
  for m = 0, nm-1 do p = errorplot(diff[*,m], ng, unc[*,m], replicate(0,nn), symbol='o', sym_size=0.75, /sym_filled, linestyle=6, $
    color=cols[m], errorbar_capsize=0, /overplot)
  p = plot(avg, ng, thick=5, color='blue', /overplot)
  t = text(0.7, 0.8, 'Average', color='blue', font_size=11, /norm)
  labels = ['1994   ','1992-3 ','1995-6 ','1997   ','1999-2000 ']
  order = [1, 0, 2, 3, 4]
  for k = 0, nm-1 do t = text(0.7, 0.76-0.04*k, labels[order[k]]+'('+names[order[k]]+')', color=cols[order[k]], font_size=11, /norm)
  fig_save, p, 'Supplementary_Figure_10'

end
