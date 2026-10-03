;+
; Supplementary Figure 5: NH midlatitude (30-60N) mean ages at 20 km (a),
; average sampled latitude (b), latitude adjustment (c) and the climatological
; seasonal latitude gradients used for the adjustment (d).
; Data: <output>/figure_data/supp_fig5.nc.
;-
pro supp_fig5

  compile_opt idl2

  d = fig_data('supp_fig5')
  yrs = d['years'] & st = d['sample_type'] & age = d['age'] & adj = d['age_lat_adj'] & unc = d['age_lat_adj_uncert']
  alt = string(d['alt'], format='(i2)')
  f = d['fit'] & fa = d['fit_lat_adj']

  p = plot(indgen(2), /nodata, /buffer, xrange=[1990,2028], xtickinterval=5, ytitle='Mean Age (years)', xtitle='Years', $
    title='a  30$\deg$N-60$\deg$N Average Mean Age at '+alt+'km Altitude', yrange=[2.05,5], font_size=11, margin=[0.1,0.1,0.03,0.08], dimensions=[700,400])
  sets = list({k:1, color:'blue', sym:'o', size:1}, {k:2, color:'purple', sym:'tu', size:2}, {k:3, color:'magenta', sym:'d', size:2})
  foreach s, sets do begin
    chk = where(st eq s.k, nchk)
    if nchk eq 0 then continue
    p = plot(yrs[chk], age[chk], color=s.color, symbol=s.sym, sym_size=s.size, linestyle=6, /overplot)
    p = errorplot(yrs[chk], adj[chk], replicate(0,nchk), unc[chk], color=s.color, symbol=s.sym, sym_size=s.size, /sym_filled, linestyle=6, $
      errorbar_capsize=0, /overplot)
  endforeach
  p = plot([1990,2026], f[1]*[1990,2026] + f[0], color='orange', thick=2, linestyle=1, /overplot)
  p = plot([1990,2026], fa[1]*[1990,2026] + fa[0], thick=2, color='orange', linestyle=2, /overplot)
  s = symbol(0.38, 0.24, 'tu', sym_size=2, sym_color='purple', /sym_filled, /norm)
  t = text(0.4, 0.23, 'In situ balloon', color='purple', font_size=11, /norm)
  s = symbol(0.38, 0.2, 'o', sym_size=1.5, sym_color='blue', /sym_filled, /norm)
  t = text(0.4, 0.19, 'In situ aircraft', color='blue', font_size=11, /norm)
  s = symbol(0.38, 0.16, 'd', sym_size=2, sym_color='magenta', /sym_filled, /norm)
  t = text(0.4, 0.15, 'AirCore', color='magenta', font_size=11, /norm)
  s = symbol(0.35, 0.835, 'o', /norm)
  s = symbol(0.35, 0.795, 'o', /sym_filled, /norm)
  t = text(0.37, 0.82, '= Orig', font_size=11, /norm)
  t = text(0.37, 0.78, '= Lat Adj', font_size=11, /norm)
  p = plot([2007,2009], [4.7,4.7], color='orange', linestyle=1, thick=2, /overplot)
  p = plot([2007,2009], [4.55,4.55], color='orange', linestyle=2, thick=2, /overplot)
  t = text(0.53, 0.82, ' = Fit ('+string(10.*f[1],format='(f5.2)')+'yrs/dec, RMSE='+string(d['rmse'],format='(f4.2)')+'yrs)', font_size=11, /norm)
  t = text(0.53, 0.78, ' = Fit ('+string(10.*fa[1],format='(f5.2)')+'yrs/dec, RMSE='+string(d['rmse_lat_adj'],format='(f4.2)')+'yrs)', font_size=11, /norm)
  fig_save, p, 'Supplementary_Figure_5a'

  p = plot(indgen(2), /nodata, /buffer, xrange=[1990,2028], xtickinterval=5, yrange=[27,63], ytitle='Avg sampled latitude', xtitle='Years', $
    title='b  Time series of average sampled latitudes at '+alt+'km Altitude', font_size=11, margin=[0.1,0.1,0.03,0.08], dimensions=[700,400])
  p = plot([1990,2028], [45,45], linestyle=2, /overplot)
  p = plot(yrs, d['avg_sampled_lat'], color='green', symbol='o', /sym_filled, linestyle=6, /overplot)
  fig_save, p, 'Supplementary_Figure_5b'

  p = plot(indgen(2), /nodata, /buffer, xrange=[1990,2028], xtickinterval=5, ytitle='Mean Age Adjustment (years)', xtitle='Years', $
    title='c  Mean Age Adjustment at '+alt+'km Altitude', font_size=11, margin=[0.1,0.1,0.03,0.08], dimensions=[700,400])
  p = plot([1990,2028], [0,0], linestyle=2, /overplot)
  p = plot(yrs, adj-age, color='green', symbol='o', /sym_filled, linestyle=6, /overplot)
  fig_save, p, 'Supplementary_Figure_5c'

  c = d['clim_age_rel_45n'] & lat = d['clim_lat']
  p = plot(indgen(2), /nodata, /buffer, ytitle='years', xrange=[20,70], xminor=1, yrange=[-1.2,1.2], xtitle='Latitude', $
    title='d  Climatological Seasonal Mean Age Relative to 45$\deg$N  Alt='+string(d['clim_alt'],format='(i2)')+'km', font_size=11, $
    margin=[0.1,0.1,0.03,0.08], dimensions=[700,400])
  p = plot([20,70], [0,0], linestyle=2, /overplot)
  p = plot([45,45], [-1.5,1.5], linestyle=2, /overplot)
  cols = ['blue','green','orange','red'] & seas = ['DJF','MAM','JJA','SON']
  for s = 0, 3 do begin
    p = plot(lat, c[*,s], thick=3, color=cols[s], /overplot)
    t = text(0.17, 0.81-0.04*s, seas[s], color=cols[s], font_size=11, /norm)
  endfor
  fig_save, p, 'Supplementary_Figure_5d'

end
