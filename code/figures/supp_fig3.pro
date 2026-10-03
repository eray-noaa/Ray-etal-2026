;+
; Supplementary Figure 3: N2O equivalent latitude example for AirCore data.
; Data: <output>/figure_data/supp_fig3.nc.
;-
pro supp_fig3

  compile_opt idl2

  d = fig_data('supp_fig3')
  ab = d['alt_bin']
  p = plot(indgen(2), /nodata, /buffer, yrange=[1.01,0.77], ytitle='Normalized N$_2$O', xrange=[0,90], xtitle='Latitude', $
    title='N$_2$O  Alt='+string(mean(ab),format='(f5.2)')+'km', font_size=11, margin=[0.11,0.1,0.04,0.08], dimensions=[700,400])
  p = plot(d['swoosh_lat'], d['swoosh_n2o_norm_seasonal'], color='green', thick=2, /overplot)
  p = plot(d['equiv_lat_merra2'], d['n2o_norm'], symbol='o', /sym_filled, color='blue', sym_size=0.75, linestyle=6, /overplot)
  p = plot(d['lat'], d['n2o_norm'], symbol='o', /sym_filled, color='sky blue', sym_size=0.75, linestyle=6, /overplot)
  p = plot(d['equiv_lat_n2o'], d['n2o_norm'], symbol='o', /sym_filled, color='red', sym_size=0.75, linestyle=6, /overplot)
  t = text(0.15, 0.83, 'Lat', font_size=11, color='sky blue', /norm)
  t = text(0.15, 0.78, 'Elat from MERRA2', font_size=11, color='blue', /norm)
  t = text(0.15, 0.73, 'SWOOSH seas avg', font_size=11, color='green', /norm)
  t = text(0.15, 0.68, 'Elat from SWOOSH', font_size=11, color='red', /norm)
  fig_save, p, 'Supplementary_Figure_3'

end
