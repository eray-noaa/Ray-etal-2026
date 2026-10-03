;+
; Extended Data Figure 3: correlations of normalized CH4 vs N2O (a),
; N2O vs CFC-12 (b) and CH4 vs CFC-12 (c).
; Data: <output>/figure_data/ed_fig2_3.nc.
;-
pro ed_fig3

  compile_opt idl2

  d = fig_data('ed_fig2_3')
  panels = list({x:'n2o_norm', y:'ch4_norm', xl:'Normalized N$_2$O', yl:'Normalized CH$_4$', xr:[1.1,0], yr:[1.1,0], skip:0, p:'a', t:'CH$_4$ vs. N$_2$O'}, $
    {x:'f12_norm', y:'n2o_norm', xl:'Normalized CFC-12', yl:'Normalized N$_2$O', xr:[1.04,-0.03], yr:[1.01,0.05], skip:2, p:'b', t:'N$_2$O vs. CFC-12'}, $
    {x:'f12_norm', y:'ch4_norm', xl:'Normalized CFC-12', yl:'Normalized CH$_4$', xr:[1.1,0], yr:[1.1,0], skip:2, p:'c', t:'CH$_4$ vs. CFC-12'})
  foreach pn, panels do begin
    p = plot(indgen(2), /nodata, /buffer, yrange=pn.yr, xrange=pn.xr, ytitle=pn.yl, xtitle=pn.xl, title=pn.p+'  '+pn.t, font_size=11, $
      position=[0.1,0.1,0.97,0.9], dimensions=[600,500])
    p = plot([1.1,0], [1.1,0], linestyle=2, /overplot)
    if pn.p eq 'a' then p = plot(d['norm_grid'], d['aircore_ch4_norm_on_n2o'], symbol='o', color='green', /sym_filled, /overplot)
    foreach reg, ['tropics','vortex','midlat'] do begin
      col = (reg eq 'tropics') ? 'red' : ((reg eq 'vortex') ? 'purple' : 'blue')
      xx = d['oms_'+reg+'_'+pn.x] & yy = d['oms_'+reg+'_'+pn.y]
      p = plot(xx[*,3], yy[*,3], thick=3, color=col, /overplot)
    endforeach
    fig_flask_points, d, pn.x, pn.y, pn.x+'_q', pn.y+'_q', symbol='o', sym_size=0.75, skip_profiles=pn.skip
    t = text(0.15, 0.82, 'Solid lines = In situ balloon', font_size=11, /norm)
    t = text(0.15, 0.78, 'Tropics', color='red', font_size=11, /norm)
    t = text(0.27, 0.78, 'Midlat', color='blue', font_size=11, /norm)
    t = text(0.37, 0.78, 'Vortex', color='purple', font_size=11, /norm)
    t = text(0.15, 0.73, 'Symbols = Flask balloon', color='lime green', font_size=11, /norm)
    fig_save, p, 'Extended_Data_Figure_3'+pn.p
  endforeach

end
