;+
; Extended Data Figure 2: profiles of normalized N2O (a), CH4 (b), CFC-12 (c)
; and mean age (d) from OMS in situ balloon and cryo-flask balloon data.
; Data: <output>/figure_data/ed_fig2_3.nc.
;-
pro ed_fig2

  compile_opt idl2

  d = fig_data('ed_fig2_3')
  alt = d['alt_grid'] & alt2 = d['alt_grid2'] & nz2 = n_elements(alt2)
  panels = list({k:'n2o_norm', lab:'Normalized N$_2$O', xr:[1.03,0], p:'a'}, {k:'ch4_norm', lab:'Normalized CH$_4$', xr:[1.05,0.1], p:'b'}, $
    {k:'f12_norm', lab:'Normalized CFC-12', xr:[1.03,0], p:'c'}, {k:'age', lab:'Mean Age (years)', xr:[0,7.5], p:'d'})
  foreach pn, panels do begin
    p = plot(indgen(2), /nodata, /buffer, yrange=[14,36], xrange=pn.xr, ytitle='Altitude (km)', xtitle=pn.lab, font_size=11, $
      title=pn.p+'  '+((pn.k eq 'age') ? 'Mean Age' : pn.lab)+' Profiles', margin=[0.09,0.09,0.04,0.08], dimensions=[550,500])
    tr = d['oms_tropics_'+pn.k] & ml = d['oms_midlat_'+pn.k] & vx = d['oms_vortex_'+pn.k]
    if pn.k ne 'f12_norm' then p = plot(tr[*,3], alt, thick=3, color='red', /overplot)
    p = plot(ml[*,0], alt, thick=1, color='sky blue', /overplot)
    p = plot(ml[*,3], alt, thick=3, color='blue', /overplot)
    p = plot(ml[*,6], alt, thick=1, color='sky blue', /overplot)
    p = plot(vx[*,0], alt, thick=1, color='orchid', /overplot)
    p = plot(vx[*,3], alt, thick=(pn.k eq 'age') ? 2 : 3, color='purple', /overplot)
    p = plot(vx[*,6], alt, thick=1, color='orchid', /overplot)
    if pn.k eq 'f12_norm' then p = plot(tr[*,3], alt, thick=3, color='red', /overplot)
    if pn.k ne 'age' then begin
      fig_flask_points, d, pn.k, 'alt', pn.k+'_q', 'alt', symbol='o', sym_size=0.75
      fa = d['flask_avg_'+pn.k]
      p = errorplot(fa[*,0], alt2, fa[*,1], replicate(0,nz2), symbol='s', /sym_filled, thick=3, color='green', linestyle=2, $
        errorbar_capsize=0, /overplot)
    endif else begin
      fig_flask_points, d, 'sf6_age', 'alt', 'sf6_age_q', 'alt', symbol='o', sym_size=0.75
      p = plot(d['flask_co2_age_q'], d['flask_alt'], color='lime green', /sym_filled, symbol='tu', sym_size=1, linestyle=6, /overplot)
      fc = d['flask_avg_co2_age'] & fs = d['flask_avg_sf6_age']
      p = errorplot(fc[*,0], alt2, fc[*,1], replicate(0,nz2), symbol='tu', /sym_filled, sym_size=1.5, thick=3, color='green', linestyle=2, $
        errorbar_capsize=0, /overplot)
      p = errorplot(fs[*,0], alt2, fs[*,1], replicate(0,nz2), symbol='o', /sym_filled, thick=3, color='green', linestyle=1, $
        errorbar_capsize=0, /overplot)
    endelse
    t = text(0.13, 0.85, 'Solid lines = In situ balloon', font_size=11, /norm)
    t = text(0.13, 0.81, 'Tropics', color='red', font_size=11, /norm)
    t = text(0.25, 0.81, 'Midlat', color='blue', font_size=11, /norm)
    t = text(0.35, 0.81, 'Vortex', color='purple', font_size=11, /norm)
    if pn.k eq 'age' then begin
      t = text(0.65, 0.17, 'Flask balloon ', color='lime green', font_size=11, /norm)
      t = text(0.56, 0.13, 'Circles=SF$_6$, Triangles=CO$_2$', color='lime green', font_size=11, /norm)
    endif else t = text(0.13, 0.76, 'Symbols = Flask balloon', color='lime green', font_size=11, /norm)
    fig_save, p, 'Extended_Data_Figure_2'+pn.p
  endforeach

end
