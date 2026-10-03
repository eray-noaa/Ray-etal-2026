;+
; Supplementary Figure 4: 1990s gridded (a) mean age and (b) normalized N2O vs
; equivalent latitude and altitude; latitude distributions at 20 km of
; (c) mean age and (d) N2O compared with CCMI-2022 refD1.
; Data: <output>/figure_data/supp_fig4.nc.
;-
pro supp_fig4_grid, d, key, scale, offset, cbtitle, ticknames, tickvalues, title, name

  compile_opt idl2

  g = d[key] & lat = d['equiv_lat'] & alt = d['alt'] & tp = d['tropopause_alt_on_grid']
  nyg = n_elements(lat) & nz = n_elements(alt)
  loadct, 39, rgb_table=rgb, /silent
  p = plot(indgen(2), /nodata, /buffer, xrange=[-90,90], xtickinterval=20, xminor=1, yrange=[9,25], position=[0.1,0.1,0.89,0.9], $
    xtitle='Equivalent Latitude', ytitle='Altitude (km)', font_size=11, xticklen=0.04, yticklen=0.03, dimensions=[700,450], title=title)
  for z = 0, nz-1 do begin
    ; Boxes are drawn as three overlapping squares; full size where >= 3 seasons were sampled, small where 1-2.
    sel = list()
    ii = where(finite(g[0,0:-3,z]) and g[2,0:-3,z] ge 3 and tp[0:-3]-1.5 le alt[z], nii)
    if nii gt 0 then sel.add, {ii:ii, size:1.6, dx:0.6}
    ii = where(finite(g[0,0:-2,z]) and g[2,0:-2,z] ge 1 and g[2,0:-2,z] le 2 and tp[0:-2]-1.5 le alt[z], nii)
    if nii gt 0 then sel.add, {ii:ii, size:0.8, dx:0.5}
    ii = where(finite(g[0,-2:-1,z]) and g[2,-2:-1,z] ge 1, nii) + nyg-2
    if nii gt 0 then sel.add, {ii:ii, size:0.8, dx:0.5}
    if alt[z] le 21 then begin
      ii = where(finite(g[0,-2:-1,z]) and g[2,-2:-1,z] ge 3, nii) + nyg-2
      if nii gt 0 then sel.add, {ii:ii, size:1.6, dx:0.5}
    endif
    foreach s, sel do begin
      vcols = ((reform(g[0,s.ii,z]) - offset)*scale > 0.) < 240.
      foreach dx, [0., -s.dx, s.dx] do p = plot(lat[s.ii]+dx, replicate(alt[z],n_elements(s.ii)), symbol='s', /sym_filled, $
        linestyle=6, rgb_table=39, vert_colors=vcols, sym_size=s.size, /overplot)
    endforeach
  endfor
  p = plot(d['merra2_lat'], d['merra2_tropopause_alt'], thick=2, color='magenta', /overplot)
  cb = colorbar(position=[0.905,0.2,0.925,0.75], tickname=ticknames, tickvalues=tickvalues, rgb_table=rgb, border=1, orientation=1, $
    /textpos, /tickdir, font_size=10, title=cbtitle)
  fig_save, p, name

end

pro supp_fig4

  compile_opt idl2

  d = fig_data('supp_fig4')
  supp_fig4_grid, d, 'age_1990s', 46., 0., 'years', string(findgen(12)*0.5, format='(f3.1)'), findgen(12)*23., $
    'a  Mean Ages 1990s', 'Supplementary_Figure_4a'
  supp_fig4_grid, d, 'n2o_norm_1990s', 416., 0.4, 'Normalized N$_2$O', string(findgen(7)*0.1+0.4, format='(f3.1)'), findgen(7)*41.6, $
    'b  Normalized N$_2$O 1990s', 'Supplementary_Figure_4b'

  lat = d['equiv_lat'] & lc = d['ccmi_lat']
  ca = d['ccmi_age_20km_1990s'] & cn = d['ccmi_n2o_norm_20km_1990s'] & nm = (size(ca,/dim))[1]
  ga = d['age_1990s'] & gn = d['n2o_norm_1990s']
  iz = (where(d['alt'] eq d['alt_20km']))[0]

  p = plot(indgen(2), /nodata, /buffer, xrange=[-90,90], xtickinterval=20, xtitle='Latitude', yrange=[0,6], ytitle='Mean Age (years)', $
    xticklen=0.02, yticklen=0.02, title='c  Mean age at 20 km  1990s', margin=[0.08,0.1,0.04,0.08], font_size=11, dimensions=[700,400])
  for m = 0, nm-1 do p = plot(lc, ca[*,m], color='lime green', thick=1, /overplot)
  p = plot(lat[0:-2], reform(ga[0,0:-2,iz]), color='blue', thick=3, /overplot)
  t = text(0.12, 0.28, 'In situ', color='blue', font_size=11, /norm)
  t = text(0.12, 0.23, 'CCMI-2022', color='lime green', font_size=11, /norm)
  fig_save, p, 'Supplementary_Figure_4c'

  p = plot(indgen(2), /nodata, /buffer, xrange=[-90,90], xtickinterval=20, xtitle='Latitude', yrange=[0.4,1], ytitle='Normalized N$_2$O', $
    xticklen=0.02, yticklen=0.02, title='d  N$_2$O at 20 km  1990s', margin=[0.09,0.1,0.03,0.08], font_size=11, dimensions=[700,400])
  for m = 3, nm-1 do if m ne 15 then p = plot(lc, cn[*,m], color='lime green', thick=1, /overplot)
  p = plot(lat[0:-2], reform(gn[0,0:-2,iz]), color='blue', thick=3, /overplot)
  t = text(0.12, 0.82, 'In situ', color='blue', font_size=11, /norm)
  t = text(0.12, 0.77, 'CCMI-2022', color='lime green', font_size=11, /norm)
  fig_save, p, 'Supplementary_Figure_4d'

end
