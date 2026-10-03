;+
; Supplementary Figure 2: NH mean age trends vs normalized N2O (a) and CH4 (b)
; for the simple and age-convolution normalizations.
; Data: <output>/figure_data/supp_fig2.nc. Trends are plotted per decade.
;-
pro supp_fig2

  compile_opt idl2

  d = fig_data('supp_fig2')
  ng = d['norm_grid'] & nn = n_elements(ng)
  panels = list({tr:'n2o', lab:'Normalized N$_2$O', yr:[1.02,0.2], xr:[-0.31,0.11], i0:6, p:'a'}, $
    {tr:'ch4', lab:'Normalized CH$_4$', yr:[1.01,0.45], xr:[-0.3,0.15], i0:5, p:'b'})
  foreach pn, panels do begin
    i = lindgen(nn-pn.i0-1) + pn.i0
    p = plot(indgen(2), /nodata, /buffer, xrange=pn.xr, yrange=pn.yr, ytitle=pn.lab, xtitle='Years/Decade', title=pn.p+'  NH Mean Age Trends', $
      font_size=11, margin=[0.13,0.09,0.12,0.07], dimensions=[450,500])
    p = plot([0,0], [1,0.1], linestyle=2, /overplot)
    foreach nm, ['simple','age_conv'], k do begin
      tr = d['trend_on_'+pn.tr+'_'+nm] & sg = d['trend_on_'+pn.tr+'_'+nm+'_sigma']
      col = (k eq 0) ? 'blue' : 'sky blue'
      p = errorplot(10.*tr[i], ng[i], 10.*sg[i], replicate(!values.f_nan,n_elements(i)), symbol='o', sym_size=1.5, /sym_filled, color=col, $
        errorbar_capsize=0, /overplot)
      t = text(0.19, 0.86-0.04*k, (k eq 0) ? 'Simple norm' : 'Age conv norm', color=col, font_size=10, /norm)
    endforeach
    fig_save, p, 'Supplementary_Figure_2'+pn.p
  endforeach

end
