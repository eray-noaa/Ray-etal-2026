;+
; Supplementary Figure 1: mean age vs normalized N2O (a) and CH4 (b) for the
; simple and age-convolution normalizations, 1990s and 2020s in situ data.
; Data: <output>/figure_data/supp_fig1.nc.
;-
pro supp_fig1

  compile_opt idl2

  d = fig_data('supp_fig1')
  ng = d['norm_grid'] & i = lindgen(n_elements(ng)-13) + 12   ; bins plotted in the paper
  cols = ['purple','pink','red','orange']
  labs = ['1990s simple norm','1990s age conv norm','2020s simple norm','2020s age conv norm']

  p = plot(indgen(2), /nodata, /buffer, yrange=[1,0], xrange=[0,7.3], ytitle='Normalized N$_2$O', xtitle='Mean Age (years)', $
    title='a  Mean Age vs. N$_2$O', font_size=11)
  keys = ['n2o_age_1990s_simple','n2o_age_1990s_age_conv','n2o_age_2020s_simple','n2o_age_2020s_age_conv']
  unc = ['n2o_age_1990s_uncert','n2o_age_1990s_uncert','n2o_age_2020s_uncert','n2o_age_2020s_uncert']
  for k = 0, 3 do begin
    a = d[keys[k]] & u = d[unc[k]]
    p = errorplot(a[i], ng[i], u[i], replicate(0,n_elements(i)), symbol='o', color=cols[k], sym_size=0.75, /sym_filled, linestyle=6, $
      errorbar_capsize=0, /overplot)
    t = text(0.18, 0.8-0.04*k, labs[k], color=cols[k], font_size=11, /norm)
  endfor
  fig_save, p, 'Supplementary_Figure_1a'

  p = plot(indgen(2), /nodata, /buffer, yrange=[1,0.3], xrange=[0,7], ytitle='Normalized CH$_4$', xtitle='Mean Age (years)', $
    title='b  Mean Age vs. CH$_4$', font_size=11)
  keys = ['ch4_age_1990s_simple','ch4_age_1990s_age_conv','ch4_age_2020s_simple','ch4_age_2020s_age_conv']
  for k = 0, 3 do begin
    a = d[keys[k]]
    p = errorplot(a[i,0], ng[i], a[i,1], replicate(0,n_elements(i)), symbol='o', color=cols[k], sym_size=0.75, /sym_filled, linestyle=6, $
      errorbar_capsize=0, /overplot)
    t = text(0.18, 0.8-0.04*k, labs[k], color=cols[k], font_size=11, /norm)
  endfor
  fig_save, p, 'Supplementary_Figure_1b'

end
