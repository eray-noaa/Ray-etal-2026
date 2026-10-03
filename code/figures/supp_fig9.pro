;+
; Supplementary Figure 9: surface CO2 time series and mean ages for the
; 6 August 1994 ASHOE-MAESA ER-2 flight. Data: <output>/figure_data/supp_fig9.nc.
;-
pro supp_fig9

  compile_opt idl2

  d = fig_data('supp_fig9')
  t0 = d['flight_date_float']
  p = plot(indgen(2), /nodata, /buffer, xrange=[1988.5,1995.5], xtickinterval=1, xminor=3, yrange=[350.2,360], ytitle='CO$_2$ (ppmv)', $
    xtitle='Year', font_size=11, title='Surface and ASHOE  '+d['flight_date']+' ER-2 CO$_2$', dimensions=[700,400])
  p = plot(d['years_surface'], d['co2_mlo_smo'], thick=2, color='lime green', /overplot)
  p = plot(d['years_surface'], d['co2_mlo_smo_smooth'], thick=2, color='green', /overplot)
  p = plot(d['steady_growth_time'], d['steady_growth_co2'], thick=2, color='black', /overplot)
  p = plot(t0 - d['flight_co2_age_archive'], d['flight_co2_archive_adj'], symbol='o', color='blue', /overplot)
  p = plot(t0 - d['flight_co2_lag_age'], d['flight_co2_adj'], symbol='o', color='red', linestyle=6, /overplot)
  n = n_elements(d['flight_co2_adj'])
  p = plot(replicate(t0,n), d['flight_co2_adj'], symbol='o', sym_size=0.75, linestyle=6, color='purple', /overplot)
  p = plot(replicate(t0,n), d['flight_co2_archive_adj'], symbol='o', sym_size=0.75, linestyle=6, color='purple', /overplot)
  t = text(0.16, 0.8, 'MLO-SMO Avg w/2 mon delay', color='lime green', font_size=11, /norm)
  t = text(0.16, 0.76, 'MLO-SMO Avg smooth', color='green', font_size=11, /norm)
  t = text(0.16, 0.72, 'MLO-SMO Avg const growth rate', color='black', font_size=11, /norm)
  t = text(0.16, 0.68, 'Adj CO$_2$ and lag age', color='red', font_size=11, /norm)
  t = text(0.16, 0.64, 'Archive Adj CO$_2$ and age', color='blue', font_size=11, /norm)
  t = text(0.16, 0.6, 'Flight data', color='purple', font_size=11, /norm)
  fig_save, p, 'Supplementary_Figure_9'

end
