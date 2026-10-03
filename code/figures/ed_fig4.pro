;+
; Extended Data Figure 4: (a) sampling altitudes by platform and (b) 24-35 km
; average mean age time series. Data: <output>/figure_data/ed_fig4.nc.
;-
pro ed_fig4

  compile_opt idl2

  d = fig_data('ed_fig4')
  yrs = d['years'] & st = d['balloon_sample_type'] & mx = d['balloon_max_alt']

  p = plot(indgen(2), /nodata, /buffer, xrange=[1970,2030], yrange=[14.5,36], ytitle='Altitude (km)', xtitle='Years', $
    title='a  Sampling Altitude', font_size=11, margin=[0.08,0.1,0.04,0.08], dimensions=[700,400])
  p = plot([1970,2030], [24,24], linestyle=2, /overplot)
  cols = ['lime green','purple','magenta']
  for k = 0, 2 do begin
    ti = where(finite(mx) and st eq k, nti)
    for i = 0, nti-1 do p = plot(replicate(yrs[ti[i]],2), [15,mx[ti[i]]], color=cols[k], thick=2, /overplot)
  endfor
  ya = d['aircraft_years'] & ma = d['aircraft_max_alt']
  ti = where(finite(ma), nti)
  for i = 0, nti-1 do p = plot(replicate(ya[ti[i]],2), [15,ma[ti[i]]], color='blue', thick=2, /overplot)
  t = text(0.93, 0.83, 'Flask', color='lime green', font_size=11, alignment=1, /norm)
  t = text(0.93, 0.79, 'In situ balloon', color='purple', font_size=11, alignment=1, /norm)
  t = text(0.93, 0.75, 'In situ aircraft', color='blue', font_size=11, alignment=1, /norm)
  t = text(0.93, 0.71, 'AirCore', color='magenta', font_size=11, alignment=1, /norm)
  fig_save, p, 'Extended_Data_Figure_4a'

  fl = d['age_24_35km_with_flask'] & nf = d['age_24_35km_no_flask']
  p = plot(indgen(2), /nodata, /buffer, xrange=[1970,2030], yrange=[3,7], ytitle='Mean Age (years)', xtitle='Years', $
    title='b  Mean Age in 24-35km Altitude Range', font_size=11, margin=[0.08,0.1,0.04,0.08], dimensions=[700,400])
  ti = where(st eq 0 or st eq 3, nti)
  p = errorplot(yrs[ti], fl[ti,0], replicate(0,nti), fl[ti,1], color='lime green', symbol='o', linestyle=6, /sym_filled, errorbar_capsize=0, /overplot)
  si = where(st eq 2, nsi)
  p = errorplot(yrs[si], nf[si,0], replicate(0,nsi), nf[si,1], color='magenta', symbol='s', /sym_filled, linestyle=6, errorbar_capsize=0, /overplot)
  si = where(st eq 1, nsi)
  p = errorplot(yrs[si], nf[si,0], replicate(0,nsi), nf[si,1], color='purple', symbol='s', /sym_filled, linestyle=6, errorbar_capsize=0, /overplot)
  f = d['fit_no_flask'] & g = d['fit_with_flask']
  p = plot([1995,2026], f[0] + f[1]*[1995,2026], thick=2, linestyle=2, color='purple', /overplot)
  p = plot([1975,2026], g[0] + g[1]*[1975,2026], thick=2, linestyle=2, color='green', /overplot)
  t = text(0.13, 0.82, 'In situ balloon', color='purple', font_size=11, /norm)
  t = text(0.13, 0.78, 'AirCore', color='magenta', font_size=11, /norm)
  t = text(0.13, 0.74, 'Flask', color='lime green', font_size=11, /norm)
  fig_save, p, 'Extended_Data_Figure_4b'

end
