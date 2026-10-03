;+
; MAKE_FIGURES
;
; Draws every figure from the NetCDF files in <output>/figure_data/ and saves
; PDF and PNG versions in <output>/figures/. A figure whose data file is
; missing is skipped with a message.
;-
pro make_figures

  compile_opt idl2

  routines = ['fig1','fig2','fig3','fig4', $
    'ed_fig1','ed_fig2','ed_fig3','ed_fig4','ed_fig5','ed_fig6', $
    'supp_fig1','supp_fig2','supp_fig3','supp_fig4','supp_fig5','supp_fig6','supp_fig7','supp_fig8','supp_fig9','supp_fig10']
  foreach r, routines do begin
    catch, err
    if err ne 0 then begin
      catch, /cancel
      print, 'Skipped ', r, ': ', !error_state.msg
      continue
    endif
    call_procedure, r
    catch, /cancel
  endforeach

end
