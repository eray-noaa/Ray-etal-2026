function convol_extra, arr, kern, NORM=norm, EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                       EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                       EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast

  compile_opt idl2, strictarrsubs

  if total(keyword_set(EDGE_TRUNCATE) + keyword_set(EDGE_MIRROR) + keyword_set(EDGE_WRAP) + $
           keyword_set(EDGE_ZERO) + keyword_set(EDGE_NAN) + keyword_set(EDGE_FCAST) + $
           keyword_set(EDGE_PFCAST)) gt 1 then $
              message, 'can only set one EDGE_* keyword'

  n_kern = n_elements(kern)

  if ~keyword_set(EDGE_FCAST) then begin
     arr_lp = convol(arr, kern, NORM=norm, EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                     EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero)
  endif else begin
     n_arr = n_elements(arr)
     if n_elements(fcast_order) eq 0 then $
        fcast_order = min([n_arr-1, max([10L, long(0.05 * n_arr)])])
     fcast_left = n_kern/2
     fcast_right = n_kern/2
     arr_fcast = make_array(n_arr+fcast_left+fcast_right, TYPE=size(arr,/TYPE))
     arr_fcast[0] = ts_fcast(arr, fcast_order, fcast_left, /BACKCAST)
     arr_fcast[fcast_left] = arr
     arr_fcast[fcast_left+n_arr] = ts_fcast(arr, fcast_order, fcast_right)
     arr_lp_fcast = convol(arr_fcast, kern, NORM=norm)
     arr_lp = arr_lp_fcast[fcast_left:n_arr+fcast_left-1]
  endelse

  if keyword_set(EDGE_PFCAST) then begin
     fcast_left = n_kern/2
     fcast_right = (n_kern-1)/2
     arr_lp = fcast_array(arr_lp, NLEFT=fcast_left, NRIGHT=fcast_right, FCAST_ORDER=fcast_order)
  endif

  if keyword_set(EDGE_NAN) then begin
     if n_kern gt 1 then $
        arr_lp[0:n_kern/2-1] = isa(arr_lp, 'float') ? !values.f_nan : !values.d_nan
     if n_kern gt 2 then $
        arr_lp[-(n_kern-1)/2:-1] = isa(arr_lp, 'float') ? !values.f_nan : !values.d_nan
  endif

  return, arr_lp
end

function lowpass_cfc_1d, arr, BOXCAR=boxcar, $
                         EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                         EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                         INTERP_TO_X_INDICES=interp_to_x_indices, $
                         LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                         SG_WIDTH=sg_width, SG_DEGREE=sg_degree, $
                         SG_DERIV_ORDER=sg_deriv_order, SG_DELTA_T=sg_delta_t, $
                         SG_NLEFT=sg_nleft, SG_NRIGHT=sg_nright, $
                         EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast

  ;; Performs a boxcar smooth (like smooth) but allows even boxcar sizes.
  ;;
  ;; ARR: input array (must be 1-dimensional)
  ;; BOXCAR: width of the boxcar smoother (defaults to 12)
  ;; EDGE_TRUNCATE: pass EDGE_TRUNCATE to convol
  ;; EDGE_MIRROR: pass EDGE_MIRROR to convol
  ;; EDGE_WRAP: pass EDGE_WRAP to convol
  ;; EDGE_ZERO: pass EDGE_ZERO to convol
  ;; EDGE_NAN: Puts NaNs at edges where box extends beyond end
  ;; INTERP_TO_X_INDICES: Interpolates arr back to input grid if boxcar is even
  ;;                      CAUTION: Should set one of the EDGE_* keywords (not enforced)
  ;; LSQUADRATIC, QUADRATIC, SPLINE: passed to interpol if interp_to_x_indices set

  compile_opt idl2, strictarrsubs

  ndim =(size(arr))[0]
  if ndim eq 0 or ndim gt 1 then message, 'arr must be a 1 dim array'

  sg_input = total(n_elements(SG_WIDTH) gt 0 or n_elements(SG_NLEFT) gt 0 or n_elements(SG_NRIGHT) gt 0 or $
                   n_elements(SG_DEGREE) gt 0 or n_elements(SG_DERIV_ORDER) gt 0 or n_elements(SG_DELTA_T) gt 0)

  boxcar_input = total(n_elements(boxcar) gt 0 or n_elements(INTERP_TO_X_INDICES) gt 0 or $
                       n_elements(LSQUADRATIC) gt 0 or n_elements(QUADRATIC) gt 0 or n_elements(SPLINE) gt 0)

  if sg_input gt 0 and boxcar_input gt 0 then message, 'Cannot input both BOXCAR and SG inputs'

  if sg_input gt 0 then begin
     arr_lp = sg_smooth(arr, sg_width, DEGREE=sg_degree, EDGE_TRUNCATE=edge_truncate, $
                        EDGE_MIRROR=edge_mirror, EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, $
                        DERIV_ORDER=sg_deriv_order, DELTA_T=sg_delta_t, $
                        NLEFT=sg_nleft, NRIGHT=sg_nright, $
                        EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast)
     return, arr_lp
  endif

  if n_elements(boxcar) eq 0 then boxcar = 12
  if boxcar eq 1 then return, arr

  if n_elements(boxcar) gt 0 then begin
     n = boxcar
     kern = replicate(1.0, n)
     arr_lp = convol_extra(arr, kern, /NORM, EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                           EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                           EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast)
     if keyword_set(INTERP_TO_X_INDICES) and 2*(n/2) eq n then begin
        n_arr = n_elements(arr)
        x = isa(arr_lp, 'float') ? findgen(n_arr) : dindgen(n_arr)
        x_lp = convol_extra(x, kern, /NORM, EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                            EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=0, $
                            EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast)
        arr_interp = interpol(arr_lp, x_lp, x, $
                              LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline)
        arr_lp =arr_interp
     endif
  endif

  return, arr_lp
end

function lowpass_cfc_2call_1d, arr, BOXCAR=boxcar, $
                               EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                               EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                               INTERP_TO_X_INDICES=interp_to_x_indices, $
                               LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                               SG_WIDTH=sg_width, SG_DEGREE=sg_degree, $
                               EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast

  compile_opt idl2, strictarrsubs

  if n_elements(SG_WIDTH) gt 0 or n_elements(SG_DEGREE) gt 0 then $
     message, '2CALL is incompatable with SG_WIDTH & SG_DEGREE'

  if n_elements(boxcar) eq 0 then boxcar = 12
  if boxcar eq 1 then return, arr
  if 2*(long(boxcar)/2) ne boxcar then message, 'not setup for odd boxcar size'
  n = boxcar

  edge_total = keyword_set(EDGE_TRUNCATE) + $
               keyword_set(EDGE_MIRROR) + $
               keyword_set(EDGE_WRAP) + $
               keyword_set(EDGE_ZERO) + $
               keyword_set(EDGE_NAN)

  if edge_total eq 0 then message, 'probably want to set of of the EDGE_* keywords', /INFO

  ;; only apply EDGE_NAN at the end, lose many points if EDGE_TRUNCATE=0 so set
  arr_lp = lowpass_cfc_1d(arr, BOXCAR=boxcar, $
                          EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                          EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=0, $
                          INTERP_TO_X_INDICES=interp_to_x_indices, $
                          LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                          EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast)

  arr_lp = [arr_lp, arr_lp[-1]]

  ;; only apply EDGE_NAN at the end, lose many points if EDGE_TRUNCATE=0 so set
  arr_lp = lowpass_cfc_1d(arr_lp, BOXCAR=boxcar, $
                          EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                          EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                          INTERP_TO_X_INDICES=interp_to_x_indices, $
                          LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                          EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast)

  arr_lp = [arr_lp[1:-1]]

  if keyword_set(EDGE_NAN) then begin
     arr_lp[0:n/2-1] = isa(arr_lp, 'float') ? !values.f_nan : !values.d_nan
     arr_lp[(-(n-1)/2-1):-1] = isa(arr_lp, 'float') ? !values.f_nan : !values.d_nan
  endif

  return, arr_lp
end

function lowpass_cfc, arr, BOXCAR=boxcar, EDGE_TRUNCATE=edge_truncate, $
                      EDGE_MIRROR=edge_mirror, EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, $
                      EDGE_NAN=edge_nan, REMOVE_EDGE=remove_edge, $
                      INTERP_TO_X_INDICES=interp_to_x_indices, $
                      LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                      TWICE_CALL=twice_call, DIM=dim, $
                      SG_WIDTH=sg_width, SG_DEGREE=sg_degree, $
                      SG_DERIV_ORDER=sg_deriv_order, SG_DELTA_T=sg_delta_t, $
                      SG_NLEFT=sg_nleft, SG_NRIGHT=sg_nright, $
                      EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, EDGE_PFCAST=edge_pfcast

  ;; Performs a boxcar smooth (like smooth) but allows even boxcar sizes.
  ;;
  ;; ARR: input array (can be n-dimensional if dim is set)
  ;; DIM: dimension to apply lowpass over
  ;; BOXCAR: width of the boxcar smoother
  ;; EDGE_TRUNCATE: pass EDGE_TRUNCATE to convol
  ;; EDGE_MIRROR: pass EDGE_MIRROR to convol
  ;; EDGE_WRAP: pass EDGE_WRAP to convol
  ;; EDGE_ZERO: pass EDGE_ZERO to convol
  ;; EDGE_NAN: Puts NaNs at edges where box extends beyond end
  ;; REMOVE_EDGE: Puts NaNs at edges where box extends beyond end (obsolete,
  ;;              use EDGE_NAN instead)
  ;; INTERP_TO_X_INDICES: Interpolates arr back to input grid if boxcar is even
  ;; LSQUADRATIC, QUADRATIC, SPLINE: passed to interpol if interp_to_x_indices set
  ;; EDGE_FCAST: use ts_fcast to predect points beyond edges
  ;; FCAST_ORDER: number of points to use in ts_fcast (default: 5% of arr
  ;; length but not less than 10, except when there are less than 10 point in arr)
  ;; EDGE_PFCAST: like EDGE_FCAST but done after smoothing

  compile_opt idl2, strictarrsubs

  ndim =(size(arr))[0]
  if ndim eq 0 then message, 'arr must be a 1 dim array'

  if keyword_set(TWICE_CALL) then func_1d = 'lowpass_cfc_2call_1d' else func_1d = 'lowpass_cfc_1d'

  if (arg_present(REMOVE_EDGE) and arg_present(EDGE_NAN)) and $
     (keyword_set(REMOVE_EDGE) ne keyword_set(EDGE_NAN)) then $
     message, 'REMOVE_EDGE and EDGE_NAN settings inconsistent'

  if keyword_set(REMOVE_EDGE) and ~keyword_set(EDGE_NAN) then edge_nan = 1

  if ndim eq 1 then begin
     if n_elements(dim) gt 0 then $
        if dim ne 1 then message, 'DIM must be 1 or unset if ARR is 1-dimensional'
     arr_lp = call_function(func_1d, arr, BOXCAR=boxcar,  $
                            EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                            EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                            INTERP_TO_X_INDICES=interp_to_x_indices, $
                            LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                            SG_WIDTH=sg_width, SG_DEGREE=sg_degree, $
                            SG_DERIV_ORDER=sg_deriv_order, SG_DELTA_T=sg_delta_t, $
                            SG_NLEFT=sg_nleft, SG_NRIGHT=sg_nright, $
                            EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, $
                            EDGE_PFCAST=edge_pfcast)
  endif else begin
     if n_elements(dim) eq 0 then $
        message, 'DIM must be set if ARR is greater than 1-dimensional'
     arr_lp = rwp_apply_func_to_dim(func_1d, arr, dim, BOXCAR=boxcar, $
                                    EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                                    EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero, EDGE_NAN=edge_nan, $
                                    INTERP_TO_X_INDICES=interp_to_x_indices, $
                                    LSQUADRATIC=lsquadratic, QUADRATIC=quadratic, SPLINE=spline, $
                                    SG_WIDTH=sg_width, SG_DEGREE=sg_degree, $
                                    SG_DERIV_ORDER=sg_deriv_order, SG_DELTA_T=sg_delta_t, $
                                    SG_NLEFT=sg_nleft, SG_NRIGHT=sg_nright, $
                                    EDGE_FCAST=edge_fcast, FCAST_ORDER=fcast_order, $
                                    EDGE_PFCAST=edge_pfcast)
  endelse

  return, arr_lp
end
