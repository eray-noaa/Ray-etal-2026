function sg_smooth_1dim, arr, width, $
                         DEGREE=degree, $
                         NLEFT=nleft_in, $
                         NRIGHT=nright_in, $
                         DERIV_ORDER=deriv_order, $
                         DELTA_T=delta_t, $
                         DOUBLE=double, $
                         EDGE_TRUNCATE=edge_truncate, $
                         EDGE_MIRROR=edge_mirror, $
                         EDGE_WRAP=edge_wrap, $
                         EDGE_ZERO=edge_zero, $
                         EDGE_NAN=edge_nan, $
                         EDGE_FCAST=edge_fcast, $
                         FCAST_ORDER=fcast_order, $
                         EDGE_PFCAST=edge_pfcast

  ;; if ORDER > 0 then need to pass in the sampling interval DELTA_T
  ;; EDGE_FCAST: use ts_fcast to predect points beyond edges
  ;; FCAST_ORDER: number of points to use in ts_fcast (default: 5% of arr
  ;; length but not less than 10, except when there are less than 10 point in arr)
  ;; EDGE_PFCAST: like EDGE_FCAST but done after smoothing

  compile_opt idl2, strictarrsubs

  if n_elements(size(arr, /DIM)) ne 1 then message, 'ARR must be 1-dimensional'
  if n_elements(nleft_in) gt 0 then nleft = nleft_in
  if n_elements(nright_in) gt 0 then nright = nright_in

  flag1 = n_elements(width) gt 0
  flag2 = n_elements(nleft) gt 0
  flag3 = n_elements(nright) gt 0

  ;print, flag1, flag2, flag3
  if flag1 + flag2 + flag3 eq 3 then $
     if width ne nleft + nright + 1 then message, 'width, nleft, & nright are incompatable'

  if flag1 + flag2 + flag3 eq 0 then $
     message, 'must input width and/or nleft & nright'

  if flag1 eq 1 and flag2 + flag3 eq 0 then begin
     nleft = floor(width/2)
     nright = floor(width/2)
  endif else if flag1 eq 1 and flag2 + flag3 eq 1 then begin
     if flag2 eq 0 then nleft = width - nright - 1
     if flag3 eq 0 then nright = width - nleft - 1
  endif

  if n_elements(degree) eq 0 then degree = 2
  if n_elements(deriv_order) eq 0 then deriv_order = 0

  if ~keyword_set(double) and n_elements(double) gt 0 and isa(arr, 'double') then $
     message, 'ARR input is type double but DOUBLE explictly set to 0'

  if total(keyword_set(EDGE_TRUNCATE) + keyword_set(EDGE_MIRROR) + keyword_set(EDGE_WRAP) + $
           keyword_set(EDGE_ZERO) + keyword_set(EDGE_NAN) + keyword_set(EDGE_FCAST) + $
           keyword_set(EDGE_PFCAST)) gt 1 then $
              message, 'can only set one EDGE_* keyword'

  ;; needed since sagvol sets double even if not input
  if keyword_set(DOUBLE) or isa(arr, 'double') then $
     sg_filter = savgol(nleft, nright, deriv_order, degree, DOUBLE=1) $
  else $
     sg_filter = savgol(nleft, nright, deriv_order, degree)

  if deriv_order gt 0 then begin
     if n_elements(delta_t) eq 0 then $
        message, 'Must pass in sampling interval DELTA_T if ORDER > 0' $
     else $
        sg_filter /= factorial(deriv_order) / (delta_t^deriv_order)
  endif

  n_filter = n_elements(sg_filter)

  if ~keyword_set(EDGE_FCAST) then begin
     arr_lp = convol(arr, sg_filter, EDGE_TRUNCATE=edge_truncate, EDGE_MIRROR=edge_mirror, $
                     EDGE_WRAP=edge_wrap, EDGE_ZERO=edge_zero)
  endif else begin
     n_arr = n_elements(arr)
     if n_elements(fcast_order) eq 0 then $
        fcast_order = min([n_arr-1, max([10L, long(0.05 * n_arr)])])
     fcast_left = n_filter/2
     fcast_right = n_filter/2
     arr_fcast = make_array(n_arr+fcast_left+fcast_right, TYPE=size(arr,/TYPE))
     arr_fcast[0] = ts_fcast(arr, fcast_order, fcast_left, /BACKCAST, DOUBLE=double)
     arr_fcast[fcast_left] = arr
     arr_fcast[fcast_left+n_arr] = ts_fcast(arr, fcast_order, fcast_right, DOUBLE=double)
     arr_lp_fcast = convol(arr_fcast, sg_filter)
     arr_lp = arr_lp_fcast[fcast_left:n_arr+fcast_left-1]
  endelse

  if keyword_set(EDGE_PFCAST) then begin
     fcast_left = n_filter/2
     fcast_right = (n_filter-1)/2
     arr_lp = fcast_array(arr_lp, NLEFT=fcast_left, NRIGHT=fcast_right, FCAST_ORDER=fcast_order)
  endif

  if keyword_set(EDGE_NAN) then begin
     if ~isa(arr, /FLOAT) then message, 'ARR must be floating point if EDGE_NAN is set'
     if n_filter gt 1 then $
        arr_lp[0:n_filter/2-1] = isa(arr_lp, 'float') ? !values.f_nan : !values.d_nan
     if n_filter gt 2 then $
        arr_lp[-(n_filter-1)/2:-1] = isa(arr_lp, 'float') ? !values.f_nan : !values.d_nan
  endif

  return, arr_lp
end

function sg_smooth, arr, width, $
                    DIM=dim, $
                    DEGREE=degree, $
                    NLEFT=nleft, $
                    NRIGHT=nright, $
                    DERIV_ORDER=deriv_order, $
                    DELTA_T=delta_t, $
                    DOUBLE=double, $
                    EDGE_TRUNCATE=edge_truncate, $
                    EDGE_MIRROR=edge_mirror, $
                    EDGE_WRAP=edge_wrap, $
                    EDGE_ZERO=edge_zero, $
                    EDGE_NAN=edge_nan, $
                    EDGE_FCAST=edge_fcast, $
                    FCAST_ORDER=fcast_order, $
                    EDGE_PFCAST=edge_pfcast

  ;; if DERIV_ORDER > 0 then need to pass in the sampling interval DELTA_T
  ;; EDGE_FCAST: use ts_fcast to predect points beyond edges
  ;; FCAST_ORDER: number of points to use in ts_fcast (default: 5% of arr
  ;; length but not less than 10, except when there are less than 10 point in arr)
  ;; EDGE_PFCAST: like EDGE_FCAST but done after smoothing

  compile_opt idl2, strictarrsubs

  ndim = n_elements(size(arr,/DIM))

  if n_elements(dim) gt 0 then begin
     if dim lt 1 or dim gt ndim then message, 'DIM is not present in ARR'
  endif

  if ndim eq 1 then begin
     arr_lp = sg_smooth_1dim(arr, width, $
                             DEGREE=degree, $
                             NLEFT=nleft, $
                             NRIGHT=nright, $
                             DERIV_ORDER=deriv_order, $
                             DELTA_T=delta_t, $
                             DOUBLE=double, $
                             EDGE_TRUNCATE=edge_truncate, $
                             EDGE_MIRROR=edge_mirror, $
                             EDGE_WRAP=edge_wrap, $
                             EDGE_ZERO=edge_zero, $
                             EDGE_NAN=edge_nan, $
                             EDGE_FCAST=edge_fcast, $
                             FCAST_ORDER=fcast_order, $
                             EDGE_PFCAST=edge_pfcast)
  endif else begin
     if n_elements(dim) eq 0 then message, 'must set DIM if ARR is greater than 1 dimensional'
     arr_lp = rwp_apply_func_to_dim('sg_smooth_1dim',arr, dim, width, $
                                    DEGREE=degree, $
                                    NLEFT=nleft, $
                                    NRIGHT=nright, $
                                    DERIV_ORDER=deriv_order, $
                                    DELTA_T=delta_t, $
                                    DOUBLE=double, $
                                    EDGE_TRUNCATE=edge_truncate, $
                                    EDGE_MIRROR=edge_mirror, $
                                    EDGE_WRAP=edge_wrap, $
                                    EDGE_ZERO=edge_zero, $
                                    EDGE_NAN=edge_nan, $
                                    EDGE_FCAST=edge_fcast, $
                                    FCAST_ORDER=fcast_order, $
                                    EDGE_PFCAST=edge_pfcast)
  endelse

  return, arr_lp
end
