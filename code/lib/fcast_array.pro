function fcast_array_1dim, arr, npoints, $
                           NLEFT=nleft_in, NRIGHT=nright_in, $
                           EXPAND_ARRAY=expand_array, $
                           DOUBLE=double, $
                           FCAST_ORDER=fcast_order

  compile_opt idl2, strictarrsubs

  n = n_elements(arr)
  if n_elements(size(arr, /DIM)) ne 1 then message, 'ARR must be 1-dimensional'
  if n_elements(nleft_in) gt 0 then nleft = nleft_in
  if n_elements(nright_in) gt 0 then nright = nright_in

  flag1 = n_elements(npoints) gt 0
  flag2 = n_elements(nleft) gt 0
  flag3 = n_elements(nright) gt 0

  ;print, flag1, flag2, flag3
  if flag1 + flag2 + flag3 eq 3 then $
     if npoints ne nleft + nright then message, 'npoints, nleft, & nright are incompatable'

  if flag1 + flag2 + flag3 eq 0 then $
     message, 'must input npoints and/or nleft & nright'

  if flag1 eq 1 and flag2 + flag3 eq 0 then begin
     if 2*(long(npoints)/2) ne long(npoints) then $
        message, 'NPOINTS must be even (consider using NLEFT and NRIGHT)'
     nleft = floor(npoints/2)
     nright = floor(npoints/2)
  endif else if flag1 eq 1 and flag2 + flag3 eq 1 then begin
     if flag2 eq 0 then nleft = npoints - nright
     if flag3 eq 0 then nright = npoints - nleft
  endif

  if ~keyword_set(double) and n_elements(double) gt 0 and isa(arr, 'double') then $
     message, 'ARR input is type double but DOUBLE explictly set to 0'

  n_arr = n_elements(arr)
  if ~keyword_set(EXPAND_ARRAY) then begin
     if nleft+nright ge n_arr then message, 'not enough points in array'
     n_arr -= nleft + nright                                ; only used for fcast_order
  endif

  if n_elements(fcast_order) eq 0 then $
     fcast_order = min([n_arr-1, max([10L, long(0.05 * n_arr)])])

  if keyword_set(EXPAND_ARRAY) then begin
     sz = size(arr)
     sz[1] += nleft + nright
     arr_fc = make_array(SIZE=sz, DOUBLE=double)
     arr_fc[nleft] = arr
     if nleft gt 0 then $
        arr_fc[0] = ts_fcast(arr, fcast_order, nleft, /BACKCAST, DOUBLE=double)
     if nright gt 0 then $
        arr_fc[n+nleft] = ts_fcast(arr, fcast_order, nright, DOUBLE=double)
  endif else begin
     if keyword_set(DOUBLE) then arr_fc = double(arr) else arr_fc = arr
     if nleft gt 0 then $
        arr_fc[0] = ts_fcast(arr[nleft:n-nright-1], fcast_order, nleft, /BACKCAST, DOUBLE=double)
     if nright gt 0 then $
        arr_fc[n-nright] = ts_fcast(arr[nleft:n-nright-1], fcast_order, nright, DOUBLE=double)
  endelse

  return, arr_fc
end

function fcast_array, arr, npoints, $
                      DIM=dim, $
                      NLEFT=nleft, NRIGHT=nright, $
                      EXPAND_ARRAY=expand_array, $
                      DOUBLE=double, $
                      FCAST_ORDER=fcast_order

  ;; FCAST_ORDER: number of points to use in ts_fcast (default: 5% of arr
  ;; length but not less than 10, except when there are less than 10 point in arr)

  compile_opt idl2, strictarrsubs

  ndim = n_elements(size(arr,/DIM))

  if n_elements(dim) gt 0 then begin
     if dim lt 1 or dim gt ndim then message, 'DIM is not present in ARR'
  endif

  if ndim eq 1 then begin
     arr_fc = fcast_array_1dim(arr, npoints, $
                               NLEFT=nleft, NRIGHT=nright, $
                               EXPAND_ARRAY=expand_array, $
                               DOUBLE=double, $
                               FCAST_ORDER=fcast_order)
  endif else begin
     if n_elements(dim) eq 0 then message, 'must set DIM if ARR is greater than 1 dimensional'
     arr_fc = rwp_apply_func_to_dim('fcast_array_1dim', arr, dim, npoints, $
                                    NLEFT=nleft, NRIGHT=nright, $
                                    EXPAND_ARRAY=expand_array, $
                                    DOUBLE=double, $
                                    FCAST_ORDER=fcast_order)
  endelse

  return, arr_fc
end
