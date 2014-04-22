function coadd_data $
   , data $
   , weight = weight $
   , use_median = use_median
  
;+
; NAME:
;
; coadd_data
;
; PURPOSE:
;
; Add together data, optionally applying weights and optionally using the
; median (with or without weights).
;
; CATEGORY:
;
; Data handling tool
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; PARSE THE DATA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  sz = size(data)
  
  if sz[0] ne 2 then begin
     message, 'Require a two-dimensional array for co-adding.'
  endif

; MAKE DEFAULT WEIGHTS AND TAKE CARE OF NOT-A-NUMBERS  
  if n_elements(weights) eq 0 then begin
     weights = finite(data)*1.0
  endif else begin
     weights *= finite(data)*1.0
  endelse

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; IF WE ARE NOT USING THE MEDIAN, DO A WEIGHTED MEAN
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
  
  if keyword_set(use_median) eq 0 then begin     
     return, total(weights * data, 2, /nan) / total(weights, 2)
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; IF WE ARE USING THE MEDIAN, RUN THROUGH THINGS STEP BY STEP
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  output = fltarr(sz[1])*!values.f_nan

  if keyword_set(use_median) then begin
     
     for i = 0, sz[1]-1 do begin
;       GET THE DATA FOR THIS ELEMENT        
        data_vec = reform(data[i,*])
        weight_vec = reform(weight[i,*])
        total_weight = total(weight_vec)
        if total_weight gt 0 then begin
           ind = sort(data_vec)
           data_vec = data_vec[ind]
           norm_weight = weight_vec[ind] / total_weight
           weight_vec[0] = 0.5*norm_weight[0]
           for j = 1, sz[2]-1 do $
              weight_vec[j] = total(norm_weight[0:j-1]) + 0.5*norm_weight[j]
           output[i] = interpol(data_vec,weight_vec,0.5)
        endif
     endfor

     return, output

  endif

end
