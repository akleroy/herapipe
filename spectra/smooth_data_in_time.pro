pro smooth_data_in_time $
   , data $
   , telescope $
   , timestamp $
   , kernel_width = kernel_width $
   , propogate_blanks = propogate_blanks
  
;+
;
; Reduction subroutine, called by other programs. Takes a data image (spectral
; channel x record), separates it into blocks according to pixel or instrument
; used, sorts these blocks by time, and then smooths each block by a kernel of
; a specified width along the "record" dimension, effectively a smoothing in
; time. Called by other programs (particularly the flagging routines) to pull
; baseline ripples and other pathologies out of otherwise noisy data.
;
; Avoids reordering the data.
;
;-

  if n_elements(kernel_width) eq 0 then $
     kernel_width = 10

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IDENTIFY UNIQUE TELESCOPE FIELDS (MAKES THIS A BIT GENERAL)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  uniq_telescope = telescope[uniq(telescope, sort(telescope))]

  n_telescope = n_elements(uniq_telescope)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER TELESCOPE ENTRIES, EXTRACT DATA FOR EACH AS A BLOCK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, n_telescope-1 do begin

     tel_ind = where(telescope eq uniq_telescope[i], tel_ct)

;    SHOULD NOT HAPPEN!
     if tel_ct eq 0 then continue

     sub_data = data[*,tel_ind]

     sub_time = timestamp[tel_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SORT BY TIME, SMOOTH THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     sort_time = sort(sub_time)

     sub_data = sub_data[*,sort_time]

     if keyword_set(propogate_blanks) then begin
        blank_im = finite(sub_data)*1.0
        blank_im = smooth(blank_im, [1,kernel_width])
     endif

     sub_data = smooth(sub_data, [1, kernel_width], /nan)

     if keyword_set(propogate_blanks) then begin
        blank_ind = where(blank_im lt 1.0, blank_ct)
        if blank_ct gt 0 then $
           sub_data[blank_ind] = !values.f_nan
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; REPLACE IN THE ORIGINAL ARRAY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     data[*,tel_ind[sort_time]] = sub_data

  endfor 

end
