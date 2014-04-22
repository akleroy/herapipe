function make_signal_mask $
   , cube $
   , rms = rms $
;  CONDITIONS FOR THE MASK
   , prior = prior $
   , sig_thresh = sig_thresh $
   , nchan = nchan $
   , min_pix = min_pix $
   , min_area = min_area $
;  ADDITIONAL PROCESSING
   , grow_xy = grow_xy $
   , grow_z = grow_z

;
; makes a byte mask of locations within the user-specified fitting window
;

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET SOME DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(nchan) eq 0 then $
     nchan = 2

  if n_elements(sig_thresh) eq 0 then $
     sig_thresh = 3.0

  if n_elements(rms) eq 0 then rms = mad(cube)

  sz = size(cube)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT THE DATA TO A SIGNIFICANCE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  rms_sz = size(rms)

; ... WE HAVE ONE NUMBER FOR THE RMS
  if rms_sz[0] eq 0 then $
     sig_cube = cube / rms

; ... WE HAVE AN RMS VECTOR (ASSUME Z)
  if rms_sz[0] eq 1 then begin
     sig_cube = cube
     for i = 0, sz[3]-1 do $
        sig_cube[*,*,i] = cube[*,*,i] / rms[i]
  endif

; ... WE HAVE AN RMS MAP (ASSUME X-Y)
  if rms_sz[0] eq 2 then begin
     sig_cube = cube
     for i = 0, sz[3]-1 do $
        sig_cube[*,*,i] = cube[*,*,i] / rms
  endif

; ... WE HAVE AN RMS CUBE
  if rms_sz[0] eq 3 then begin
     sig_cube = cube / rms
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CHECK AGAINST THE SPECIFIED CONDITION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
; IDENTIFY ALL REGIONS WHERE nchan CHANNELS ARE ABOVE sig SIGMA
  conj = sig_cube gt sig_thresh
  for i = 1, nchan-1 do $
     conj *= shift(sig_cube gt sig_thresh,0,0,i)

; SET ALL OF THE PIXELS IN THESE REGIONS TO 1 IN THE MASK
  for i = 1, nchan-1 do $
     conj += shift(conj, 0,0,-1*i)
  mask = conj ge 1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARE CONTIGUOUS REGIONS BY PIXELS OR AREA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; GET RID OF REGIONS SMALLER THAN A USER-SPECIFIED SIZE
  if n_elements(min_pix) gt 0 then begin
     reg = label_region(mask)
     for i = 1, max(reg) do $
        if (total(reg eq i) lt min_pix) then $
           mask[where(reg eq i)] = 0B     
  endif

  if n_elements(min_area) gt 0 then begin
     mask_2d = total(mask,3) gt 0
     reg = label_region(mask_2d)
     for i = 1, max(reg) do $
        if (total(reg eq i) lt min_area) then $
           mask_2d[where(reg eq i)] = 0B     
     for i = 0, (size(cube))[3]-1 do $
        mask[*,*,i] = mask[*,*,i]*mask_2d
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;  GROW THE FINAL MASK IN THE XY OR Z DIRECTIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(grow_xy) gt 0 then begin
     for i = 0, sz[3]-1 do $
        mask[*,*,i] =  exp_mask(mask[*,*,i], rad=grow_xy, /quiet)
  endif
  
  if n_elements(grow_z) gt 0 then begin
     for i = 0, grow_z do $
        mask = mask or shift(mask,0,0,i) or shift(mask,0,0,-i)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY A PRIOR, IF ONE IS SUPPLIED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(prior) gt 0 then begin
     mask *= prior
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RETURN THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  return, mask

end                             ; of make_signal_mask
