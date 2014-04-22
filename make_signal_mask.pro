pro make_signal_mask $
   , in_file $
   , prior = prior_file $
   , sig = sig $
   , chan = chan $
   , min_pix = min_pix $
   , min_area = min_area $
   , grow = grow $
   , out_file = out_file $
   , map_out_file = map_out_file $
   , rms = rms $
   , rms_mask_file = rms_mask_file

;
; makes a byte mask of locations within the user-specified fitting window
;

  if n_elements(chan) eq 0 then $
     chan = 2

  if n_elements(sig) eq 0 then $
     sig = 3.0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  cube = readfits(in_file, hdr)

  if n_elements(rms_mask_file) gt 0 then begin
     test = file_search(rms_mask_file, count=count)
     have_rms_mask = count gt 0
     if have_rms_mask then $
        rms_mask = readfits(rms_mask_file)
  endif else begin
     have_rms_mask = 0
  endelse

  test = file_search(prior_file, count=count)
  have_prior = count gt 0
  if have_prior then begin
     prior = readfits(prior_file, prior_hdr)
     if total((size(prior))[1:3] ne (size(cube))[1:3]) gt 0 then $
        message, 'Prior does not match size of cube. Stopping.'
  endif

  if have_rms_mask eq 0 and have_prior then begin
     if total(prior eq 0B and finite(cube) eq 1B) gt 0 then begin
        message, 'Using complement of prior as RMS mask.', /info
        rms_mask = prior eq 0B and finite(cube) eq 1
        have_rms_mask = 1
     endif
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF THE NOISE IS NOT SUPPLIED, ESTIMATE IT FROM THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(rms) eq 0 then begin
        rms = have_rms_mask ? mad(cube[where(rms_mask)]) : mad(cube)        
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT THE CUBE TO A SIGNIFICANCE MEASUREMENT AT EACH PIXEL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sig_cube = cube / rms

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CHECK AGAINST THE SPECIFIED CONDITION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
; IDENTIFY ALL REGIONS WHERE chan CHANNELS ARE ABOVE sig SIGMA
  conj = sig_cube gt sig
  for i = 1, chan-1 do $
     conj *= shift(sig_cube gt sig,0,0,i)

; SET ALL OF THE PIXELS IN THESE REGIONS TO 1 IN THE MASK
  for i = 1, chan-1 do $
     conj += shift(conj, 0,0,-1*i)
  mask = conj ge 1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PARE AND GROW, IF REQUESTED
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

; GROW THE FINAL MASK (DONE IN 3D NOW, COULD BE REWORKED TO 2D)
  if n_elements(grow) gt 0 then $
     mask = exp_mask(temporary(mask), rad=grow)
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY A PRIOR, IF ONE IS SUPPLIED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if have_prior then begin
     mask *= prior
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT TO A FITS FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sxaddpar, hdr, 'BUNIT', 'MASK', '1=USEABLE,0=BAD'
  sxaddpar, hdr, 'HISTORY', 'Mask based on signal identification.'

  if n_elements(out_file) gt 0 then begin
     writefits, out_file, mask, hdr
  endif

  if n_elements(map_out_file) gt 0 then begin     
     mask_2d = total(mask,3) gt 0
     hdr_2d = twod_head(hdr)
     writefits, map_out_file, mask_2d, hdr_2d
  endif


end                             ; of make_signal_mask
