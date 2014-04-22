pro mask_to_windows $
;  INPUT
   , mask_file $
;  OUTPUT FILE STEM
   , out_root = out_root $
;  OPTIONS RELATING TO ASTROMETRY AND UNITS
   , target_hdr = target_hdr $
   , hel_to_lsr = hel_to_lsr $
   , ms_to_kms = ms_to_kms $
;  OPTIONS TO PROCESS THE MASK
   , pad_vel = pad_vel $
   , exp_rad = exp_rad $
   , min_region_size = min_region_size

;+
; NAME:
;
; mask_to_windows
;
; PURPOSE:
;
; Create two .fits files that indicate part of a spectrum to blank. This
; procedure uses a mask as a template.
;
; OUTPUTS:
;
; Two fits files: 
;
; window_low.fits, window_high.fits
;
; These give the locations of the upper and lower boundaries of the
; region to be ignored while fitting.
;
; MODIFICATION HISTORY:
;
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS, DEFINITIONS, ERROR-CHECKING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; PHYSICAL CONSTANTS BATCH FILE
  @constants.bat
  
; READ THE MASK
  mask = readfits(mask_file, hdr, /silent)

; CONSTRUCT THE VELOCITY AXIS
  make_axes, hdr, vaxis=vaxis, rimg=rimg, dimg=dimg

; CONVERT M/S TO KM/S
  if keyword_set(ms_to_kms) then begin
     vaxis /= 1000.
  endif

; CONVERT HELIOCENTRIC TO LSR VELOCITY
  if keyword_set(hel_to_lsr) then begin
  
;    CONVERT FROM HELIOCENTRIC VELOCITY TO LSR
     helio2lsr, 0., vlsr, ra = mean(rimg), dec = mean(dimg)
     vfield += vlsr
  endif

; MOVE TO TARGET ASTROMETRY 
  if n_elements(target_hdr) gt 0 then begin
     mask = cube_hastrom(cube_in = mask $
                         , hdr_in = hdr $
                         , target_hdr = target_hdr)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS THE MASK AS REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; REJECT REGIONS BELOW A CERTAIN SIZE
  if n_elements(min_region_size) gt 0 then $
     reject_regions, mask, min_region_size

; GROW THE MASK
  if n_elements(exp_rad) gt 0 then $
     mask = exp_mask(mask, rad=exp_rad)
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT THE HIGHEST AND LOWEST VELOCITY AS A FUNCTION OF LOCATION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sz = size(mask)

; INITIALIZE THE WINDOWS
  mean_vel = mean(vaxis)
  window_high = fltarr(sz[1],sz[2])+mean_vel
  window_low = fltarr(sz[1],sz[2])+mean_vel

  for i = 0, sz[1]-1 do begin
     for j = 0, sz[2]-1 do begin

        spec = mask[i,j,*]        
        ind = where(spec eq 1, ct)
        if ct eq 0 then continue
                
        vel_here = vaxis[ind]
        window_high[i,j] = max(vel_here)
        window_low[i,j]  = min(vel_here)        
     endfor
  endfor

  if n_elements(pad_vel) ne 0 then begin
     window_low -= pad_vel
     window_high += pad_vel
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
  
  sxaddpar, hdr, 'BUNIT', 'KM/S', 'LSR VELOCITY'

  writefits, out_root+'_low.fits', window_low, hdr
  writefits, out_root+'_high.fits', window_high, hdr

end                             ; of mask_to_windows
