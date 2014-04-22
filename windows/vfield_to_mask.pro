pro vfield_to_mask $
   , vfield = vfield_in $ 
   , hdr_vfield = vfield_hdr_in $
   , target_hdr = target_hdr $
   , mask = mask $
   , deltav = deltav $
   , ms_to_kms = ms_to_kms $
   , lsr_to_hel = vfield_lsr_to_hel $
   , hel_to_lsr = vfield_hel_to_lsr

;+
; NAME:
;
; vfield_to_mask
;
; PURPOSE:
;
; Convert a velocity field to a byte mask on a target header.
;
; CATEGORY:
;
; Data analysis / science tool.
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
; !!!IMPORTANT!!! Works in units of km/s. Key if you are using the
; LSR/Heliocentric functionality.
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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; COPY TO AVOID EDITING THE INPUT VFIELD
  vfield = vfield_in
  vfield_hdr = vfield_hdr_in

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ALIGN THE VELOCITY FIELD AND CONVERT UNITS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; ... CONVERT FROM M/S TO KM/S
  if keyword_set(ms_to_kms) then $
     vfield  /= 1e3

; ... VFIELD LSR VELOCITY TO HELIO(BARY)CENTRIC
  if keyword_set(vfield_lsr_to_hel) then begin
     make_axes, vfield_hdr, ra=ra, da=da
     mean_ra = mean(ra)
     mean_dec = mean(da)
     helio2lsr, 0., vlsr, ra = mean_ra, dec = mean_dec
     vfield -= vlsr
  endif

; ... VFIELD HELIO(BARY)CENTRIC VELOCITY TO LSR
  if keyword_set(vfield_hel_to_lsr) then begin
     make_axes, vfield_hdr, ra=ra, da=da
     mean_ra = mean(ra)
     mean_dec = mean(da)
     helio2lsr, 0., vlsr, ra = mean_ra, dec = mean_dec
     vfield += vlsr
  endif

; ... ALIGN TO THE TARGET HEADER
  hastrom, vfield, vfield_hdr, twod_head(target_hdr) $
           , interp=2, cubic=-0.5, missing=!values.f_nan

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE (EMPTY) MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nx = sxpar(target_hdr,'NAXIS1')
  ny = sxpar(target_hdr,'NAXIS2')
  nv = sxpar(target_hdr,'NAXIS3')

  mask = bytarr(nx,ny,nv)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE VELOCITY AXIS FOR THE TARGET HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  make_axes, target_hdr, vaxis=vaxis, /vonly

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  vdiff = float(mask)*!values.f_nan

  for i = 0, nv-1 do $
     vdiff[*,*,i] = vaxis[i] - vfield

  mask = abs(vdiff) le deltav

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RETURN / CLEANUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  return

end                             ; of vfield_to_mask
