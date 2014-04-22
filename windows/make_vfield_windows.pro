pro make_vfield_windows $
;  INPUT
   , vfield_file $
   , gname = gname $
;  OUTPUT FILE STEM
   , out_root = out_root $
;  OPTIONAL PROCESSING OF THE VELOCITY FIELD
   , target_hdr = target_hdr $
   , hel_to_lsr = hel_to_lsr $
   , ms_to_kms = ms_to_kms $
   , smooth = smooth $
;  TUNING PARAMETERS
   , offset = offset $
   , window = window $
   , outer_window = outer_window $
   , transition = transition

;+
; NAME:
;
; make_vfield_windows
;
; PURPOSE:
;
; Create two .fits ("_low" and "_high") that bracket the part of a
; spectrum to be ignored.
;
; CATEGORY:
;
;
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
; Two fits files: 
;
; out_root_low.fits, out_root_high.fits
;
; These give the locations of the upper and lower boundary of the
; part of the spectrum to be avoided when fitting.
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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS, DEFINITIONS, ERROR-CHECKING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; PHYSICAL CONSTANTS BATCH FILE
  @constants.bat

; ROOT OF OUTPUT FILE NAMES
  if n_elements(out_root) eq 0 then $
     out_root = '../windows/'+gname+'_vfield_window'
  
; READ THE VELOCITY FIELD
  vfield = readfits(vfield_file, hdr, /silent)
  blank = sxpar(hdr, 'BLANK')
  blank_ind = where(vfield eq blank, blank_ct)
  if blank_ct gt 0 then vfield[blank_ind] = !values.f_nan
  
; GET INFO FROM MY DATABASE ON THIS GALAXY
  s = things_galaxies(gname)
  
; DEFAULT SOURCE WIDTH = 100 KM/S
  if n_elements(source_window) eq 0 then $
     source_window = 100.       ; KM/S
       
; DEFAULT TO NO OFFSET
  if n_elements(offset) eq 0 then $
     offset = 0.0
 
; FIX NGC 6946 (BUSTED VELOCITY AXIS IN THIS THINGS REDUCTION)
;  if median(vfield) gt 150. and gname eq 'ngc6946' then begin
;     message, 'Correcting NGC6946.', /info
;     if keyword_set(ms_to_kms) then $
;        vfield += -163.593d3 $
;     else $
;        vfield += -163.593
;  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS THE VELOCITY FIELD AS REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; SMOOTH THE VELOCITY FIELD
  if keyword_set(smooth) then $
     for i = 0, 3 do $
        vfield = median(vfield,7)

; CONVERT HELIOCENTRIC TO LSR VELOCITY
  if keyword_set(hel_to_lsr) then begin
  
     if n_elements(gname) eq 0 then begin
        message, 'Need galaxy name to convert helio -> lsr', /info
        return
     endif

;    CONVERT FROM HELIOCENTRIC VELOCITY TO LSR
     helio2lsr, 0., vlsr, ra = s.ra_deg, dec = s.dec_deg
     vfield += vlsr
  endif
  
; CONVERT M/S TO KM/S
  if keyword_set(ms_to_kms) then begin
     vfield /= 1000.
  endif
     
; MOVE TO TARGET ASTROMETRY 
  if n_elements(target_hdr) gt 0 then begin
     hastrom, vfield, hdr, target_hdr $
              , cubic=-0.5, interp=2, missing=!values.f_nan
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT THE SOURCE WINDOW AS A FUNCTION OF LOCATION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; IF THERE IS AN INNER WINDOW AND AN OUTER WINDOW, WORK OUT THE GRID TO
; TRANSITION BETWEEN THEM.
  if (n_elements(outer_source_window) gt 0) and $
     (n_elements(transition) gt 0) then begin
     make_axes, hdr, ri=ri, di=di
     deproject, ri, di, gal=s, rgrid=rgrid
     rgrid *= 3600.
     window = (rgrid le transition)*window + $
              (rgrid gt transition)*outer_window
     window = median(window,7)
  endif

; NOTE HALF THE SOURCE WIDTH (USEFUL)
  half_window = window/2.0
 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE MAPS THAT WE WILL WRITE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  edge_low = vfield + offset - half_window
  edge_high = vfield + offset + half_window
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
  
  sxaddpar, hdr, 'BUNIT', 'KM/S', 'LSR VELOCITY'

  writefits, out_root+'_low.fits', edge_low, hdr
  writefits, out_root+'_high.fits', edge_high, hdr

end                             ; of make_vfield_windows
