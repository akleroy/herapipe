pro make_flat_windows $
   , gname = gname $
   , vcenter = vcenter $
   , window = window $   
   , offset = offset $
   , out_root = out_root $
   , target_hdr = target_hdr $
   , hel_to_lsr = hel_to_lsr $
   , ms_to_kms = ms_to_kms

;+
; NAME:
;
; make_flat_windows
;
; PURPOSE:
;
; Create two .fits files that indicate part of a spectrum to
; blank. This version makes a simple set of "flat" windows, centered
; at a particular velocity with a single width.
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

; ROOT OF OUTPUT FILE NAMES
  if n_elements(out_root) eq 0 then $
     out_root = '../windows/'+gname+'_flat'

; DEFAULT TO NO OFFSET
  if n_elements(offset) eq 0 then $
     offset = 0.0
    
; GET INFO FROM MY DATABASE ON THIS GALAXY
; s = things_galaxies(gname)
  s = gal_data(gname, found=found) ; data_dir='local_root/galbase/gal_data/'
  if found eq 0 then STOP
  
  if n_elements(vcenter) eq 0 then begin
     vcenter = s.vhel_kms
     hel_to_lsr = 1
     ms_to_kms = 0
  endif
  
  if n_elements(window) eq 0 then begin
     window = (2.*s.vrot_kms) > 50.
     if finite(window) eq 0 then $
        window = 50.
  endif

  half_window = window / 2.0
       
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A DUMMY HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; IDEA IS FOR IT TO COVER THE ENTIRE PLAUSIBLE AREA OF THE GALAXY WITHOUT
; USING MANY PIXELS (I.E., DISK SPACE). SO MAKE THE CDELTS HUGE ...
  
  naxis_vec = [11, 11]
  mkhdr, hdr, 2, naxis_vec

  ; MAKE THE POSITION AXES
  sxaddpar, hdr, 'CTYPE1', 'RA---TAN'
  sxaddpar, hdr, 'CRVAL1', s.ra_deg
  sxaddpar, hdr, 'CRPIX1', 6.
  sxaddpar, hdr, 'CDELT1', -1.0*20./60.

  sxaddpar, hdr, 'CTYPE2', 'DEC--TAN'
  sxaddpar, hdr, 'CRVAL2', s.dec_deg
  sxaddpar, hdr, 'CRPIX2', 6.
  sxaddpar, hdr, 'CDELT2', 20./60.

; ADD THE J2000 EQUINOX
  sxaddpar, hdr, 'EQUINOX', 2000

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS THE VELOCITY FIELD AS REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; CONVERT HELIOCENTRIC TO LSR VELOCITY
  if keyword_set(hel_to_lsr) then begin
  
     if n_elements(gname) eq 0 then begin
        message, 'Need galaxy name to convert helio -> lsr', /info
        return
     endif

;    CONVERT FROM HELIOCENTRIC VELOCITY TO LSR
     helio2lsr, 0., vlsr, ra = s.ra_deg, dec = s.dec_deg
     vcenter += vlsr
  endif
  
; CONVERT M/S TO KM/S
  if keyword_set(ms_to_kms) then begin
     vcenter /= 1000.
  endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE SOME MAPS WITH THE RIGHT SIZE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  nx = sxpar(hdr,'NAXIS1')
  ny = sxpar(hdr,'NAXIS2')

  edge_low = fltarr(nx,ny) + vcenter + offset - half_window
  edge_high = fltarr(nx,ny) + vcenter + offset + half_window

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
  
  sxaddpar, hdr, 'BUNIT', 'KM/S', 'LSR VELOCITY'

  writefits, out_root+'_low.fits', edge_low, hdr
  writefits, out_root+'_high.fits', edge_high, hdr

end                             ; of make_flat_windows
