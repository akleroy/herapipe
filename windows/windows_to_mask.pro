pro windows_to_mask $
   , window_root $
   , target_hdr = target_hdr $
   , out_file = out_file $
   , ms_to_kms = ms_to_kms

;
; makes a byte mask of locations within the user-specified fitting
; window(s). 1 = this location appears inside one of the supplied windows, 0 =
; it does not.
;

  if n_elements(out_file) eq 0 then $
     out_file = 'temp_window_mask.fits'

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

  if keyword_set(ms_to_kms) then vaxis /= 1000.

  vres = abs(vaxis[1] - vaxis[0])

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER ALL OF THE WINDOWS SUPPLIED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  n_windows = n_elements(window_root)

  for k = 0, n_windows-1 do begin

;    READ THE WINDOW FILE  
     window_low  = readfits(window_root[k]+'_low.fits', lo_hdr, /silent)
     window_high  = readfits(window_root[k]+'_high.fits', hi_hdr, /silent)

;    ALIGN THE WINDOW MAPS TO THE TARGET HEADER
     hastrom, window_low, lo_hdr, target_hdr $
              , interp=2, cubic=-0.5, missing=!values.f_nan

     hastrom, window_high, hi_hdr, target_hdr $
              , interp=2, cubic=-0.5, missing=!values.f_nan
     
;    APPLY THE WINDOWSGO SPECTRUM-BY-SPECTRUM THROUGH THE CUBE
     for i = 0, nx-1 do begin
        for j = 0, ny-1 do begin

           mask[i,j,*] += $
              (vaxis ge window_low[i,j]) and $
              (vaxis le window_high[i,j])

        endfor
     endfor

  endfor

; COLLAPSE FROM "1 FOR EACH MASK" TO JUST 1S AND 0S
  mask = (mask ge 1)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT TO A FITS FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sxaddpar, target_hdr, 'BUNIT', 'MASK', '1=USE,0=IGNORE'
  sxaddpar, target_hdr, 'HISTORY', 'Mask based on spectral fitting windows.'

  writefits, out_file, mask, target_hdr

end                             ; of window_to_mask
