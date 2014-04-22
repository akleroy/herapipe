pro window_to_mask $
   , window_root $
   , target_hdr $
   , out_file = out_file $
   , noise = only_noise $
   , line = only_line $
   , border = only_border $
   , n_border = n_border $
   , ms_to_kms = ms_to_kms

;
; makes a byte mask of locations within the user-specified fitting window
;

  if n_elements(out_file) eq 0 then $
     out_file = 'temp_window_mask.fits'

  if n_elements(n_border) eq 0 then $
     n_border = 5

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RESTORE THE WINDOW FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  noise_low  = readfits(window_root+'_noise_low.fits', nlo_hdr, /silent)
  line_low   = readfits(window_root+'_line_low.fits', llo_hdr, /silent)
  line_high  = readfits(window_root+'_line_high.fits', lhi_hdr, /silent)
  noise_high = readfits(window_root+'_noise_high.fits', nhi_hdr, /silent)

  win_sz = size(window_map)

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

  border = n_border*vres

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ALIGN THE WINDOW MAPS TO THE TARGET HEADER
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  hastrom, noise_low, nlo_hdr, target_hdr $
           , interp=2, cubic=-0.5, missing=!values.f_nan

  hastrom, line_low, llo_hdr, target_hdr $
           , interp=2, cubic=-0.5, missing=!values.f_nan

  hastrom, line_high, lhi_hdr, target_hdr $
           , interp=2, cubic=-0.5, missing=!values.f_nan
  
  hastrom, noise_high, nhi_hdr, target_hdr $
           , interp=2, cubic=-0.5, missing=!values.f_nan

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% 
; GO SPECTRUM-BY-SPECTRUM THROUGH THE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, nx-1 do begin
     for j = 0, ny-1 do begin

        if keyword_set(only_border) then begin
           mask[i,j,*] = $
           ((vaxis ge line_low[i,j] - border) and $
            (vaxis lt line_low[i,j] + border)) or $
           ((vaxis ge line_high[i,j] - border) and $
            (vaxis lt line_high[i,j] + border))
        endif else begin
           mask[i,j,*] = $
           (vaxis ge noise_low[i,j]) and $
           (vaxis le noise_high[i,j])
        endelse        

;       CONSIDER ONLY THE FITTING REGIONS, OMITTING THE LINE
        if keyword_set(only_noise) then begin
           if keyword_set(only_border) then begin
              mask[i,j,*] = $
                 ((vaxis ge line_low[i,j] - border) and $
                  (vaxis lt line_low[i,j])) or $
                 ((vaxis ge line_high[i,j]) and $
                  (vaxis lt line_high[i,j] + border))
           endif else begin
              mask[i,j,*] = $
                 ((vaxis ge noise_low[i,j]) and $
                  (vaxis lt line_low[i,j])) or $
                 ((vaxis ge line_high[i,j]) and $
                  (vaxis le noise_high[i,j]))
           endelse
        endif

;       CONSIDER ONLY THE LINE, OMITTING THE FITTING
        if keyword_set(only_line) then begin
           if keyword_set(only_border) then begin
              mask[i,j,*] = $
                 ((vaxis ge line_low[i,j]) and $
                  (vaxis lt line_low[i,j] + border)) or $
                 ((vaxis ge line_high[i,j] - border) and $
                  (vaxis lt line_high[i,j]))
           endif else begin
              mask[i,j,*] = $
                 (vaxis ge line_low[i,j]) and $
                 (vaxis le line_high[i,j])
           endelse
        endif

     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT TO A FITS FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sxaddpar, target_hdr, 'BUNIT', 'MASK', '1=USEABLE,0=BAD'
  sxaddpar, target_hdr, 'HISTORY', 'Mask based on spectral fitting windows.'

  writefits, out_file, mask, target_hdr

end                             ; of make_window_mask
