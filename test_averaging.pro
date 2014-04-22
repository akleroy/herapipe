pro test_averaging, file $
                    , flagged=flagged


;+
; NAME:
;
; test_averaging
;
; PURPOSE:
;
;
; CATEGORY:
;
; Data reduction tool.
;
; CALLING SEQUENCE:
;
; test_averaging, "some_hera_data.fits"
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

  @define_hera_pixels.bat

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SOME SMART DEFAULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if n_elements(use_windows) eq 0 then use_windows = 1

  if n_elements(flagged) eq 0 then flagged = 0

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; READ THE DATA AND EXTRACT THE PART WE ARE INTERESTED IN
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; ... DATA FILE
  t = mrdfits(file,1,hdr)

; ... VELOCITY AXIS
  nchan = n_elements(t[0].spectrum)

  vaxis = t[0].v0 + findgen(nchan)*t[0].deltav

; ... KEEP ONLY FLAGGED/UNFLAGGED, REFERENCE-SUBTRACTED, ON-SORUCE DATA
  if keyword_set(flagged) eq 0 then begin
     ind = where(t.flagged eq 0 and t.on eq 1 and $
                 t.ref_subtracted eq 1, ct)
  endif else begin
     ind = where(t.flagged eq 1 and t.on eq 1 and $
                 t.ref_subtracted eq 1, ct)
  endelse

; ... CRASH IF EMPTY
  if ct eq 0 then return

  t = t[ind]  

; REDUCE PIXEL LIST TO PIXELS WITH DATA
  keep_pix = bytarr(npix)+1
  for i = 0, npix-1 do begin

     if total(t.telescop eq pixel_list[i]) eq 0 then $
        keep_pix[i] = 0
     
  endfor

  pixel_list = pixel_list[where(keep_pix)]
  npix = n_elements(pixel_list)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; READ THE SPECTRAL FITTING WINDOWS INTO A MASK
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  blank_mask = (t.spectrum eq t.spectrum)*0B

  if keyword_set(use_windows) then begin

;       EXTRACT THE WINDOWS AS AN ARRAY OF LO-HI x NSPEC 
        win_arr = extract_windows(t)
        win_sz = size(win_arr)

;       BLANK THE DATA INSIDE THE WINDOWS
        vaxis_arr = vaxis # (fltarr(ct)+1.)
        for m = 0, win_sz[1]-1, 2 do begin
           win_arr_lo = (fltarr(nchan) + 1.) # win_arr[m,*]
           win_arr_hi = (fltarr(nchan) + 1.) # win_arr[m+1,*]
           blank = where(vaxis_arr ge win_arr_lo and $
                         vaxis_arr le win_arr_hi, blank_ct)
           if blank_ct gt 0 then blank_mask[blank] = 1B
        endfor
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; RUN THE TEST
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; CURRENT STATUS OF BROWSING
  pixel = 0

  n_to_avg = floor(10.^(alog10(n_elements(t)-1)*findgen(101)/100.))

  mad_ra = fltarr(n_elements(n_to_avg))*!values.f_nan
  
  im = t.spectrum
  if blank_ct gt 0 then im[where(blank_mask)] = !values.f_nan

  mad_ra[0] = mad(im[*,0])
  for i = 1L, n_elements(n_to_avg)-1 do $
     mad_ra[i] = $
     mad(total(im[*,0:n_to_avg[i]],2,/nan) $
         /total(finite(im[*,0:n_to_avg[i]])*1.0,2))
  
  plot, n_to_avg, mad_ra, /ylo, ps=1, /xlo, xrange=[1, max(n_to_avg)]
  ind = where(finite(alog10(n_to_avg)) and finite(alog10(mad_ra)) and n_to_avg ge 10)
  sixlin, alog10(n_to_avg[ind]), alog10(mad_ra[ind]), a, sa, b, sb
  oplot, n_to_avg, mad_ra[0] / sqrt(n_to_avg), lines=2, color=getcolor('red')
  oplot, 10^(findgen(100)-50.), 10.^(b[2]*(findgen(100)-50.) + a[2]) $
         , color=getcolor('blue')
  print, a, b

  stop
  
end
