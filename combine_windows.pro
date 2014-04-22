pro combine_windows $
   , in_root_1 = in_root_1 $
   , in_root_2 = in_root_2 $
   , out_root = out_root $
   , line = line

;+
; NAME:
;
; combine windows
;
; PURPOSE:
;
; Combine a pair of spectral fitting windows into a new set of
; windows. Combination can by either by "or" or "and."
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
; Two sets of four fits files:
;
; in_root_1_noise_low.fits
; in_root_1_line_low.fits
; in_root_1_line_high.fits
; in_root_1_noise_high.fits
; (and the same set for in_root_2)
;
; These give the locations of the upper and lower boundaries of the fitting
; region (the "noise" files) and the upper and lower boundary of the likely
; location of the line (the "line" files).
;
; OPTIONAL INPUTS:
;
; None
;
; KEYWORD PARAMETERS:
;
; USE_OR - either this or USE_AND must be set. If this is set, perform the
;          combination using the "or" operator (i.e., take the narrower
;          window).
;
; USE_AND - either this or USE_OR must be set. If this is set, perform the
;           combination using the "and" operator (i.e., take the wider
;           window)
;
; N.B., in either case, the program cares about the line window and then the
; width of the noise window around the line (rather than operating on the
; noise window itself).
;
; OUTPUTS:
;
; Four fits files: 
;
; out_root_noise_low.fits, 
; out_root_line_low.fits, 
; out_root_line_high.fits,
; out_root_noise_high.fits
;
; (see inputs)
;
; OPTIONAL OUTPUTS:
;
; None
;
; COMMON BLOCKS:
;
; None
;
; SIDE EFFECTS:
;
; None
;
; RESTRICTIONS:
;
; None
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

  on_error, 2

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE TWO INPUT WINDOWS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  in_1_noise_low = readfits(in_root_1+'_noise_low.fits', hdr_1)
  in_1_line_low = readfits(in_root_1+'_line_low.fits', hdr_1)
  in_1_line_high = readfits(in_root_1+'_line_high.fits', hdr_1)
  in_1_noise_high = readfits(in_root_1+'_noise_high.fits', hdr_1)

  in_2_noise_low = readfits(in_root_2+'_noise_low.fits', hdr_2)
  in_2_line_low = readfits(in_root_2+'_line_low.fits', hdr_2)
  in_2_line_high = readfits(in_root_2+'_line_high.fits', hdr_2)
  in_2_noise_high = readfits(in_root_2+'_noise_high.fits', hdr_2)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ALIGN THE SECOND SET OF MAPS TO THE ASTROMETRY OF THE FIRST SET
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  hastrom, in_2_noise_low, hdr_2, in_2_noise_low, hdr_2_temp, hdr_1 $
           , missing=!values.f_nan, ngrid=11, interp=2, cubic=-0.5

  hastrom, in_2_line_low, hdr_2, in_2_line_low, hdr_2_temp, hdr_1 $
           , missing=!values.f_nan, ngrid=11, interp=2, cubic=-0.5

  hastrom, in_2_line_high, hdr_2, in_2_line_high, hdr_2_temp, hdr_1 $
           , missing=!values.f_nan, ngrid=11, interp=2, cubic=-0.5

  hastrom, in_2_noise_high, hdr_2, in_2_noise_high, hdr_2_temp, hdr_1 $
           , missing=!values.f_nan, ngrid=11, interp=2, cubic=-0.5

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIRST DO THE COMBINATION ON THE NOISE WINDOW
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  out_noise_low = in_1_noise_low*0.0

  out_noise_high = in_1_noise_low*0.0

  finite_noise_1 = in_1_noise_low ne in_1_noise_high and $
                   finite(in_1_noise_low) and $
                   finite(in_1_noise_high)

  finite_noise_2 = in_2_noise_low ne in_2_noise_high and $
                   finite(in_2_noise_low) and $
                   finite(in_2_noise_high)

  if keyword_set(use_and) then begin
     use_ind = where(finite_noise_1 eq 1 and finite_noise_2 eq 1, use_ct)
     if use_ct gt 0 then begin
        out_noise_low[use_ind] = (in_1_noise_low > in_2_noise_low)[use_ind]
        out_noise_high[use_ind] = (in_1_noise_high < in_2_noise_high)[use_ind]
     endif
  endif else begin
     only_1 = where(finite_noise_1 eq 1 and finite_noise_2 eq 0, use_1_ct)
     if use_1_ct gt 0 then begin
        out_noise_low[only_1] = in_1_noise_low[only_1]
        out_noise_high[only_1] = in_1_noise_high[only_1]
     endif

     only_2 = where(finite_noise_1 eq 0 and finite_noise_2 eq 1, use_2_ct)
     if use_2_ct gt 0 then begin
        out_noise_low[only_2] = in_2_noise_low[only_2]
        out_noise_high[only_2] = in_2_noise_high[only_2]
     endif

     both = where(finite_noise_1 eq 1 and finite_noise_2 eq 1, use_both_ct)
     if use_both_ct gt 0 then begin
        out_noise_low[both] = (in_1_noise_low < in_2_noise_low)[both]
        out_noise_high[both] = (in_1_noise_high > in_2_noise_high)[both]
     endif
  endelse
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOW DO THE SAME COMBINATION ON THE LINE WINDOW
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  out_line_low = in_1_line_low*0.0

  out_line_high = in_1_line_low*0.0
  
  finite_line_1 = in_1_line_low ne in_1_line_high and $
                  finite(in_1_line_low) and $
                  finite(in_1_line_high)

  finite_line_2 = in_2_line_low ne in_2_line_high and $
                  finite(in_2_line_low) and $
                  finite(in_2_line_high)

  if keyword_set(use_and) then begin
     use_ind = where(finite_line_1 eq 1 and finite_line_2 eq 1, use_ct)
     if use_ct gt 0 then begin
        out_line_low[use_ind] = (in_1_line_low > in_2_line_low)[use_ind]
        out_line_high[use_ind] = (in_1_line_high < in_2_line_high)[use_ind]
     endif
  endif else begin
     only_1 = where(finite_line_1 eq 1 and finite_line_2 eq 0, use_1_ct)
     if use_1_ct gt 0 then begin
        out_line_low[only_1] = in_1_line_low[only_1]
        out_line_high[only_1] = in_1_line_high[only_1]
     endif

     only_2 = where(finite_line_1 eq 0 and finite_line_2 eq 1, use_2_ct)
     if use_2_ct gt 0 then begin
        out_line_low[only_2] = in_2_line_low[only_2]
        out_line_high[only_2] = in_2_line_high[only_2]
     endif

     both = where(finite_line_1 eq 1 and finite_line_2 eq 1, use_both_ct)
     if use_both_ct gt 0 then begin
        out_line_low[both] = (in_1_line_low < in_2_line_low)[both]
        out_line_high[both] = (in_1_line_high > in_2_line_high)[both]
     endif
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE THE OUTPUT WINDOW
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  hdr = hdr_1

  writefits, out_root+'_noise_low.fits', out_noise_low, hdr
  writefits, out_root+'_line_low.fits', out_line_low, hdr
  writefits, out_root+'_line_high.fits', out_line_high, hdr
  writefits, out_root+'_noise_high.fits', out_noise_high, hdr

end                             ; of combine_windows
