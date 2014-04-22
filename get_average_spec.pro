pro get_average_spec $
   , data=data $                ; THE DATA (# SPECTRA, # CHAN)
   , ra=ra $                    ; R.A. FOR EACH SPECTRUM
   , dec=dec $                  ; DEC FOR EACH SPECTRUM
   , rms=rms $                  ; RMS FOR EACH SPECTRUM
   , weight=weight $            ; WEIGHT FOR EACH SPECTRUM
   , window_mask=window_mask $  ; A 3D MASK SHOWING REGIONS WHERE BASELINES HAVE BEEN FIT
   , mask = mask $              ; A 2D MASK SHOWING THE LOCATION OF BRIGHT EMISSION
   , hdr_mask = hdr_mask $      ; A HEADER CONTAINING THE ASTROMETRY OF THE MASK
   , ref_cube = ref_cube $      ; THE REFERENCE CUBE, CONTAINING "TRUE" VALUES
   , hdr_cube = hdr_cube $      ; A HEADER CONTAINING ASTROMETRY FOR THE REFERENCE CUBE
   , avg_spec = bright_spec $   ; OUTPUT: AVERAGE SPECTRUM
   , unc_spec = unc_spec $      ; OUTPUT: UNCERTAINTY IN THE AVERAGE SPECTRUM
   , sum = sum $                ; OUTPUT: SUM OF DATA OVER COMPARISON REGION
   , unc_in_sum = unc_in_sum $  ; OUTPUT: UNCERTAINTY IN SUM OF DATA
   , ref_spec = ref_spec $      ; OUTPUT: AVERAGE SPECTRUM OF REFERENCE DATA
   , ref_unc_spec = ref_unc_spec $ ; OUTPUT:  UNCERTAINTY IN THE REFERENCE SPECTRUM
   , ref_sum = ref_sum $           ; OUTPUT: SUM OF REFERENCE DATA OVER COMPARISON REGION
   , ref_unc_sum = ref_unc_in_sum $ ; OUTPUT: UNCERTAINTY IN SUM OF REFERENCE DATA
   , show = show

; GET THE SIZE OF THE DATA 
  sz = size(data)
  n_chan = sz[2]
  n_spec = sz[1]

; GET THE ASTROMETRY FOR THE MASK AND THE CUBE
  extast, hdr_mask, astrom_mask
  extast, hdr_cube, astrom_cube

; CONVERT THE RA AND DEC OF EACH SPECTURM INTO A PIXEL
  ad2xy, ra, dec, astrom_mask, x_pix_mask, y_pix_mask
  ad2xy, ra, dec, astrom_cube, x_pix_cube, y_pix_cube
  
; CHECK IF THESE PIXELS CORRESPOND TO VALID LOCATIONS IN THE MASK
  bright = where(mask[round(x_pix_mask), round(y_pix_mask)] eq 1, bright_ct)
  
; SHOW THE MASK
  if keyword_set(show) then begin
     loadct, 0, /silent
     reversect
     disp, mask, /sq, xtitle='!6Pixel', ytitle='Pixel' $
           , title='!6Mask of Bright Emission (3 at 3!7r!6)' $
           , color=255
     oplot, x_pix_mask, y_pix_mask, ps=3, color=getcolor('blue')
     if bright_ct gt 0 then $
        oplot, x_pix_mask[bright], y_pix_mask[bright] $
               , ps=1, color=getcolor('red') $
               , symsize = 0.5
  endif

; WORK OUT THE TOTAL WEIGHTS
  if bright_ct gt 0 then $
     sum_of_weight = total(weight[bright],/nan)  

; GET THE AVERAGE SPECTRUM OF THESE DATA
  bright_spec = fltarr(n_chan)
  unc_spec_sq = fltarr(n_chan)

  for i = 0, bright_ct-1 do begin
     spec = reform(data[bright[i],*]*weight[bright[i]]/sum_of_weight)
     unc = rms[bright[i]]*weight[bright[i]]/sum_of_weight
     x = round(x_pix_cube[bright[i]])
     y = round(y_pix_cube[bright[i]])
     ind = where(window_mask[x,y,*] eq 1, ct)
     if ct eq 0 then continue
     bright_spec[ind] = bright_spec[ind] + spec[ind]
     unc_spec_sq[ind] = unc_spec_sq[ind] + unc^2     
  endfor
  unc_spec = sqrt(unc_spec_sq)

; GET THE NOISE IN THE REFERENCE CUBE
  ref_rms = mad(ref_cube[where(window_mask)])

; NOW GO GET THE SAME SPECTRUM FROM THE REFERENCE CUBE
  ref_spec = fltarr(n_chan)
  ref_unc_spec_sq = fltarr(n_chan)

  for i = 0, bright_ct-1 do begin
     spec = (ref_cube[x,y,*])*weight[bright[i]]/sum_of_weight
     unc = ref_rms*weight[bright[i]]/sum_of_weight
     x = round(x_pix_cube[bright[i]])
     y = round(y_pix_cube[bright[i]])
     ind = where(window_mask[x,y,*] eq 1, ct)
     if ct eq 0 then continue
     ref_spec[ind] = ref_spec[ind] + spec[ind]
     ref_unc_spec_sq[ind] = ref_unc_spec_sq[ind] + unc^2
  endfor
  ref_unc_spec = sqrt(ref_unc_spec_sq)

; WORK OUT THE SUM AND UNCERTAINTY IN THE SUM FOR THE DATA
  sum = total(bright_spec)
  unc_in_sum = sqrt(total(unc_spec^2))
  ref_sum = total(ref_spec)
  ref_unc_in_sum = sqrt(total(ref_unc_spec^2))

end
