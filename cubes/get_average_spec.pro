pro get_average_spec $
   , data=data $
   , mask_2d = mask_2d $
   , hdr_mask_2d = hdr_mask_2d $
   , mask_3d = mask_3d $
   , ref_cube = ref_cube $
   , hdr_cube = hdr_cube $
   , avg_spec = bright_spec $
   , unc_spec = unc_spec $
   , sum = sum $          
   , unc_sum = unc_sum $
   , ref_spec = ref_spec $    
   , ref_unc_spec = ref_unc_spec $
   , ref_sum = ref_sum $          
   , ref_unc_sum = ref_unc_in_sum $
   , show = show

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT SIZES AND POSITIONS OF THING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; GET THE SIZE OF THE DATA 
  sz = size(data)
  n_chan = n_elements(data[0].spectrum)
  n_spec = sz[1]

; GET THE ASTROMETRY FOR THE MASK AND THE CUBE
  extast, hdr_mask_2d, astrom_mask_2d
  extast, hdr_cube, astrom_cube

; CONVERT THE RA AND DEC OF EACH SPECTURM INTO A PIXEL
  ad2xy, data.ra_deg, data.dec_deg, astrom_mask_2d, x_pix_mask, y_pix_mask
  ad2xy, data.ra_deg, data.dec_deg, astrom_cube, x_pix_cube, y_pix_cube
  
; CHECK IF THESE PIXELS CORRESPOND TO VALID LOCATIONS IN THE MASK
  bright = where(mask_2d[round(x_pix_mask), round(y_pix_mask)] eq 1, bright_ct)
  
; GET A VECTOR OF WEIGHTS AND RMS VALUES
  rms = data.fit_rms
  weight = 1/rms^2

; SHOW THE MASK
  if keyword_set(show) then begin
     loadct, 0, /silent
     reversect
     disp, mask_2d, /sq, xtitle='!6Pixel', ytitle='Pixel' $
           , title='!6Mask of Bright Emission (3 at 3!7r!6)' $
           , color=255
     oplot, [x_pix_mask], [y_pix_mask], ps=3, color=getcolor('blue')
     if bright_ct gt 0 then $
        oplot, [x_pix_mask[bright]], [y_pix_mask[bright]] $
               , ps=1, color=getcolor('red') $
               , symsize = 0.5
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ADD UP THE DATA THAT LIE IN THE BRIGHT REGION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; WORK OUT THE TOTAL WEIGHTS
  if bright_ct gt 0 then $
     sum_of_weight = total(weight[bright],/nan)  

; GET THE AVERAGE SPECTRUM OF THESE DATA
  bright_spec = fltarr(n_chan)
  unc_spec_sq = fltarr(n_chan)

; LOOP OVER ALL SPECTRA THAT LIE IN THE BRIGHT REGION
  for i = 0, bright_ct-1 do begin
     spec = reform(data[bright[i]].spectrum*weight[bright[i]]/sum_of_weight)
     unc = rms[bright[i]]*weight[bright[i]]/sum_of_weight
     x = round(x_pix_cube[bright[i]])
     y = round(y_pix_cube[bright[i]])
     ind = where(mask_3d[x,y,*] eq 1, ct)
     if ct eq 0 then continue
     bright_spec[ind] = bright_spec[ind] + spec[ind]
     unc_spec_sq[ind] = unc_spec_sq[ind] + unc^2     
  endfor
  unc_spec = sqrt(unc_spec_sq)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PERFORM AN IDENTICAL OPERATION ON THE REFERENCE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; GET THE NOISE IN THE REFERENCE CUBE
  ref_rms = mad(ref_cube[where(mask_3d)])

; NOW GO GET THE SAME SPECTRUM FROM THE REFERENCE CUBE
  ref_spec = fltarr(n_chan)
  ref_unc_spec_sq = fltarr(n_chan)

  for i = 0, bright_ct-1 do begin
     spec = (ref_cube[x,y,*])*weight[bright[i]]/sum_of_weight
     unc = ref_rms*weight[bright[i]]/sum_of_weight
     x = round(x_pix_cube[bright[i]])
     y = round(y_pix_cube[bright[i]])
     ind = where(mask_3d[x,y,*] eq 1, ct)
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
