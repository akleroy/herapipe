pro make_iram_noise_cube $
   , cube_file = cube_file $
   , cube_in = cube $
   , out_file = output_file $
   , mask_file = mask_file $
   , mask_in = mask_in $
   , show = show $
   , rms_cube = noise_cube

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; DATA CUBE
  if n_elements(cube) eq 0 and n_elements(cube_file) gt 0 then $
     cube = readfits(cube_file, hdr)  

  sz = size(cube)

; MASK
; GET THE MASK, 1=FIT, 0=IGNORE
  if n_elements(mask_in) eq n_elements(cube) then begin
;    ... INVERT THE USER SUPPLIED MASK
     mask = (mask_in eq 0) and finite(cube)
  endif else begin
     if n_elements(mask_file) eq 0 then $
        mask = finite(cube) $
     else begin
;    READ IN THE MASK
        mask = readfits(mask_file, mask_hdr) 
        
;    ... "INVERT" THE MASK TO GET WHERE ITS OKAY TO FIT
        mask = (mask eq 0) and finite(cube)
     endelse
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MEASURE SIGMA(X,Y) 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; INITIALIZE OUTPUT
  sigma_xy = fltarr(sz[1],sz[2])*!values.f_nan

; IDENTIFY A REASONABLE OUTLIER CRITERIA
; (HONESTLY 3 SIGMA REJECTION WOULD PROBABLY BE FINE)
  sig_false = bisection(3., 'erf0' $
                        , erftarg = (1d0-(5d-1)/sz[3]), /double)
  
; LOOP OVER PLANES AND FIT
  for i = 0, sz[1]-1 do begin
     for j = 0, sz[2]-1 do begin
        noise_ind = where(mask[i,j,*], noise_ct)

        if noise_ct lt 30 then $
           continue

        data = (cube[i,j,*])[noise_ind]

;       DEFAULT ESTIMATE IS M.A.D. BASED RMS OF THE DATA
        sigma_xy[i,j] = mad(data)
        

        if keyword_set(iterate) then begin
           
;          FIRST ESTIMATE (FROM ONLY THE NEGATIVES)
           neg_ind = where(data lt 0, neg_ct)
           if (neg_ct lt 25) then continue
           sigma = mad([data[neg_ind], -1.*data[neg_ind]])

;          FINAL ESTIMATE (ALL DATA WITH OUTLIER REJECTION)
           use_ind = where(data lt sig_false*sigma)
           sigma_xy[i,j] = mad(data[use_ind])
        endif

     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MEASURE SIGMA(Z)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; INITIALIZE OUTPUT
  sigma_z = fltarr(sz[3])*!values.f_nan

; LOOP OVER PLANES AND FIT
  for i = 0, sz[3]-1 do begin
     noise_ind = where(mask[*,*,i], noise_ct)

     if noise_ct lt 50 then $
        continue

     data = (cube[*,*,i])[noise_ind]

;    DEFAULT ESTIMATE IS M.A.D. BASED RMS OF THE DATA
     sigma_z[i] = mad(data)

     if keyword_set(iterate) then begin
;    IDENTIFY A REASONABLE OUTLIER CRITERIA
;    (HONESTLY 3 SIGMA REJECTION WOULD PROBABLY BE FINE)
        sig_false = bisection(3., 'erf0' $
                              , erftarg = (1d0-(5d-1)/noise_ct), /double)
        
;    FIRST ESTIMATE (FROM ONLY THE NEGATIVES)
        neg_ind = where(data lt 0, neg_ct)
        if (neg_ct lt 25) then continue
        sigma = mad([data[neg_ind], -1.*data[neg_ind]])

;    FINAL ESTIMATE (ALL DATA WITH OUTLIER REJECTION)
        use_ind = where(data lt sig_false*sigma)
        sigma_z[i] = mad(data[use_ind])
     endif

  endfor
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMBINE THEM INTO A NOISE CUBE, INTERPOLATING OVER EMPTIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; SMOOTH THE RMS UNCERTAINTY MAP
  nan_ind = where(finite(sigma_xy) eq 0, nan_ct)
  sigma_xy = median(sigma_xy, 13) 
  if nan_ct gt 0 then sigma_xy[nan_ind] = !values.f_nan

; SMOOTH THE RMS CHANNEL PROFILE  
  nan_ind = where(finite(sigma_z) eq 0, nan_ct)
  sigma_z = median(sigma_z, 13)
  if nan_ct gt 0 then sigma_z[nan_ind] = !values.f_nan

; FIT A POLYNOMIAL TO THE RMS VS. CHANNEL
  chan = findgen(sz[3])
  fit_ind = where(finite(sigma_z))
  coeffs = robust_poly_fit(chan[fit_ind], sigma_z[fit_ind], 2, yfit)

; SHOW OUR RESULTS
  if keyword_set(show) then begin
     
     !p.multi=[0,2,1]

     plot, chan, sigma_z, ps=1 $
           , xtitle='Channel', ytitle='RMS (T!uA!n!d*!n) [K]'
     oplot, chan[fit_ind], yfit, color=getcolor('red'), thick=3
     oplot, median(chan)*[1,1], median(sigma_xy)*[1,1] $
            , psy=7, symsize=5, color=getcolor('blue'), thick=10

     disp, median(sigma_xy,13)

     !p.multi=0
     
  endif

; MAKE A NOISE CUBE: USE THE POLYNOMIAL FIT TO THE NOISE AND SCALE THE
; SIGMA(X,Y) MAP TO HAVE THIS MEDIAN RMS AT EACH CHANNEL

  noise_cube = cube*!values.f_nan
  med_sigma_xy = median(sigma_xy)

  for i = 0, sz[3]-1 do begin
     
     fit_rms = coeffs[0] + (i*1.0)*coeffs[1] + (i*1.0)^2*coeffs[2]

     noise_cube[*,*,i] = fit_rms * sigma_xy / med_sigma_xy

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(output_file) gt 0 then begin

     sxaddpar, hdr, 'HISTORY', 'Now an RMS noise cube.'
     sxaddpar, hdr, 'HISTORY', 'Fit RMS(Z) and RMS(X,Y) to make noise cube.'
     
     writefits, output_file, noise_cube, hdr

  endif

end                             ; of make_noise_cube
