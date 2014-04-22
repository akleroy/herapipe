pro reject_spectra $
   , data = data $
   , hdr = hdr $
   , thresh_in_sigma = thresh_in_sigma $
   , longest_sign = longest_sign $
   , rejected_spectra = rejected_spectra $
   , decimate = decimate $
   , kern = kern $
   , show = show $
   , quiet = quiet

;+
; NAME:
;
; reject_spectra
;
; PURPOSE:
;
; Reject noisy spectra from an IRAM .30m structure. Removes the offending
; spectra from the main structure (optionally passing them out as their own
; structure).
;
; CATEGORY:
;
; IRAM 30m data reduction script.
;
; CALLING SEQUENCE:
;
; reject_spectra, data, hdr $
;               , thresh_in_sigma = thresh_in_sigma $
;               , rejected_spectra = rejected_spectra $
;               , quiet = quiet
;
; INPUTS:
;
; data 
;
; OPTIONAL INPUTS:
;
; thresh_in_sigma - threshold to reject spectra. Units are sigma (measured
;                   via the median absolute deviation) from the median ratio
;                   of measured/theoretical noise.
;
; decimate        - (warning! you probably don't want to do this) reject the
;                   noisiest 10% of the data, regardless of Tsys.
;
; KEYWORD PARAMETERS:
;
; show - plot a histogram illustrating the rejection
; quiet - don't print any messages
;
; OUTPUTS:
;
; Input DATA is altered by the program.
;
; OPTIONAL OUTPUTS:
;
; Optionally, pass out a block of the rejected spectra as rejected_spectra.
;
; COMMON BLOCKS:
;
; None.
;
; SIDE EFFECTS:
;
; DATA is altered by the program.
;
; RESTRICTIONS:
;
; None.
;
; PROCEDURE:
;
; None.
;
; EXAMPLE:
;
; None.
;
; MODIFICATION HISTORY:
;
; 29mar09 - spun off from main baseline fitting routine. leroy@mpia.de
;
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; SET OF SMOOTHING KERNELS TO USE
  if n_elements(kern) eq 0 then $
     kern = [0, 25]        

; THRESHOLD FOR TSYS-BASED REJECTION IN UNITS OF SIGMA
  if n_elements(thresh_in_sigma) eq 0 then $
     thresh_in_sigma = 2.0
  
; THRESHOLDS FOR "RIPPLE-BASED" REJECTION; THESE NUMBERS ARE THE ACCEPTABLE
; LENGTH OF THE 1ST, 2ND, AND 3RD LONGEST SINGLE-SIGN-CONTIGUOUS REGIONS.

  if n_elements(longest_sign) eq 0 then $
     longest_sign = 20

  if n_elements(longest_sign_2) eq 0 then $
     longest_sign_2 = 15

  if n_elements(longest_sign_3) eq 0 then $
     longest_sign_3 = 13
  
; FREQUENCY STEP USED IN TSYS-CALCULATION
  if n_elements(hdr) eq 0 then begin
     message, 'Warning! Using fiducial frequency step (2MHz).', /info
     dfreq = 2.0d6
  endif else begin
     dfreq = sxpar(hdr,'CDELT1')
  endelse
  
; INIALIZE COUNTER
  rejected_so_far = 0

; BUILD THE VELOCITY AXIS
  nchan = n_elements(data[0].spectrum)
  crval = sxpar(hdr,'VELO-LSR')
  crpix = sxpar(hdr,'CRPIX1')
  cdelt = sxpar(hdr,'DELTAV')
  v = findgen(nchan)
  vdif = v - (crpix-1.0)
  vaxis = vdif * cdelt + crval
  if abs(cdelt) gt 100. then $
     vaxis /= 1e3    

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% 
; IF THE "DECIMATE" FLAG IS THROWN, START BY SORTING THE DATA ACCORDING TO
; THEIR RMS AND THEN REJECTING THE NOISIEST 10%. THIS ISN'T USUALLY ALL THAT
; HELPFUL.
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(decimate) then begin
     ind = sort(data.fit_rms)
     ndata = n_elements(data)
     kept_ind =ind[0:(0.9*ndata)-1]
     rejected_ind = ind[(0.9*ndata):*]
     rejected_spectra = data[rejected_ind]
     rejected_so_far = n_elements(rejected_spectra)
     data = data[kept_ind]     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GET RID OF EMPTY SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  spec_sum = total(data.spectrum,1,/nan)

  kept_ind = where(spec_sum ne 0., kept_ct)
  rejected_ind = where(spec_sum eq 0, rejected_ct)
  
  if rejected_ct gt 0 then begin
     
     if rejected_so_far eq 0 then $
        rejected_spectra = data[rejected_ind] $
     else $
        rejected_spectra = [rejected_spectra, data[rejected_ind]]

     rejected_so_far += rejected_ct
  endif

  if kept_ct eq 0 then begin
     data = -1
     goto, report
  endif

  data = data[kept_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FLAG DATA WITH DRAMATICALLY DIFFERENT UPPER AND LOWER SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

      delta = (data.mean_low - data.mean_high)
      rms_delta = mad(delta)
      
      kept_ind = where(abs(delta) lt 5.*rms_delta, kept_ct)
      rejected_ind = where(abs(delta) ge 5.*rms_delta, rejected_ct)

      print, 'delta ', rejected_ct, kept_ct

      if rejected_ct gt 0 then begin
         
         if rejected_so_far eq 0 then $
            rejected_spectra = data[rejected_ind] $
         else $
            rejected_spectra = [rejected_spectra, data[rejected_ind]]

         rejected_so_far += rejected_ct
      endif

      if kept_ct eq 0 then begin
         data = -1
         goto, report
      endif

      data = data[kept_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% 
; REJECT INDIVIDUAL SPECTRA BASED ON HAVING SYSTEMATIC REGIONS WITH A FIXED
; SIGN (+ OR -). THE IDEA IS TO DETECT RIPPLES ABOUT THE BASELINE FIT.
;
; THEN REPEAT THE EXERCISE BASED ON SMOOTHED VERSIONS OF THE DATA.
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

   @define_hera_pixels.bat     

   for m = 0, n_elements(kern)-1 do begin
      
      flagged_bad = bytarr(n_elements(data))
      
;  LOOP OVER INDIVIDUAL PIXELS
      for j = 0, npix-1 do begin
         ind = where(data.telescop eq pixel_list[j], ct)
         if ct eq 0 then continue

;        ... NO POINT IF THE SMOOTHING KERNEL IS TOO BIG
         if ct lt 3.*kern[m] then continue

;     SMOOTH ADJACENT SPECTRA TOGETHER
         im = data[ind].spectrum

;     BLANK REGIONS OUTSIDE THE FIT
         for k = 0L, ct-1 do begin
            blank_ind = where(vaxis gt data[ind[k]].windows[0] and $
                              vaxis lt data[ind[k]].windows[1], blank_ct)
            if blank_ct gt 0 then im[blank_ind,k] = !values.f_nan
            blank_ind = where(vaxis gt data[ind[k]].windows[2] and $
                              vaxis lt data[ind[k]].windows[3], blank_ct)
            if blank_ct gt 0 then im[blank_ind,k] = !values.f_nan
            
            if n_elements(data[ind[k]].windows) eq 4 then continue
            
            blank_ind = where(vaxis gt data[ind[k]].windows[4] and $
                              vaxis lt data[ind[k]].windows[5], blank_ct)
            if blank_ct gt 0 then im[blank_ind,k] = !values.f_nan        
         endfor
                        
;        NOW SMOOTH THE SPECTRA IN TIME TO BRING OUT LOW-LYING PATHOLOGIES
         im = smooth(im, [1, kern[m]], /nan, /edge_truncate)         

;        INITIALIZE OUTPUT
         longest_1 = lonarr(ct)
         longest_2 = lonarr(ct)
         longest_3 = lonarr(ct)
         
;        FIND THE LONGEST ONE-SIGN REGION
         for k = 0L, ct-1 do begin
            spec = reform(im[*,k])
            length = length_one_sign(spec)
            longest_1[k] = length[0]
            longest_2[k] = length[1]
            longest_3[k] = length[2]             
         endfor

;        REMOVE DATA OVER THE LIMIT
         temp_keep = (longest_1 lt longest_sign and $
                      longest_2 lt longest_sign_2 and $
                      longest_3 lt longest_sign_3)

         keep = temp_keep
         for k = -kern[m]/2, kern[m]/2 do begin
            keep *= shift(temp_keep,k)
         endfor

;        NOTE WHERE DATA HAVE BEEN FLAGGED BAD         
         flag_these_bad = where(keep eq 0, flag_bad_ct)
         if flag_bad_ct gt 0 then $
            flagged_bad[ind[flag_these_bad]] = 1B

         if total(keep) gt 0 then begin
            !p.multi=[0,2,2]
            loadct,5,/silent
            scale = median(abs(im))*2.
            disp, im[*,where(keep)], min=-1*scale, max=scale $
                  , title=pixel_list[j]+' Keeping'       
            if total(keep eq 0) gt 0 then $
            disp, im[*,where(keep eq 0)] $
                  , min=-1*scale, max=scale $
                  , title=pixel_list[j]+' Tossing'         
            if total(keep eq 0) gt 0 then $
            disp, im[*,where(keep eq 0)] gt 0 $
                  , title=pixel_list[j]+' Tossing'         
            if total(keep eq 0) gt 0 then $
            disp, im[*,where(keep eq 0)] lt 0 $
                  , title=pixel_list[j]+' Tossing'         
         endif

      endfor

;  GET RID OF THE DATA REJECTED THIS GO-ROUND
      kept_ind = where(flagged_bad eq 0, kept_ct)
      rejected_ind = where(flagged_bad $
                           , rejected_ct)

      print, 'ripple test '+str(m), rejected_ct, kept_ct

      if rejected_ct gt 0 then begin
         
         if rejected_so_far eq 0 then $
            rejected_spectra = data[rejected_ind] $
         else $
            rejected_spectra = [rejected_spectra, data[rejected_ind]]
         
         rejected_so_far += rejected_ct
      endif
      
      if kept_ct eq 0 then begin
         data = -1
         goto, report
      endif
      
      data = data[kept_ind]

   endfor
   
            
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT THE THEORETICAL NOISE AND COMPARE TO REAL RATIO
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; RADIOMETER FORMULA
      theoretical_noise = data.tsys/sqrt(dfreq*data.obstime)

; THE RATIO OF MEASURED TO THEORETICAL NOISE
      noise_ratio = data.fit_rms/theoretical_noise

; MEDIAN AND ONE SIGMA OF THIS RATIO
      median_ratio = median(noise_ratio)
      sig_ratio = mad(noise_ratio)
      
; THRESHOLD FOR REJECCTION
      thresh = median_ratio + thresh_in_sigma*sig_ratio

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT THE HISTOGRAM/REJECTION THRESHOLD (IF REQUESTED)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

      if keyword_set(show) then begin
         dummy = fitnoise(noise_ratio, nbins=500, /show $
                          , xtitle='Measured/Theoretical Noise')
         oplot, thresh*[1,1], [-1e6,1e6], lines=2, color=getcolor('red')
      endif     
      
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; REJECT DATA WITH NOISE ABOVE THE THRESHOLD
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; THESE DATA ARE TOO NOISY, A LIKELY SIGNAL OF PATHOLOGICAL BASELINES
      rejected_ind = where(noise_ratio ge thresh,rejected_ct)

      if rejected_ct gt 0 then begin

         if rejected_so_far eq 0 then $
            rejected_spectra = data[rejected_ind] $
         else $
            rejected_spectra = [rejected_spectra, data[rejected_ind]]

         rejected_so_far += rejected_ct

      endif else $
         rejected_spectra = -1

; THESE DATA ARE OKAY
      kept_ind = where(noise_ratio lt thresh,kept_ct)
      if kept_ct gt 0 then $
         data = data[kept_ind] $
      else $
         data = -1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ALSO REJECT DATA WITH VERY LOW (LIKELY PATHOLOGICAL) RMS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; LOW RMS
; (RECALCULATE THE RATIO)
      noise_ratio = data.fit_rms/theoretical_noise
      hard_lower_limit = 0.1
      low_ind = where(noise_ratio lt hard_lower_limit, low_ct)
      if low_ct gt 0 then begin

         if rejected_so_far eq 0 then $
            rejected_spectra = data[rejected_ind] $
         else $
            rejected_spectra = [rejected_spectra, data[rejected_ind]]

         rejected_so_far += rejected_ct

      endif

; THESE DATA ARE OKAY
      kept_ind = where(noise_ratio gt hard_lower_limit, kept_ct)
      if kept_ct gt 0 then $
         data = data[kept_ind] $
      else $
         data = -1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF NOT TOLD TO KEEP QUIET, REPORT SOME BASIC RESULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

      report:

      if keyword_set(quiet) eq 0 then begin
         message, 'Report on data rejection:', /info
         message, 'Total spectra : '+str(n_elements(data)+rejected_so_far) $
                  , /info
         message, 'Total rejected: '+str(rejected_so_far), /info
         
         expected = round( $
                    p(thresh_in_sigma)*(rejected_ct+kept_ct) + $
                    0.015*(rejected_ct+kept_ct) $ ; TWIDDLE FOR LONGEST == 20
                         )
         message, 'Expected to reject: '+str(expected), /info
         
         message, 'Median ratio of measured/theoretical noise = ' $
                  + str(median_ratio), /info
         
         message, 'Percentage rejected by pixel:', /info
     @define_hera_pixels.bat
     for i = 0, npix-1 do begin
        kept = total(data.telescop eq pixel_list[i])
        zapped = total(rejected_spectra.telescop eq pixel_list[i])
        message, '... '+pixel_list[i]+' : '+ $
                 string(zapped/(zapped+kept)*100.) + ' %' $
                 , /info
     endfor

  endif

   end
