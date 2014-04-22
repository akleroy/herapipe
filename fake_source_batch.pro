pro fake_source_batch $
;  INPUT
   , list_file $                        ; INPUT DATA FILE
   , in_ext = in_ext $                  ; EXTENSION FOR INPUT FILES
   , out_ext = out_ext $
 ;  DATA FLAGGING
   , skip_day = skip_day $              ; DAYS TO OMIT
   , bad_data_file = bad_data_file_in $ ; LIST OF DAY PIXEL SCAN TO SKIP
;  FAKE SOURCE
   , x_ctr = x_ctr $
   , y_ctr = y_ctr $
   , v_ctr = v_ctr $
   , fwhm = fwhm $
   , v_fwhm = v_fwhm $
   , i_peak = i_peak $
   , flux = flux $
   , eff = eff $
;  CONTROL FEEDBACK
   , quiet=quiet $                      ; FLAG: SUPPRESS OUTPUT
   , show=show                          ; FLAG: SHOW OUTPUT ON SCREEN
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(in_ext) eq 0 then in_ext = ''

  if n_elements(skip_day) eq 0 then skip_day = ''  

  if n_elements(out_ext) eq 0 then out_ext = '_fake'

  if out_ext eq '' then begin
     message, 'WARNING! It looks like you are trying to replace the original data.', /info
     message, '*DO* *NOT* *DO* *THIS*.'
     stop
     return
  endif

  if n_elements(bad_data_file_in) eq 0 then bad_data_file_in = ''

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF MAPS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, day, map_name, format='A,X,A'
  day = strcompress(day, /remove_all)
  map_name = strcompress(map_name, /remove_all)
  nmaps = n_elements(map_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN *ALL* DATA (MEMORY INEFFICIENT, BUT SIMPLE AND FAST)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, nmaps-1 do begin

;    CHECK IF WE ARE SKIPPING THIS DAY
     if total(strcompress(skip_day,/remove_all) eq (day[i])[0]) gt 0 then $
        continue

;    READ THE DATA
     infile = '../'+day[i]+'/calib/'+map_name[i]+in_ext+'.fits'
     data = mrdfits(infile,1,hdr,/silent)

;    REMOVE BAD DATA
     remove_bad_data $
        , data $
        , day[i] $
        , bad_data_file = bad_data_file_in $
        , all_bad = all_bad

;    MAKE SURE THERE IS SOMETHING LEFT
     if all_bad then $
        continue
     
;    THE RIGHT ASCENSION AND DECLINATIONS
;    (COMPLICATED BY CONVENTION CHANGE)
     if sxpar(hdr,'CRVAL2') eq 0.0 then begin
        ra = double(data.crval2) + $
                  double(data.cdelt2/cos(!dtor*data.crval3))
        dec = double(data.crval3 + data.cdelt3)
     endif else begin
        ra = double(sxpar(hdr,'CRVAL2')) + $
                  double(data.cdelt2/cos(!dtor*sxpar(hdr,'CRVAL3')))
        dec = double(sxpar(hdr,'CRVAL3') + data.cdelt3)
     endelse

;    MAKE THE VELOCITY AXIS
     crval = sxpar(hdr,'VELO-LSR')
     crpix = sxpar(hdr,'CRPIX1')
     cdelt = sxpar(hdr,'DELTAV')
     
     chan = findgen(sxpar(hdr,'MAXIS1'))
     chan_offset = chan - (crpix-1.0)
     vaxis = chan_offset * cdelt + crval
     if abs(cdelt) gt 100. then $
        vaxis /= 1e3          
  
;    PULL OUT THE SPECTRA
     spec_data = data.spectrum

     add_fake_source $
        , data = spec_data $
        , ra = ra $
        , dec = dec $
        , vaxis = vaxis $
        , source_x = x_ctr $
        , source_y = y_ctr $
        , source_vctr = v_ctr $
        , source_fwhm = fwhm $
        , source_vfwhm = v_fwhm $
        , source_peak = i_peak $
        , source_flux = flux $
        , efficiency = eff

;    REPLACE THE SPECTRA
     data.spectrum = spec_data

;    WRITE THE DATA WE WILL KEEP TO DISK
     outfile = '../'+day[i]+'/calib/'+map_name[i]+out_ext+'.fits'
     if outfile eq infile then begin
        message, 'WARNING! It looks like you are trying to replace the original data.', /info
        message, '*DO* *NOT* *DO* *THIS*.'
        stop
        return
     endif

     spawn, 'rm '+outfile       ; NEEDED TO KEEP MWRFITS FROM APPENDING
     mwrfits, data, outfile, hdr

  endfor


end                             ; of fake_source_batch
