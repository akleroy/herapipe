pro spec_batch $
   , list_file $
   , in_ext = in_ext $
   , out_ext = out_ext $
   , rej_ext = rej_ext $
   , skip_day = skip_day $
   , bad_data = bad_data_file_in $
   , window_root = window_root $
   , zap_fourier = zap_fourier $
   , bad_fft_chan = bad_fft_chan $
   , reject_sigma = reject_sigma $
   , degree = degree $
   , no_rejection = no_rejection $
   , decimate = decimate $
   , show = show $
   , report_dir=report_dir $
   , force_v0 = force_v0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(skip_day) eq 0 then skip_day = ''

  if n_elements(in_ext) eq 0 then in_ext = '_refsub'

  if n_elements(out_ext) eq 0 then out_ext = '_base'

  if n_elements(rej_ext) eq 0 then rej_ext = '_rejected'
  
  if n_elements(degree) eq 0 then degree = 1

  if n_elements(bad_data_file_in) eq 0 then bad_data_file_in = ''

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF MAPS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, day, map_name, ref_subtracted, format='A,X,A,A'
  day = strcompress(day, /remove_all)
  map_name = strcompress(map_name, /remove_all)
  nmaps = n_elements(map_name)
  ref_subtracted = strupcase(strcompress(ref_subtracted))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER MAPS AND FIT BASELINES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, nmaps-1 do begin

;    CHECK IF WE ARE SKIPPING THIS DAY
     if total(strcompress(skip_day,/remove_all) eq (day[i])[0]) gt 0 then $
        continue

;    READ THE DATA
     infile = '../'+day[i]+'/calib/'+map_name[i]+in_ext+'.fits'
     data = mrdfits(infile,1,hdr, /silent)

;    CONVERT ANY OLD-FORMAT TELESCOP ENTRIES TO NEW ONES
     update_telescop_field, data

;    KLUGE TO FORCE THE ZERO-POINT VELOCITY (NEEDED FOR OLDER DATA)
     if n_elements(force_v0) ne 0 then $
        sxaddpar, hdr, 'VELO-LSR', force_v0, 'SET BY HAND, BE CAREFUL!'     

;    REMOVE BAD DATA (SPECIFIED BY THE USER)
     remove_bad_data $
        , data $
        , day[i] $
        , bad_data_file = bad_data_file_in $
        , all_bad = all_bad

;    MAKE SURE THERE IS SOMETHING LEFT
     if all_bad then $
        continue
         
;    CONSTRUCT A FOURIER POWER SPECTRUM AND REMOVE BAD FOURIER COMPONENTS
     if keyword_set(zap_fourier) gt 0 then begin
        fourier_analyze $
           , data = data $
           , hdr = hdr $
           , zap = bad_fft_chan $
           , show = show
     endif

;    SUBTRACT A POLYNOMIAL BASELINE FROM EACH SPECTRUM
     hera_base_sub $
        , data = data $
        , hdr = hdr $
        , degree = degree $
        , window_root = window_root $
        , show = show       

     if keyword_set(no_rejection) eq 0 then begin

;    REJECT UNEXPECTEDLY NOISY SPECTRA
        reject_spectra $
           , data = data $
           , hdr = hdr $
           , thresh = reject_sigma $
           , rejected = rejected $
           , decimate = decimate $
           , show = show

;    WRITE OUT THE DATA THAT THE PROGRAM HAS REJECTED ... BUT NOT REFERENCE
;    MEASUREMENTS OR DAY/DATA THAT THE USER HAS ASKED TO OMIT
        outfile = '../'+day[i]+'/base_sub/'+map_name[i]+rej_ext+'.fits'
        spawn, 'rm '+outfile    ; NEEDED TO KEEP MWRFITS FROM APPENDING
        mwrfits, rejected, outfile, hdr
        
     endif    

;    WRITE THE DATA WE WILL KEEP TO DISK
     outfile = '../'+day[i]+'/base_sub/'+map_name[i]+out_ext+'.fits'
     spawn, 'rm '+outfile       ; NEEDED TO KEEP MWRFITS FROM APPENDING
     mwrfits, data, outfile, hdr

  endfor

end                             ; of spec_batch
