pro ref_batch $
   , list_file $
   , out_ext = out_ext $
   , skip_day = skip_day $
   , bad_data = bad_data_file_in $
   , window_root = window_root $
   , force_v0 = force_v0 $
   , min_ref = min_ref_time $
   , max_ref = max_ref_time $   
   , ref_from_mask_only = ref_from_mask_only $
   , ref_mask_file=ref_mask_file $
   , show = show

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(skip_day) eq 0 then skip_day = ''

  if n_elements(out_ext) eq 0 then out_ext = '_refsub'

  if n_elements(ref_ext) eq 0 then ref_ext = '_ref'
  
  if n_elements(bad_data_file_in) eq 0 then bad_data_file_in = ''

  if n_elements(min_ref_time) eq 0 then min_ref_time = 5 

  if n_elements(max_ref_time) eq 0 then max_ref_time = 5 

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
     infile = '../'+day[i]+'/calib/'+map_name[i]+'.fits'
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
         
;    CHECK IF THERE HAS BEEN A REFERENCE SUBTRACTION
     if ref_subtracted[i] eq 'Y' then $
        already_ref_subtracted = 1B $
     else if ref_subtracted[i] eq 'N' then $
        already_ref_subtracted = 0B $
     else $
        already_ref_subtracted = is_ref_subtracted(data = data, hdr = hdr)     

;    IF REFERENCES STILL NEED TO BE BUILT AND SUBTRACTED, DO SO
     if already_ref_subtracted eq 0 then begin
        identify_references, data=data $
                             , hdr=hdr $
                             , /rule_of_thumb $
                             , tref=min_ref_size $
                             , mask=ref_mask_file $
                             , /add_field     
        make_and_subtract_reference, data = data $
                                     , hdr = hdr $
                                     , ref_data = ref_data $
                                     , min_ref = min_ref_size $
                                     , max_ref = max_ref_size
     endif          

;    SHOW THE LOCATIONS OF THE ON/OFF MEASUREMENTS
     if keyword_set(show) then begin
        plot, data.cdelt2, data.cdelt3, ps=3, ystyle=16
        if n_elements(ref_data) gt 0 then $
           oplot, ref_data.cdelt2, ref_data.cdelt3, ps=3, color=getcolor('red')
     endif    

;    WRITE THE DATA
     outfile = '../'+day[i]+'/calib/'+map_name[i]+out_ext+'.fits'
     spawn, 'rm '+outfile       ; NEEDED TO KEEP MWRFITS FROM APPENDING
     mwrfits, data, outfile, hdr

     if n_elements(ref_data) gt 0 then begin
;       WRITE THE REFERENCES
        ref_outfile = '../'+day[i]+'/calib/'+map_name[i]+ref_ext+'.fits'
        spawn, 'rm '+ref_outfile ; NEEDED TO KEEP MWRFITS FROM APPENDING
        mwrfits, ref_data, ref_outfile, hdr
     endif
     
  endfor

end                             ; of spec_batch
