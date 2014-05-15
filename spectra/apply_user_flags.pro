pro apply_user_flags $
   , list_file $
   , tag = tag $
   , bad_data_file = bad_data_file
  
;+
;
; Accepts a text file with bad data identified by scan, channel, day,
; and pixel and applies this flagging to a data set. Does so by
; setting the "FLAGGED" field to "1" and adding "U" to the
; "WHY_FLAGGED" string in the data structure.
;
; Not really intended for close user interaction.
;
;-
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE FLAGGING FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  all_bad = 0

  bad_data_file = file_search(bad_data_file, count=ct)
  if ct eq 0 then begin
     message, "Did not find a data flagging file. Skipping this step.", /info
     return
  endif

  readcol, bad_data_file[0], format='A,A,A,A,F,F' $
           , bad_data, bad_pixel $
           , bad_scan_lo, bad_scan_hi $
           , bad_chan_lo, bad_chan_hi $
           , count = num_bad_line $
           , comment="#"
  if num_bad_line gt 0 then begin
     bad_data = strcompress(bad_data, /remove_all)
     bad_pixel = strcompress(bad_pixel, /remove_all)
     bad_scan_lo = strcompress(bad_scan_lo, /remove_all)
     bad_scan_hi = strcompress(bad_scan_hi, /remove_all)
  endif else begin
     message, "Bad data file appears empty. Skipping this step.", /info
     return
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND FLAG DATA FOR EACH
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin
    
;    READ THE DATA
     indir = '../spectra/'
     infile = indir+working_name[i]+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
      
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)

; SOME QUICK CALCULATIONS
     n_start = n_elements(data)
     n_chan = n_elements(data[0].spectrum)

; A 1 OR 0 INDICATING WHETHER DATA ARE TO BE DISCARDED OR SPECIFIC CHANNELS
; ARE TO BE FLAGGED
     discard = (bad_chan_lo le 0) and (bad_chan_hi ge n_chan)

; GET CHANNELS TO ACTUALLY INDEX WITH
     flag_lo = bad_chan_lo > 0
     flag_hi = bad_chan_hi < (n_chan-1)

; GET A LIST OF FLAGGING FOR THE DAY BEING CONSIDERED
     day_ind = where(working_name[i] eq bad_data or bad_data eq 'all', day_ct)

; LOOP OVER THOSE FLAGS
     for j = 0L, day_ct-1 do begin

        here = day_ind[j]

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    CASE WHERE ALL PIXELS *AND* ALL SCANS ARE FLAGGED
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

        if (bad_pixel[here] eq 'all') and $
           (bad_scan_lo[here] eq 'all') then begin
;       ... TRIVIAL CASE IF WE ARE DISCARDING
           if discard[here] then begin
              data.flagged = data.flagged*0 + 1
              data.why_flagged += 'H '
              all_bad = 1B
              return
           endif else begin
;       ... ELSE FLAG BAD CHANNELS (SET THEM TO NANS)
              data[*,flag_lo[here]:flag_hi[here]] = !values.f_nan
           endelse
        endif

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    CASE WHERE ALL PIXELS ARE FLAGGED
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

        if bad_pixel[here] eq 'all' then begin
           
           flag = where(data.scan ge bad_scan_lo[here] and $
                        data.scan le bad_scan_hi[here], flag_ct)

;       ... CASE 1: THROW AWAY THE DATA
           if flag_ct gt 0 then begin
              if discard[here] then begin        
                 data[flag].flagged = 1B
                 data[flag].why_flagged += 'U'
              endif else begin
;       ... CASE 2: FLAG SPECIFIC CHANNELS (SETTING THEM TO NANS)
                 data[flag,flag_lo[here]:flag_hi[here]] = !values.f_nan
              endelse
           endif

        endif

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    CASE WHERE SPECIFIC PIXELS ARE FLAGGED
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

        if bad_pixel[here] ne 'all' then begin
;       ... ALL SCANS ARE FLAGGED
           if bad_scan_lo[here] eq 'all' then begin           
              
              flag = where(data.telescop eq bad_pixel[here], flag_ct)

;       ... CASE 1: THROW AWAY THE DATA
              if flag_ct gt 0 then begin
                 if discard[here] then begin        
                    data[flag].flagged = 1B
                    data[flag].why_flagged += 'U'
                 endif else begin
;       ... CASE 2: FLAG SPECIFIC CHANNELS (SETTING THEM TO NANS)
                    data[flag,flag_lo[here]:flag_hi[here]] = !values.f_nan
                 endelse
              endif

           endif else begin
;       ... SPECIFIC SCANS ARE FLAGGED
              flag = where((data.telescop eq bad_pixel[here]) and $
                           (data.scan ge long(bad_scan_lo[here])) and $
                           (data.scan le long(bad_scan_hi[here])), flag_ct)

;       ... CASE 1: THROW AWAY THE DATA
              if flag_ct gt 0 then begin
                 if discard[here] then begin        
                    data[flag].flagged = 1B
                    data[flag].why_flagged += 'U'
                 endif else begin
;       ... CASE 2: FLAG SPECIFIC CHANNELS (SETTING THEM TO NANS)
                    data[flag,flag_lo[here]:flag_hi[here]] = !values.f_nan
                 endelse
              endif
           endelse
        endif
        
     endfor

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    WRITE OUT THE DATA TO THE SAME FILE WE READ FROM
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of apply_user_flags
