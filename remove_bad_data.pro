pro remove_bad_data $
   , data $
   , today $
   , bad_data_file = bad_data_file_in $
   , all_bad = all_bad
  

; READ THE FLAGGING FILE
  bad_data_file = file_search(bad_data_file_in, count=ct)
  if ct gt 0 then begin
     readcol, bad_data_file[0], format='A,A,A,F,F' $
              , bad_day, bad_pixel, bad_scan $
              , bad_chan_lo, bad_chan_hi $
              , count = num_bad_line
     if num_bad_line gt 0 then begin
        bad_day = strcompress(bad_day, /remove_all)
        bad_pixel = strcompress(bad_pixel, /remove_all)
        bad_scan = strcompress(bad_scan, /remove_all)
     endif
  endif else begin
     all_bad = 0
     return
  endelse

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
  day_ind = where(today eq bad_day or bad_day eq 'all', day_ct)

; LOOP OVER THOSE FLAGS
  for i = 0L, day_ct-1 do begin

     here = day_ind[i]

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    CASE WHERE ALL PIXELS *AND* ALL SCANS ARE FLAGGED
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

     if (bad_pixel[here] eq 'all') and $
        (bad_scan[here] eq 'all') then begin
;       ... TRIVIAL CASE IF WE ARE DISCARDING
        if discard[here] then begin
           data = -1
           all_bad = 1B
           return
;       ... ELSE FLAG BAD CHANNELS (SET THEM TO NANS)
        endif else begin
           data[*,flag_lo[here]:flag_hi[here]] = !values.f_nan
        endelse
     endif

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    CASE WHERE ALL PIXELS ARE FLAGGED
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

     if bad_pixel[here] eq 'all' then begin

;       ... CASE 1: THROW AWAY THE DATA
        if discard[here] then begin
           keep = where(data.scan ne bad_scan[here], keep_ct)
           if keep_ct eq 0 then begin
              data = -1
              all_bad = 1B
              return
           endif else begin
              data = data[keep]
           endelse
        endif else begin
;       ... CASE 2: FLAG SPECIFIC CHANNELS (SETTING THEM TO NANS)
           flag = where(data.scan eq bad_scan[here], flag_ct)
           if flag_ct gt 0 then begin
              data[flag,flag_lo[here]:flag_hi[here]] = !values.f_nan
           endif
        endelse
     endif

;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
;    CASE WHERE SPECIFIC PIXELS ARE FLAGGED
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

     if bad_pixel[here] ne 'all' then begin
;       ... ALL SCANS ARE FLAGGED
        if bad_scan[here] eq 'all' then begin           
;          ... AND WE ARE DISCARDING DATA
           if discard[here] then begin
              keep = where(data.telescop ne bad_pixel[here], keep_ct)
              if keep_ct eq 0 then begin
                 data = -1
                 all_bad = 1B
                 return
              endif else begin
                 data = data[keep]
              endelse           
           endif else begin
;          ... AND WE ARE FLAGGING DATA
              flag = where(data.telescop eq bad_pixel[here], flag_ct)
              if flag_ct gt 0 then begin
                 data[flag,flag_lo[here]:flag_hi[here]] = !values.f_nan
              endif
           endelse           
        endif else begin
;       ... SPECIFIC SCANS ARE FLAGGED
           if discard[here] then begin
;          ... AND WE ARE DISCARDING DATA
              keep = where((data.telescop ne bad_pixel[here]) or $
                           (data.scan ne long(bad_scan[here])), keep_ct)
              if keep_ct eq 0 then begin
                 data = -1
                 all_bad = 1B
                 return
              endif else begin
                 data = data[keep]
              endelse           
           endif else begin
;          ... AND WE ARE FLAGGING DATA
              flag = where((data.telescop eq bad_pixel[here]) and $
                           (data.scan eq long(bad_scan[here])), flag_ct)
              if flag_ct gt 0 then begin
                 data[flag,flag_lo[here]:flag_hi[here]] = !values.f_nan
              endif
           endelse                      
        endelse
     endif
     
  endfor

  n_stop = n_elements(data)
  n_reject = n_start - n_stop

  all_bad = 0B

end
