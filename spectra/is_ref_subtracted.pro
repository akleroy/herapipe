function is_ref_subtracted $
   , list_file $
   , tag = tag
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file $
           , working_name, already_ref_sub $
           , format='X,A,A' $
           , comment="#", /silent
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)
  already_ref_sub = strupcase(strcompress(already_ref_sub, /remove_all))

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND INITIALIZE A NEW STRUCTURE FILE FOR EACH
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; INITIALIZE OUR GUESS - THAT ALL DATA ARE ALREADY REF-SUBTRACTED
  all_already_subtracted = 1B

  for i = 0, ndata-1 do begin

;    FLAG INDICATING THAT WE NEED TO REWRITE THESE DATA
     need_to_write = 0

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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIRST HAVE A LOOK AT THE FIELDS, THEY MAY MAKE THIS EASY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    IF THE STRUCTURE SAYS THEY'RE ALL REFERENCE SUBTRACTED, TRUST IT
     if total(data.ref_subtracted eq 0) eq 0 then $
        continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CHECK FOR THE EXISTENCE OF BLANK DATA IN EACH SCAN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    IDENTIFY BLANK SPECTRA
     is_blank = bytarr(sz[1])
     for j = 0L, sz[1]-1 do $
        is_blank[j] = max(data[j].spectrum lt -999.)

;    IDENTIFY INDIVIDUAL SCANS 
     uniq_scans = $
        (data.scan)[uniq(data.scan, sort(data.scan))]
     n_scans = n_elements(uniq_scans)

;    CHECK IF EACH SCAN HAS A BLANK SPECTRA
;    (THIS CAN SIGNIFY A REFERENCE SUBRACTION)
     has_blank = bytarr(n_scans)
     for j = 0L, n_scans-1 do begin
        ind = where((data.scan) eq uniq_scans[j],ct)
        has_blank[j] = total(is_blank[ind]) ne 0
     endfor
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CHECK THE TYPICAL AMPLITUDE OF THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     median_amp = median(abs(data.spectrum))    

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RETURN OUR BEST GUESS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%    

;    IF THE USER SAYS THAT THEY'RE REFERENCE SUBTRACTED, TRUST THEM
     if already_ref_sub[i] eq 'Y' then begin
        data.ref_subtracted = 1
        blank_ind = where(is_blank, blank_ct)
        if blank_ct gt 0 then $
           data[is_blank].on = 0
        need_to_write = 1
     endif else begin

        if (median_amp lt 10.) and total(has_blank eq 0) eq 0 then begin
;       ... HAVE BLANK SPECTRA AND LOW AMPLITUDE = REF SUBTRACTION
           data.ref_subtracted = 1        
           blank_ind = where(is_blank, blank_ct)
           if blank_ct gt 0 then $
              data[is_blank].on = 0
           need_to_write = 1
        endif else begin
;       ... HIGH AMPLITUDE AND NO BLANK SPECTRA = NO REF SUBTRACTION
           if (median_amp gt 10.) and total(has_blank) eq 0 then begin
              data.ref_subtracted = 0
              all_already_subtracted = 0
              need_to_write = 0
           endif

;       ... HIGH AMPLITUDE AND SOME BLANK SPECTRA = NO REF SUBTRACTION + WARNING
           if (median_amp gt 10.) and total(has_blank) ne 0 then begin
              data.ref_subtracted = 0
              all_already_subtracted = 0
              need_to_write = 0
              message $
                 , 'Look at data file '+infile, /info
              message $
                 , 'These data have high amplitude but some blank (reference?) spectra.', /info
              message $
                 , 'I will reference subtract them, but you may want to check this.', /info
           endif

;       ... LOW AMPLITUDE BUT NO BLANKS / AMBIGUOUS (USER NEEDS TO KNOW)
           if (median_amp lt 10.) then begin
              message $
                 , 'Look at data file '+infile, /info
              message $
                 , 'These data have low amplitude and may already be reference-subtracted.', /info
              message $
                 , 'Please check/guess and specify Y or N in orig_data.txt.', /info
           endif
        endelse

     endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF WE HAVE UPDATED OUR STRUCTURE, WE NEED TO REWRITE THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if need_to_write then begin
        outfile = infile
        spawn, 'rm '+outfile
        mwrfits, data, outfile, hdr
     endif

  endfor

  return, all_already_subtracted

end                             ; of is_ref_subtracted

