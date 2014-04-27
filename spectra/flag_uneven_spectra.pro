pro flag_uneven_spectra $
   , list_file $
   , tag = tag $
   , smooth = smooth $
   , blank = blank $
   , show = show $
   , report = report $
   , narrow = narrow $
   , wide = wide

; ROUGHLY A 5-SIGMA CUT BASED ON STATTING ~20 DATA SETS
  abs_bad_narrow = 1.5       
  narrow_sigma_cut = 3.

; ROUGHLY A 7-SIGMA CUT BASED ON STATTING ~20 DATA SETS
  abs_bad_wide = 0.5       
  wide_sigma_cut = 3.    

  if n_elements(narrow) eq 0 then $
     narrow = 1

  if n_elements(wide) eq 0 then $
     wide = 1

  if narrow eq 0 and wide eq 0 then $
     return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND EXTRACT RELEVANT DATA FOR EACH FILE
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

;    LOOK AT UNFLAGGED, REF-SUBTRACTED DATA
     fit_ind = where(data.on eq 1 and $
                     data.ref_subtracted eq 1 and $
                     data.flagged eq 0, fit_ct)

     if fit_ct eq 0 then continue

     im = data[fit_ind].spectrum

;    BLANK THE FITTING WINDOWS, IF REQUESTED
     if keyword_set(blank) then begin
        win_arr = extract_windows(data[fit_ind])
        win_sz = size(win_arr)
        
        for k = 0L, fit_ct-1 do begin
           for m = 0, win_sz[1]-1, 2 do begin
              blank = where(vaxis ge win_arr[m,k] and $
                            vaxis le win_arr[m+1,k], blank_ct)
              if blank_ct gt 0 then im[blank,k] = !values.f_nan
           endfor
        endfor
     endif

;    SMOOTH THE DATA IN TIME, IF REQUESTED
     if n_elements(smooth) gt 0 then begin
        if smooth gt 1 then begin
           smooth_data_in_time $
              , im $
              , data[fit_ind].telescop $
              , data[fit_ind].ut $
              , kernel_width = smooth
        endif
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% 
;
; CONSTRUCT SEVERAL "DELTA"-PARAMETERS; THESE ARE THE ABSOLUTE VALUE OF THE
; DIFFERENCE BETWEEN DIFFERENT PARTS OF THE SPECTRUM. IT'S A USEFUL WAY TO
; FIND PLATFORMING AND PATHOLOGIES. THERE'S TWO WAYS TO FLAG A BAD DELTA: ON
; STATISTICS (I.E. ALL 3-SIGMA OUTLIERS) AND ON ABSOLUTE VALUE (I.E., A 1
; KELVIN JUMP IS BAD-CAPITAL-B). THIS IS PRETTY REDUNDANT AS THE ABSOLUTE CUT
; IS ALMOST ALWAYS ABOVE THE 3SIGMA CUTOFF. STILL USEFUL TO HAVE IN CASE OF A
; REALLY PATHOLOGICAL DATA SET
;
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     total_bad = 0

     if keyword_set(narrow) then begin
;    ... FIRST OPERATE ON THE MIDDLE OF THE SPECTRUM
;    (THIS GIVES US A CHANCE TO CATCH PLATFORMED DATA)
        
;    CALCULATE DELTA     
        lo = median(im[225:230,*],dim=1)
        hi = median(im[235:240,*],dim=1)
        delta = (hi - lo)

;    SAVE IN STRUCTURE
        data[fit_ind].delta_narrow = delta

;    TOSS OUT THOSE WITH BAD DELTAS ON THE ABSOLUTE CUT
        bad = where(abs(delta) gt abs_bad_narrow, bad_ct)
        if bad_ct gt 0 then begin
           data[fit_ind[bad]].flagged = 1
           data[fit_ind[bad]].why_flagged = 'UAN '
        endif
        total_bad += bad_ct

;    NOW ALSO FLAG 3SIGMA OUTLIER FOR THIS DATA SET     
        mad_delta = mad(delta)
        
        bad = where(abs(delta) gt narrow_sigma_cut*mad(delta), bad_ct)
        if bad_ct gt 0 then begin
           data[fit_ind[bad]].flagged = 1
           data[fit_ind[bad]].why_flagged = 'USN '
        endif
        total_bad += bad_ct
     endif

     if keyword_set(wide) then begin
;    ... NOW GET A GRADIENT ACROSS A LARGE CHUNK OF THE SPECTRUM
;   (THIS FLAGS GENERAL PATHOLOGIES)
        
;   CALCULATE DELTA
        lo = median(im[100:150,*],dim=1)
        hi = median(im[310:360,*],dim=1)
        delta = (hi - lo)

;   SAVE IN STRUCTURE
        data[fit_ind].delta_wide = delta
        
;   TOSS OUT THOSE WITH BAD DELTAS ON THE ABSOLUTE CUT
        bad = where(abs(delta) gt abs_bad_wide and $
                    data[fit_ind].flagged eq 0, bad_ct)
        if bad_ct gt 0 then begin
           data[fit_ind[bad]].flagged = 1
           data[fit_ind[bad]].why_flagged = 'UAW '
        endif
        total_bad += bad_ct

;   NOW ALSO FLAG 3SIGMA OUTLIER FOR THIS DATA SET     
        mad_delta = mad(delta)
        
        bad = where(abs(delta) gt wide_sigma_cut*mad(delta) and $
                    data[fit_ind].flagged eq 0, bad_ct)
        if bad_ct gt 0 then begin
           data[fit_ind[bad]].flagged = 1
           data[fit_ind[bad]].why_flagged = 'USW '
        endif
        total_bad += bad_ct
     endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOTE IN THE HEADER AND PRINT TO SCREEN THE STATS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     message, 'Flagged '+str(total_bad)+' spectra for uneven-ness', /info
     message, '... this is '+sigfig(total_bad*1./fit_ct*100., 3)+' %', /info

     sxaddpar, hdr, 'HISTORY', 'FLAGGING FROM DELTA: ' + $
               sigfig(total_bad/fit_ct*100., 3) + ' %'    

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUT THE DATA TO THE SAME FILE
; (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of flag_pathological_spectra
