pro noise_report $
   , list_file $
   , tag = tag $
   , blank = blank

  @define_hera_pixels.bat

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)
 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND EXTRACT RELEVANT DATA FOR EACH FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
  first_plot = 1

  get_lun, lun
  report_file = '../reports/noise_report'+tag+'.txt'
  openw, lun, report_file

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
     
;    LOOK AT UNFLAGGED, REF-SUBTRACTED DATA
     fit_ind = where(data.on eq 1 and $
                     data.ref_subtracted eq 1 and $
                     data.flagged eq 0, fit_ct)

;    TRAP PATHOLOGICAL CASE
     if fit_ct eq 0 then continue

;    DATA TO IMAGE
     im = data[fit_ind].spectrum

;    NUMBER OF ELEMENTS IN ONE SPECTRUM
     n_chan = n_elements(data[0].spectrum)
  
;    AN ARRAY OF CHANNEL NUMBERS
     chan = findgen(n_chan)

;    THE VELOCITY AXIS
     vaxis = chan*data[0].deltav + data[0].v0

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

;    GET THE SIZE OF THE DATA
     sz = size(im)

;    KEEP ONLY AN EVEN POWER OF TWO
     pow2 = long(floor(alog10(sz[2])*1.0/alog10(2.)))
     im = im[*,0:2L^pow2-1L]

     noise_ra = fltarr(pow2+1)*!values.f_nan

;    SUCCESSIVELY REBIN THE DATA, REMEASURING NOISE
     for j = 0, pow2 do begin
        
        if j gt 0 then begin
           sz = size(im)
           im = rebin(im, sz[1], sz[2]/2)
        endif

        sz = size(im)
        
        noise_ra[j] = mad(im)
        
     endfor

;    FIT A POWER LAW TO THE NOISE
     sixlin, alog10(2.^(findgen(pow2+1))), alog10(noise_ra) $
             , fit_a, sig_a, fit_b, sig_b
     
     printf,lun, 'Data set'+working_name[i]
     printf,lun, '... base noise: '+sigfig(1e3*noise_ra[0],3)+' mK (TA*)'
     printf,lun, '... averaging index: '+sigfig(fit_b[2],2)
     printf,lun, '... usable spectra: '+str(fit_ct)
     printf,lun, '... flagged spectra: '+$
            str(total(data.on eq 1 and data.flagged eq 1))
     
;    PLOT
     if first_plot eq 1 then begin
        !p.multi=[0,0,0]
        loadct, 0
        reversect
        circle
        plot, 2^(findgen(pow2+1)), noise_ra, ps=-8, /xlo, /ylo $
              , xtitle='!6Spectra Averaged Together' $
              , ytitle='!6Noise [K, T!dA!n!u*!n]'
        first_plot = 0
     endif else begin
        oplot, 2^(findgen(pow2+1)), noise_ra, ps=-8
     endelse

  endfor

; SAVE THE IMAGE
  im = tvrd(true=1)
  write_jpeg, '../reports/'+'noise_vs_scale'+tag+'.jpeg', im, true=1

; PRINT THE REPORT
  close, lun
  free_lun, lun
  spawn, 'cat '+report_file

end                             ; of noise_report
