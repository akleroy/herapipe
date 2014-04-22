pro flag_noisy_spectra $
   , list_file $
   , tag = tag $
   , blank = blank $
   , smooth = smooth $
   , show = show $
   , report = report $
   , relative_noise_only = relative_noise_only $
   , abs_noise_cute = abs_noise_cut

; AN ABSOLUTE CUTOFF AT ABOUT 3SIGMA FOR A SIGNIFICANT SUBAMPLE
  if n_elements(abs_noise_cut) eq 0 then $
     abs_noise_cut = 1.0

; A LESS STRINGENT STATISTICAL CUTOFF
  sigma_cut = 3. 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND INITIALIZE A NEW STRUCTURE FILE FOR EACH
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

;    FREQUENCY STEP USED IN TSYS-CALCULATION
     if n_elements(hdr) eq 0 then begin
        message, 'Warning! Using fiducial frequency step (2MHz).', /info
        dfreq = 2.0d6
     endif else begin
        dfreq = sxpar(hdr,'CDELT1')
     endelse

;    NUMBER OF ELEMENTS IN ONE SPECTRUM
     n_chan = n_elements(data[0].spectrum)
     
;    AN ARRAY OF CHANNEL NUMBERS
     chan = findgen(n_chan)

;    THE VELOCITY AXIS
     vaxis = chan*data[0].deltav + data[0].v0

;    LOOK AT ON SOURCE, REF-SUBTRACTED DATA
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
              , kernel_width = smooth $
              , /propogate_blanks
        endif
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MEASURE THE NOISE, COMPARE TO THEORY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     total_bad = 0

;    MEASURE THE NOISE IN EACH SPECTRUM
     measured_noise = fltarr(fit_ct)
     for k = 0L, fit_ct-1 do $
        measured_noise[k] = mad(im[*,k])

;    SAVE THIS NUMBER IF WE AREN'T SMOOTHING
     if n_elements(smooth) eq 0 then begin
        data[fit_ind].fit_rms = measured_noise
     endif else begin
        if smooth le 1 then begin
           data[fit_ind].fit_rms = measured_noise
        endif
     endelse

;    WORK OUT THE THEORETICAL NOISE FOR EACH SPECTRUM  
     if n_elements(smooth) eq 0 then begin
        theory_noise = $
           (data[fit_ind].tsys/sqrt(dfreq*data[fit_ind].obstime))
     endif else begin
        theory_noise = $
           (data[fit_ind].tsys/sqrt(dfreq*data[fit_ind].obstime*smooth))
     endelse

;    TAKE THE RATIO OF THE TWO
     rat = measured_noise / theory_noise
     
     if keyword_set(relative_noise_only) eq 0 then begin
;       FLAG ON ABSOLUTE CUT
        bad = where(rat gt abs_noise_cut, bad_ct)
        if bad_ct gt 0 then begin
           data[fit_ind[bad]].flagged = 1
           data[fit_ind[bad]].why_flagged = 'NA '
        endif
        total_bad += bad_ct
     endif

;    FLAG ON SIGMA CUT
     sigma_cut_val = sigma_cut*mad(rat) + median(rat)
     bad = where(rat gt sigma_cut_val and $
                 data[fit_ind].flagged ne 1, bad_ct)
     if bad_ct gt 0 then begin
        data[fit_ind[bad]].flagged = 1
        data[fit_ind[bad]].why_flagged = 'NS '
     endif
     total_bad += bad_ct
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOTE IN THE HEADER AND PRINT TO SCREEN THE STATS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    message, 'Flagged '+str(total_bad)+' spectra for noise', /info
    message, '... this is '+sigfig(total_bad*1./fit_ct*100., 3)+' %', /info

    sxaddpar, hdr, 'HISTORY', 'FLAGGING FROM NOISE: ' + $
              sigfig(total_bad/fit_ct*100., 3) + ' %'    
    
    if keyword_set(report) and keyword_set(show) then begin
       
       loadct, 0, /silent
       reversect
       @define_hera_pixels.bat
       !p.multi=[0,5,4]
       for j = 0, npix-1 do begin

          ind = where(data[fit_ind].telescop eq pixel_list[j], pix_ct)

          if pix_ct eq 0 then begin
             loadct, 0, /silent
             reversect
             plot, findgen(10), title=pixel_list[j]
             oplot, 10.-1.*findgen(10)
             continue
          endif

          loadct, 0, /silent
          reversect
          fasthist, rat[ind] < 1.2, title=pixel_list[j], charsize=1.5 $
                    , /poly

          oplot, abs_noise_cut*[1.,1.], [-1e6,1e6], color=getcolor('red')

          oplot, sigma_cut_val*[1.,1.], [-1e6,1e6], color=getcolor('magenta')

          im = tvrd(true=1)
          write_jpeg, '../reports/flag_noise_'+working_name[i]+'.jpeg' $
                      , im, true=1

       endfor

    endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUT THE DATA TO THE SAME FILE
; (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr

  endfor

end                             ; of flag_noisy_spectra
