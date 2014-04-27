pro flag_ripply_spectra $
   , list_file $
   , tag = tag $
   , blank = blank $
   , smooth = smooth $
   , show = show $
   , report = report $
   , fts = fts

; ABSOLUTE CUTOFF - ABOUT 5% BASED ON A SUBSTANTIAL SWATH OF DATA
  abs_ripple_cut = 20           

; A (USUALLY) LESS STRINGENT STATISTICAL CUT
  sigma_cut = 5. 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
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

;    NUMBER OF ELEMENTS IN ONE SPECTRUM
     n_chan = n_elements(data[0].spectrum)
  
;    AN ARRAY OF CHANNEL NUMBERS
     chan = findgen(n_chan)

;    THE VELOCITY AXIS
     vaxis = chan*data[0].deltav + data[0].v0

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
; FIGURE OUT THE LONGEST SINGLE-SIGN REGION FOR EACH SPECTRUM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     total_bad = 0

;    CALCULATE THE LONGEST RIPPLE IN EACH SPECTRUM
     longest_ripple = intarr(fit_ct)
     for k = 0L, fit_ct-1 do $
        longest_ripple[k] = length_one_sign(im[*,k], /longest)

;    SAVE THIS NUMBER IF WE AREN'T SMOOTHING
     if n_elements(smooth) eq 0 then begin
        data[fit_ind].longest_ripple = longest_ripple
     endif else begin
        if smooth le 1 then begin
           data[fit_ind].longest_ripple = longest_ripple
        endif
     endelse
            
;    FLAG ON ABSOLUTE CUT
     bad = where(longest_ripple gt abs_ripple_cut, bad_ct)
     if bad_ct gt 0 then begin
        data[fit_ind[bad]].flagged = 1
        data[fit_ind[bad]].why_flagged += 'RA '
     endif
     total_bad += bad_ct

;    FLAG ON SIGMA CUT
     sigma_cut_val = sigma_cut*mad(longest_ripple) + median(longest_ripple)
     bad = where(longest_ripple gt sigma_cut_val and $
                 data[fit_ind].flagged ne 1, bad_ct)
     if bad_ct gt 0 then begin
        data[fit_ind[bad]].flagged = 1
        data[fit_ind[bad]].why_flagged += 'RS '
     endif
     total_bad += bad_ct
         
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOTE IN THE HEADER AND PRINT TO SCREEN THE STATS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

    message, 'Flagged '+str(total_bad)+' spectra for ripples', /info
    message, '... this is '+sigfig(total_bad*1./fit_ct*100., 3)+' %', /info

    sxaddpar, hdr, 'HISTORY', 'FLAGGING FROM RIPPLES: ' + $
              sigfig(total_bad*1./fit_ct*100., 3) + ' %'    

    if keyword_set(report) and keyword_set(show) then begin
       
       loadct, 0, /silent
       reversect
       @define_hera_pixels.pro
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

;         AUTOHIST HATES INTS
          loadct, 0, /silent
          reversect
          fasthist, longest_ripple[ind]*1.0, title=pixel_list[j], charsize=1.5

          oplot, abs_ripple_cut*[1.,1.], [-1e6,1e6], color=getcolor('red')

          oplot, sigma_cut_val*[1.,1.], [-1e6,1e6], color=getcolor('magenta')

          im = tvrd(true=1)
          write_jpeg, '../reports/flag_ripple_'+working_name[i]+'.jpeg' $
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

end                             ; of flag_ripply_spectra
