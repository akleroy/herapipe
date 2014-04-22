pro gain_batch, list_file $
                , in_ext = in_ext $
                , skip_day = skip_day $
                , bad_data = bad_data_file_in $                
                , bright_map = bright_map_file $
                , ref_cube_file = ref_cube_file $
                , window_mask_file = window_mask_file $
                , gain_ascii_file = gain_ascii_file $
                , report_dir = report_dir $
                , show = show

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(skip_day) eq 0 then skip_day = ''

  if n_elements(in_ext) eq 0 then in_ext = '_base_simple'

  if n_elements(report_dir) eq 0 then report_dir = '../reports/'

  pixel_list = '30M-'+ [ $
               'W01-1H01','W01-1H02','W01-1H03', $
               'W01-1H04','W01-1H05','W01-1H06', $
               'W01-1H07','W01-1H08','W01-1H09', $
               'W02-2H01','W02-2H02','W02-2H03', $
               'W02-2H04','W02-2H05','W02-2H06', $
               'W02-2H07','W02-2H08','W02-2H09']

  npix = n_elements(pixel_list)
  
  if n_elements(bad_data_file_in) eq 0 then bad_data_file_in = ''

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE REFERENCE CUBE AND THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  mask = readfits(bright_map_file, hdr_mask)

  ref_cube = readfits(ref_cube_file, hdr_cube)

  window_mask = readfits(window_mask_file)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF MAPS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, day, map_name, format='A,X,A'
  day = strcompress(day, /remove_all)
  map_name = strcompress(map_name, /remove_all)
  nmaps = n_elements(map_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A LIST OF UNIQUE DAYS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sorted_days = day[sort(day)]
  day_list = day[uniq(sorted_days)]
  n_days = n_elements(day_list)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE AN OUTPUT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  get_lun, u
  openw, u, gain_ascii_file
  printf,u,'day telescop gain 1sig_uncertainty'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DAYS AND SOLVE FOR GAINS, WRITING THEM TO A TEXT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, n_days-1 do begin

;    CHECK IF WE ARE SKIPPING THIS DAY
     if total(day_list[i] eq skip_day) ne 0 then continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... ASSEMBLE ALL DATA FOR THIS DAY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     files_for_today = where(day_list[i] eq day, file_ct)
     if file_ct eq 0 then continue

     for j = 0, file_ct-1 do begin
        infile = '../'+day[files_for_today[j]]+'/base_sub/'+map_name[files_for_today[j]]+in_ext+'.fits'
        this_data = mrdfits(infile,1,hdr,/silent)
        

;    REMOVE BAD DATA (SPECIFIED BY THE USER)
     remove_bad_data $
        , this_data $
        , day[i] $
        , bad_data_file = bad_data_file_in $
        , all_bad = all_bad

;       MAKE SURE THERE IS SOMETHING LEFT
        if all_bad then $
           continue

;       THE RIGHT ASCENSION AND DECLINATIONS       
        this_ra = sxpar(hdr,'CRVAL2') + $
                  this_data.cdelt2/cos(!dtor*sxpar(hdr,'CRVAL3'))
        this_dec = sxpar(hdr,'CRVAL3') + this_data.cdelt3
        
        this_rms = this_data.fit_rms

        this_weight = 1./this_rms^2

;       STORE THINGS IN BIG FANCY ARRAYS
        if n_elements(telescop) eq 0 then telescop = this_data.telescop else $
           telescop = [telescop, this_data.telescop]
        if n_elements(spectra) eq 0 then spectra = this_data.spectrum else $
           spectra = [[spectra], [this_data.spectrum]]
        if n_elements(ra) eq 0 then ra = this_ra else ra = [ra, this_ra]     
        if n_elements(dec) eq 0 then dec = this_dec else dec = [dec, this_dec]
        if n_elements(rms) eq 0 then rms = this_rms else rms = [rms, this_rms]
        if n_elements(weight) eq 0 then weight = this_weight else $
           weight = [weight, this_weight]
     endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... LOOP OVER PIXELS FOR THIS DAY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    INITIALIZE A .PS FILE
     if keyword_set(show) then begin
        psfile = report_dir+'gain_'+day_list[i]+'.ps'
        ps, /ps, /def, /color, xsize=7, ysize=10, file=psfile
        !p.multi=[0,2,3]
     endif

     gain_val = fltarr(npix)*!values.f_nan
     gain_unc = fltarr(npix)*!values.f_nan

     for j = 0, npix-1 do begin

;       GET DATA FOR THIS PIXEL
        pix_ind = where(telescop eq pixel_list[j], pix_ct)        
        if pix_ct eq 0 then continue

;       PULL A FEW FIELDS OF INTEREST OUT OF THE TABLE
        this_data   = transpose(spectra[*,pix_ind])
        this_ra     = ra[pix_ind]
        this_dec    = dec[pix_ind]
        this_rms    = rms[pix_ind]
        this_weight = weight[pix_ind]
        
;       GET THE AVERAGE SPECTRUM FROM THESE DATA AND THE REFERENCE CUBE
        get_average_spec $
           , data=this_data $
           , ra=this_ra $
           , dec=this_dec $
           , rms = this_rms $
           , weight=this_weight $
           , window_mask = window_mask $
           , mask = mask $
           , hdr_mask = hdr_mask $
           , ref_cube = ref_cube $
           , hdr_cube = hdr_cube $
           , avg_spec = bright_spec $
           , unc_spec = unc_bright_spec $
           , ref_spec = ref_spec $                       
           , ref_unc_spec = unc_ref_spec $
           , sum = sum $
           , unc_in_sum = unc_sum $
           , ref_sum = ref_sum $
           , ref_unc_sum = unc_ref_sum $
           , show=show

;       PROMINENTLY LABEL THE PIXEL
        loadct, 0, /silent
        reversect
        legend, /top, /left, box=1, clear=1 $                
                , textcolor=255 $
                , lines=[-99,-99] $
                , ['Day: '+day_list[i],'Pixel: '+pixel_list[j]]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SOLVE FOR THE GAIN IN THIS PIXEL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       IDENTIFY A COMPARISON REGION IN THE SPECTRUM
        ind = where(ref_spec ne 0, ct)

        if ct eq 0 then begin
           plot, [1], [1], /nodata, xstyle = 5, ystyle = 5, xrange=[0.,1.], yrange=[0.,1]
           xyouts, [0.5], [0.5], align=0.5, 'No overlap with bright region.'           
           continue
        endif

;       EXTRACT MEASUREMENTS OVER THE COMPARISON REGION
        x = ref_spec[ind]
        xe = unc_ref_spec[ind]
        y = bright_spec[ind]
        ye = unc_bright_spec[ind]

;       FIND THE GAIN THAT RELATES THE TWO WITH THE LOWEST CHI-SQUARED
        nrat = 1201.
        ratra = findgen(nrat)/200.
        chisq = fltarr(nrat)
        for k = 0, nrat-1 do $
           chisq[k] = total((y - x*ratra[k])^2/((xe^2+ye^2)))
        chisq_min = min(chisq,minind)        
        red_chisq_min = chisq_min / (ct-1.)

        if (minind eq 0) or (minind eq nrat-1) then begin
           message, 'Failed a gain solution. Setting it to 1.', /info
           gain_val[j] = 1.0
           gain_unc[j] = !values.f_nan
        endif else begin           
;       THE GAIN IS THE INVERSE OF THIS NUMBER
           min_rat = ratra[minind]
           gain_val[j] = 1./min_rat
           
;       ERROR ESTIMATE BASED ON CHI-SQUARED (1SIGMA ERROR = PERTURBATION
;       REQUIRED TO SHIFT CHI-SQUARED BY 1)
           delta_rat = ratra - min_rat
           delta_chisq = chisq - chisq_min
           low_ind = where(delta_rat lt 0, low_ct)
           if low_ct gt 0 then $
              sig_minus = interpol(delta_rat[low_ind], delta_chisq[low_ind], 1.0) $
           else sig_minus = !values.f_nan
           high_ind = where(delta_rat gt 0, high_ct)
           if high_ct gt 0 then $
              sig_plus = interpol(delta_rat[high_ind], delta_chisq[high_ind], 1.0) $
           else sig_plus = !values.f_nan
           sig = 0.5*(-1.0*sig_minus + sig_plus)
           gain_unc[j] = sig / min_rat  * gain_val[j]
        endelse

;       PLOT
        if keyword_set(show) then begin
           loadct, 0, /silent
           reversect
           nchan = n_elements(spectra[*,0])
           chan = findgen(nchan)

           circle, /fill
           ploterror, x, y, xe, ye $
                      , xtitle='!6Reference Cube Brightness' $
                      , ytitle='!6Pixel Brightness' $
                      , color=255, errcolor=255, psym=8 $
                      , title=day_list[i]+'  '+pixel_list[j]
           xfid = findgen(100)-10.
           oplot, xfid, xfid, color=255
           oplot, xfid, xfid/gain_val[j], lines=2, color=255
           oplot, xfid, xfid/(gain_val[j]-gain_unc[j]), lines=1, color=255
           oplot, xfid, xfid/(gain_val[j]+gain_unc[j]), lines=1, color=255

           gain_val_str = finite(gain_val[j]) ? sigfig(gain_val[j],3) : 'NaN'
           gain_unc_str = finite(gain_unc[j]) ? sigfig(gain_unc[j],2) : 'NaN'
           loadct, 0, /silent           
           legend, /top, /left $
                   , box=0, clear=0 $
                   , textcolor=255 $
                   , lines=[-99,-99] $
                   , ['!6red. !7v!6!u2!n='+sigfig(red_chisq_min,2) $
                      , 'Gain='+gain_val_str+'!9+!6'+gain_unc_str]
        endif

     endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RESCALE TO HAVE A MEDIAN GAIN OF 1.0
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;  message, 'Rescaling gains so that median is 1.0', /info
;  med_gain = median(gain_val,/even)
;  gain_val /= med_gain

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A SUMMARY PLOT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  !p.multi=[0,1,2]
  loadct,0, /silent
  reversect
  circle
  
  ind = where(finite(gain_val) and finite(gain_unc), ct)

  ploterror, findgen(npix)+1, gain_val, gain_unc, ps=8 $
             , xtitle='!6Pixel Number (1H01 -> 2H09)' $
             , ytitle='!6Gain', yrange=[0.,2.0] $
             , color=255, errcolor=255
  oplot, [-1e6,1e6], [1,1], lines=1, color=255

  if keyword_set(show) then begin
     ps, /x
     spawn, 'gv '+psfile+' &'
  endif

;  im = tvrd(true=1)
;  write_jpeg, '../reports/'+day_list[i]+'_gains.jpeg', im, true=1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT OUTPUT TO THE TEXT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for j = 0, npix-1 do begin
     if finite(gain_val[j]) eq 0 then continue

     printf, u, day_list[i] + ' ' + pixel_list[j] + $
             ' '+string(gain_val[j],format='(F6.3)') + $
             ' ' +string(gain_unc[j],format='(F6.3)')
  endfor
     
  endfor

  close, u

end

