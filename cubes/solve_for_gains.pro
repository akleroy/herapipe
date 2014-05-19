pro solve_for_gains $
   , list_file $
   , working_dir = working_dir $
   , tag = tag $
   , prev_cube = prev_cube $
   , prev_mask_2d = prev_mask_2d $
   , prev_mask_3d = prev_mask_3d $
   , gain_ascii_file = gain_ascii_file $
   , report = report $
   , show = show $
   , fts = fts

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  @define_hera_pixels.pro

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE PREVIOUS DATA CUBE AND THE MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  cube_file = file_search(working_dir+"cubes/"+prev_cube, count=cube_ct)

  mask_2d_file = file_search(working_dir+"cubes/"+prev_mask_2d, count=mask_2d_ct)

  mask_3d_file = file_search(working_dir+"cubes/"+prev_mask_3d, count=mask_3d_ct)

  if (cube_ct + mask_2d_ct + mask_3d_ct) ne 3 then begin
     message, 'Need a previous mask, map, and cube to solve for the gains.'
  endif else begin
     cube = readfits(cube_file, hdr_cube)
     mask_2d = readfits(mask_2d_file, hdr_mask_2d)
     mask_3d = readfits(mask_3d_file, hdr_mask_3d)
  endelse

; FIX NEEDED TO CLEAN UP THE MASK SO THAT IT DOESN'T INCLUDE NOT-A-NUMBERS
; (THIS SOMETIMES HAPPENS WITH BRIGHT STUFF NEAR THE EDGE OF THE MAP - IT
; SHOULDN'T, BUT IT DOES)

  finite_cube_map = total(finite(cube),3)
  bad_ind = where(mask_2d eq 1 and finite_cube_map eq 0, bad_ct)
  if bad_ct gt 0 then mask_2d[bad_ind] = 0

  bad_ind = where(mask_3d eq 1 and finite(cube) eq 0, bad_ct)
  if bad_ct gt 0 then mask_3d[bad_ind] = 0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE AN OUTPUT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  get_lun, u
  openw, u, gain_ascii_file
  printf,u,'dataset telescop gain 1sig_uncertainty'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND SOLVE FOR GAINS IN EACH
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... ASSEMBLE ALL DATA FOR THIS DAY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
;    READ THE DATA
     indir = working_dir+'spectra/'
     infile = indir+working_name[i]+'_'+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)

;    MEASURE THE SIZE OF THE DATA
     sz = size(data)

;    CHECK THAT THERE IS SOMETHING TO USE
     if total(data.flagged eq 0) eq 0 then continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... LOOP OVER PIXELS FOR THIS DAY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    INITIALIZE A .PS FILE
     if keyword_set(report) then begin
        psfile = working_dir+'reports/gain_'+working_name[i]+'.ps'
        ps, /ps, /def, /color, xsize=7, ysize=10, file=psfile
     endif

     if keyword_set(show) then !p.multi=[0,2,3]

;    INITIALIZE THE GAINS
     gain_val = fltarr(npix)*!values.f_nan
     gain_unc = fltarr(npix)*!values.f_nan

     for j = 0, npix-1 do begin

;       GET DATA FOR THIS PIXEL
        pix_ind = where(data.telescop eq pixel_list[j] and $
                        data.flagged eq 0 and $
                        data.on eq 1 and $
                        data.ref_subtracted eq 1, pix_ct)        
        if pix_ct eq 0 then continue
        
;       GET THE AVERAGE SPECTRUM FROM THESE DATA AND THE REFERENCE CUBE
        get_average_spec $
           , data=data[pix_ind] $
           , mask_2d = mask_2d $
           , hdr_mask_2d = hdr_mask_2d $
           , mask_3d = mask_3d $
           , ref_cube = cube $
           , hdr_cube = hdr_cube $
           , avg_spec = bright_spec $
           , unc_spec = unc_bright_spec $
           , ref_spec = ref_spec $                       
           , ref_unc_spec = unc_ref_spec $
           , sum = sum $
           , unc_sum = unc_sum $
           , ref_sum = ref_sum $
           , ref_unc_sum = unc_ref_sum $
           , show=show
        
;       PROMINENTLY LABEL THE PIXEL
        if keyword_set(show) then begin
           loadct, 0, /silent
           reversect
           legend, /top, /left, box=1, clear=1 $                
                   , textcolor=255 $
                   , lines=[-99,-99] $
                   , ['Data: '+working_name[i],'Pixel: '+pixel_list[j]]
        endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SOLVE FOR THE GAIN IN THIS PIXEL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       IDENTIFY A COMPARISON REGION IN THE SPECTRUM
        ind = where(ref_spec ne 0, ct)

;       CATCH THE CASE OF NO OVERLAP
        if ct eq 0 then begin
           if keyword_set(show) then begin
              plot, [1], [1] $
                    , /nodata, xstyle = 5, ystyle = 5 $
                    , xrange=[0.,1.], yrange=[0.,1]
              xyouts, [0.5], [0.5] $
                      , align=0.5, 'No overlap with bright region.'           
              continue
           endif
        endif
        
;       EXTRACT MEASUREMENTS OVER THE COMPARISON REGION
        x = ref_spec[ind]
        xe = unc_ref_spec[ind]
        y = bright_spec[ind]
        ye = unc_bright_spec[ind]

;       FIND THE GAIN THAT RELATES THE TWO WITH THE LOWEST CHI-SQUARED
;       ... VERY KLUGY, DIRECT SEARCH 0 TO 6
        nrat = 1201.
        ratra = findgen(nrat)/200.
        chisq = fltarr(nrat)
        for k = 0, nrat-1 do $
           chisq[k] = total((y - x*ratra[k])^2/((xe^2+ye^2)))
        chisq_min = min(chisq,minind)
        red_chisq_min = chisq_min / (ct-1.)

        if (minind eq 0) or (minind eq nrat-1) then begin
           message, 'Failed a gain solution. Leaving it unset.', /info
           gain_val[j] = !values.f_nan
           gain_unc[j] = !values.f_nan
        endif else begin           
;          THE GAIN IS THE INVERSE OF THIS NUMBER
           min_rat = ratra[minind]
           gain_val[j] = 1./min_rat
           
;          UNCERTAINTY ESTIMATE BASED ON CHI-SQUARED 
;          (1SIGMA = PERTURBATION REQUIRED TO SHIFT CHI-SQUARED BY 1)
           delta_rat = ratra - min_rat
           delta_chisq = chisq - chisq_min
           low_ind = where(delta_rat lt 0, low_ct)
           if low_ct gt 0 then $
              sig_minus = $
              interpol(delta_rat[low_ind], delta_chisq[low_ind], 1.0) $
           else sig_minus = !values.f_nan
           high_ind = where(delta_rat gt 0, high_ct)
           if high_ct gt 0 then $
              sig_plus = $
              interpol(delta_rat[high_ind], delta_chisq[high_ind], 1.0) $
           else sig_plus = !values.f_nan
           sig = 0.5*(-1.0*sig_minus + sig_plus)
           gain_unc[j] = sig / min_rat  * gain_val[j]
        endelse

;       PLOT
        if keyword_set(show) then begin
           loadct, 0, /silent
           reversect
           nchan = n_elements(data[0].spectrum)
           chan = findgen(nchan)

           circle, /fill
           ploterror, x, y, xe, ye $
                      , xtitle='!6Reference Cube Brightness' $
                      , ytitle='!6Pixel Brightness' $
                      , color=255, errcolor=255, psym=8 $
                      , title=working_name[i]+'  '+pixel_list[j]
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
                   , textcolor=getcolor('black') $
                   , lines=[-99,-99] $
                   , ['!6red. !7v!6!u2!n='+sigfig(red_chisq_min,2) $
                      , 'Gain='+gain_val_str+'!9+!6'+gain_unc_str]
        endif

     endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE A SUMMARY PLOT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     if keyword_set(show) then begin
        !p.multi=[0,1,2]
        loadct,0, /silent
        reversect
        circle
        
        ind = where(finite(gain_val) and finite(gain_unc), ct)
        
        if ct eq 0 then $
           stop

        ploterror, findgen(npix)+1, gain_val, gain_unc, ps=8 $
                   , xtitle='!6Pixel Number (1H01 -> 2H09)' $
                   , ytitle='!6Gain', yrange=[0.,2.0] $
                   , color=255, errcolor=255
        oplot, [-1e6,1e6], [1,1], lines=1, color=255
        
        if keyword_set(report) then begin
           ps, /x
           spawn, 'gv '+psfile+' &'
        endif
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRINT OUTPUT TO THE TEXT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     for j = 0, npix-1 do begin
        if finite(gain_val[j]) eq 0 then continue
        
        printf, u, working_name[i] + ' ' + pixel_list[j] + $
                ' '+string(gain_val[j],format='(F6.3)') + $
                ' ' +string(gain_unc[j],format='(F6.3)')
     endfor
     
  endfor

  close, u

end

