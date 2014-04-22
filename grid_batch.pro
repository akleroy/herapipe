pro grid_batch $
;  INPUT
   , list_file $                        ; INPUT DATA FILE
   , in_ext = in_ext $                  ; EXTENSION FOR INPUT FILES
;  DATA FLAGGING
   , skip_day = skip_day $              ; DAYS TO OMIT
   , bad_data_file = bad_data_file_in $ ; LIST OF DAY PIXEL SCAN TO SKIP
;  GRID DEFINITION
   , galaxy = gname $                   ; USED TO LOOK UP GALAXY CENTER
   , xctr=xctr $                        ; HDR: FORCE/RETURN RA CENTER
   , yctr=yctr $                        ; HDR: FORCE/RETURN DEC CENTER
   , xsize=xsize $                      ; HDR: FORCE/RETURN CUBE X-SIZE
   , ysize=ysize $                      ; HDR: FORCE/RETURN CUBE Y-SIZE
   , pix_scale=pix_scale $              ; HDR: FORCE/RETURN PIXEL SCALE
;  DATA GROUPING
   , grid_together = together $         ; FLAG: GRID ALL THESE DATA INTO ONE CUBE?
   , grid_by_coverage = by_coverage $   ; FLAG: GRID BY INDIVIDUAL GALAXY COVERAGE?
   , coverage_list = coverage_list $    ; FILENAME: num day orig_file scan angle
   , split_by_pol = split_by_pol $      ; FLAG: SPLIT DATA BY POLARIZATION
   , split_by_angle = split_by_angle $  ; FLAG: SPLIT DATA BY POLARIZATION
;  OUTPUT DEFINITION
   , cube_list_out = cube_list_out $    ; STRING: A FILE TO WHICH WE WILL WRITE THE NAMES OF THE CUBES WE MAKE
   , together_out = together_out_root $ ; STRING: NAME OF A BIG, COMBINED CUBE
   , tag = tag $                        ; STRING APPENDED TO FILE NAMES (e.g., "_iter1")
;  METHODOLOGY
   , grid_by_median=use_median $        ; FLAG: USE MEDIAN-BASED TECHNIQUE TO GRID?
   , apply_gain = apply_gain $          ; FLAG: APPLY THE PIXEL IMAGE GAINS?
   , gain_file = gain_file $            ; FILENAME: day telescop gain uncertainty
   , gauss_kern = gauss_kern $          ; FLAG: USE THE GAUSSIAN KERNEL (IF NOT USE MANGUM KERNEL)
   , gauss_fwhm = gauss_fwhm $          ; FWHM OF GAUSSIAN KERNEL TO USE
   , use_mean_ctr = use_mean_ctr $      ; FLAG: USE THE MEAN CENTER OF THE MAP (NOT THE GALAXY CENTER)
;  CONTROL FEEDBACK
   , quiet=quiet $                      ; FLAG: SUPPRESS OUTPUT
   , show=show                          ; FLAG: SHOW OUTPUT ON SCREEN
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(in_ext) eq 0 then in_ext = ''

  if n_elements(skip_day) eq 0 then skip_day = ''  

  if n_elements(pix_scale) eq 0 then pix_scale = 2.0/3600.

  if n_elements(cube_list_out) eq 0 then cube_list_out = 'grid_batch_record.txt'

  if n_elements(bad_data_file_in) eq 0 then bad_data_file_in = ''

  pad_pix = 10                  ; pixel padding when making the cube

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF MAPS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, day, map_name, format='A,X,A'
  day = strcompress(day, /remove_all)
  map_name = strcompress(map_name, /remove_all)
  nmaps = n_elements(map_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN *ALL* DATA (MEMORY INEFFICIENT, BUT SIMPLE AND FAST)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, nmaps-1 do begin

;    CHECK IF WE ARE SKIPPING THIS DAY
     if total(strcompress(skip_day,/remove_all) eq (day[i])[0]) gt 0 then $
        continue

;    READ THE DATA
     infile = '../'+day[i]+'/base_sub/'+map_name[i]+in_ext+'.fits'
     this_data = mrdfits(infile,1,hdr,/silent)

;    REMOVE BAD DATA
     remove_bad_data $
        , this_data $
        , day[i] $
        , bad_data_file = bad_data_file_in $
        , all_bad = all_bad

;    MAKE SURE THERE IS SOMETHING LEFT
     if all_bad then $
        continue

;    SAVE THE FIRST HEADER
     if n_elements(first_hdr) eq 0 then $
        first_hdr = hdr

;    THE SYSTEM TEMPERATURE AND THE RMS ABOUT THE FIT (TO USE AS WEIGHTS)
     this_tsys = this_data.tsys
     if total(tag_names(this_data) eq 'FIT_RMS') ne 0 then $
        this_fit_rms = this_data.fit_rms $
     else $
        this_fit_rms = this_data.tsys*!values.f_nan

;    A VECTOR NOTING THE FILE THAT THE DATA CAME FROM
     this_orig_file = replicate(map_name[i],n_elements(this_data))

;    A VECTOR NOTING THE DAY (POSSIBLY REDUNDANT BUT CONVENIENT)
     this_orig_day = replicate(day[i], n_elements(this_data))

;    A VECTOR NOTING THE POLARIZATION
     this_pol = (strpos(this_data.telescop,'W01') ne -1)*1 + $
                (strpos(this_data.telescop,'W02') ne -1)*2
                
;    A VECTOR OF SCAN NUMBERS
     this_scan = this_data.scan

;    A VECTOR OF PIXEL IDENTIFIERS
     this_telescop = this_data.telescop

;    THE SPECTRA
     this_spec = transpose(this_data.spectrum)
     
;    THE RIGHT ASCENSION AND DECLINATIONS
;    (COMPLICATED BY CONVENTION CHANGE)
     if sxpar(hdr,'CRVAL2') eq 0.0 then begin
        this_ra = double(this_data.crval2) + $
                  double(this_data.cdelt2/cos(!dtor*this_data.crval3))
        this_dec = double(this_data.crval3 + this_data.cdelt3)
     endif else begin
        this_ra = double(sxpar(hdr,'CRVAL2')) + $
                  double(this_data.cdelt2/cos(!dtor*sxpar(hdr,'CRVAL3')))
        this_dec = double(sxpar(hdr,'CRVAL3') + this_data.cdelt3)
     endelse
     
;    CHECK THAT THE VELOCITY AXIS IS THE SAME FOR EACH DATA SET
     crval = sxpar(hdr,'VELO-LSR')
     crpix = sxpar(hdr,'CRPIX1')
     cdelt = sxpar(hdr,'DELTAV')
     
     chan = findgen(sxpar(hdr,'MAXIS1'))
     chan_offset = chan - (crpix-1.0)
     this_vaxis = chan_offset * cdelt + crval
     if abs(cdelt) gt 100. then $
        this_vaxis /= 1e3               
     
     if (i eq 0) then $
        ref_vaxis = this_vaxis $
     else if total(ref_vaxis ne this_vaxis) gt 0 then $
        message, 'Velocity axis mismatch. Stopping for inspection.'    

;    CONCATENATE DATA INTO GIANT ARRAYS    
     if n_elements(data) eq 0 then begin
        data = this_spec
        ra = this_ra
        dec = this_dec
        tsys = this_tsys
        fit_rms = this_fit_rms
        scan = this_scan
        telescop = this_telescop
        pol = this_pol
        orig_file = this_orig_file
        orig_day = this_orig_day
     endif else begin
        data = [data, this_spec]
        ra = [ra, this_ra]
        dec = [dec, this_dec]
        tsys = [tsys, this_tsys]
        fit_rms = [fit_rms, this_fit_rms]
        scan = [scan, this_scan]
        this_telescop = [telescop, this_telescop]
        pol = [pol, this_pol]
        orig_file = [orig_file, this_orig_file]
        orig_day = [orig_day, this_orig_day]
     endelse

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CHOOSE THE RMS ABOUT THE FIT OR TSYS AS WEIGHTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if total(finite(fit_rms) eq 0) gt 0 then begin
     message, 'RMS about the fit not available for all spectra. Using Tsys-based weighting.', /info
     weight = 1./tsys^2 
  endif else begin 
     weight = 1./fit_rms^2
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT THE REQUIRED SIZE AND EXTENT FOR THE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; WORK OUT THE EXTENT AND AVERAGE RA AND DEC
  min_ra = min(ra)
  max_ra = max(ra)

  min_dec = min(dec)
  max_dec = max(dec)

  delta_ra = (max_ra - min_ra)*cos(!dtor*(max_dec+min_dec)*0.5)*3600.
  delta_dec = (max_dec - min_dec)*3600.

  mean_ra = mean(ra)
  mean_dec = mean(dec)

  if keyword_set(quiet) eq 0 then begin
     print, 'Total extent in RA (arcseconds):', delta_ra
     print, 'Total extent in DEC (arcseconds):',  delta_dec
     print, 'Mean RA (degrees):', mean_ra
     print, 'Mean DEC (degrees):',  mean_dec
  endif

; WORK OUT THE CENTER OF THE MAP
  if (n_elements(xctr) eq 0) or (n_elements(yctr) eq 0) then begin
     s = things_galaxies(gname)
     if keyword_set(use_mean_ctr) or (size(s))[0] eq 0 then begin
        if keyword_set(quiet) eq 0 then $
           message, 'Using mean RA and DEC as grid center', /info
        xctr = mean_ra
        yctr = mean_dec
     endif else begin
        message, 'Using known galaxy center as grid center.', /info
        xctr = s.ra_deg
        yctr = s.dec_deg
     endelse
  endif

; FIGURE OUT THE SIZE OF OUR GRID
  if (n_elements(xsize) eq 0) or (n_elements(ysize) eq 0) then begin

;    THE SIZES REQUIRED TO GET ALL OF THE POINTINGS IN THE MAP
     required_x = 2.0 * $
                  (abs((max_ra - xctr)*cos(!dtor*yctr)) > $
                   abs((xctr - min_ra)*cos(!dtor*yctr))) $
                  / pix_scale

     required_y = 2.0 * (abs(max_dec - yctr) > abs(yctr - min_dec)) $
                  / pix_scale

;    ADD PADDING TO ENSURE THE CONVOLUTION DOESN'T RUN OVER THE EDGE
     xsize = ceil(required_x + pad_pix)
     ysize = ceil(required_y + pad_pix)

     if keyword_set(quiet) eq 0 then begin
        message $
           , 'I think a '+str(xsize)+' x '+str(ysize)+' cube is required.' $
           , /info
     endif
  endif

; TURN THIS INTO A BASIC WCS-COMPLIANT HEADER
  target_hdr = make_iram_cube_header( $
               tab_hdr = first_hdr $
               , pix_scale = pix_scale $
               , xsize = xsize $
               , ysize = ysize $
               , xctr = xctr $
               , yctr = yctr)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF REQUESTED, APPLY IMAGE GAIN CORRECTIONS TO THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(apply_gain) then begin
     readcol, gain_file, format='A,A,F,F' $
              , gain_day, gain_pixel, gain_val, gain_unc $
              , count = nlines
     
;    HARDCODE SOME PRUNING
     gain_cap = 1.5
     bad_gains = where(gain_val gt gain_cap or gain_val lt 1./gain_cap, bad_ct)
     if bad_ct gt 0 then begin
        message, 'Some gains outside reasonable range. Capping these at edge of range.', /info
        gain_val[bad_gains] = ((gain_val[bad_gains] > 1./gain_cap) < gain_cap)
     endif

;    APPLY THE GAINS
     for i = 0, nlines-1 do begin
        ind = where(telescop eq gain_pixel[i] and orig_day eq gain_day[i], ct)
        if (ct gt 0) then $
           data[ind,*] = data[ind,*] * gain_val[i]
     endfor
          
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LABEL THE DATA ACCORDING TO WHAT MAP WE THINK IT WILL GO INTO
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  n_spec = (size(data))[1]
  out_map_root = replicate('',n_spec)
  
; IF GRIDDING EVERYTHING TOGETHER, THE ROOT IS SPECIFIED BY THE USER, THEN
; (MAYBE) MODIFIED BY POLARIZATION AND ANGLE STRINGS
  if keyword_set(together) then begin
     message, 'Gridding all data together.', /info

     out_map_root += together_out_root
  endif

; LACKING ANY OTHER DIRECTION, GRID THE INDIVIDUAL INPUT FILES, APPENDING A
; DEFAULT "_CUBE" EXTENSION
  if keyword_set(together) eq 0 and keyword_set(by_coverage) eq 0 then begin
     message, 'Gridding data using same groupings as input spectra.', /info

     out_map_root += '../'+orig_day+'/cubes/'+orig_file + '_cube'
  endif
  
; IF GRIDDING BY COVERAGE *OR SPLITTING BY ANGLE*, USE THE COVERAGE LIST
  if keyword_set(by_coverage) or keyword_set(split_by_angle) then begin
;    INITIALIZE AN ANGLE VECTOR
     angle = scan*0L
     readcol, coverage_list, format='A,A,A,I,I' $
              , cov_name, cov_day, cov_in_file, cov_angle, cov_scan $
              , count = nlines
     
     if keyword_set(by_coverage) then $
        message, 'Gridding data into individual coverages.', /info

     for i = 0L, nlines-1 do begin
        ind = where(orig_day eq cov_day[i] and scan eq cov_scan[i], ct)
        if ct gt 0 then begin
           if keyword_set(by_coverage) then $
              out_map_root[ind] += '../'+orig_day[ind]+'/cubes/'+cov_name[i]
           angle[ind] = cov_angle[i]
        endif
     endfor

     if keyword_set(by_coverage) eq 0 then begin
        uniq_angles = angle[uniq(angle, sort(angle))]
        for i = 0L, n_elements(uniq_angles)-1 do begin
           this_angle = uniq_angles[i]
           angle_str = '_ang'+strcompress(string(long(this_angle)),/remove_all)
           ind = where(angle eq this_angle)           
           out_map_root[ind] += angle_str
        endfor
     endif

  endif

  if keyword_set(split_by_pol) then begin        
     message, 'Splitting data by polarization.', /info

     pol_1_ind = where(pol eq 1, pol_1_ct)
     if pol_1_ct gt 0 then out_map_root[pol_1_ind] += '_pol1'

     pol_2_ind = where(pol eq 2, pol_2_ct)
     if pol_2_ct gt 0 then out_map_root[pol_2_ind] += '_pol2'

     no_pol_ind = where(pol ne 1 and pol ne 2, no_pol_ct)
     if no_pol_ct gt 0 then begin
        out_map_root[no_pol_ind] = ''
        message, 'Some data with no polarization. Stopping.'
     endif
     
  endif

  if n_elements(tag) gt 0 then $
     out_map_root += tag

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BREAK THE DATA DOWN AS REQUESTED AND GRID IT INTO CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  uniq_out_map_roots = out_map_root[uniq(out_map_root, sort(out_map_root))]
  n_cubes = n_elements(uniq_out_map_roots)
  message, 'Gridding '+strcompress(string(n_cubes))+' cubes.', /info

  get_lun, u
  openw, u, cube_list_out

  for i = 0, n_cubes-1 do begin
     ind = where(out_map_root eq uniq_out_map_roots[i])

     if keyword_set(quiet) eq 0 then spawn, 'date'          

     this_hdr = target_hdr
     if keyword_set(by_coverage) or keyword_set(split_by_angle) then $
        sxaddpar, this_hdr, 'SCANANGL', median(angle[ind]), 'MEDIAN, ROUNDED SCAN ANGLE'

     grid_otf, data=data[ind,*] $
               , ra=ra[ind] $
               , dec=dec[ind] $
               , weight=weight[ind] $
               , target = this_hdr $
               , out_root = uniq_out_map_roots[i] $
               , median = use_median $
               , gauss_kern = gauss_kern $
               , gauss_fwhm = gauss_fwhm $
               , show = show

     pol_str = (keyword_set(split_by_pol) eq 0) ? 'both' : strcompress(string(long(median(pol[ind]))))

     if (keyword_set(by_coverage) eq 0) and (keyword_set(split_by_angle) eq 0) then $
        angle_str = 'both' $
     else $
        angle_str = strcompress(string(long(median(angle[ind]))))

     printf, u, uniq_out_map_roots[i]+' '+pol_str+' '+angle_str
  endfor

  close, u

end                             ; of grid_batch
