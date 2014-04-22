pro split_into_maps $
   , list_file $                       ; INPUT DATA FILE
   , skip_day = skip_day $             ; DAYS TO OMIT
   , bad_data_file = bad_data_file_in $ ; LIST OF DAY PIXEL SCAN TO SKIP
   , output_file = output_file $       ; TEXT FILE WRITTEN BY THE PROGRAM
   , out_root = out_root $             ; STEM OF THE NAME TO GIVE EACH COVERAGE (e.g., "map")
   , first_out = first_out $           ; STARTING NUMBER USE TO LABEL COVERAGES
   , n_coverage = n_coverage $         ; NUMBER OF COVERAGES FOUND
   , map_start_file = map_start_file $ ; A FILE LISTING DAY+KNOWN SCAN STARTS
   , div_angle = div_angle $           ; FORCE/REPORT DIVIDING ANGLE TO SEPARATE SCAN DIRECTIONS
   , ref_subscan = ref_subscan $       ; A LIST OF SUBSCANS KNOWN TO BE REFERENCES
   , pixel = ref_pixel $               ; THE PIXEL TO BE USED IN THE CALCULATION
   , show=show $                       ; FLAG: SHOW OUTPUT ON SCREEN
   , jpeg_name=jpeg_name               

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SET DEFAULTS, CHECK INPUT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
  
; DEFINE THE HERA PIXELS
  pixel_list = '30M-'+ [ $
               'W01-1H01','W01-1H02','W01-1H03', $
               'W01-1H04','W01-1H05','W01-1H06', $
               'W01-1H07','W01-1H08','W01-1H09', $
               'W02-2H01','W02-2H02','W02-2H03', $
               'W02-2H04','W02-2H05','W02-2H06', $
               'W02-2H07','W02-2H08','W02-2H09']
  
  npix = n_elements(pixel_list)

  if n_elements(first_out) eq 0 then first_out = 1

  if n_elements(skip_day) eq 0 then skip_day = ''  

  if n_elements(output_file) eq 0 then output_file = 'coverage_list.txt'

  if n_elements(out_root) eq 0 then out_root = 'map_'

  if n_elements(bad_data_file_in) eq 0 then bad_data_file_in = ''

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF MAPS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, day, map_name, format='A,X,A'
  day = strcompress(day, /remove_all)
  map_name = strcompress(map_name, /remove_all)
  nmaps = n_elements(map_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF THE USER SPECIFIED A LIST OF MAP STARTS, READ THESE IN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, day, map_name, format='A,X,A'
  day = strcompress(day, /remove_all)
  map_name = strcompress(map_name, /remove_all)
  nmaps = n_elements(map_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; OPEN THE OUTPUT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  get_lun, u
  openw, u, output_file

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER INPUT DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for j = 0, nmaps-1 do begin

;    INITIALIZE OUR MAP COUNTER (OR TICK IT BY ONE IF WE ARE CHANGING FILES/DAYS)
     current_map = (j eq 0) ? first_out : current_map+1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIRST  READ IN THE DATA AND NARROW IT DOWN TO A SINGLE PIXEL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    CHECK IF WE ARE SKIPPING THIS DAY
     if total(strcompress(skip_day,/remove_all) eq (day[j])[0]) gt 0 then $
        continue

;    READ THE DATA
     infile = '../'+day[j]+'/calib/'+map_name[j]+'.fits'
     data = mrdfits(infile,1,hdr,/silent)

;    CONVERT ANY OLD-FORMAT TELESCOP ENTRIES TO NEW ONES
     update_telescop_field, data


;    REMOVE BAD DATA (SPECIFIED BY THE USER)
     remove_bad_data $
        , data $
        , day[j] $
        , bad_data_file = bad_data_file_in $
        , all_bad = all_bad

;    MAKE SURE THERE IS SOMETHING LEFT
     if all_bad then $
        continue

;    REMOVE REFERENCE MEASUREMENTS FROM THE DATA
     sz = size(data)

     is_ref = bytarr(sz[1])
     for i = 0, n_elements(ref_subscan)-1 do begin
        ind = where(data.subscan eq ref_subscan[i],ct)
        if ct gt 0 then is_ref[i] = 1B
     endfor
     
     if total(is_ref) eq sz[1] then begin
        message, 'All data flagged as reference. Breaking out of loop...', /info
        continue
     endif
     
;    IF NOT FORCED, THE REFERENCE PIXEL IS THE MOST COMMON PIXEL
     define_ref = 0B
     if n_elements(ref_pixel) eq 0 then $
        define_ref = 1B
     if n_elements(ref_pixel) ne 0 then $
        if total(data.telescop eq ref_pixel) eq 0 then $
           define_ref = 1B
              
     if define_ref then begin
        count = lonarr(npix)
        for i = 0, npix-1 do $
           count[i] = total(data.telescop eq pixel_list[i])
        dummy = max(count, maxind)
        ref_pixel = pixel_list[maxind]
     endif

;    REMOVE REFERENCES AND PIXELS OTHER THAN THE REFERENCE PIXEL
     ind = where((data.telescop eq ref_pixel) and (is_ref eq 0B))
     data = data[ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MEASURE THE SLEW ANGLE FROM EACH POINT TO THE NEXT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    MEASURE STEP BETWEEN SCANS IN R.A.
     delta_ra = data.cdelt2 - shift(data.cdelt2,1)
     delta_ra[0] = delta_ra[1]
     
;    MEASURE STEP BETWEEN SCANS IN DEC
     delta_dec = data.cdelt3 - shift(data.cdelt3,1)
     delta_dec[0] = delta_dec[1]

;    AVOID 180 DEGREE DEGENERACY: 

;    ... IF BOTH NEGATIVE, MAKE THEM POSITIVE
     ind = where(delta_ra le 0.0 and delta_dec le 0.0, ct)
     if ct gt 0 then begin
        delta_ra[ind] *= -1.0
        delta_dec[ind] *= -1.0
     endif

;    ... IF ONE NEGATIVE, MAKE THAT ONE THE RA OFFSET
     ind = where(delta_ra ge 0.0 and delta_dec lt 0.0, ct)
     if ct gt 0 then begin
        delta_ra[ind] *= -1.0
        delta_dec[ind] *= -1.0
     endif
     
;    THE SLEWING ANGLE FOR EACH POINT
     slew_angle = atan(delta_ra, delta_dec)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NOW MEASURE A MEDIAN SCAN ANGLE FOR EACH SCAN (OUR BASIC UNIT)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    IDENTIFY INDIVIDUAL SCANS 
     uniq_scans = $
        (data.scan)[uniq(data.scan, sort(data.scan))]
     n_scans = n_elements(uniq_scans)
     
;    INITIALIZE OUTPUT
     scan_angle = fltarr(n_scans)*!values.f_nan

;    LOOP OVER SCANS
     for i = 0, n_scans-1 do begin
        ind = where((data.scan) eq uniq_scans[i],ct)
        scan_angle[i] = median(slew_angle[ind])
     endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; GUESS AT A DIVIDING LINE BETWEEN ORTHOGONAL COVERAGES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     if n_elements(div_angle) eq 0 then begin
        angle_in_deg = round(scan_angle/!dtor)
        
        uniq_angles = angle_in_deg[uniq(angle_in_deg, sort(angle_in_deg))]
        div_angle = $
           mean(1.0*uniq_angles)*!dtor
     endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; IDENTIFY SCANS THAT BEGIN MAPS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     start_of_map = bytarr(n_scans)
     
     start_of_map[0] = 1B
     for i = 1, n_scans-1 do begin
;    IF WE SWITCHED DIRECTIONS, THIS IS THE START OF A NEW MAP
        if (scan_angle[i-1] gt div_angle) ne (scan_angle[i] gt div_angle) then $
           start_of_map[i] = 1B

;    IF THE USER SPECIFIED THIS AS THE START OF A MAP, IT IS...
        if n_elements(map_start_user) gt 0 then begin
           if total(user_map_start eq uniq_scans[i]) gt 0 then  $
              start_of_map[i] = 1B
        endif
     endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; STEP THROUGH SCANS AND BUILD A LIST OF MAPS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     map_number = lonarr(n_scans)

     map_number[0] = current_map
     for i = 1, n_scans-1 do begin
        if (start_of_map[i]) then current_map +=1
        map_number[i] = current_map
     endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; 
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     for i = 0, n_scans-1 do begin
        this_map_name = out_root+strcompress(string(long(map_number[i])),/remove_all)
        scan_angle_str = string(long(round(scan_angle[i]/!dtor)))
        scan_num_str = string(long(uniq_scans[i]))
        printf, u, this_map_name+' '+day[j]+' '+map_name[j]+' '+scan_angle_str+scan_num_str
     endfor

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SHOW WHAT WE DID
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     if keyword_set(show) then begin
        loadct, 0
        reversect
        circle
        plot, uniq_scans, scan_angle/!dtor, ps=8 $
              , xtitle='!6Scan Number', ytitle='Median Slew Angle' $
              , title='FILE: '+day[j]+' '+map_name[j]
        oplot, uniq_scans, div_angle/!dtor*(1.0+uniq_scans*0.0)
        ind = where(start_of_map,ct)
        for i = 0, ct-1 do $
           oplot, uniq_scans[ind[i]]*[1,1], [-1e6,1e6], lines=1, thick=2
        if n_elements(jpeg_name) gt 0 then begin
           im = tvrd(true=1)
           write_jpeg, jpeg_name, im, true=1
        endif else begin
;           ch = get_kbrd(1)
        endelse
     endif

  endfor                        ; OF THE LOOP OVER INPUT FILES

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CLOSE THE OUTPUT FILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  close, u

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; NOTE THE NUMBER OF COVERAGES IDENTIFIED
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  n_coverage = current_map - first_out + 1

end                             ; of split_into_maps
