pro measure_scan_angle $
   , list_file $
   , tag = tag

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE THE HERA PIXELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  @define_hera_pixels.bat

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND CALCULATE COORDINATES FOR EACH
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

;    IDENTIFY UNIQUE SCANS
     uniq_scans = $
        (data.scan)[uniq(data.scan, sort(data.scan))]
     n_scans = n_elements(uniq_scans)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER SCANS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     for j = 0, n_scans-1 do begin

        scan_ind = where(data.scan eq uniq_scans[j], scan_ct)

        ref_pix = data[scan_ind[0]].telescop

        sub_ind = where(data.scan eq uniq_scans[j] and $
                        data.telescop eq ref_pix, sub_ct)

;       TRAP THE CASE OF ONLY A SINGLE DUMP IN THE SCAN
        if sub_ct eq 1 then $
           continue

        sub_data = data[sub_ind]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MEASURE THE SLEW ANGLE FROM EACH POINT TO THE NEXT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       MEASURE STEP BETWEEN SCANS IN R.A.
        delta_ra = sub_data.cdelt2 - shift(sub_data.cdelt2,1)
        delta_ra[0] = delta_ra[1]
        
;       MEASURE STEP BETWEEN SCANS IN DEC
        delta_dec = sub_data.cdelt3 - shift(sub_data.cdelt3,1)
        delta_dec[0] = delta_dec[1]
        
;       AVOID 180 DEGREE DEGENERACY: 
        
;       ... IF BOTH NEGATIVE, MAKE THEM POSITIVE
        ind = where(delta_ra le 0.0 and delta_dec le 0.0, ct)
        if ct gt 0 then begin
           delta_ra[ind] *= -1.0
           delta_dec[ind] *= -1.0
        endif

;       ... IF ONE NEGATIVE, MAKE THAT ONE THE RA OFFSET
        ind = where(delta_ra ge 0.0 and delta_dec lt 0.0, ct)
        if ct gt 0 then begin
           delta_ra[ind] *= -1.0
           delta_dec[ind] *= -1.0
        endif
        
;       THE SLEWING ANGLE FOR EACH POINT
        slew_angle = atan(delta_ra, delta_dec)

;       THE MEDIAN SLEW ANGLE OVER THE SCAN
        scan_slew_angle = median(slew_angle)       

        data[scan_ind].scan_ang = scan_slew_angle/!dtor
     endfor

;    WRITE OUT THE DATA
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor


end                             ; of measure_scan_angle
