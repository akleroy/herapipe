pro fill_in_coords $
   , list_file $
   , tag = tag $
   , force_v0 = force_v0
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  working_name = strcompress(working_name, /remove_all)
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
     
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)

;    BUILD SKY COORDINATES (COMPLICATED BY CONVENTION CHANGE)
     if sxpar(hdr,'CRVAL3') eq 0.0 then begin
        dec = double(data.crval3 + data.cdelt3)
     endif else begin
        dec = double(sxpar(hdr,'CRVAL3') + data.cdelt3)
     endelse

     median_dec = median(dec)

     if sxpar(hdr,'CRVAL2') eq 0.0 then begin
        ra = double(data.crval2) + $
             double(data.cdelt2/cos(!dtor*median_dec))
     endif else begin
        ra = double(sxpar(hdr,'CRVAL2')) + $
             double(data.cdelt2/cos(!dtor*median_dec))        
     endelse

     data.ra_deg = ra
     data.dec_deg = dec

;    KLUGE TO FORCE THE ZERO-POINT VELOCITY (NEEDED FOR OLDER DATA)
     if n_elements(force_v0) ne 0 then $
        sxaddpar, hdr, 'VELO-LSR', force_v0, 'SET BY HAND, BE CAREFUL!'

;    PARSE THE HEADER TO CALCULATE THE VELOCITY AXIS  
     crval = sxpar(hdr,'VELO-LSR')
     crpix = sxpar(hdr,'CRPIX1')
     cdelt = sxpar(hdr,'DELTAV')
     
     v = findgen(n_elements(data[0].spectrum))
     vdif = v - (crpix-1.0)
     vaxis = vdif * cdelt + crval
     if abs(cdelt) gt 100. then $
        vaxis /= 1e3    
;    ... UNITS OF VAXIS SHOULD NOW BE V_LSR, KM/S

;    PUT IN VELOCITY INFO (KM/S)
     data.v0 = vaxis[0]
     data.deltav = vaxis[1] - vaxis[0]
     
;    WRITE OUT THE DATA TO THE SAME FILE
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of fill_in_coords
