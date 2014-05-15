function build_header $
   , list_file $
   , tag = tag $
   , galaxy = galaxy $
   , use_mean_ctr = use_mean_ctr $
   , pix_scale=pix_scale $
   , split_by_angle = split_by_angle $
   , uniq_angles = uniq_angles

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(pix_scale) eq 0 then pix_scale = 2.0/3600.

  pad_pix = 10                  ; pixel padding when making the cube

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

;    SAVE THE FIRST HEADER
     if i eq 0 then first_hdr = hdr

     ind = where(data.flagged eq 0 and $
                 data.on eq 1 and $
                 data.ref_subtracted eq 1, ct)

;    CHECK THAT WE'VE STILL GOT SOMETHING
     if ct eq 0 then continue    

;    BUILD US A GIANT HONKIN' ARRAY OF SPECTRA POSITIONS!
     if n_elements(ra) eq 0 then begin
        ra = data[ind].ra_deg
        dec = data[ind].dec_deg
        ang = round(data[ind].scan_ang)
     endif else begin
        ra = [ra, data[ind].ra_deg]
        dec = [dec, data[ind].dec_deg]
        ang = [ang, round(data[ind].scan_ang)]
     endelse

  endfor     

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT THE REQUIRED SIZE AND EXTENT FOR THE CUB
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

  print, 'Total extent in RA (arcseconds):', delta_ra
  print, 'Total extent in DEC (arcseconds):',  delta_dec
  print, 'Mean RA (degrees):', mean_ra
  print, 'Mean DEC (degrees):',  mean_dec

; WORK OUT THE CENTER OF THE MAP
  if (n_elements(xctr) eq 0) or (n_elements(yctr) eq 0) then begin
     if n_elements(galaxy) eq 0 then $
        s = -1 $
     else $
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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; TURN THIS INTO A BASIC WCS-COMPLIANT HEADER AND RETURN IT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  target_header = make_iram_cube_header( $
           tab_hdr = first_hdr $
           , pix_scale = pix_scale $
           , xsize = xsize $
           , ysize = ysize $
           , xctr = xctr $
           , yctr = yctr)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF WE WANT TO PLAIT, WORK OUT THE SCAN ANGLES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(split_by_angle) then begin
     uniq_angles = ang[uniq(ang, sort(ang))]
     n_ang = n_elements(uniq_angles)     
     message, 'I found '+strcompress(string(n_ang),/remove)+$
              ' unique scan angles.', /info     
  endif

  return, target_header

end                             ; of build_header
