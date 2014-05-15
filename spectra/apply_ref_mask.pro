pro apply_ref_mask $
   , list_file $
   , tag = tag $
   , working_dir = working_dir $
   , ref_mask_file=ref_mask_file

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CHECK EXISTENCE OF MASK, READ IT, EXTRACT ASTROMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(ref_mask_file) eq 0 then $
     mask_found = 0 $
  else $
     dummy = file_search(working_dir+ref_mask_file $
                         , count=mask_found)

  if mask_found eq 0 then begin
     message, "No reference mask found. Skipping this step.", /info
     return
  endif
  
  mask = readfits(working_dir+ref_mask_file $
                  , mask_hdr, /silent)

; ... NOTE THE DIMENSIONS
  sz_mask = size(mask)

; ... EXTRACT AN ASTROMETRY STRUCTURE
  extast,mask_hdr,astr

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA SETS AND IDENTIFY REFERENCES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin
         
;    READ THE DATA
     indir = working_data+'spectra/'
     infile = indir+working_name[i]+'_'+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
      
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)
     
;    FIND THE POSITIONS OF THE CURRENT DATA IN THE MASK
     ad2xy,data.ra_deg,data.dec_deg,astr,x,y

;    CHECK WHICH COORDS ARE INSIDE THE IMAGE
     in_mask_image = $
        (x ge 0) and (y ge 0) and $
        (x lt sz_mask[1]) and (y lt sz_mask[2])
     
;    FLAG DATA OUTSIDE THE MASK IMAGE AS REFERENCES
     out_of_image = where(in_mask_image eq 0, out_ct)
     if (out_ct) gt 0 then $
        data[out_of_image].on = 0

;    FOR DATA INSIDE THE MAP, REFERENCES ARE SPECTRA WHERE MASK = 0
     in_image = where(in_mask_image eq 1, in_ct)
     if in_ct gt 0 then begin

;       ... WHERE THE MASK IS ZERO, WE ARE "OFF SOURCE"
        mask_is_zero = $
           where(mask[x[in_image], y[in_image]] eq 0, zero_ct)
        if zero_ct gt 0 then $
           data[in_image[mask_is_zero]].on = 0
        
;       ... WHERE THE MASK IS ON, WE ARE "ON SOURCE"
        mask_is_one = $
           where(mask[x[in_image], y[in_image]] eq 1, one_ct)
        if one_ct gt 0 then $
           data[in_image[mask_is_one]].on = 1

     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%    

;    WRITE THE DATA (FIRST DELETE TO KEEP FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
     
  endfor

end                             ; of apply_ref_mask
