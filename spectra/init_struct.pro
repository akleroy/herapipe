pro init_struct $
   , list_file $
   , tag = tag $
   , working_dir = working_dir

;+
;
; Reduction subroutine:
;
; Essentially a "filler" to make sure that the 30m table read in from class
; has all of the fields that will be used later in the pipeline.
;
; Loops over a list of files and performs the operation on each.
;
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, orig_name, working_name $
           , format='A,A', /silent $
           , comment="#"
  orig_name = strcompress(orig_name, /remove_all)
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(orig_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND INITIALIZE A NEW STRUCTURE FILE FOR EACH
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin
     
;    READ THE DATA, CHECK THAT WE HAVE A FILE
     infile = orig_name[i]
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(orig_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
     
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)
     n_spec = n_elements(data[0].spectrum)

;    MAKE A COPY OF THE DATA STRUCTURE THAT HAS NEW FIELDS
     new_data = $
        replicate( $
        create_struct(data[0] $
                      , 'RA_DEG', 0.0d $
                      , 'DEC_DEG', 0.0d $
                      , 'V0', 0.0d $
                      , 'DELTAV', 0.0d $
                      , 'SCAN_ANG', 0.0d $
                      , 'SUBSCAN_START_UT', 0.0d $
                      , 'SUBSCAN_STOP_UT', 0.0d $
                      , 'ON', 1 $
                      , 'REF_SUBTRACTED', 0 $
                      , 'TREF_BEFORE', 0.0 $
                      , 'TIME_TO_REF_BEFORE', 0.0 $
                      , 'TREF_AFTER', 0.0 $
                      , 'TIME_TO_REF_AFTER', 0.0 $
                      , 'FLAGGED', 0 $
                      , 'WHY_FLAGGED', '' $
                      , 'BASE_DEGREE', 0 $
                      , 'BASE_COEFFS', '' $
                      , 'NWINDOWS', 0 $
                      , 'WINDOWS', '' $
                      , 'FIT_RMS', 0.0*!values.f_nan $
                      , 'DELTA_WIDE', 0.0*!values.f_nan $
                      , 'DELTA_NARROW', 0.0*!values.f_nan $
                      , 'LONGEST_RIPPLE', 0.0*!values.f_nan $
                      , 'RAW_SPECTRUM', fltarr(n_spec)*!values.f_nan $
                      , 'REF_SPECTRUM', fltarr(n_spec)*!values.f_nan $
                     ), sz[1])
     
;    COPY THE VALUES FROM THE OLD STRUCTURE TO THE NEW STRUCTURE
     struct_assign, data, new_data, /nozero

;    NOW OVERWRITE THE OLD DATA WITH THE NEW DATA
     data = new_data

;    COPY THE DATA INTO THE RAW SPECTRUM FIELD (SPECTRUM WILL CHANGE)
     data.raw_spectrum = data.spectrum

;    ADD A SUBSCAN FIELD IF WE DON'T HAVE ONE (VARIES)
     fields = tag_names(data)
     if total(strupcase(fields) eq 'SUBSCAN') eq 0 then begin
        new_data = $
           replicate(create_struct(data[0], 'SUBSCAN', 1), sz[1])        
        struct_assign, data, new_data     
        data = new_data
     endif

;    CONVERT ANY OLD-FORMAT TELESCOP ENTRIES TO NEW ONES
     update_telescop_field, data

;    WRITE OUT THE DATA 
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outdir = working_dir+'spectra/'
     outfile = outdir+working_name[i]+'_'+tag+'.processed.fits'
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of init_struct
