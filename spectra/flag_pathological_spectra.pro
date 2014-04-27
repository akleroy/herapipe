pro flag_pathological_spectra $
   , list_file $
   , tag = tag $
   , allow_high_tsys = allow_high_tsys

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

  if keyword_set(allow_high_tsys) then $
     high_tsys = 1300. $
  else $
     high_tsys = 700.

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

;    CHECK FOR BLANK SPECTRA
;    (ARE THERE OTHER SINGLE-VALUES THAT I HAVE TO WORRY ABOUT?)
     blank = where(min(data.spectrum,dim=1,/nan) eq $
                   max(data.spectrum,dim=1,/nan), blank_ct)
     if blank_ct gt 0 then begin
        data[blank].flagged = 1
        data[blank].why_flagged += 'P '
     endif

;    CHECK FOR NON-FINITE SPECTRA
     not_finite_ind = where(total(finite(data.spectrum),1) eq 0 $
                            , not_finite_ct)
     if not_finite_ct gt 0 then begin
        data[not_finite_ind].flagged = 1
        data[not_finite_ind].why_flagged += 'P '
     endif

;    CUT ON TSYS BEING TOO HIGH
     high_tsys_ind = where(data.tsys gt high_tsys, high_tsys_ct)
     if high_tsys_ct gt 0 then begin
        data[high_tsys_ind].flagged = 1
        data[high_tsys_ind].why_flagged += 'P '
     endif

;    CUT ON TSYS BEING TOO LOW
     low_tsys_ind = where(data.tsys le 0., low_tsys_ct)
     if low_tsys_ct gt 0 then begin
        data[low_tsys_ind].flagged = 1
        data[low_tsys_ind].why_flagged += 'P '
     endif

;    WRITE OUT THE DATA TO THE SAME FILE
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of flag_pathological_spectra
