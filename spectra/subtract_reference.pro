pro subtract_reference $
   , list_file $
   , tag = tag $
   , working_dir = working_dir

;+
;
; Reduction subroutine:
;
; From a structure that already contains a reference (OFF) spectrum, replace
; the working spectrum (SPECTRUM) with ON-OFF, where ON is the RAW
; spectrum. Essentially resets the working spectrum. Primarily useful to
; iterate baseline fitting or fourier flagging without having to relabel and
; rebuild the off spectra.
;
;-
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA SETS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for m = 0, ndata-1 do begin

;    READ THE DATA
     indir = working_dir+'spectra/'
     infile = indir+working_name[m]+'_'+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+ $
                 string(working_name[m])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
     
;    CHECK IF WE ARE ALREADY REFERENCE SUBTRACTED
     need_sub = where(data.on eq 1 and $
                      data.ref_subtracted eq 0 and $
                      data.flagged eq 0, need_sub_ct)
     if need_sub_ct eq 0 then continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; REPLACE THE WORKING SPECTRUM FOR "ON" DATA WITH ON-OFF
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     data[need_sub].spectrum = data[need_sub].raw_spectrum - $
                               data[need_sub].ref_spectrum
     
     data[need_sub].ref_subtracted = 1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUT THE DATA
; (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr

  endfor

  
end
