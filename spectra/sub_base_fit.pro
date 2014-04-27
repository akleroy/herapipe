pro sub_base_fit $
   , list_file $
   , tag = tag $
   , show = show
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
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
          
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)

;    NUMBER OF ELEMENTS IN ONE SPECTRUM
     n_chan = n_elements(data[0].spectrum)
     
;    AN ARRAY OF CHANNEL NUMBERS
     chan = findgen(n_chan)

;    THE VELOCITY AXIS
     vaxis = chan*data[0].deltav + data[0].v0

;    IDENTIFY ON-SOURCE, UNFLAGGED DATA
     fit_ind = where(data.on eq 1 and $
                     data.flagged eq 0 and $
                     data.ref_subtracted eq 1, fit_ct)

     if keyword_set(show) then $
        !p.multi=[0,2,2]

;    LOOP OVER SPECTRA
     for j = 0L, fit_ct-1 do begin
        coeffs = dblarr(data[fit_ind[j]].base_degree+1)
        if data[fit_ind[j]].base_coeffs eq '' then continue           
        reads, data[fit_ind[j]].base_coeffs, coeffs
        base = double(vaxis*0.0)
        for k = 0, data[fit_ind[j]].base_degree do $
           base += coeffs[k]*double(vaxis)^k     

        if keyword_set(show) then begin 
           if j mod 500 eq 0 then begin
              plot, vaxis, data[fit_ind[j]].spectrum, ps=1   
              oplot, vaxis, base, color=getcolor('red'), thick=3
;           ch = get_kbrd(1)
           endif
        endif
      
        data[fit_ind[j]].spectrum -= base
     endfor     

     !p.multi=0

;    WRITE OUT THE DATA
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of init_struct
