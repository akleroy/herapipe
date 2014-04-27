pro measure_subscan_timing $
   , list_file $
   , tag = tag $
   , fts = fts
  
  @define_hera_pixels.pro

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND CALCULATE TIMING FOR EACH SUBSCAN 
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
     
;    IDENTIFY UNIQUE SUBSCANS
     uniq_subscans = $
        (data.subscan)[uniq(data.subscan, sort(data.subscan))]
     n_subscans = n_elements(uniq_subscans)

;    LOOP OVER ALL SCAN/SUBSCAN COMBOS AND NOTE THE START/STOP TIME

     for j = 0, n_scans-1 do begin

        for k = 0, n_subscans-1 do begin

           this_ss = $
              where((data.scan eq uniq_scans[j]) and $
                    (data.subscan eq uniq_subscans[k]), ss_ct)
           
           this_ss_one_pix = $
              where((data.scan eq uniq_scans[j]) and $
                    (data.subscan eq uniq_subscans[k]) and $
                    (data.telescop eq pixel_list[0]), one_pix_ct)
           
           if ss_ct eq 0 then continue

           if one_pix_ct eq 0 then begin
              this_ss_all_pix = $
                 where((data.scan eq uniq_scans[j]) and $
                       (data.subscan eq uniq_subscans[k]))
              new_ref_pix = data[this_ss_all_pix[0]].telescop
              this_ss_one_pix = $
                 where((data.scan eq uniq_scans[j]) and $
                       (data.subscan eq uniq_subscans[k]) and $
                       (data.telescop eq new_ref_pix), one_pix_ct)
           endif 

           if one_pix_ct eq 0 then $
              stop

;          GRAB THE TIME STAMP AND SORT
           subscan_ut = (data.ut)[this_ss_one_pix]
           subscan_ut = subscan_ut[sort(subscan_ut)]

;          ... FUDGE THIS TO AVOID MUCH PAIN: USE THE 2ND AND NEXT TO LAST
;          SCANS RATHER THAN 1ST AND LAST. WHEN WE BREAK IT'S ON THE
;          FIRST AND LAST.


;          ONLY WORKS IF WE HAVE AT LEAST 3 ELEMENTS
           if one_pix_ct lt 3 then begin
              subscan_start = subscan_ut[0]
              subscan_stop = subscan_ut[one_pix_ct-1]
           endif else begin
              subscan_start = subscan_ut[2]        ; min(data[this_ss].ut)
              subscan_stop = subscan_ut[one_pix_ct-2] ; max(data[this_ss].ut)
           endelse

           data[this_ss].subscan_start_ut = subscan_start
           data[this_ss].subscan_stop_ut = subscan_stop
           
        endfor

     endfor

;    WRITE OUT THE DATA
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of measure_subscan_timing
