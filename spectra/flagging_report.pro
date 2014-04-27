pro flagging_report $
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
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)
 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND EXTRACT RELEVANT DATA FOR EACH FILE
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
     
     fname= '../reports/flagging_'+working_name[i]+'.txt'
     openw, 1,  fname

     for j = 0, npix-1 do begin

;       NUMBER OF KEPT AND REJECTED SPECTRA
        good = total(data.telescop eq pixel_list[j] and $
                     data.flagged eq 0 and $
                     data.on eq 1 and $
                     data.ref_subtracted eq 1)

        bad = total(data.telescop eq pixel_list[j] and $
                     data.flagged eq 1 and $
                     data.on eq 1 and $
                     data.ref_subtracted eq 1)

        total_spec = 1.0*good + 1.0*bad

;       REPORT ON THIS PIXEL
        printf, 1, 'Pixel ', pixel_list[j]

;       CATCH THE PATHOLOGICAL CASE
        if good eq 0 then begin
           printf,1, '... no useable data for pixel.'
           continue
        endif

;       ELSE DETAIL THE REJECTION
        rej_line = '... rejected '+sigfig(bad*1.0/total_spec*100.,3)+' % ('
        
        uneven = total(data.telescop eq pixel_list[j] and $
                        data.flagged eq 1 and $
                        (strpos(data.why_flagged, 'US') ne -1 or $
                         strpos(data.why_flagged, 'UA') ne -1) and $
                        data.on eq 1 and $
                        data.ref_subtracted eq 1)
        rej_line += 'uneven '+sigfig(uneven*1.0/total_spec*100,3)+' %, '

        ripples = total(data.telescop eq pixel_list[j] and $
                        data.flagged eq 1 and $
                        (strpos(data.why_flagged, 'RS') ne -1 or $
                         strpos(data.why_flagged, 'RA') ne -1) and $
                        data.on eq 1 and $
                        data.ref_subtracted eq 1)
        rej_line += 'ripples '+sigfig(ripples*1.0/total_spec*100,3)+' %, '

        noisy = total(data.telescop eq pixel_list[j] and $
                      data.flagged eq 1 and $
                      (strpos(data.why_flagged, 'NS') ne -1 or $
                       strpos(data.why_flagged, 'NA') ne -1) and $
                      data.on eq 1 and $
                      data.ref_subtracted eq 1)
        rej_line += 'noisy '+sigfig(noisy*1.0/total_spec*100,3)+' %, '

        path = total(data.telescop eq pixel_list[j] and $
                     data.flagged eq 1 and $
                     (strpos(data.why_flagged, 'P') ne -1) and $
                     data.on eq 1 and $
                     data.ref_subtracted eq 1)
        rej_line += 'pathological '+sigfig(path*1.0/total_spec*100,3)+' %, '


        hand = total(data.telescop eq pixel_list[j] and $
                     data.flagged eq 1 and $
                     (strpos(data.why_flagged, 'H') ne -1) and $
                     data.on eq 1 and $
                     data.ref_subtracted eq 1)
        rej_line += 'hand '+sigfig(hand*1.0/total_spec*100,3)+' %) '

        printf,1,rej_line

     endfor

     close, 1

     spawn, 'cat '+fname

  endfor

end                             ; of flagging_report
