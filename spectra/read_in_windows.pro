pro read_in_windows $
   , list_file $
   , tag = tag $
   , window_root = window_root

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE WINDOW MAP AND EXTRACT ASTROMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  win_low = $
     readfits(window_root+'_low.fits', window_hdr, /silent)
  win_high = $
     readfits(window_root+'_high.fits', /silent)

  extast, window_hdr, astrom

  sz_win = size(win_low)

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
      
;    PLACE THE DATA IN THE MASK
     ad2xy, data.ra_deg, data.dec_deg, astrom, win_x, win_y     

     if total(win_x lt 0) or total(win_x ge sz_win[1]) or $
        total(win_y lt 0) or total(win_y ge sz_win[2]) then begin
        message, 'WARNING! Some spectra outside the fitting window.', /info
        win_x = win_x > 0
        win_x = win_x < (sz_win[1]-1)
        win_y = win_y > 0
        win_y = win_y < (sz_win[2]-1)
     endif

;    NOTE THIS WINDOW PAIR FOR EACH LOCAL SPECTRA ...

;    ... NOW PRINT THE CURRENT WINDOW PAIR, SPACE DELIMITED
     window_str = string(win_low[win_x,win_y],format='(F8.1)') + $
                  ' ' + $
                  string(win_high[win_x,win_y],format='(F8.1)') + $
                  ' '
     data.nwindows += 1
     data.windows += window_str

;    WRITE OUT THE DATA
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr
  endfor

end                             ; of init_struct
