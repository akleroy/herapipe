pro show_ref $
   , list_file $
   , tag = tag $
   , pause = do_pause $
   , delay = delay $
   , report = report

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#", /silent
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA SETS AND IDENTIFY REFERENCES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin

;    FLAG TELLING US TO WRITE
     need_to_write = 0

;    SHOW THE LOCATIONS OF THE ON/OFF MEASUREMENTS
     loadct, 0, /silent
     reversect
     
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

;    PLOT SPECTRA AS DOTS
     plot, data.cdelt2, data.cdelt3, ps=3, ystyle=16 $
           , title='ON/OFF location for '+data_name[i]
     off_ind = where(data.on eq 0, off_ct)
     if off_ct gt 0 then $
        oplot, data[off_ind].cdelt2, data[off_ind].cdelt3 $
               , ps=1, color=getcolor('red')               
     
     circle, /fill
     legend, /bottom, /right, box=0, symsize=1.5*[1,1] $
             , color=[getcolor('red'),getcolor('black')] $
             , ['OFF Spectrum', 'ON Spectrum'] $
             , textcolor=getcolor('black') $
             , psym=[8,8]
     

;    AFTER PLOTTING WAIT, STOP, AND SAVE AS REQUESTED...
     
     if keyword_set(pause) then begin
        message, 'Pausing for inspection. Please hit a key to continue.', /info
        ch = get_kbrd(1)
     endif

     if keyword_set(delay) then begin
        wait, 1.0
     endif

     if keyword_set(report) then begin
        im = tvrd(true=1)
        write_jpeg, '../reports/onoff_'+working_name[i]+'.jpeg', im, true=1
     endif

  endfor

end                             ; of show_ref

