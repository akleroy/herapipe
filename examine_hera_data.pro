pro examine_hera_data, file $
                       , flagged=flagged $
                       , use_windows=use_windows $
                       , sm=sm $
                       , median=use_med $
                       , pos = pos $
                       , rad = rad $
                       , fts = fts


;+
; NAME:
;
; examine_hera_data
;
; PURPOSE:
;
; An sort-of-interactive tool to browse through a HERA spectra data
; structure. Intended to be used to identify and flag pathological pixels,
; scans, etc. and to get an idea of whether an increase in baseline order is
; needed.
;
; CATEGORY:
;
; Data reduction tool.
;
; CALLING SEQUENCE:
;
; examine_hera_data, "some_hera_data.fits", sm=[25], /blank
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-

  @define_hera_pixels.pro

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SOME SMART DEFAULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if n_elements(sm) eq 0 then sm = 25

  if n_elements(use_windows) eq 0 then use_windows = 1

  if n_elements(flagged) eq 0 then flagged = 0

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; READ THE DATA AND EXTRACT THE PART WE ARE INTERESTED IN
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; ... DATA FILE
  t = mrdfits(file,1,hdr)

; ... VELOCITY AXIS
  nchan = n_elements(t[0].spectrum)

  vaxis = t[0].v0 + findgen(nchan)*t[0].deltav

; ... KEEP ONLY FLAGGED/UNFLAGGED, REFERENCE-SUBTRACTED, ON-SORUCE DATA
  if keyword_set(flagged) eq 0 then begin
     ind = where(t.flagged eq 0 and t.on eq 1 and $
                 t.ref_subtracted eq 1, ct)
  endif else begin
     ind = where(t.flagged eq 1 and t.on eq 1 and $
                 t.ref_subtracted eq 1, ct)
  endelse

; ... CRASH IF EMPTY
  if ct eq 0 then return

  t = t[ind]  

; ... FOCUS ON A SPECIFIC POSITION
  if n_elements(pos) gt 0 then begin
     dist = sphdist(t.ra_deg, t.dec_deg, pos[0], pos[1], /deg)
     if n_elements(rad) eq 0 then $
        rad = 15.
     ind = where(dist lt rad*1.0/3600., ct)
     if ct eq 0 then return
     t = t[ind]
  endif

; REDUCE PIXEL LIST TO PIXELS WITH DATA
  keep_pix = bytarr(npix)+1
  for i = 0, npix-1 do begin

     if total(t.telescop eq pixel_list[i]) eq 0 then $
        keep_pix[i] = 0
     
  endfor

  pixel_list = pixel_list[where(keep_pix)]
  npix = n_elements(pixel_list)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; READ THE SPECTRAL FITTING WINDOWS INTO A MASK
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  blank_mask = (t.spectrum eq t.spectrum)*0B

  if keyword_set(use_windows) then begin

;       EXTRACT THE WINDOWS AS AN ARRAY OF LO-HI x NSPEC 
        win_arr = extract_windows(t)
        win_sz = size(win_arr)

;       BLANK THE DATA INSIDE THE WINDOWS
        vaxis_arr = vaxis # (fltarr(ct)+1.)
        for m = 0, win_sz[1]-1, 2 do begin
           win_arr_lo = (fltarr(nchan) + 1.) # win_arr[m,*]
           win_arr_hi = (fltarr(nchan) + 1.) # win_arr[m+1,*]
           blank = where(vaxis_arr ge win_arr_lo and $
                         vaxis_arr le win_arr_hi, blank_ct)
           if blank_ct gt 0 then blank_mask[blank] = 1B
        endfor
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; BROWSE DATA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; A NEW WINDOW
  window, xsize = 1000, ysize=1000

; CURRENT STATUS OF BROWSING
  all_at_once = 1B
  start_index = 0L
  index_increment = 1000
  pixel = 0
  blank = use_windows

; INITIALIZE USER INPUT  
  inp = ""
  refresh = 1

; INPUT LOOP
  while strupcase(inp) ne "Q" do begin

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    GET INPUT FROM THE MOUSE
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     if strupcase(inp) eq "M" then begin
        cursor, x, y
        print, "Channel ", x, " Spectrum ", y
     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    SWAP TO/FROM AN ALL-AT-ONCE VIEW
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     if strupcase(inp) eq "A" then begin
        if all_at_once eq 1 then begin
           all_at_once = 0
           start_index = 0
           refresh = 1
        endif else begin
           all_at_once = 1
           start_index = 0
           refresh = 1
        endelse
     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    INCREMENT OR DECREMENT THE SCROLL
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     if inp eq "+" then begin
        if start_index lt (data_ct - index_increment - 1) then begin
           start_index = (start_index + index_increment) < $
                         (data_ct - index_increment - 1)
           refresh = 1
        endif
     endif

     if inp eq "-" then begin
        if start_index gt 0 then begin
           start_index = (start_index - index_increment) > 0
           refresh = 1
        endif
     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    INCREMENT OR DECREMENT PIXEL BEING VIEWED
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     if strupcase(inp) eq "P" then begin
        if pixel gt 0 then begin
           pixel -= 1
           start_index = 0
           refresh = 1
        endif
     endif

     if strupcase(inp) eq "N" then begin
        if pixel lt npix-1 then begin
           pixel += 1
           start_index = 0
           refresh = 1
        endif
     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    REFRESH DSPLAY IF NECESSARY     
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

     if keyword_set(refresh) then begin

;       ... GET DATA FOR RELEVANT PIXEL
        data_ind = where(t.telescop eq pixel_list[pixel], data_ct)        
        im = t[data_ind].spectrum

;       ... BLANK THE DATA IF NEEDED
        if keyword_set(blank) then begin
           blank_ind = where(blank_mask[*,data_ind], blank_ct)
           if blank_ct gt 0 then $
              im[blank_ind] = !values.f_nan
        endif

;       ... BUILD TRIVIAL AXES
        sz = size(im)
        xaxis = findgen(sz[1])
        yaxis = findgen(sz[2])

;       ... EXTRACT A SUB-IMAGE IF THAT'S WHAT WE'RE DOING
        if all_at_once eq 0 then begin
           y_lo = start_index
           y_hi = (start_index + index_increment) < (data_ct-1)
           im = im[*,y_lo:y_hi]
           yaxis = yaxis[y_lo:y_hi]
        endif

;       ... GET THE INTEGRATED SPECTRUM       
        int_spec = total(im,2,/nan)
        if keyword_set(use_med) then int_spec = median(im,dim=2)

        erase

;       ... DISPLAY
        !p.multi=[1,2,2]
        loadct, 5, /silent
        disp, im, xaxis, yaxis $
              , min=-0.5, max=0.5 $
              , title=pixel_list[pixel]+' '+file

        !p.multi=[3,2,2]
        loadct, 5, /silent
        plot, int_spec

        !p.multi=[0,2,1]
        loadct, 5, /silent
        disp, im gt 0, xaxis, yaxis $
              , min=0, max=1 $
              , title=pixel_list[pixel]+' '+file $
              , /noerase

        refresh = 0

     endif

;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;    GET USER INPUT
;    -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
     print, "Options:"
     print, "(A) - Toggle all spectra / subset of spectra"
     print, "(C) - Click mouse for info on target spectrum"
     print, "(N) - Next pixel"
     print, "(P) - Previous pixel"
     print, "(+) - Scroll up in subset mode"
     print, "(-) - Scroll down in subset mode"
     print, "(Q) - Quit"

     inp = get_kbrd(1)

  endwhile
  
end
