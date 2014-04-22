pro hera_base_fit $
   , list_file $
   , tag = tag $
   , degree = degree $
   , edge = edge $
   , smooth_in_time = smooth_in_time $
   , smooth_in_vel = smooth_in_vel $
   , infile = infile $
   , keep_blank=keep_blank $
   , show=show $
   , fts=fts

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET SOME DEFAULT VALUES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(skip_day) eq 0 then skip_day = ''

; DEFINE THE HERA PIXELS
@define_hera_pixels.pro

; THE DEFAULT DEGREE OF POLYNOMIAL FIT
  if n_elements(degree) eq 0 then $
     degree = 1

; FLAG TELLING THE SPECTRUM WHETHER OR NOT TO SMOOTH IN VELOCITY
  if n_elements(smooth_in_vel) eq 0 then $
     smooth_in_vel = 1

; FLAG TELLING THE SPECTRUM WHETHER OR NOT TO SMOOTH IN TIME
  if n_elements(smooth_in_time) eq 0 then $
     smooth_in_time = 1

; NUMBER OF EDGE CHANNELS TO AUTOMATICALLY BLANK
  if n_elements(edge) eq 0 then $
     edge = 11

; DEFINE A "LARGE TIME STEP" IN SECONDS (WE CHECK THE TIME GAP BETWEEN
; SUCCESSIVE DATA DUMPS AND DATA BEFORE/AFTER SUCH A STEP CANNOT BE
; SMOOTHED TOGETHER WHEN BUILDING A BASELINE)
  large_t_step = 21.

; TIME SMOOTHING KERNEL (IN SECONDS)
  time_smooth_sec = 5.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA FILES AND THEN LOOP OVER PIXELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; LOOP OVER DATA FILES
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
 
;    SPECIFY THE DEGREE OF ALL SPECTRA TO BE THE SAME (FOR NOW)
     data.base_degree = degree

;    NUMBER OF ELEMENTS IN ONE SPECTRUM
     n_chan = n_elements(data[0].spectrum)
     
;    AN ARRAY OF CHANNEL NUMBERS
     chan = findgen(n_chan)

;    THE VELOCITY AXIS
     vaxis = chan*data[0].deltav + data[0].v0
     
;    LOOP OVER PIXELS
     for j = 0, npix-1 do begin
        counter, j+1, npix, 'Pixel '

        pix_ind = where(data.telescop eq pixel_list[j], pix_ct)
        
;       IF WE HAVE NO DATA HERE, PROCEED TO THE NEXT PIXEL
        if pix_ct eq 0 then $
           continue

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY ANY DESIRED SMOOTHING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       ELSE EXTRACT THESE DATA TO PLAY WITH
;       ... PROCEDURE IS BLANK ALL DATA NOT INVOLVED IN THE FIT
;       ... THEN ONLY THE FIT DEGREE/COEFFICIENTS WILL BE KEPT
        pix_data = data[pix_ind].spectrum

;       BLANK EDGES TO USER SPECIFICATIONS
        pix_data[0:edge-1,*] = !values.f_nan
        pix_data[(n_chan-edge):*,*] = !values.f_nan

;       BLANK OFF-SOURCE DATA
        off_ind = where(data[pix_ind].on eq 0 or $
                        data[pix_ind].ref_subtracted eq 0 or $
                        data[pix_ind].flagged eq 1, off_ct)
        if off_ct gt 0 then $
           pix_data[*,off_ind] = !values.f_nan

;       EXTRACT THE WINDOWS AS AN ARRAY OF LO-HI x NSPEC 
        win_arr = extract_windows(data[pix_ind])
;        win_arr = win_arr_big[pix_ind,*]
        win_sz = size(win_arr)

;       BLANK THE DATA INSIDE THE WINDOWS
        vaxis_arr = vaxis # (fltarr(pix_ct)+1.)
        for m = 0, win_sz[1]-1, 2 do begin
           win_arr_lo = (fltarr(n_chan) + 1.) # win_arr[m,*]
           win_arr_hi = (fltarr(n_chan) + 1.) # win_arr[m+1,*]
           blank = where(vaxis_arr ge win_arr_lo and $
                         vaxis_arr le win_arr_hi, blank_ct)
           if blank_ct gt 0 then pix_data[blank] = !values.f_nan
        endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY ANY DESIRED SMOOTHING, KEEPING TRACK OF WEIGHTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       INITIALIZE A WEIGHT IMAGE
        blank = where(finite(pix_data) eq 0, blank_ct)
        wt_image = 1.0*finite(pix_data)

;       SMOOTH IN VELOCITY (IF REQUESTED)
        if keyword_set(smooth_in_vel) then begin

;          ... A MEDIAN FILTER TO DESPIKE
           for jj = 0L, pix_ct-1 do $
              pix_data[*,jj] = median(pix_data[*,jj], 5)
           wt_image = smooth(wt_image, [5,1])
           
        endif
        
        if blank_ct gt 0 then wt_image[blank] = 0
        if blank_ct gt 0 then pix_data[blank] = !values.f_nan

;       SMOOTH IN TIME (IF REQUESTED)
        if keyword_set(smooth_in_time) then begin
           
;          TRANSLATE THE SMOOTHING KERNEL IN SECONDS INTO CHANNELS
           time_kern = ceil(time_smooth_sec / median(data[pix_ind].obstime))

;          SORT DATA BY UT (SHUFFLE IT INSIDE THE STRUCTURE)
           sort_ind = sort(data[pix_ind].ut)
           data[pix_ind] = data[pix_ind[sort_ind]]

           pix_data = smooth(pix_data, [1,time_kern], /nan)
           wt_image = smooth(wt_image*1.0, [1,time_kern], /nan)

        endif

        not_one = where(wt_image lt 1.0, not_one_ct)
        if not_one_ct gt 0 then begin
           wt_image[not_one] = 0
           pix_data[not_one] = !values.f_nan
        endif

        if blank_ct gt 0 then wt_image[blank] = 0
        if blank_ct gt 0 then pix_data[blank] = !values.f_nan

;       THE WEIGHT IMAGE IS NOW KIND OF PROPORTIONAL TO THE NUMBER OF
;       SORT-OF-INDEPENDENT CHANNELS/DUMPS AVERAGED INTO THE
;       POINT. WE'D EXPECT THE RMS IN THE POINT TO BE TWIDDLE
;       THE SQUARE ROOT OF THIS, SO TAKE ROOT WEIGHT IMAGE.
        wt_image = sqrt(wt_image)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FIT AND RECORD THE RESULTS IN THE STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%        

;      USE FOUR PANELS TO VIEW FITS
        if keyword_set(show) then $
           !p.multi=[0,2,2]
        
        for k = 0L, pix_ct-1 do begin                         
           if total(finite(pix_data[*,k])) eq 0 then $
              continue

;          PLOT EVERY 500TH FIT IF REQUESTED
           show_this = (keyword_set(show) and (k mod 1000) eq 0) ? 1B : 0B
           
           if total(reform(pix_data[*,k]),/nan) eq 0 then stop

;          FIT A BASELINE (EXTRAS GET PASSED TO PLOT)
              coeffs = $
                 base_fit(xaxis = vaxis $
                          , yaxis = reform(pix_data[*,k]) $
                          , flag = reform(wt_image[*,k]) $
                          , degree = data[pix_ind[k]].base_degree $
                          , show = show_this $
                          , orig_y = data[pix_ind[k]].spectrum $
                          , status = status)
           
;          IF THE FIT FAILED THEN LEAVE THE COEFFICIENTS BLANK
              if status eq -1 then $
                 stop
              
;          STORE THE FIT COEFS
           data[pix_ind[k]].base_coeffs = ''
           for m = 0, data[pix_ind[k]].base_degree do $              
              data[pix_ind[k]].base_coeffs += string(coeffs[m]) + ' '

           if strcompress(data[pix_ind[k]].base_coeffs,/remove_all) eq '' $
              then stop

        endfor                  ; LOOP OVER SPECTRA

;    UNSET PANELS
     if keyword_set(show) then $
        !p.multi=[0,1,1]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SAVE OUR RESULTS (FIRST IN THE BIG DATA STRUCTURE, THEN WRITE)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%        
           
     endfor                     ; LOOP OVER PIXELS

;    WRITE PROCESSED DATA TO DISK
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
        outfile = infile
        spawn, 'rm '+outfile
        mwrfits, data, outfile, hdr

     endfor
     

     
  end
