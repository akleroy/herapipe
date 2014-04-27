pro hera_base_sub $
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
           , format='X,A', /silent $
           , comment="#"
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

;       REFERENCE-SUBTRACTED, UNFLAGGED, ON SOURCE DATA FOR THIS PIXEL
;        pix_ind = where(data.on eq 1 and $
;                        data.flagged eq 0 and $
;                        data.ref_subtracted eq 1 and $
;                        data.telescop eq pixel_list[j], pix_ct)

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
        win_sz = size(win_arr)

;       BLANK THE DATA INSIDE THE WINDOWS
        for k = 0, pix_ct-1 do begin
           for m = 0, win_sz[1]-1, 2 do begin
              blank = where(vaxis ge win_arr[m,k] and $
                            vaxis le win_arr[m+1,k], blank_ct)
              if blank_ct gt 0 then $
                 pix_data[blank,k] = !values.f_nan
           endfor
        endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY ANY DESIRED SMOOTHING, KEEPING TRACK OF WEIGHTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       INITIALIZE A WEIGHT IMAGE
        blank = where(finite(pix_data) eq 0, blank_ct)
        wt_image = 1.0*finite(pix_data)

;       SMOOTH IN VELOCITY (IF REQUESTED)
        if keyword_set(smooth_in_vel) then begin

           pix_data = smooth(pix_data, [11,1], /nan)
           wt_image = smooth(wt_image, [11,1])

;          ... FIRST A MEDIAN FILTER TO DESPIKE
;           pix_data = median(pix_data, 5)
;           wt_image = smooth(wt_image, [5,1])

;          ... THEN THREE ITERATIONS OF A 7-PIXEL HANNING KERNEL
;           for k = 0, pix_ct-1 do begin
;              spec = pix_data[*,k]
;              wt   = wt_image[*,k]
;              for q = 0, 3 do begin
;                 spec = hans(7,spec,/nan)
;                 wt   = hans(7,wt,/nan)
;              endfor
;              pix_data[*,k] = spec
;              wt_image[*,k] = wt
;           endfor
        endif
        
        if blank_ct gt 0 then wt_image[blank] = 0
        if blank_ct gt 0 then pix_data[blank] = !values.f_nan

;       SMOOTH IN TIME (IF REQUESTED)
        if keyword_set(smooth_in_vel) then begin
           
;          TRANSLATE THE SMOOTHING KERNEL IN SECONDS INTO CHANNELS
           time_kern = ceil(time_smooth_sec / median(data[pix_ind].obstime))

;          SORT DATA BY UT (SHUFFLE IT INSIDE THE STRUCTURE)
           sort_ind = sort(data[pix_ind].ut)
           data[pix_ind] = data[pix_ind[sort_ind]]

;          GET THE TIME STEP BETWEEN SUCCESSIVE DATA DUMPS
           delta_t = data[pix_ind].ut - shift(data[pix_ind].ut,1)
           delta_t[0] = 0.
           
           obs_block = lonarr(pix_ct)
           
;          ... DEFINE BLOCKS BY TIME DISCONTINUITIES
           for k = 1L, pix_ct-1 do $
              obs_block[k] = obs_block[k-1] + (delta_t[k] gt large_t_step)
           
;          ... NOTE THE NUMBER OF BLOCKS
           n_blocks = max(obs_block)+1

           for k = 0L, n_blocks-1 do begin
              counter, k, n_blocks, 'Out of'

;             ... COPY THE DATA
              wt_copy = wt_image
              data_copy = pix_data

;             ... BLANK DATA FROM OTHER TIME BLOCKS
              other_block = where(obs_block ne k)
              wt_copy[*,other_block] = 0.0
              data_copy[*,other_block] = !values.f_nan

;             ... SMOOTH
              data_copy = smooth(data_copy, [1,time_kern], /nan)
;              for m = 0, n_chan-1 do begin   
;                 tvec = reform(data_copy[m,*])
;                 tvec = median(tvec, time_kern)
;                 data_copy[m,*] = tvec
;              endfor
;              wt_copy = smooth(wt_copy, [1,time_kern], /nan)
              
;             ... COPY THESE DATA BACK INTO THE REAL ARRAYS
              this_block = where(obs_block eq k)
              wt_image[*,this_block] = wt_copy[*,this_block]
              pix_data[*,this_block] = data_copy[*,this_block]
           endfor

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
        
        for k = 0, pix_ct-1 do begin                         
        
;          PLOT EVERY 500TH FIT IF REQUESTED
           show_this = (keyword_set(show) and (k mod 1000) eq 0) ? 1B : 0B
           
;          FIT A BASELINE (EXTRAS GET PASSED TO PLOT)
           coeffs = $
              base_fit(xaxis = vaxis $
                       , yaxis = reform(pix_data[*,k]) $
                       , weight = reform(wt_image[*,k]) $
                       , degree = data[pix_ind[k]].base_degree $
                       , show = show_this $
                       , orig_y = data[pix_ind[k]].spectrum $
                       , fit = fit $
                       , status = status)
           
;          IF THE FIT FAILED THEN LEAVE THE COEFFICIENTS BLANK
           if status eq -1 then $
              continue

;          STORE THE FIT COEFS
           data[pix_ind].base_coeffs = ''
           for m = 0, data[pix_ind[k]].base_degree do $              
              data[pix_ind].base_coeffs += string(coeffs[m]) + ' '

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
