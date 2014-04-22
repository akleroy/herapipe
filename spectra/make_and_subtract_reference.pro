pro make_and_subtract_reference $
   , list_file $
   , tag = tag $
   , equal_weight = equal_weight $
   , median = use_median $
   , min_ref = min_ref_size $
   , max_ref = max_ref_size
  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SET DEFAULTS, CHECK INPUT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; DEFINE THE HERA PIXELS
@define_hera_pixels.bat
  
; DEFINE A "LARGE TIME STEP" IN SECONDS (WE CHECK THE TIME GAP BETWEEN
; SUCCESSIVE DATA DUMPS AND DATA BEFORE/AFTER SUCH A STEP CANNOT BE
; PART OF THE SAME REFERENCE MEASUREMENTS)
  large_t_step = 21.

; DEFINE THE MAXIMUM REFERENCE TIME
; (ONLY THIS MUCH REFERENCE DATA IN EACH DIRECTION WILL BE USED)
  if n_elements(max_ref_size) eq 0 then $
     max_ref_size = 20.

; ... AND A CORRESPONDING MINIMUM REFERENCE TIME
; (BLOCKS BELOW THIS SIZE CANNOT MAKE A REFERENCE)
  if n_elements(min_ref_size) eq 0 then $
     min_ref_size = 2.
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA SETS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for m = 0, ndata-1 do begin

;    READ THE DATA
     indir = '../spectra/'
     infile = indir+working_name[m]+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[m])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
      
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)

;    NUMBER OF ELEMENTS IN ONE SPECTRUM
     nchan = n_elements(data[0].spectrum)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER PIXELS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     for i = 0L, npix - 1 do begin    
        
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... FIRST COLLECT THE DATA FOR THIS PIXEL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    FIND THE DATA FOR THIS PIXEL
        pix_ind = where((data.telescop) eq pixel_list[i] $
                        ,pix_ct)        
        if pix_ct eq 0 then continue
        
;    EXTRACT THE KEY INFO HERE: SPECTRA, REFERENCE FLAG, CALIBRATED FLAG
        pix_data = data[pix_ind]

;    SORT DATA BY UT
        sort_ind = sort(pix_data.ut)
        pix_data = pix_data[sort_ind]
        
;    GET THE TIME STEP BETWEEN SUCCESSIVE DATA DUMPS
        delta_t = pix_data.ut - shift(pix_data.ut,1)
        delta_t[0] = 0.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... NOW BREAK THE DATA INTO INDIVIDUAL "BLOCKS"
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        obs_block = lonarr(pix_ct)

;    ... DEFINE BLOCKS BY TIME DISCONTINUITIES
        for j = 1L, pix_ct-1 do $
           obs_block[j] = obs_block[j-1] + (delta_t[j] gt large_t_step)

;    ... NOTE THE NUMBER OF BLOCKS
        n_blocks = max(obs_block)+1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... NOW LOOK AT EACH "TIME BLOCK" IN TURN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    COMBINE REFERENCE AND OBSERVING BLOCK FLAGS TO ID REFERENCE SPECTRA
        for j = 0L, n_blocks-1 do begin
           
;       IDENTIFY THE "ON" AND "OFF" SPECTRA FOR THIS OBSERVING BLOCK
           off_this_block = (obs_block eq j) and (pix_data.on eq 0)
           on_this_block  = (obs_block eq j) and (pix_data.on eq 1)

;       LABEL REFERENCE REGIONS
           reg = label_region(off_this_block)        

;       REJECT REFERENCES SMALLER THAN THE MINIMUM SIZE
           for k = 1L, max(reg) do begin
              this_reg = where(reg eq k)
              total_time_off = total(pix_data[this_reg].obstime)
;          ... REMOVE THESE FROM THE OFF-SPECTRUM LIST IF THE TIMING
;          IS TOO SHORT (BUT N.B. THAT THEY DO NOT BECOME "ON" SPECTRA)
              if total_time_off lt min_ref_size then begin
                 off_this_block[this_reg] = 0B
              endif
           endfor

;       RELABEL REFERENCE REGIONS IN CASE THEY HAVE CHANGED
           reg = label_region(off_this_block)

;       ***CAUTION*** WE ARE RELYING ON "LABEL_REGION" TO FOLLOW A
;       PARTICULAR PATTERN HERE (REGIONS INCREASING FROM LEFT TO RIGHT
;       ACROSS THE ARRAY).

           for k = 1L, max(reg)-1 do begin
;          INDICES OF THIS REFERENCE AND THE NEXT ONE
              start_ind = where(reg eq k, start_ct)  
              stop_ind = where(reg eq (k+1), stop_ct)
              
;          KEEP ONLY DATA UP TO THE MAXIMUM SIZE
;          (WE ASSUME ROUGHLY CONSTANT TIME PER SPECTRUM HERE, COULD
;          GET BURNED BY ODD DATA)
              if (total(pix_data[start_ind].obstime) gt max_ref_size) then begin
                 t_obs = median(pix_data[stop_ind].obstime)
                 n_obs = floor(max_ref_size / t_obs)
                 start_ind = start_ind[(start_ct - n_obs):*]
              endif

              if (total(pix_data[stop_ind].obstime) gt max_ref_size) then begin
                 t_obs = median(pix_data[stop_ind].obstime)
                 n_obs = floor(max_ref_size / t_obs)
                 stop_ind = stop_ind[0:n_obs-1]
              endif

;          CLEAN UP THE INFO ON THE OFF MEASUREMENTS
              start_ct = n_elements(start_ind)
              stop_ct = n_elements(stop_ind)
              
;          BUILD THE REFERENCES BEFORE AND AFTER THE SCAN
              if keyword_set(use_median) then begin
                 start_ref_spec = $
                    median(pix_data[start_ind].spectrum,dim=2)
                 stop_ref_spec = $
                    median(pix_data[stop_ind].spectrum,dim=2)
              endif else begin
                 start_ref_spec = $
                    total(pix_data[start_ind].spectrum,2,/nan) / (1.0*start_ct)
                 stop_ref_spec = $
                    total(pix_data[stop_ind].spectrum,2,/nan) / (1.0*stop_ct)
              endelse

;          NOTE THE MEAN CLOCK TIME (UT) AND INTEGRATED OBSERVING TIME OF EACH REFERENCE
              start_ref_time = $
                 1.0*mean(pix_data[start_ind].ut,/nan)           
              start_int_time = $
                 1.0*total(pix_data[start_ind].obstime,/nan)           
              stop_ref_time = $
                 1.0*mean(pix_data[stop_ind].ut,/nan)
              stop_int_time = $
                 1.0*total(pix_data[stop_ind].obstime,/nan)

;          FIND ALL DATA "ON" SOURCE IN BETWEEN THE TWO
              on_ind = where((pix_data.ut gt start_ref_time) and $
                             (pix_data.ut lt stop_ref_time) and $
                             on_this_block, on_ct)
              
;          AKL - NOT POSITIVE HOW TO TREAT THE "OFF" DATA HERE. SHOULD
;                WE REFERENCE SUBTRACT THEM?

              if on_ct gt 0 then begin
;             ... CASE 1: EQUAL WEIGHT TO EACH "REFERENCE" DATA DUMP
                 if keyword_set(equal_weight) then begin
                    ref_spectra = $
                       (start_ref_spec*sqrt(start_int_time) + $
                        stop_ref_spec*sqrt(stop_int_time)) / $
                       (sqrt(start_int_time)+sqrt(stop_int_time))
                 endif else begin
;             ... CASE 2: TIME-WEIGHTED REFERENCE

;                MAKE ARRAYS OUT OF THE UT OF EACH SPECT AND THE
;                REFERENCE SPECTRA.
                    tobs_array = (fltarr(nchan)+1.) # pix_data[on_ind].ut
                    start_spec_array = start_ref_spec # (fltarr(on_ct)+1.)
                    stop_spec_array = stop_ref_spec # (fltarr(on_ct)+1.)

                    delta_time = (stop_ref_time - start_ref_time)
                    ref_spectra = $
                       (tobs_array - start_ref_time) / delta_time * start_spec_array + $
                       (stop_ref_time - tobs_array) / delta_time * stop_spec_array
                 endelse
                 pix_data[on_ind].spectrum = (pix_data[on_ind].spectrum - ref_spectra)
                 pix_data[on_ind].ref_subtracted = 1
                 pix_data[on_ind].tref_before = start_int_time
                 pix_data[on_ind].tref_after = stop_int_time
                 pix_data[on_ind].time_to_ref_before = pix_data[on_ind].ut - start_ref_time
                 pix_data[on_ind].time_to_ref_after = stop_ref_time - pix_data[on_ind].ut
              endif
           endfor
        endfor  
        
        data[pix_ind] = pix_data
        
     endfor
     
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$     
;    WRITE OUT THE DATA
;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
;    &$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr

  endfor

  
end
