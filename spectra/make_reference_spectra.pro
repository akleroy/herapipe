pro make_reference_spectra $
   , list_file $
   , tag = tag $
   , equal_weight = equal_weight $
   , sliding_window = sliding_window $
   , median = use_median $
   , min_ref = min_ref_size $
   , max_ref = max_ref_size $
   , relaxed = relaxed

;+
;
; Reduction subroutine:
;
; From a structure where spectra have already been labeled as either "ON" or
; "OFF," construct a reference spectrum from average of nearby OFF spectra for
; each ON position. Save that spectrum and the original (raw) spectrum inside
; the data structure.
;
;-

  
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
  
; ALSO DEFINE THE LONGEST TIME THAT WE WILL ALLOW A SPECTRUM TO BE SEPARATED
; FROM ANY OFF.
  large_t_separation = 90.

; SET THE TUNING PARAMETERS TO A MORE RELAXED SETTING, USEFUL FOR A COUPLE
; EARLY CASES WHERE THE OBSERVING STRATEGY WASN'T SO GREAT. PLEASE NOTE THAT
; THIS IS *REALLY* NOT THE PREFERRED APPROACH. SOMETIMES,
; UNFORTUNATELY, YOU SUBTRACT THE REFERENCE YOU HAVE
  if keyword_set(relaxed) then begin
     min_ref_size = 2.0
     large_t_step = 90.
     large_t_separation = 120.
  endif

; ANOTHER WAY TO DO THINGS (FOR HARD CASES)
  if keyword_set(sliding_window) then begin
     window_size = 100.
     min_ref_size = 5.0
  endif

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
        message, 'File not found '+ $
                 string(working_name[m])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
     
;    CHECK IF WE ARE ALREADY REFERENCE SUBTRACTED
     need_sub = total(data.on eq 1 and $
                      data.ref_subtracted eq 0 and $
                      data.flagged eq 0)
     if need_sub eq 0 then continue

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
; IF WE'RE USING A SLIDING WINDOW, DO THAT NOW
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

        if keyword_set(sliding_window) then begin
           on_ind = where(pix_data.on, on_ct)
           off_data = pix_data[where(pix_data.on eq 0)]
           
           for j = 0, on_ct-1 do begin
              off_for_this_spec = $
                 (abs(off_data.ut - pix_data[on_ind[j]].ut) lt $
                  window_size)
              
              off_ind = $
                 where(off_for_this_spec, off_ct)
              
              if off_ct eq 0 then $
                 continue
              
              if off_ct eq 1 then $
                 ref_spectrum = off_data[off_ind].raw_spectrum $
              else $
                 ref_spectrum = $
                 total(off_data[off_ind].raw_spectrum,2,/nan)/ (1.0*off_ct)
              
;             STORE THE OFF SPECTRUM
              pix_data[on_ind[j]].ref_spectrum = ref_spectrum

;             SAVE THE REFERENCE-SUBTRACTED SPECTRUM AS THE CURRENT ONE
              pix_data[on_ind[j]].spectrum = $
                 (pix_data[on_ind[j]].raw_spectrum - ref_spectrum)

;             NOTE THAT THESE DATA HAVE BEEN REFERENCE SUBTRACTED
              pix_data[on_ind].ref_subtracted = 1

;             STORE SOME DETAILS ABOUT THE REFERENCE SPECTRUM
              pix_data[on_ind[j]].tref_before = $
                 total(pix_data[off_ind].obstime)
              pix_data[on_ind[j]].tref_after = !values.f_nan
              pix_data[on_ind[j]].time_to_ref_before = $
                 window_size
              pix_data[on_ind[j]].time_to_ref_after = $
                 window_size
           endfor

           data[pix_ind] = pix_data    

           continue

        endif

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

;       DIAGNOSTIC PLOT COMMENTED OUT
;           bind = where(on_this_block or off_this_block)
;           plot, on_this_block[bind], yrange=[-0.1, 1.1]
;           oplot, off_this_block[bind], color=getcolor('red')
;           ch = get_kbrd(1)

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
                 print, 'Zapping short OFF block.'
              endif
           endfor

;       RELABEL REFERENCE REGIONS IN CASE THEY HAVE CHANGED
           reg = label_region(off_this_block)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... DEAL WITH TWO SPECIAL CASES: "ON" BEFORE AND AFTER ALL OFFS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;          THE FIRST OFF-SOURCE REGION IN THIS OBSERVING BLOCK
           first_off = where(reg eq 1, first_off_ct)
           if (total(pix_data[first_off].obstime) gt max_ref_size) then begin
              t_obs = median(pix_data[first_off].obstime)
              n_obs = floor(max_ref_size / t_obs)
              first_off = first_off[0:n_obs-1]
           endif
           first_off_ct = n_elements(first_off)
           
           first_ref_time = $
              1.0*mean(pix_data[first_off].ut,/nan)           
           first_int_time = $
              1.0*total(pix_data[first_off].obstime,/nan)           
           first_ref_spec = $
              total(pix_data[first_off].raw_spectrum,2,/nan) / $
              (1.0*first_off_ct)

           on_before_first_off = $
              where(on_this_block and (pix_data.ut lt first_ref_time) $
                    , on_before_first_off_ct)
           
           if (on_before_first_off_ct gt 0) then begin
;             STORE THE OFF SPECTRUM (INDEXING DOES WORK THIS WAY)
              pix_data[on_before_first_off].ref_spectrum = first_ref_spec

;             SAVE THE REFERENCE-SUBTRACTED SPECTRUM AS THE CURRENT ONE
              pix_data[on_before_first_off].spectrum = $
                 (pix_data[on_before_first_off].raw_spectrum - $
                  pix_data[on_before_first_off].ref_spectrum)
              
;             NOTE THAT THESE DATA HAVE BEEN REFERENCE SUBTRACTED
              pix_data[on_before_first_off].ref_subtracted = 1

;             STORE SOME DETAILS ABOUT THE REFERENCE SPECTRUM
              pix_data[on_before_first_off].tref_before = 0.0
              pix_data[on_before_first_off].tref_after = first_int_time
              pix_data[on_before_first_off].time_to_ref_before = !values.f_nan
              pix_data[on_before_first_off].time_to_ref_after = $
                 first_ref_time - pix_data[on_before_first_off].ut
           endif

;          THE LAST OFF-SOURCE REGION IN THIS OBSERVING BLOCK
           last_off = where(reg eq max(reg), last_off_ct)
           if (total(pix_data[last_off].obstime) gt max_ref_size) then begin
              t_obs = median(pix_data[last_off].obstime)
              n_obs = floor(max_ref_size / t_obs)
              last_off = last_off[(last_off_ct - n_obs):*]
           endif
           last_off_ct = n_elements(last_off)

           last_ref_time = $
              1.0*mean(pix_data[last_off].ut,/nan)
           last_int_time = $
              1.0*total(pix_data[last_off].obstime,/nan)
           if last_off_ct eq 1 then $
              last_ref_spec = pix_data[last_off].raw_spectrum $
           else $
              last_ref_spec = $
              total(pix_data[last_off].raw_spectrum,2,/nan) / (1.0*last_off_ct)
           
           on_after_last_off = $
              where(on_this_block and pix_data.ut gt last_ref_time $
                    , on_after_last_off_ct)
           
           if (on_after_last_off_ct gt 0) then begin
;             STORE THE OFF SPECTRUM
              pix_data[on_after_last_off].ref_spectrum = last_ref_spec

;             SAVE THE REFERENCE-SUBTRACTED SPECTRUM AS THE CURRENT ONE
              pix_data[on_after_last_off].spectrum = $
                 (pix_data[on_after_last_off].raw_spectrum - $
                  pix_data[on_after_last_off].ref_spectrum)

;             NOTE THAT THESE DATA HAVE BEEN REFERENCE SUBTRACTED
              pix_data[on_after_last_off].ref_subtracted = 1

;             STORE SOME DETAILS ABOUT THE REFERENCE SPECTRUM
              pix_data[on_after_last_off].tref_before = last_int_time
              pix_data[on_after_last_off].tref_after = 0.0
              pix_data[on_after_last_off].time_to_ref_before = $
                 pix_data[on_after_last_off].ut - last_ref_time
              pix_data[on_after_last_off].time_to_ref_after = !values.f_nan
           endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... NOW LOOP OVER DISTINCT PAIRS OF OFF REGIONS INSIDE THE BLOCK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;       ***CAUTION*** WE ARE RELYING ON "LABEL_REGION" TO FOLLOW A PARTICULAR
;       PATTERN HERE (REGIONS INCREASING FROM LEFT TO RIGHT ACROSS THE
;       ARRAY). IT DOES AS OF THE WRITING OF THIS PROGRAM BUT THE BEHAVIOR
;       COULD CONCEIVABLY CHANGE.

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
                    median(pix_data[start_ind].raw_spectrum,dim=2)
                 stop_ref_spec = $
                    median(pix_data[stop_ind].raw_spectrum,dim=2)
              endif else begin
                 start_ref_spec = $
                    total(pix_data[start_ind].raw_spectrum,2,/nan) / (1.0*start_ct)
                 stop_ref_spec = $
                    total(pix_data[stop_ind].raw_spectrum,2,/nan) / (1.0*stop_ct)
              endelse

;          NOTE THE MEAN CLOCK TIME (UT) AND INTEGRATED OBSERVING TIME OF EACH
;          REFERENCE
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
;                WE REFERENCE-SUBTRACT THEM? FOR NOW, NO

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE REFERENCE SPECTRUM FOR EACH ON SPECTRUM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

              if on_ct gt 0 then begin

;             ... CASE 1: EQUAL WEIGHT TO EACH "REFERENCE" DATA DUMP
                 if keyword_set(equal_weight) then begin

;                   ... WEIGHT BY INTEGRATION TIME
                    weight_to_start = (start_int_time*1.0) 
                    weight_to_stop = (stop_int_time*1.0)

                 endif else begin
;             ... CASE 2: TIME-WEIGHTED REFERENCE

;                MAKE ARRAYS OUT OF THE UT OF EACH SPECT AND THE
;                REFERENCE SPECTRA.
                    tobs_array = (fltarr(nchan)+1.) # pix_data[on_ind].ut

                    delta_time = (stop_ref_time - start_ref_time)
                    weight_to_start = $
                       (stop_ref_time - tobs_array) / delta_time
                    weight_to_stop = $
                       (tobs_array - start_ref_time) / delta_time
                 endelse

;               COMBINE THE STARTING AND STOPPING SPECTRUM
                 start_spec_array = start_ref_spec # (fltarr(on_ct)+1.)
                 stop_spec_array = stop_ref_spec # (fltarr(on_ct)+1.)

                 total_weight = weight_to_start + weight_to_stop
                 weight_to_start = weight_to_start / total_weight
                 weight_to_stop = weight_to_stop / total_weight
                 
                 if total(abs((weight_to_start + weight_to_stop) - 1.0) gt 1d-6) gt 0 then $
                    message, 'Weighting problem!'
                 
                 ref_spectra = $
                    weight_to_start * start_spec_array + $
                    weight_to_stop * stop_spec_array

;                STORE THE OFF SPECTRUM
                 pix_data[on_ind].ref_spectrum = ref_spectra

;                SAVE THE REFERENCE-SUBTRACTED SPECTRUM AS THE CURRENT ONE
                 pix_data[on_ind].spectrum = $
                    (pix_data[on_ind].raw_spectrum - ref_spectra)

;                NOTE THAT THESE DATA HAVE BEEN REFERENCE SUBTRACTED
                 pix_data[on_ind].ref_subtracted = 1

;                STORE SOME DETAILS ABOUT THE REFERENCE SPECTRUM
                 pix_data[on_ind].tref_before = start_int_time
                 pix_data[on_ind].tref_after = stop_int_time
                 pix_data[on_ind].time_to_ref_before = $
                    pix_data[on_ind].ut - start_ref_time
                 pix_data[on_ind].time_to_ref_after = $
                    stop_ref_time - pix_data[on_ind].ut

              endif
           endfor
        endfor  
        
        data[pix_ind] = pix_data
        
     endfor
   
     missed = where(data.on and data.ref_subtracted eq 0, missed_ct)
     total = where(data.on, total_ct)
     message $
        , 'Missed '+str(missed_ct)+' of '+ $
        str(total_ct)+' spectra during reference subtraction.', /info

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% 
; WRITE OUT THE DATA
; (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% 

     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr

  endfor

  
end
