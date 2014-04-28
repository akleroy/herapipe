;+
; NAME:
;
; spectra_pipeline
;
; PURPOSE:
;
; Reduce HERA data taken in OTF mapping mode into spectra ready to be
; analyzed and gridded.
;
; CATEGORY:
;
; one of two data reduction master scripts
;
; CALLING SEQUENCE:
;
; spectra_pipeline
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
; 
; OUTPUTS:
;
; baseline subtracted, flagged, spectra
;
; OPTIONAL OUTPUTS:
;
; visualization (/show) and reports (/report)
;
; RESTRICTIONS: In addition to the IDL routines that make up the pipeline, it
; requires a specific directory structure:
;
; top_level/
; code/ cubes/ masks/ other_data/ reports/ spectra/
;
; The pipeline is assumed to be running out of the "code" directory.
;
; The pipeline also *REQUIRES* a file called "orig_data.txt" (the name can be
; changed via input), a space-delimited text file with these columns:
;
; "original_data_file" "working_file_name" "ref. sub. flag" "ref. sub. rules"
;
; Optional text input files allow the user to flag data as bad using some
; combination of day, pixel, and scan number ("bad_data.txt")
;
; PROCEDURE:
;
; 1) initialize an output structure to hold the processed spectra
; 2) reference-subtract the data
; 3) prune problematic fourier components
; 4) fit baselines to each spectrum
; 5) flag problematic spectra
;
; MODIFICATION HISTORY:
;
; spun out of single pipeline, streamlined - Oct 2009 leroy@mpia.de
;
; attempt update so that you never need to rerun the whole thing; only restart
; from the part of the pipeline you want to alter - Dec 2009 leroy@mpia.de
;
; at same time, update data flagging to produce a clear set of reports and
; store ON and OFF spectra as independent fields, increases size of data x3 -
; Dec 2009 leroy@mpia.de
;
;-

pro spectra_pipeline $
;  GENERAL INPUT/OUTPUT
   , identifier = tag $
   , orig_data_file = orig_data_file $
   , bad_data_file = bad_data_file $
;  TUNING PARAMETERS::: REFERENCE SUBTRACTION
   , equal_weight_ref = equal_weight_ref $
   , ref_mask_file = ref_mask_file $
   , min_ref = min_ref_time $
   , max_ref = max_ref_time $
   , relaxed_referencing = relaxed_referencing $
   , sliding_window = sliding_window $
;  TUNING PARAMETERS::: FOURIER PRUNING
   , bad_fft_chan = bad_fft_chan $
;  TUNING PARAMETERS::: BASELINE FITTING
   , window_root = window_root $
   , degree = degree $          
;  TUNING PARAMETERS::: DATA FLAGGING
   , allow_high_tsys = allow_high_tsys $
   , skip_uneven = skip_uneven $
   , smooth_for_flagging = smooth_for_flagging $
   , narrow_only = narrow_only $
   , relative_noise_only = relative_noise_only $
;  SUPPRESS PARTS OF THE PIPELINE
   , no_ref = no_ref $
   , no_fourier = no_fourier $
   , no_baseline = no_baseline $
   , no_rejection = no_rejection $
;  OUTPUT SWITCHES
   , show=show $
   , report=report $
;  ROUTE TO SPECIFIC PARTS OF THE PIPELINE
   , goto_init = goto_init $
   , goto_ref = goto_ref $
   , goto_sub = goto_sub $
   , goto_fft = goto_fft $
   , goto_fit = goto_fit $
   , goto_rej = goto_rej $
   , just = just $
;  ARE WE USING THE FTS?
   , fts = fts

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; DEFAULTS AND ERROR CHECKING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%

; ---------------------------------------------------------------------
; Input files
; ---------------------------------------------------------------------

; A FILE LISTING THE ORIGINAL .30M DATA
; FORMAT: "day" ".30m input file" "output file root (leave off '.fits')"
  if n_elements(orig_data_file) eq 0 then $
     orig_data_file = 'orig_data.txt'

; CHECK THAT THE FILE IS THERE...
  test = file_search(orig_data_file, count = ct)
  if ct ne 1 then begin
     message, 'Cannot find file '+orig_data_file+'.', /info
     message, 'A list of files to process is required.', /info 
     message, 'Specify via ORIG_DATA_FILE.'
  endif

; FILE LISTING DATA KNOWN TO BE BAD
; FORMAT: "day" "pixel" "scan #" (first line ignored)
  if n_elements(bad_data_file) eq 0 then $
     bad_data_file = 'bad_data.txt'

; AN IDENTIFIER FOR THIS PIPELINE RUN
  if n_elements(tag) eq 0 then tag = ''

; ---------------------------------------------------------------------
; Directory structure
; ---------------------------------------------------------------------

; OUTPUT DIRECTORY FOR REPORTS/IMAGES
  report_dir = '../reports/'       
; CHECK THAT THE DIRECTORY IS THERE...
  test = file_search(report_dir, count = ct)
  if ct ne 1 then begin
     message, 'The report directory '+report_dir+' is missing.', /info
     message, 'Please create it.'
  endif
  
; ---------------------------------------------------------------------
; Reference subtraction
; ---------------------------------------------------------------------
  
; DEFAULT TO SUBTRACTING REFERENCE MEASUREMENTS
  sub_reference = 1
  if keyword_set(no_ref) then begin
     sub_reference = 0
     message, 'Suppressing reference subtraction.', /info
  endif

  if n_elements(equal_weight_ref) eq 0 then begin
     equal_weight_ref = 0
  endif

; ---------------------------------------------------------------------
; Fourier pruning
; ---------------------------------------------------------------------

; DEFAULT TO PRUNING PROBLEMATIC FOURIER MODES
  zap_fourier = 1
  if keyword_set(no_fourier) then begin
     zap_fourier = 0
     message, 'Suppressing Fourier pruning.', /info
  endif

; THE LOCATION (IN THE FFT) OF BAD FOURIER CHANNELS/STANDING WAVES
  if n_elements(bad_fft_chan) eq 0 then $
     bad_fft_chan = [60, 135]

; ---------------------------------------------------------------------
; Baseline subtraction
; ---------------------------------------------------------------------

; DEFAULT TO SUBTRACTING BASELINE FROM THE DATA
  sub_baseline = 1
  if keyword_set(no_baseline) then begin
     sub_baseline = 0
     message, 'Suppressing baseline subtraction.', /info
  endif
  
; DEFAULT TO A LINEAR BASELINE FIT
  if n_elements(degree) eq 0 then degree = 1

; ---------------------------------------------------------------------
; Data flagging
; ---------------------------------------------------------------------

; DEFAULT TO FLAGGING DATA THAT SHOW PATHOLOGIES
  flag_data = 1
  if keyword_set(no_rejection) then begin
     flag_data = 0
     message, 'Suppressing data flagging.', /info
  endif

; SMOOTHING FOR DATA FLAGGING (OPTIONAL)
  if n_elements(smooth_for_flagging) eq 0 then begin
     smooth_for_flagging = 1
  endif

; ALLOW ONLY THE NARROW PART OF THE UNEVENNESS FLAGGING
  if n_elements(narrow_only) eq 0 then begin
     narrow_only = 0
  endif

; ... END DEFAULTS
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; SKIP TO THE USER-SPECIFIED PART OF THE PIPELINE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Using the reviled goto, jump to the requested part of the script. If the
; JUST keyword was also set by the user, we will do one part of the reduction
; and then exit. Otherwise, we step unit the end of one of the iterations.
;
  if keyword_set(goto_init) then goto, init
  if keyword_set(goto_ref) then goto, ref
  if keyword_set(goto_sub) then goto, sub
  if keyword_set(goto_fft) then goto, fft
  if keyword_set(goto_fit) then goto, fit
  if keyword_set(goto_rej) then goto, rej

; ... END GOTO RE-ROUTING STATEMENT
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; INITIALIZE OUTPUT STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; The data come in as a .fits table. Take this structure and add several new
; fields to hold the results of the reference subtraction, data flagging, and
; baseline fit. Then call several subroutines to fill out basic aspects of the
; data. The main slowness here is the constant reading and writing of
; files. This could be streamlined significantly with a bit of effort.

  init:

  message, 'Initializing output file.', /info
  message, 'Time stamp: ', /info
  spawn, 'date'

; ... MAKE SURE THE STRUCTURE HAS THE NEEDED FIELDS
  init_struct $
     , orig_data_file $
     , tag = tag

; ... FILL IN POSITION AND VELOCITY COORDINATES
  fill_in_coords $
     , orig_data_file $
     , tag = tag $
     , force_v0 = force_v0

; ... APPLY ANY USER-SUPPLIED FLAGGING
  apply_user_flags $
     , orig_data_file $
     , tag = tag $
     , bad_data_file = bad_data_file

; ... FLAG PATHOLOGICAL SPECTRA (BLANK, HIGH TSYS)
  flag_pathological_spectra $
     , orig_data_file $
     , tag = tag $
     , allow_high_tsys = allow_high_tsys

; ... ESTIMATE THE POSITION ANGLE OF EACH SCAN
  measure_scan_angle $
     , orig_data_file $
     , tag = tag $
     , fts = fts

; ... MEASURE THE TIMING OF INDIVIDUAL SUBSCANS
  measure_subscan_timing $
     , orig_data_file $
     , tag = tag 

  if keyword_set(just) then return

; ... END INITIALIZATION AND POPULATION OF OUTPUT STRUCTURE

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; CONSTRUCT REFERENCES FOR INDIVIDUAL SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Identify a subset of spectra as OFF source measurements and use
; these to build reference (blank sky) spectra
;

  ref:
  message, 'Building reference spectra.',/info
  message, 'Time stamp: ', /info
  spawn, 'date'

  if sub_reference then begin
;    ... RESET THE REFERENCES
     reset_references $
        , orig_data_file $
        , tag = tag

;    ... FIRST CHECK WHETHER A SUBTRACTION IS NEEDED     
     all_already_subtracted = $
        is_ref_subtracted(orig_data_file, tag = tag)

;    ... TURN OFF THE SUBTRACTION IF ALL DATA ALREADY PROCESSED
     if all_already_subtracted then sub_reference = 0
  endif

  if sub_reference then begin
     
;    ... FIRST, APPLY A REFERENCE MASK IF WE HAVE ONE
     apply_ref_mask $
        , orig_data_file $
        , tag = tag $
        , ref_mask = ref_mask_file

;    ... THEN APPLY ANY USER-SPECIFIED RULES FOR REFERENCES
     apply_user_ref $
        , orig_data_file $
        , tag = tag

;    ... ILLUSTRATE THE REFERENCES
     if keyword_set(show) then begin
        show_ref $
           , orig_data_file $
           , tag = tag $
           , report = report $
           , /pause
     endif

;    ... NOW ASSEMBLE THE OFF DATA INTO SPECTRA
     make_reference_spectra $
        , orig_data_file $
        , tag = tag $
        , equal_weight = equal_weight_ref $
        , relaxed = relaxed_referencing $
        , sliding_window = sliding_window $
        , fts=fts
  endif

  if keyword_set(just) then return

; ... END REFERENCE CONSTRUCTION

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; SUBTRACT REFERENCES FROM INDIVIDUAL SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Subtract the reference spectra that we just created. This is a simple step
; that effectively resets the second half of the spectra pipeline. If the
; fitting windows, rejection scheme, or FFT pruning change the pipeline can be
; reset at this step, avoiding the need to reconstruct references.
;

  sub:
  message, 'Subtracting reference spectra.',/info
  message, 'Time stamp: ', /info
  spawn, 'date'

; SUBTRACT THE REFERENCE SPECTRUM FROM THE 
  subtract_reference $
     , orig_data_file $
     , tag = tag

  if keyword_set(just) then return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; FFT PRUNING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; The 30m is known to have several standing waves induced by the
; optical path. This step removes these by hand, replacing the
; affected fourier components with new components that have a random
; phase and an amplitude interpolated from neighboring channels.
;
  
  fft:
  message, 'Pruning bad Fourier components.', /info
  message, 'Time stamp: ', /info
  spawn, 'date'

  if zap_fourier then begin

;    FOURIER EDITING
     fourier_prune $
        , orig_data_file $
        , tag = tag $
        , bad_channels =  bad_fft_chan $
        , show = show $
        , report = report $
        , channel_prof = kern_in $
        , /monte

  endif
        
  if keyword_set(just) then return

; ... END FOURIER PRUNING

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; BASELINE SUBTRACTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Fit a polynomial baseline to each spectrum after optionally
; smoothing the data in time and velocity.
;

  fit:
  message, 'Fitting baselines.',/info
  message, 'Time stamp: ', /info
  spawn, 'date'

; FIRST READ IN THE SPECTRAL BASELINE WINDOWS
  reset_windows $
     , orig_data_file $
     , tag = tag

  for i = 0, n_elements(window_root)-1 do begin
     read_in_windows $
        , orig_data_file $
        , tag = tag $
        , window_root = window_root[i]
  endfor

; NOW FIT AND THEN SUBTRACT A BASELINE FROM EACH SPECTRUM
  if sub_baseline then begin

     hera_base_fit $
        , orig_data_file $
        , tag = tag $
        , degree = degree $
        , show = show $
        , smooth_in_time = 0 $
        , /smooth_in_vel $
        , fts=fts
    
     sub_base_fit $
        , orig_data_file $
        , tag = tag $
        , show = show
     
  endif

  if keyword_set(just) then return

; ... END BASELINE FITTING

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; FLAG PROBLEMATIC SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; Despite the steps taken so far, pathological spectra (i.e., those
; with poorly behaved baselines due to instabilities in the receivers
; or atmosphere) can persist to this step and cause problems moving
; forward. Here we use several ad hoc tests to identify and flag these
; data.
;

  rej:
  message, 'Flagging problematic data.', /info
  message, 'Time stamp: ', /info
  spawn, 'date'

  if (flag_data eq 1B) then begin

;    RESET THE FLAGGING
     reset_flags $
        , orig_data_file $
        , tag = tag

; ... APPLY ANY USER-SUPPLIED FLAGGING
     apply_user_flags $
        , orig_data_file $
        , tag = tag $
        , bad_data_file = bad_data_file

;    FLAG PATHOLOGICAL SPECTRA (BLANK, HIGH TSYS)
     flag_pathological_spectra $
        , orig_data_file $
        , tag = tag $
        , allow_high_tsys = allow_high_tsys

;    FLAG UNEVEN SPECTRA
     if keyword_set(skip_uneven) eq 0 then begin
        for i = 0, n_elements(smooth_for_flagging)-1 do begin
           flag_uneven_spectra $
              , orig_data_file $
              , tag = tag $
              , show = show $
              , report = report $
              , smooth = smooth_for_flagging[i] $
              , wide = (narrow_only eq 0)
        endfor
     endif
    
;    FLAG RIPPLY SPECTRA (BASED ON VELOCITY-CONTIGUOUS SINGLE-SIGN REGIONS)
     for i = 0, n_elements(smooth_for_flagging)-1 do begin
        flag_ripply_spectra $
           , orig_data_file $
           , tag = tag $
           , show = show $
           , report = report $
           , smooth = smooth_for_flagging[i] $
           , /blank $
           , fts = fts
     endfor

;    FLAG NOISY SPECTRA (BASED ON RMS ABOUT THE BASELINE FIT)
     for i = 0, n_elements(smooth_for_flagging)-1 do begin
        flag_noisy_spectra $
           , orig_data_file $
           , tag = tag $
           , show = show $
           , report = report $
           , smooth = smooth_for_flagging[i] $
           , relative_noise_only = relative_noise_only $
           , /blank $
           , fts = fts
     endfor

;    MAKE A QUICK REPORT ON THE DATA FLAGGING
     flagging_report $
        , orig_data_file $
        , tag = tag $
        , fts = fts

;    MEASURE THE NOISE AND ASSESS ITS AVERAGING PROPERTIES
     noise_report $
        , orig_data_file $
        , tag = tag $
        , /blank $
        , fts=fts

  endif

  if keyword_set(just) then return 

; ... END DATA FLAGGING

  message, 'Time stamp: ', /info
  spawn, 'date'

end                             ; of HERACLES "SPECTRA" pipeline

; ... now you should move to the "CUBE" portion of the pipeline
