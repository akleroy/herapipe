;+
; NAME:
;
; heracles_pipeline
;
; PURPOSE:
;
; Reduce HERA spectra taken in OTF mapping mode into science-ready data cubes.
;
; CATEGORY:
;
; data reduction master script
;
; CALLING SEQUENCE:
;
; heracles_pipeline (then whatever options you want)
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
; baseline subtracted spectra, cubes
;
; OPTIONAL OUTPUTS:
;
; reports
;
; RESTRICTIONS: In addition to the IDL routines that make up the pipeline, it
; requires a specific directory structure:
;
; top_level/
; code/ reports/ output/ windows/ DDmmmYY/ DDmmmYY/ DDmmmYY/
;
; where DDmmmYY/ is a "day" of observing (e.g., 31oct07) and contains the
; subdirectories:
;
; calib/ base_sub/ cubes/
;
; The pipeline is assumed to be running out of the "code" directory.
;
; The pipeline also *REQUIRES* a file called "orig_data.txt" (the name can be
; changed via input), a space-delimited text file with these columns:
;
; "day" "original_30m_file" "original_fits_file" "reference subtracted flag"
;
; It will convert the .30m files to the specified .fits files and then use this
; file as list of the data sets being reduced.
;
; Optional text input files allow the user to flag data as bad using some
; combination of day, pixel, and scan number ("bad_data.txt")
;
; PROCEDURE: The broad procedure in its most complicated form is:
;
; 1) dump the data from CLASS to FITS
; 2) reference subtract the data
; 3) fit baselines to each spectrum
; 4) sort the original spectra by the direction they scan across the galaxy
; 5) use the original reduction and these spectra to calibrate pixel gains
; 6) grid data sharing a common scan direction into cubes
; 8) combine different scan directions in the FFT domain
; 9) identify regions of bright emission inside the cube
; 
; Not all of these steps can be done during the first reduction. The approach
; is designed to be somewhat iterative. Each step is documented more
; thoroughly inside this script and in the routines related to that step.
;
; EXAMPLE:
;
; 
;
; MODIFICATION HISTORY:
;
;-

pro heracles_pipeline $
;  BASIC INPUT/OUTPUT FROM THIS PIPELINE RUN
   , gname = gname $            ; GALAXY NAME
   , tag = tag $                ; IDENTIFIER FOR THIS REDUCTION RUN
   , orig_data_file = orig_data_file $ ; TEXT FILE POINTING TO ORIGINAL DATA
   , bad_data_file = bad_data_file $   ; TEXT FILE IDENTIFYING BAD DATA
   , window_root = window_root $       ; LOCATION OF SPECTRAL BASELINE WINDOWS
   , skip_day = skip_day $             ; LIST OF 'DAYS' TO SKIP
   , map_start_file = map_start_file $ ; USER-SUPPLIED START OF MAPS
;  INPUT FROM A PREVIOUS PIPELINE RUN
   , ref_mask_file = ref_mask_file $  ; MASK SHOWING WHAT IS 'ON' SOURCE
   , prev_bright_map = prev_bright_map_in $ ; PREVIOUS MAP - USED FOR GAINS
   , prev_cube = prev_cube_in $ ; PREVIOUS CUBE - USED FOR GAINS
   , prev_mask = prev_mask_in $ ; PREVIOUS MASK - USED FOR GAINS
;  TUNING PARAMETERS
   , degree = degree $          ; DEGREE OF POLYNOMIAL FIT
   , bad_fft_chan = bad_fft_chan $ ; FFT CHANNELS TO BE INTERPOLATED ACROSS
   , resid_reject_sigma = resid_reject_sigma $ ; TSYS REJECTION THRESHOLD
   , coverage_thresh = coverage_thresh $       ; COVERAGE/MAX TO KEEP
   , force_v0 = force_v0 $                     ; HACK FOR M51 ONLY?
   , min_ref = min_ref_time $                  ; MINIMUM TIME FOR A REFERENCE
   , max_ref = max_ref_time $                  ; MAXIMUM TIME FOR A REFERENCE
;  OUTPUT SWITCHES
   , show=show $
;  ROUTE TO SPECIFIC PARTS OF THE PIPELINE
   , goto_ref = goto_ref $
   , goto_fit = goto_fit $
   , goto_cal = goto_cal $
   , goto_split = goto_split $
   , goto_grid = goto_grid $
   , goto_winmask = goto_winmask $
   , goto_coadd = goto_coadd $
   , goto_mask = goto_mask $
   , goto_plait = goto_plait $
   , goto_report = goto_report $
   , just = just $              ; ONLY DO ONE PART OF THE REDUCTION
;  SWITCHES INDICATING WHAT TO DO
   , no_ref = no_ref $
   , no_rejection = no_rejection $
   , decimate = decimate $
   , zap_fourier = zap_fourier $
   , one_big_grid = one_big_grid $
   , split_by_pol = split_by_pol $
   , split_by_ang = split_by_ang $
   , gauss_kern = gauss_kern $
   , gauss_fwhm = gauss_fwhm $
   , grid_by_median = grid_by_median $
   , use_mean_ctr = use_mean_ctr $
   , coverage_coadd = coverage_coadd $
   , apodize = apodize $
   , cal_pixel_gain = cal_pixel_gain $
   , apply_gain = apply_gain $
   , plait = plait


; EXIT THE PROGRAM IF THERE IS AN ERROR
  ON_ERROR, 2

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; DEFAULTS AND ERROR CHECKING ON USER INPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%

; ---------------------------------------------------------------------
; Basic input/output
; ---------------------------------------------------------------------

; GALAXY NAME
  if n_elements(gname) eq 0 then begin
     message, 'Field GNAME (unique galaxy identifier ) required.'
  endif

; DEFAULT REDUCTION IDENTIFIER
  default_tag = ''

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

; A LIST OF DAYS TO SKIP IN THE CURRENT REDUCTION
; FORMAT: string array of form DDmmmYY, e.g., 01apr08
  if n_elements(skip_day) eq 0 then $
     skip_day = ''                     

; A LIST OF SCANS KNOWN TO START MAPS
; format is "day" "scan number"
  if n_elements(map_start_file) eq 0 then $
     map_start_file = 'map_start_list.txt'

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

; OUTPUT DIRECTORY FOR MASKS, CUBES, ETC.
  output_dir = '../output/'
; CHECK THAT THE DIRECTORY IS THERE...
  test = file_search(output_dir, count = ct)
  if ct ne 1 then begin
     message, 'The output directory '+output_dir+' is missing.', /info 
     message, 'Please create it.'
     return
  endif

; WINDOW WORKING DIRECTORY
  win_dir = '../windows/'       
; CHECK THAT THE DIRECTORY IS THERE...
  test = file_search(win_dir, count = ct)
  if ct ne 1 then begin
     message, 'The spectral windows directory '+win_dir+' is missing.', /info
     message, 'Please create it.'
     return
  endif

; ---------------------------------------------------------------------
; Input based on a previous reduction
; ---------------------------------------------------------------------

; NAME OF A 2D MASK SHOWING BRIGHT SIGNAL IDENTIFIED FROM THE DATA CUBE
  if n_elements(prev_bright_map_in) eq 0 then $
     prev_bright_map = 'prelim/'+gname+'_bright_map_prelim.fits' $
  else $
     prev_bright_map = prev_bright_map_in

; NAME OF THE DATA CUBE PRODUCED DURING THE FIRST ROUND OF REDUCTION
  if n_elements(prev_cube_in) eq 0 then $
     prev_cube = 'prelim/'+gname+'_prelim.fits' $
  else $
     prev_cube = prev_cube_in

; ... ASSUME THE BRIGHT MAP IS INSIDE THE OUTPUT DIRECTORY
  prev_bright_map = output_dir+prev_bright_map

; ... CHECK THAT IT IS THERE
  test = file_search(prev_bright_map, count = ct)
  if ct eq 0 then begin
     message, 'No "bright emission map" from previous reduction found.', /info
     message, 'Cannot calibrate pixel gains. Turning that part of the pipeline OFF.', /info
     cal_pixel_gain = 0
  endif

  prev_cube = output_dir + prev_cube
  test = file_search(prev_cube, count = ct)
  if ct eq 0 then begin
     message, 'No "reference cube" from previous reduction found.', /info
     message, 'Cannot calibrate pixel gains. Turning that part of the pipeline OFF.', /info
     cal_pixel_gain = 0
  endif
  
; A MASK MATCHED TO THIS EARLIER DATA CUBE, USED TO IDENTIFY WHICH PART OF THE
; SPECTRUM HAS BEEN FIT
  if n_elements(prev_mask_in) eq 0 then $
     prev_mask = 'prelim/'+gname+'_prelim_line_mask.fits' $
  else prev_mask = prev_mask_in
  prev_mask = output_dir + prev_mask
  test = file_search(prev_mask, count = ct)
  if ct eq 0 then begin
     message, 'No mask from previous reduction found.', /info
     message, 'Cannot calibrate pixel gains. Turning that part of the pipeline OFF.', /info
     cal_pixel_gain = 0
  endif

; ---------------------------------------------------------------------
; What to do (these are flags, where 1 == do this, 0 == don't do this)
; ---------------------------------------------------------------------
  
; FLAGS APPROPRIATE FOR A FIRST PASS WITH NOTHING KNOWN
  if keyword_set(preliminary) then begin
     message, 'Attempting a preliminary reduction.', /info
     if n_elements(degree) eq 0 then degree = 1
     default_decimate       = 0B
     default_zap_fourier    = 0B
     default_one_big_grid   = 1B
     default_split_by_cover = 0B
     default_split_by_pol   = 0B
     default_split_by_ang   = 1B
     default_gauss_kern     = 0B
     default_coverage_coadd = 1B
     default_apodize        = 0B
     default_cal_pixel_gain = 0B
     default_apply_gain     = 0B
     default_plait          = 1B
     default_tag            = 'prelim'
  endif
  
; SET THE FLAGS, ALLOWING USER INPUT TOP PRIORITY
  if n_elements(no_ref) eq 0 then $
     no_ref = default_no_ref
  if n_elements(no_rejection) eq 0 then $
     no_rejection = default_no_rejection
  if n_elements(decimate) eq 0 then $
     decimate = default_decimate
  if n_elements(zap_fourier) eq 0 then $
     zap_fourier = default_zap_fourier
  if n_elements(one_big_grid) eq 0 then $$
     one_big_grid = default_one_big_grid
  if n_elements(split_by_cover) eq 0 then $
     split_by_cover = default_split_by_cover
  if n_elements(split_by_pol) eq 0 then $
     split_by_pol = default_split_by_pol
  if n_elements(split_by_ang) eq 0 then $
     split_by_ang = default_split_by_ang
  if n_elements(gauss_kern) eq 0 then $
     gauss_kern = default_gauss_kern
  if n_elements(coverage_coadd) eq 0 then $
     coverage_coadd = default_coverage_coadd
  if n_elements(apodize) eq 0 then $
     apodize = default_apodize
  if n_elements(cal_pixel_gain) eq 0 then $
     cal_pixel_gain = default_cal_pixel_gain
  if n_elements(apply_gain) eq 0 then $
     apply_gain = default_apply_gain
  if n_elements(plait) eq 0 then begin
     plait = default_plait
     if split_by_ang eq 0 then begin
        message, 'Turning on splitting by angle (needed to PLAIT).', /info
     endif
  endif

; ---------------------------------------------------------------------
; Figure out the identifier for this reduction
; ---------------------------------------------------------------------

; SET THE IDENTIFIER FOR THE REDUCTION
  if n_elements(tag) eq 0 then $
     tag = default_tag
; SET THE OUTPUT SUBDIRECTORY TO USE THE TAG NAME
  output_dir += tag+'/'
  test = file_search(output_dir, count = ct)
  if ct ne 1 then begin
     message, 'The output subdirectory '+output_dir+' does not exist. Creating it.', /info 
     spawn, 'mkdir '+output_dir
  endif

; ---------------------------------------------------------------------
; Tuning parameters
; ---------------------------------------------------------------------

; DEGREE OF THE POLYNOMIAL USED TO FIT THE BASELINE
  if n_elements(degree) eq 0 then $
     degree = 1

; MAGNITUDE OF HIGH DEVIANTS ABOUT EXPECTED NOISE TO FLAG AS BAD
  if n_elements(resid_reject_sigma) eq 0 then $
     resid_reject_sigma = 2.0

; THE LOCATION (IN THE FFT) OF BAD FOURIER CHANNELS/STANDING WAVES
  if n_elements(bad_fft_chan) eq 0 then $
     bad_fft_chan = [135] ; [60, 135] ?
  
; FRACTION OF MEDIAN COVERAGE NEEDED TO BE INCLUDED IN THE CUBE WHEN COADDING
  if n_elements(coverage_thresh) eq 0 then $
     coverage_thresh = 0.5

; ---------------------------------------------------------------------
; Extensions and roots for files output
; ---------------------------------------------------------------------

; A MASK SHOWING THE BASELINE FITTING WINDOWS AND LINE ON THE FINAL GRID
  window_mask_ext  = '_window_mask'

; A MASK SHOWING ONLY THE LINE (NOT THE FITTING REGIONS)
  window_line_ext  = '_line_mask'

; A MASK SHOWING SEVERAL CHANNELS ON EACH SIDE OF THE LINE WINDOW
  window_line_border_ext  = '_line_border_mask'

; A MASK SHOWING ONLY THE FITTING REGION (EXPECTED TO BE SIGNAL FREE)
  window_noise_ext = '_noise_mask'

; A MASK SHOWING SEVERAL CHANNELS ON THE EDGE OF THE FITTING WINDOW
  window_noise_border_ext = '_noise_border_mask'

; EXTENSIONS FOR FILES HOLDING REFERENCE-SUBTRACTED SPECTRA
  ref_sub_spectra_ext = '_refsub'

; EXTENSIONS FOR FILES HOLDING BASELINE-SUBTRACTED SPECTRA
  base_sub_spectra_ext = '_base_'+tag

; THE STEM OF CUBES CONTAINING ONLY ONE COVERAGE
  coverage_cube_root = 'map_'+tag

; STEM OF CUBES CONTAINING DATA AVERAGED OVER ALL DAYS
  final_cube_root = output_dir+gname+'_'+tag       

; 3D BYTE MASK OF SIGNAL
  signal_mask_ext = '_signal_mask' 

; 3D BYTE MASK OF BRIGHT SIGNAL
  bright_mask_ext = '_bright_mask'

; 2D BYTE MASK OF BRIGHT SIGNAL
  bright_map_ext = '_bright_map' 

; BYTE MASK OF SIGNAL FOUND IN A SMOOTHED VERSION OF THE CUBE
  smooth_mask_ext = '_smooth_mask' 

; AN ASCII FILE TO CONTAIN THE GAINS
  gain_ascii_file = gname+'_'+tag+'_pixel_gains.txt'

; IF NO GAIN FILE IS SUPPLIED BY THE USER, USE THE GAINS FROM THIS REDUCTION
; FORMAT: "day" "pixel name" "gain" "uncertainty"  
  if n_elements(gain_to_apply) eq 0 then $
     gain_to_apply = gain_ascii_file $
  else begin
     test = file_search(gain_to_apply, count=ct)
     if ct eq 0 then begin
        message, 'File '+gain_to_apply+' containing gains not found.', /info
        message, 'Using solutions from this pipeline iteration instead.', /info
        gain_to_apply = gain_ascii_file
     endif
  endelse

; ---------------------------------------------------------------------
; File names
; ---------------------------------------------------------------------

; A SCRATCH FILE USED BY THE GRIDDING AND COADDING ROUTINES
  cube_list_file = 'grid_batch_output_'+tag+'.txt'

; A SCRATCH FILE THAT INDICATES HOW SCANS GROUP INTO GALAXY COVERAGES
; (SHOULD NOT BE REDUCTION-SPECIFIC, SO WE LEAVE THE TAG OUT ... IS THIS OKAY?)
  coverage_list_file = 'coverage_list.txt'

; MASKS RELATED TO THE SPECTRAL BASELINE WINDOWS (DEFINITIONS ABOVE)
  window_mask_file = output_dir+gname+'_'+tag+window_mask_ext+'.fits'

  window_line_file = output_dir+gname+'_'+tag+window_line_ext+'.fits'

  window_line_border_file = $
     output_dir+gname+'_'+tag+window_line_border_ext+'.fits'

  window_noise_file = output_dir+gname+'_'+tag+window_noise_ext+'.fits'

  window_noise_border_file = $
     output_dir+gname+'_'+tag+window_noise_border_ext+'.fits'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; SKIP TO THE USER-SPECIFIED PART OF THE PIPELINE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Using the reviled goto, jump to the requested part of the script. If the
; JUST keyword was also set by the user, we will do one part of the reduction
; and then exit. Otherwise, we step unit the end of one of the iterations.
;
  if keyword_set(goto_dump) then goto, dump
  if keyword_set(goto_ref) then goto, ref
  if keyword_set(goto_fit) then goto, fit
  if keyword_set(goto_cal) then goto, cal
  if keyword_set(goto_split) then goto, split
  if keyword_set(goto_grid) then goto, grid
  if keyword_set(goto_winmask) then goto, winmask
  if keyword_set(goto_coadd) then goto, coadd
  if keyword_set(goto_plait) then goto, plait
  if keyword_set(goto_mask) then goto, mask
  if keyword_set(goto_report) then goto, report

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE A .CLASS SCRIPT TO DUMP THE DATA TO .FITS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; Write a simple CLASS script using the provided file names. When the user
; executes this script, CLASS will write .30m files out as .fits tables, with
; fields containing the time, system temperature, spectrum, etc. for each
; file. IDL can read these files in as structures, making for each analysis.
;

  dump:

  message, 'Writing class script.',/info
  write_class_fits_dump, orig_data_file

; PROMPT THE USER TO QUIT AND RUN THE CLASS SCRIPT  
  message, 'You needs to run the class script to dump .fits by hand.', /info
  message, 'Start CLASS and type "@write_to_fits"', /info
  message, 'Or just copy this to the IDL prompt:', /info
  message, "spawn, 'class @write_to_fits'", /info
  message, 'Pausing. Type .c to continue', /info
  
  if keyword_set(just) then return

  stop

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; REFERENCE-SUBTRACT INDIVIDUAL SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Subtract a reference (OFF, blank sky) measurement from each spectrum.
;

  ref:
  message, 'Building and subtracting reference spectra.',/info

  if keyword_set(no_ref) then begin
;    THE BASIC SPECTRA ARE ALREADY REFERENCE-SUBTRACTED
     ref_sub_spectra_ext = ''
  endif else begin
     ref_batch, orig_data_file $
        , skip_day = skip_day $
        , bad_data = bad_data_file $
        , ref_mask_file = ref_mask_file $
        , ref_from_mask_only = ref_from_mask_only $
;       CONSTRAIN THE REFERENCE MEASUREMENT TIME
        , min_ref = min_ref_time $
        , max_ref = max_ref_time $
;       EXTENSION APPENDED TO OUTPUT FILES
        , out_ext = ref_sub_spectra_ext $
;       KLUGE TO FORCE THE LSR VELOCITY OFFSET
        , force_v0 = force_v0 $
        , show=show
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; PROCESS INDIVIDUAL SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; This step takes the calibrated (output from CLASS) spectra and processes
; them until they are ready to be gridded. In addition to the treatment from
; the first iteration, we have the option of zapping parts of the FFT of each
; spectrum, replacing it with random-phase data with magnitudes interpolated
; from nearby channels.
;

  fit:
  message, 'Fitting baselines.',/info

  spec_batch $
;    INPUT
     , orig_data_file $
     , in_ext = ref_sub_spectra_ext $
;    DEGREE OF FIT
     , degree = degree $
;    DATA FLAGGING
     , skip_day = skip_day $
     , bad_data = bad_data_file $
;    BASELINE FITTING WINDOWS
     , window_root = window_root $
;    FOURIER EDITING
     , zap_fourier = zap_fourier $
     , bad_fft_chan = bad_fft_chan $
;    PATHOLOGICAL SPECTRUM REJECTION
     , no_rejection = no_rejection $
     , decimate = decimate $
     , reject = resid_reject_sigma $
;    EXTENSION APPENDED TO OUTPUT FILES
     , out_ext = base_sub_spectra_ext $
;    KLUGE TO FORCE THE LSR VELOCITY OFFSER
     , force_v0 = force_v0 $
     , show=show
     
  if keyword_set(just) then return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALIBRATE PIXEL GAINS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; This step can only be done after a preliminary reduction. Take bright
; regions in the resulting data cube and compare the spectra derived from
; combining all data (i.e., the resulting cube) to spectra measured by each
; individual pixel over matched areas. Summing over the comparison area,
; derive gain-style corrections to each pixel. The nominal motivation here is
; that the image gain is known to vary among pixels (worse for HERA2 than
; HERA1). This should apparently scatter from tuning to tuning, so we treat
; whole days here, deriving a correction for each pixel on each day.
;

  cal:

  if (cal_pixel_gain eq 1B) then begin

     gain_batch $
;       INPUT: BASELINE-SUBTRACTED SPECTRA
        , orig_data_file $
        , in_ext = base_sub_spectra_ext $
;       DATA FLAGGING
        , skip_day = skip_day $
        , bad_data = bad_data_file $                
;       INPUT: RESULTS OF A PREVIOUS REDUCTION
        , bright_map = prev_bright_map $
        , ref_cube_file = prev_cube $
;       INPUT: A MASK HOLDING THE BASELINE FITTING WINDOWS
        , window_mask_file = prev_mask $
;       OUTPUT: A TEXT FILE
        , gain_ascii_file = gain_ascii_file $
        , show = show     

  endif

  if keyword_set(just) then return 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; IDENTIFY INDIVIDUAL COVERAGES OF THE GALAXY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; In the second iteration, we want to separate coverages of the galaxy. At a
; minimum, we need to separate data observed using different scan
; directions. Ideally, we would like to split each coverage of the galaxy into
; its own file.
;
 
  split:
  
  message, 'Splitting the data into individual coverages.',/info
  
  split_into_maps $
;    INPUT
     , orig_data_file $  
     , map_start_file = map_start_file $
;    SUBSCANS KNOWN TO BE REFERENCES
     , ref_subscan = ref_subscan $
;    DATA FLAGGING                       
     , skip_day = skip_day $               
     , bad_data = bad_data_file $          
     , first_out = 1 $                     
;    FORCE DIVIDING SCAN ANGLE OR REFERENCE PIXEL
     , div_angle = div_angle $
     , pixel = ref_pixel $    
;    OUTPUT AND NAMES
     , output_file = coverage_list_file $ 
     , out_root = coverage_cube_root $
     , n_coverage = n_coverage $           
     , show=show $                         
     , jpeg_name=jpeg_name               

  if keyword_set(just) then return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; GRID BASELINE-SUBTRACTED SPECTRA INTO A DATA CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; This step creates our first data cube. The default approach is to grid each
; day / input file separately and separate them polarization, though this can
; be changed by the user by flipping various flags. The kernel used in the
; gridding may be either a Gaussian with FWHM 1/3 that of the telescope or a
; damped Bessel function.
;

  grid:
  message, 'Gridding OTF maps.',/info

  grid_batch $
;    INPUT
     , orig_data_file $             
     , in_ext = base_sub_spectra_ext $
     , coverage_list = coverage_list_file $
;    DATA FLAGGING
     , bad_data = bad_data_file $
     , skip_day = skip_day $
;    CORRECTIONS TO THE IMAGE GAINS
     , apply_gain = apply_gain $
     , gain_file = gain_to_apply $
;    SPECIFY DATA GROUPINGS TO GRID
     , grid_together = one_big_grid $
     , split_by_pol = split_by_pol $ 
     , split_by_angle = split_by_ang $
     , grid_by_coverage = split_by_cover $
;    SPECIFY GRIDDING METHODOLOGY    
     , gauss_kern = gauss_kern $     
     , gauss_fwhm = gauss_fwhm $
     , grid_by_median = grid_by_median $
     , use_mean_ctr = use_mean_ctr $
;    SPECIFY GALAXY (USED FOR GRID CENTER)
     , galaxy = gname $
;    SPECIFY OUTPUT FILE IF COMBINING DATA
     , together_out = final_cube_root $
;    A TEXT FILE LISTING THE CUBES CREATED
     , cube_list_out = cube_list_file $
     , show=show

  if keyword_set(just) then return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; MAKE MASKS OF THE BASELINE WINDOWS ON THE NEW GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%

  winmask:
  message, 'Making masks of the baseline fitting windows', /info

; IF WE HAVE A WINDOW FILE, MAKE A MASK BASED ON THIS
  if (window_root ne '') then begin
;    GET A HEADER
     readcol, cube_list_file, format='A', cube_list
     target_hdr = headfits(cube_list[0]+'.fits')

     window_to_mask $
        , window_root $
        , target_hdr $
        , out_file = window_mask_file $
        , /ms_to_kms

     window_to_mask $
        , window_root $
        , target_hdr $
        , out_file = window_line_file $
        , /ms_to_kms $
        , /line

     window_to_mask $
        , window_root $
        , target_hdr $
        , out_file = window_line_border_file $
        , /ms_to_kms $
        , /line $
        , /border $
        , n_border = 10

     window_to_mask $
        , window_root $
        , target_hdr $
        , out_file = window_noise_file $
        , /ms_to_kms $
        , /noise

     window_to_mask $
        , window_root $
        , target_hdr $
        , out_file = window_noise_border_file $
        , /ms_to_kms $
        , /noise $
        , /border $
        , n_border = 10
  endif    

  if keyword_set(just) then return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; IF WE HAVE SEVERAL CUBES, COADD THEM INTO A SINGLE CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%

  coadd:
  message, 'Coadding data cubes.',/info

  readcol, cube_list_file, format='A,X,A', cube_list, scan_angle_str
     
; IF WE ARE PLAITING, THE OUTPUT OF THIS STEP NEEDS TO BE CUBES BROKEN DOWN
; BY SCAN ANGLE. OTHERWISE WE GRID EVERYTHING AT ONCE.
  if keyword_set(plait) then begin
     scan_angle = long(scan_angle_str)
     uniq_angles = scan_angle[uniq(scan_angle,sort(scan_angle))]
     n_angles = n_elements(uniq_angles)
     angle_str = '_ang'+strcompress(long(uniq_angles),/remove_all)
  endif else begin
     n_angles = 1
     angle_str = ''
  endelse
  
  for i = 0, n_angles-1 do begin        
;    GET THE DATA FOR THIS SCAN ANGLE
     if keyword_set(plait) then begin
        ind = where(scan_angle eq uniq_angles[i])
        this_cube_list = cube_list[ind]
     endif else $
        this_cube_list = cube_list

;    ADJUST THE OUTPUT FILE NAME        
     this_cube_root= final_cube_root+angle_str[i]
     
     coadd_cube_stack $
;       INPUT
        , this_cube_list $
;       OUTPUT
        , out_root = this_cube_root $
;       LOCATION OF THE NOISE
        , noise_mask = window_noise_file $
;       METHODOLOGY
        , coverage_thresh = coverage_thresh $
        , apodize = apodize $
        , use_coverage = coverage_coadd
  endfor

  if keyword_set(just) then return
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLAIT TWO CUBES WITH DIFFERENT SCAN ANGLES TO GET A FINAL CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  plait:

  if keyword_set(plait) then begin
     message, 'Plaiting OTF maps.',/info

     readcol, cube_list_file, format='A,X,A', cube_list, scan_angle_str
     
     scan_angle = long(scan_angle_str)
     uniq_angles = scan_angle[uniq(scan_angle,sort(scan_angle))]
     n_angles = n_elements(uniq_angles)
     angle_str = '_ang'+strcompress(long(uniq_angles),/remove_all)
     
     if n_angles ne 2 then $
        message, 'Expect EXACTLY TWO angles to plait.', /info

     plait_two_cubes $
        , file_one = final_cube_root+angle_str[0]+'.fits' $
        , angle_one = uniq_angles[0] $
        , file_two = final_cube_root+angle_str[1]+'.fits' $
        , angle_two = uniq_angles[1] $
        , out_file = final_cube_root+'.fits' $
        , mask_file = window_mask_file $
        , show = show $
        , /test

  endif

  if keyword_set(just) then return  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY SIGNAL MASKING TO THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; 
; Make several byte masks that can be used to identify regions of interest in
; the data. The two most important of these are the "window mask," which
; indicates regions over which the spectral baseline fit is expected to be
; valid and the "bright 2d mask" which indicates in a 2d map what regions are
; bright enough that we might attempt the pixel gain calibration. The other
; maps identify signal in the data data cubes and so are interesting from a
; science (not reduction) point of view.
;

  mask:
  message, 'Masking final cube.', /info

  cube_file = final_cube_root+'.fits'
  masked_cube_file = final_cube_root+'.masked.fits'
  smooth_cube_file = output_dir+gname+'_30as_'+tag+'.fits'
  smooth_mask_file = output_dir+gname+smooth_mask_ext+'_'+tag+'.fits'

; MAKE A BASIC SIGNAL-IDENTIFICATION MASK
  make_signal_mask $
     , cube_file $
     , sig = 1.5 $
     , chan = 3 $
     , prior = window_mask_file $
     , out_file = output_dir+gname+signal_mask_ext+'_'+tag+'.fits' $
     , rms_mask_file = window_noise_file

; SMOOTH ...
  conv_with_gauss $
     , infile = cube_file $
     , outfile = smooth_cube_file $
     , target_fwhm = 30. $
     , quiet = quiet $
     , /grow_nan $
     , /cube

; ... AND MASK
  make_signal_mask $
     , smooth_cube_file $
     , sig = 3 $
     , chan = 2 $
     , prior = window_mask_file $
     , out_file = smooth_mask_file $
     , rms_mask_file = window_noise_file

; MAKE A MORE STRICT SIGNAL-IDENTIFICATION MASK (BOTH 2 AND 3D)
  make_signal_mask $
     , cube_file $
     , sig = 3 $
     , chan = 3 $
     , min_area = 15 $
     , grow = 3 $
     , prior = smooth_mask_file $
     , out_file = output_dir+gname+bright_mask_ext+'_'+tag+'.fits' $
     , map_out_file = output_dir+gname+bright_map_ext+'_'+tag+'.fits' $
     , rms_mask_file = window_noise_file
  
; APPLY THE WINDOW MASK TO THE DATA
  cube = readfits(cube_file, hdr)
  mask = readfits(window_mask_file)
  bad = where(mask eq 0, badct)
  if badct gt 0 then cube[bad] = !values.f_nan
  writefits, masked_cube_file, cube, hdr

  if keyword_set(just) then return

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPARE SPECTRA FROM DIFFERENT DAYS IN THE BRIGHT REGIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  report:

  if keyword_set(one_big_grid) eq 0 then begin
     readcol, cube_list_file, format='A,X,X', cube_list
     cube_list += '.fits'
   
     compare_cube_stack $
        , cube_list $
        , ref_cube = final_cube_root+'.fits' $
        , bright_map = bright_map_file $
        , window_mask = window_line_file $
        , psfile = report_dir+'cube_comparison_'+tag+'.ps'
  endif

  if keyword_set(just) then return  

end                             ; of HERACLES pipeline
