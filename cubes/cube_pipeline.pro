;+
; NAME:
;
; cube_pipeline
;
; PURPOSE:
;
; Take reduce HERA spectra and turn them into science-quality data cubes.
;
; CATEGORY:
;
; one of two data reduction master scripts
;
; CALLING SEQUENCE:
;
; cube_pipeline
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
; data cubes
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
; PROCEDURE:
;
; 1) (optionally) solve for the gain of each pixel + day
; 2) grid data according to scan direction
; 3) use windows files to blank the cube and make masks
; 3) (optionally) baseline subtract the cubes
; 4) plait the two cubes together
; 5) measure the noise properties of the cube
; 6) make a series of masks
;
; MODIFICATION HISTORY:
;
; spun out of single pipeline, streamlined - Oct 2009 leroy@mpia.de
;
;-

pro cube_pipeline $
;  GENERAL INPUT/OUTPUT
   , identifier = tag $
   , output_identifier = out_tag $
   , orig_data_file = orig_data_file $
   , out_root = out_root $
;  TUNING PARAMETERS:: GAIN SOLUTION
   , cal_pixel_gain = cal_pixel_gain $
   , prev_cube = prev_cube $
   , prev_mask_2d = prev_mask_2d $
   , prev_mask_3d = prev_mask_3d $
;  TUNING PARAMETERS:: GRIDDING
   , force_header = user_header $
   , gal_for_ctr = gname $
   , use_mean_ctr = use_mean_ctr $
   , pix_scale = pix_scale $
   , gauss_kern = gauss_kern $
   , gauss_fwhm = gauss_fwhm $
   , grid_by_median = grid_by_median $
   , coverage_thresh = coverage_thresh $   
;  TUNING PARAMETERS:: MASKING BASED ON FITTING WINDOWS
   , blank_window_root = blank_window_root $
   , mask_window_root = mask_window_root $  
;  TUNING PARAMETERS:: FITTING AFTER GRIDDING
   , refit = refit $
   , degree = degree $          
;  TUNING PARAMETERS:: PLAITING
   , plait = plait $
;   , scale_plait = plait_scale $
;  TUNING PARAMETERS:: MASKING
;  OUTPUT SWITCHES
   , show=show $
   , report = report $
;  ROUTE TO SPECIFIC PARTS OF THE PIPELINE
   , goto_cal  = goto_cal $
   , goto_astro = goto_astro $
   , goto_grid  = goto_grid $
   , goto_winmask   = goto_winmask $
   , goto_clean = goto_clean $
   , goto_refit = goto_refit $
   , goto_plait = goto_plait $
   , goto_noise = goto_noise $
   , goto_mask  = goto_mask $
   , just = just

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; DEFAULTS AND ERROR CHECKING ON USER INPUT
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

; AN IDENTIFIER FOR THE INPUT 
  if n_elements(tag) eq 0 then tag = ''

; AN IDENTIFIER FOR THE OUTPUT
  if n_elements(out_tag) eq 0 then out_tag = tag

; ---------------------------------------------------------------------
; Output files
; ---------------------------------------------------------------------

; OUTPUT DIRECTORY FOR REPORTS/IMAGES
  output_dir = '../cubes/'
  if out_tag ne '' then $
     output_dir += out_tag+'/'

; CHECK THAT THE DIRECTORY IS THERE...
  test = file_search(output_dir, count = ct)
  if ct ne 1 then begin
     message, 'The output directory '+output_dir+' is missing.', /info
     message, 'I am creating it.', /info
     spawn, 'mkdir '+output_dir
  endif

; A DEFAULT CUBE NAME
  if n_elements(out_root) eq 0 then $
     out_root = 'cube'

; ---------------------------------------------------------------------
; Gain Calibration
; ---------------------------------------------------------------------

; DEFAULT TO NO GAIN CALIBRATION
  if n_elements(cal_pixel_gain) eq 0 then cal_pixel_gain = 0

  if keyword_set(cal_pixel_gain) then begin

;    AN ASCII FILE TO CONTAIN THE GAINS
     gain_ascii_file = 'pixel_gains'+out_tag+'.txt'

;    CHECK IF WE HAVE ALL THE FILES THAT WE NEED
     missing_files = $
        (n_elements(prev_cube) eq 0) or $
        (n_elements(prev_mask_2d) eq 0) or $
        (n_elements(prev_mask_3d) eq 0)
     
     if missing_files then begin
        message, 'Need a previous cube, map, and mask to '+$
                 'solve for pixel gains.', /info
        message, 'Skipping this step.', /info
        cal_pixel_gain = 0
     endif
     
  endif
  
; ---------------------------------------------------------------------
; Gridding
; ---------------------------------------------------------------------

; FRACTION OF MEDIAN COVERAGE TO BE INCLUDED IN CUBE
  if n_elements(coverage_thresh) eq 0 then $
     coverage_thresh = 0.25

; TYPE OF CONVOLUTION KERNEL
  if n_elements(gauss_kern) eq 0 then $
     gauss_kern = 1

; SIZE OF CONVOLUTION KERNEL
  if n_elements(gauss_fwhm) eq 0 then $
     gauss_fwhm = 8.1/3600.

; ---------------------------------------------------------------------
; Spectral Fitting Windows
; ---------------------------------------------------------------------

; ---------------------------------------------------------------------
; Refitting Baselines
; ---------------------------------------------------------------------

; DEFAULT TO NOT DOING THIS STEP
  if n_elements(refit) eq 0 then $
     refit = 0

; DEGREE OF THE POLYNOMIAL USED TO FIT THE BASELINE
  if n_elements(degree) eq 0 then $
     degree = 5

; ---------------------------------------------------------------------
; Plaiting
; ---------------------------------------------------------------------

; TURN ON PLAITING BY DEFAULT, WE LIKE PLAITING
  if n_elements(plait) eq 0 then $
     plait = 1

; ---------------------------------------------------------------------
; Measuring Noise
; ---------------------------------------------------------------------

; ---------------------------------------------------------------------
; Masking
; ---------------------------------------------------------------------


; ... END OF DEFAULTS

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; SKIP TO THE USER-SPECIFIED PART OF THE PIPELINE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Using the reviled goto, jump to the requested part of the script. If the
; JUST keyword was also set by the user, we will do one part of the reduction
; and then exit. Otherwise, we step unit the end of one of the iterations.
;
  if keyword_set(goto_cal) then goto, cal
  if keyword_set(goto_astro) then goto, astro
  if keyword_set(goto_grid) then goto, grid
  if keyword_set(goto_winmask) then goto, winmask
  if keyword_set(goto_clean) then goto, clean
  if keyword_set(goto_refit) then goto, refit
  if keyword_set(goto_plait) then goto, plait
  if keyword_set(goto_noise) then goto, noise
  if keyword_set(goto_mask) then goto, mask

; ... END GOTO RE-ROUTING STATEMENT
  
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
; HERA1). Not knowing the image gain introduces an uncertainty into the
; chopper wheel calibration, which means the flux scale is uncertain. This
; should apparently scatter from tuning to tuning, so we treat whole days
; here, deriving a correction for each pixel on each day. K.S. suggests that
; order 20% effects may be expected but factor of two would be alarming.
;

  cal:

  if (cal_pixel_gain eq 1B) then begin

     message, 'Attempting pixel gain calibration.', /info

     solve_for_gains $
        , orig_data_file $
        , tag = tag $
        , prev_cube = prev_cube $
        , prev_mask_2d = prev_mask_2d $
        , prev_mask_3d = prev_mask_3d $
        , gain_ascii_file = gain_ascii_file $
        , report = report $
        , show = show     

  endif

  if keyword_set(just) then return 

; ... END OF CALIBRATING PIXEL GAINS

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; FIGURE OUT THE ASTROMETRY AND LIST OF CUBES THAT WE WILL MAKE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; This step creates our first data cube. The default approach is to grid each
; day / input file separately and separate them polarization, though this can
; be changed by the user by flipping various flags. The kernel used in the
; gridding may be either a Gaussian with FWHM 1/3 that of the telescope or a
; damped Bessel function.
;

  astro:
  message, 'Building target astrometry.', /info

; MAKE A HEADER BASED ON THE AVAILABLE DATA
  target_hdr = build_header( $
               orig_data_file $
               , tag = tag $
               , galaxy = galaxy $
               , use_mean_ctr = use_mean_ctr $
               , pix_scale = pix_scale $
               , split_by_angle = plait $
               , uniq_ang = uniq_ang)

  if n_elements(user_header) gt 0 then $
     target_hdr = user_header

; MAKE A LIST OF CUBES THAT WE INTEND TO MAKE
  if keyword_set(plait) then $
     cube_list = $
     output_dir + out_root + $
     'ang' + strcompress(indgen(n_elements(uniq_ang))+1 $
                         , /remove_all) $
     + out_tag $
  else $
     cube_list = output_dir + out_root + out_tag

; SAVE THE HEADER AND CUBE LIST TO AN IDL FILE
  save, target_hdr, uniq_ang, plait, cube_list $
        , file = 'cube_info'+out_tag+'.idl'

  if keyword_set(just) then return 

; ... END OF BUILDING ASTROMETRY AND CUBE LIST

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

  restore, 'cube_info'+out_tag+'.idl', /v
  
  if n_elements(cube_list) eq 0 then begin
     message, 'Cannot find cube information.'
  endif

  for i = 0, n_elements(cube_list)-1 do begin

     if keyword_set(plait) then $
        scan_ang = uniq_ang[i] $
     else $
        scan_ang = !values.f_nan

;    GRID DATA EITHER ALL TOGETHER
     grid_wrapper $
        , orig_data_file $
        , tag = tag $
        , outfile = cube_list[i]+'_raw' $
        , target_hdr = target_hdr $
        , apply_gain = apply_gain $
        , gain_file = gain_ascii_file $
        , split_by_angle = plait $
        , scan_angle = scan_ang $
        , grid_by_median=grid_by_median $
        , gauss_kern = gauss_kern $  
        , gauss_fwhm = gauss_fwhm $  
        , show=show                  
     
  endfor

  if keyword_set(just) then return

; ... END OF GRIDDING

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
; MAKE MASKS OF THE BASELINE WINDOWS ON THE NEW GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%%&%
;
; Take the windows used to blank parts of the spectrum during basline fittings
; and convert them into three dimensional byte masks on the same grid as the
; cube. We can then use these to refit baselines and blank the parts of the
; cube that have not been baseline subtracted.
;

  winmask:
  message, 'Making masks of the baseline fitting windows', /info

  restore, 'cube_info'+out_tag+'.idl', /v

  windows_to_mask $
     , blank_window_root $
     , target_hdr = target_hdr $
     , out_file = output_dir + 'mask_blank' + out_tag + '.fits' $
     , /ms_to_kms

  windows_to_mask $
     , mask_window_root $
     , target_hdr = target_hdr $
     , out_file = output_dir + 'mask_line' + out_tag + '.fits' $
     , /ms_to_kms

  if keyword_set(just) then return

; ... END OF APPLYING SPECTRAL FITTING WINDOWS

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CLEAN UP THE CUBES (BLANK, SUPPRESS LOW COVERAGE, APODIZE)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; Use the window masks and coverage cubes to clean up the data cubes that we
; just gridded. Blank regions that have not had baselines subtracted, remove
; regions with little or no data, and apodize to smooth out edges and fill in
; missing (if requested).
;

  clean:
  message, 'Cleaning up the data cube(s)', /info

  restore, 'cube_info'+out_tag+'.idl', /v
  
  for i = 0, n_elements(cube_list)-1 do begin

     clean_up_cube $
        , cube_list[i]+'_raw.fits' $
        , out_file = cube_list[i]+'.fits' $
        , coverage_cube = cube_list[i]+'_raw.coverage.fits' $
        , coverage_thresh = coverage_thresh $
        , blank_mask = output_dir + 'mask_blank' + out_tag + '.fits' $
        , apodize=apodize $
        , apod_kern = apod_kern     

  endfor

  if keyword_set(just) then return


; ... END OF CLEANUP

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; REFIT BASELINES, NOW TO THE GRIDDED CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; Optionally fit new baselines to the cubes that we have just
; constructed. Everything up to this point is pretty much linear, so this is
; perfectly rigorous. Gives a chance to get at subtler features that might
; have been missed by the fits to individual spectra.
;
  
  refit:
  message, '(Re)fitting baselines to the data cube(s)', /info

  restore, 'cube_info'+out_tag+'.idl', /v

  for i = 0, n_elements(cube_list)-1 do begin

     fit_base_to_cube $
        , cube_list[i]+'.fits' $
        , out_file = cube_list[i]+'.fits' $
        , degree = degree $
        , mask = output_dir + 'mask_line' + out_tag + '.fits' $
        , show = show

  endfor

  if keyword_set(just) then return

; ... END OF REFITTING BASELINES
  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLAIT TWO CUBES WITH DIFFERENT SCAN ANGLES TO GET A FINAL CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; If we have gridded based on scan angle then combine the two cubes in fourier
; space, suppressing large scale features in the scan direction. At the cost
; of sqrt(2) higher noise in low-frequency features we gain significant
; artifact suppression.
;

  plait:

  if keyword_set(plait) then begin
     message, 'Plaiting OTF maps.',/info
     
     restore, 'cube_info'+out_tag+'.idl', /v
     
     if n_elements(cube_list) eq 0 then $
        message, 'Cannot find cube information.'
     
     n_angles = n_elements(uniq_ang)     
     if n_angles ne 2 then $
        message, 'Expect EXACTLY TWO angles to plait.'

     plait_two_cubes $
        , file_one = cube_list[0]+'.fits' $
        , angle_one = uniq_ang[0] $
        , file_two = cube_list[1]+'.fits' $
        , angle_two = uniq_ang[1] $
        , out_file = output_dir + out_root + out_tag + '.fits' $
        , show = show $
        , /test

  endif

  if keyword_set(just) then return  

; ... END OF PLAITING

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MEASURE THE NOISE PROPERTIES OF EACH CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; Measure the x,y,z dependence of noise in the resulting cube. 
;

  noise:
  
  make_noise_cube $
     , cube_file = output_dir + out_root + out_tag + '.fits' $
     , out_file = output_dir + 'noise' + out_tag + '.fits' $
     , mask_file = output_dir + 'mask_line' + out_tag + '.fits' $
     , show=show

  if keyword_set(just) then return  

; ... END OF MEASURING NOISE PROPERTIES

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; APPLY SIGNAL MASKING TO THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; 
; Make several byte masks that can be used to identify regions of interest in
; the data.
;

  mask:
  message, 'Masking cube.', /info

  heracles_masking $
     , output_dir + out_root + out_tag + '.fits' $
     , prior = output_dir + 'mask_line' + out_tag + '.fits' $
     , mask_file = output_dir + out_root + out_tag + '_mask.fits' $
     , bright_map_file = output_dir + out_root + out_tag + '_bright_map.fits' $
     , show=show
  
  if keyword_set(just) then return

; ... END OF SIGNAL-ORIENTED MASKING

end                             ; of CUBE portion of HERACLES pipeline
