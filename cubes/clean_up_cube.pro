pro clean_up_cube $
   , input_file $
   , out_file = output_file $
   , coverage_cube = coverage_file $
   , coverage_thresh = coverage_thresh $
   , blank_mask = blank_mask_file $
   , apodize=apodize $
   , apod_kern = apod_kern 
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  cube = readfits(input_file, hdr)

; MEASURE THE SIZE OF THE CUBE
  sz = size(cube)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE REGIONS INDICATED IN THE "BLANK MASK" INTO NOT-A-NUMBERS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; CHECK FOR THE EXISTENCE OF THE FILE
  dummy = file_search(blank_mask_file, count=blank_mask_ct)

  if blank_mask_ct eq 1 then begin

;    READ THE MASK AND CHECK ITS SIZE AGAINST THE CUBE
     blank_mask = readfits(blank_mask_file, blank_mask_hdr)

     blank_sz = size(blank_mask)

     if (sz[1] ne blank_sz[1]) or (sz[2] ne blank_sz[2]) or  $
        (sz[3] ne blank_sz[3]) then $
           message, 'Mask and cube have mismatched sizes.'
     
;    BLANK THE CUBE ANYWHERE THAT THE BLANK MASK IS 1
     blank_ind = where(blank_mask eq 1, blank_ct)
     
     if blank_ct gt 0 then $
        cube[blank_ind] = !values.f_nan

     sxaddpar, hdr, 'HISTORY', 'Applied blanking mask.'

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BLANK REGIONS WITH LITTLE COVERAGE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%


; CHECK FOR THE EXISTENCE OF THE FILE
  dummy = file_search(coverage_file, count=cov_ct)

  if n_elements(coverage_thresh) gt 0 and cov_ct eq 1 then begin

;    READ THE COVERAGE FILE AND CHECK ITS SIZE AGAINST THE CUBE
     coverage = readfits(coverage_file, coverage_hdr)

     cov_sz = size(coverage)
     
     if (sz[1] ne cov_sz[1]) or (sz[2] ne cov_sz[2]) or  $
        (sz[3] ne cov_sz[3]) then $
           message, 'Mask and cube have mismatched sizes.'
     
;    CALCULATE THE MEDIAN COVERAGE WHERE WE HAVE DATA
     fin_ind = where(finite(cube) and coverage ne 0, fin_ct)

     if fin_ct eq 0 then $
        message, 'No regions with finite data and non-zero covearge.'

     med_coverage = median(coverage[fin_ind])

;    FLAG DATA WITH LESS THAN coverage_thresh OF THIS COVERAGE
     low_coverage = where(coverage lt med_coverage*coverage_thresh $
                          , low_ct)

     if low_ct gt 0 then $
        cube[low_coverage] = !values.f_nan

     sxaddpar, hdr, 'HISTORY', 'Applied coverage threshold.'

  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; USE A SMOOTHED VERSION OF THE CUBE TO TAPER EDGES AND FILL HOLES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(apodize) then begin

;    DEFAULT TO A KERNEL EQUAL IN SIZE TO THE BEAM
     if n_elements(apod_kern) eq 0 then begin
        apod_kern = sxpar(hdr, 'BMAJ') * 3600.        
     endif

     apod_kern_pix = $
        apod_kern / (abs(sxpar(hdr, 'CDELT1'))*3600.)

;    MAKE A COPY OF THE CUBE SMOOTHED WITH THIS KERNEL
     conv_with_gauss $
        , in_data = cube $
        , in_hdr = hdr $
        , out_data = apod_cube $
        , out_hdr = apod_hdr $
        , /cube $
        , /quiet

;    FIND LOCATIONS WITH NON-FINITE DATA CLOSE TO FINITE DATA
     finite_mask = finite(cube)

     for i = 0, sz[3]-1 do begin
        finite_mask[*,*,i] = $
           exp_mask(finite_mask[*,*,i], rad = ceil(apod_ker_pix / 2))
     endfor

     apod_targets = where(finite(cube) eq 0 and $
                          finite_mask eq 1, apod_ct)

     if apod_ct gt 0 then $
        cube[apod_targets] = apod_cube[apod_targets]

;    REPLACE MISSING DATA WITHIN A FWHM OF EXISTING WITH THIS SMOOTH DATA

     sxaddpar, hdr, 'HISTORY' $
               , 'Smoothed edges and empties with apodized cube.'

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE THE RESULTING CUBE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  writefits, output_file, cube, hdr

end

