pro heracles_masking $
   , cube_file $
   , prior = prior_file_in $
   , mask_file = mask_file_out $
   , bright_map_file = bright_map_file_out $
   , show = show

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  cube = readfits(cube_file, hdr)

  if n_elements(prior_file_in) then begin
     prior = readfits(prior_file_in, prior_hdr)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SMOOTH AND MASK, USING THE LINE MASK AS A PRIOR
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; SMOOTH THE CUBE TO 30 ARCSECOND RESOLUTION
  conv_with_gauss $
     , in_data = cube $
     , in_hdr = hdr $
     , out_data = smooth_cube $
     , out_hdr = smooth_hdr $
     , target_fwhm = 30. $
     , /cube

; GET THE RMS BY POSITION OF THE CUBE
  make_iram_noise_cube $
     , cube_in = smooth_cube $
     , mask_in = prior $
     , rms_cube = smooth_rms_cube $
     , show = show

; MAKE THE MASK
  mask = make_signal_mask( $
         smooth_cube $
         , rms = smooth_rms_cube $
         , prior = prior $
         , sig_thresh = 3. $
         , nchan = 2 $
         , grow_xy = 7 $
         , grow_z = 4)
  
; WRITE TO DISK
  if n_elements(mask_file_out) gt 0 then begin
     
     mask_hdr = hdr
     sxaddpar, mask_hdr, 'BUNIT','MASK','1=SIGNAL,0=NONE'
     sxaddpar, mask_hdr, 'History', 'Masked using smooth-and-mask approach'
     writefits, mask_file_out, mask, mask_hdr

  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOK FOR BRIGHT EMISSION, USING THE SMOOTH MASK AS A PRIOR
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; GET THE RMS BY POSITION OF THE CUBE
  make_iram_noise_cube $
     , cube_in = cube $
     , mask_in = mask $
     , rms_cube = rms_cube $
     , show = show

; MAKE THE MASK
  bright_mask = make_signal_mask( $
                cube $
                , rms = rms_cube $
                , prior = mask $
                , sig_thresh = 3. $
                , nchan = 3 $
                , min_area = 15 $
                , grow_xy = 3 $
                , grow_z = 3)
  
  bright_map = total(bright_mask,3) gt 0

  if keyword_set(show) then $
     disp, bright_map

; WRITE TO DISK
  if n_elements(bright_map_file_out) gt 0 then begin

     bright_map_hdr = twod_head(hdr)
     sxaddpar, bright_map_hdr, 'BUNIT','MASK','1=SIGNAL,0=NONE'
     sxaddpar, bright_map_hdr, 'History', 'Masked identifying bright regions.'

     writefits, bright_map_file_out, bright_map, bright_map_hdr

  endif

end
