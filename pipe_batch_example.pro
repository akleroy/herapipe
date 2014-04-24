pro pipe_batch, spec=spec, cube=cube, win=win, all=all

  if keyword_set(win) or keyword_set(all) then begin

;     make_flat_windows $
;        , out_root = '../masks/windows/ngc4303_flat' $
;        , gname = 'ngc4303' $
;        , window = 300.

;     make_flat_windows $
;        , gname = 'ngc4303' $
;        , out_root = '../masks/windows/ngc4303_loweredge' $
;        , /hel_to_lsr $
;        , /ms_to_kms $
;        , window = 1000. $
;        , offset = -900.

;     make_flat_windows $
;        , gname = 'ngc4303' $
;        , out_root = '../masks/windows/ngc4303_upperedge' $
;        , /hel_to_lsr $
;        , /ms_to_kms $
;        , window = 1000. $
;        , offset = +900.

;     mask_to_windows $
;        , '../masks/windows/ngc4303_mask.fits' $
;        , out_root = '../masks/windows/ngc4303_line' $
;        , /ms_to_kms

  endif

  if keyword_set(spec) or keyword_set(all) then begin

     spectra_pipeline $
        , identifier = '_beta' $
        , orig_data_file = 'orig_data.txt' $
;        , bad_data_file = 'bad_data.txt' $
;        , ref_mask_file = '../masks/reference/ngc4303_on_mask.fits' $
;        , window_root = $ 
;        ['../masks/windows/ngc4303_loweredge' $
;         , '../masks/windows/ngc4303_line' $
;         , '../masks/windows/ngc4303_upperedge'] $
        , degree = 3 $     
        , /show $
        , /report $
        , smooth=[1,10]
     
  endif

  if keyword_set(cube) or keyword_set(all) then begin

     cube_pipeline $
        , out_root = 'ngc4303' $
        , identifier = '_beta' $
        , orig_data_file = 'orig_data.txt' $
        , blank_window_root = $ 
        ['../masks/windows/ngc4303_loweredge' $
         , '../masks/windows/ngc4303_upperedge'] $
        , mask_window_root = $
        ['../masks/windows/ngc4303_line'] $
        , degree = 3 $     
        , cal=0 $
        , prev_cube = '../cubes/_beta/ngc4303_beta.fits' $
        , prev_mask_2d =  $
        '../cubes/_beta/ngc4303_beta_bright_map.fits' $
        , prev_mask_3d = '../cubes/_beta/ngc4303_beta_mask.fits' $
        , /show $         
        , /report

  endif

end
