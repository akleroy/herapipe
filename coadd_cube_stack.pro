pro coadd_cube_stack $
;  INPUTS
   , cube_list_in $
   , noise_mask = location_of_noise $
   , align = align $ 
   , target_hdr = target_hdr $
;  OUTPUTS
   , out_root = out_root $
;  METHODOLOGY
   , plane_by_plane = plane_by_plane $
   , use_disk = use_disk $
   , use_median = use_median $
   , use_coverage = use_coverage $
   , equal_weight = equal_weight $
   , coverage_thresh = thresh_coverage $
   , apodize = apodize $
   , apod_kern_fwhm = apod_kern_fwhm
  
;+
; NAME:
;
; coadd_cube_stack
;
; PURPOSE:
;
; Add together a stack of cubes on the target grid (default is to use the
; header for the first cube). Currently REQUIRES that they share a matched
; velocity axis.
;
; Optionally, weight by another set of (coverage) cubes or by the RMS of the
; cubes. Use either the median or mean to add. Reject elements with less than
; a specified level of coverage.
;
; CATEGORY:
;
; Reduction and science tool.
;
; CALLING SEQUENCE:
;
; coadd_cube_stack, cube_list_in, out_root = out_root $
;   , /align, target_hdr = target_hdr $
;   , noise_mask = noise_mask $
;   , /use_median, /use_coverage $
;   , coverage_thresh = 0.1 $
;   , /apodize, apod_kern_fwhm = 20./3600.
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; DEFAULTS AND ERROR CHECKING
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; READ THE INPUTPUT LIST
  cube_list_in = strcompress(cube_list_in,/remove_all)
  nf = n_elements(cube_list_in)  
  
; VERIFY THAT THE CUBES EXIST
  exists = bytarr(nf)
  for i = 0, nf-1 do begin
     check = file_search(cube_list_in[i]+'.fits',count = ct)
     exists[i] = ct eq 1
     if ct eq 0 then $
        message, 'Warning - file '+cube_list_in[i]+' is missing',/info
  endfor
  use = where(exists, nf)
  cube_list = cube_list_in[use]

; CHECK WHETHER A MASK SHOWING THE LOCATION OF NOISE HAS BEEN SUPPLIED
  if n_elements(location_of_noise) eq 0 then begin
     count = 0
     have_noise_mask = 0
  endif else begin
     noise_mask_file = file_search(location_of_noise, count=count)
     have_noise_mask = count gt 0
  endelse

  if have_noise_mask then begin
     noise_mask = readfits(noise_mask_file)    
     if total(noise_mask) eq 0 then $
        have_noise_mask = 0
  endif

; GIVE A DEFAULT OUTPUT  
  if n_elements(out_root) eq 0 then out_root = 'combined'

; DEFAULT TO REQUIRING 50% OF MEDIAN COVERAGE (TOO HIGH?)
  if n_elements(thresh_coverage) eq 0 then thresh_coverage = 0.50

; GET THE TARGET HEADER
  if n_elements(target_hdr) eq 0 or keyword_set(align) eq 0 then begin
     target_hdr = headfits(cube_list[0]+'.fits')
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; INITIALIZE OUTPUTS (EMPTY COVERAGE AND DATA CUBES)
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  coadd = fltarr(sxpar(target_hdr,'NAXIS1') $
                 , sxpar(target_hdr,'NAXIS2') $
                 , sxpar(target_hdr,'NAXIS3'))*!values.f_nan
  sz = size(coadd)  
  coverage = fltarr(sz[1],sz[2],sz[3])
  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; HANDLE THE ALL-AT-ONCE CASE...
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... FIRST READ ALL THE DATA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if keyword_set(plane_by_plane) eq 0 then begin
     data_stack = fltarr(n_elements(coadd),nf)
     for i = 0L, nf-1 do begin
        cube = readfits(cube_list[i]+'.fits', this_hdr, /silent)
        if keyword_set(align) then begin        
           cube = cube_hastrom(cube_in = cube $
                               , hdr_in = this_hdr $
                               , target = target_hdr $
                               , /cubic, missing=!values.f_nan)
        endif
        data_stack[*,i] = reform(cube)     
     endfor
     
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... THEN READ THE COVERAGE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     if keyword_set(use_coverage) then begin
        coverage_stack = fltarr(n_elements(coadd),nf)
        for i = 0L, nf-1 do begin
           this_coverage = readfits(cube_list[i]+'.coverage.fits', this_hdr, /silent)
           if keyword_set(align) then begin        
              this_coverage = cube_hastrom(cube_in = this_coverage $
                                           , hdr_in = this_hdr $
                                           , target = target_hdr $
                                           , /cubic, missing=!values.f_nan)
           endif
           coverage_stack[*,i] = reform(this_coverage)
        endfor
        coverage_stack *= finite(data_stack)*1.0
     endif else begin
        coverage_stack = finite(data_stack)*1.0
     endelse
     
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... GET A NOISE MEASUREMENT FOR EACH CUBE TO USE AS A WEIGHT
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
     
     rms_vec = fltarr(nf)*!values.f_nan
     
     if keyword_set(equal_weight) eq 0 then begin
        for i = 0L, nf-1 do begin
           if have_noise_mask then begin
              rms_vec[i] = mad((data_stack[*,i])[where(noise_mask)])
           endif else begin
              rms_vec[i] = mad((data_stack)[*,i])
           endelse
           coverage_stack[*,i] *= 1.0/(rms_vec[i])^2
        endfor
     endif
     
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... SUM IT UP AT EACH POSITION TO GET THE FINAL COVERAGE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     if nf gt 1 then $
        coverage[*,*,*] = reform(total(coverage_stack,2)) $
     else $
        coverage[*,*,*] = reform(coverage_stack)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... IDENTIFY WHERE TO COADD
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     median_coverage = $
        median(coverage[where(coverage ne 0)])
     coadd_here = coverage / median_coverage ge thresh_coverage
     where_to_coadd = where(coadd_here, coadd_ct)
     
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... COADD
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     if nf eq 1 then begin
        coadd[where_to_coadd] = data_stack[where_to_coadd]
     endif

     if coadd_ct gt 0 and nf gt 1 then begin
        coadd[where_to_coadd] = $
           coadd_data(data_stack[where_to_coadd,*] $
                      , weight = coverage_stack[where_to_coadd,*] $
                      , use_median = use_median)
     endif
  endif

; ... FROM HERE WILL SKIP TO END (BYPASSING IF STATEMENT)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; HANDLE THE PLANE BY PLANE CASE ...
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
  
  if keyword_set(plane_by_plane) then begin

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... READ THE CUBES ONE AT A TIME TO GET THE RMS ...
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
     rms_vec = fltarr(nf)*!values.f_nan

     if keyword_set(use_disk) then begin
        scratch_dir = 'coadd_cube_stack_scratch/'
        spawn, 'mkdir '+ scratch_dir
     endif

     if keyword_set(equal_weight) eq 0 then begin        
        for i = 0, nf-1 do begin
           counter, i+1, nf, 'Pre-work on cube '
           data = readfits(cube_list[i]+'.fits',this_hdr,/silent)
           if keyword_set(align) and $
              (have_noise_mask or keyword_set(use_disk)) then begin        
              data = cube_hastrom(cube_in = data $
                                  , hdr_in = this_hdr $
                                  , target = target_hdr $
                                  , /cubic, missing=!values.f_nan)
           endif

           if keyword_set(use_coverage) and keyword_set(use_disk) then begin
              this_coverage = readfits(cube_list[i]+'.coverage.fits',this_hdr,/silent)
              if keyword_set(align) then begin        
                 this_coverage = cube_hastrom(cube_in = this_coverage $
                                              , hdr_in = this_hdr $
                                              , target = target_hdr $
                                              , /cubic, missing=!values.f_nan)
              endif
           endif
           
           if have_noise_mask then begin
              rms_vec[i] = mad(data[where(noise_mask)])
           endif else begin
              rms_vec[i] = mad(data)
           endelse

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... AND, IF REQUESTED, TO WRITE OUT PLANE-BY-PLANE .FITS FILES TO REPLACE
; THE CUBES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
           if keyword_set(use_disk) then begin
              for j = 0, sz[3]-1 do begin
                 plane = data[*,*,j]
                 plane_name = scratch_dir + $
                    'cube_'+strcompress(i,/remove_all)+'_plane_'+ $
                    strcompress(j,/remove_all)+'.fits'
                 writefits, plane_name, plane, twod_head(this_hdr)
                 if keyword_set(use_coverage) then begin
                    this_coverage_plane = this_coverage[*,*,j]
                    plane_name = scratch_dir + $
                       'coverage_'+strcompress(i,/remove_all)+'_plane_'+ $
                       strcompress(j,/remove_all)+'.fits'
                    writefits, plane_name, this_coverage_plane, twod_head(this_hdr)                 
                 endif
              endfor
           endif
        endfor
     endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... NOW LOOP OVER PLANES
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
     for k = 0, sz[3]-1 do begin
        counter, k+1, sz[3], 'Plane'

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... READ THE DATA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
        
        data_stack = fltarr(sz[1]*sz[2],nf)
        for i = 0L, nf-1 do begin
           if keyword_set(use_disk) then begin
              plane_name = scratch_dir + $
                 'cube_'+strcompress(i,/remove_all)+'_plane_'+ $
                 strcompress(k,/remove_all)+'.fits'
              data = readfits(plane_name, /silent)
           endif else begin
              data = (readfits(cube_list[i]+'.fits', this_hdr, /silent))[*,*,k]
              this_hdr = twod_head(this_hdr)
              if keyword_set(align) then begin        
                 hastrom, data, this_hdr, twod_head(target_hdr) $
                          , interp=2, cubic=-0.5, missing=!values.f_nan
              endif
           endelse
           data_stack[*,i] = reform(data)     
        endfor 

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... READ THE DATA
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

        if keyword_set(use_coverage) then begin
           coverage_stack = fltarr(sz[1]*sz[2],nf)
           for i = 0L, nf-1 do begin
              if keyword_set(use_disk) then begin
                 plane_name = scratch_dir + $
                    'coverage_'+strcompress(i,/remove_all)+'_plane_'+ $
                    strcompress(k,/remove_all)+'.fits'
                 this_coverage = readfits(plane_name,/silent)
              endif else begin
                 this_coverage = $
                    (readfits(cube_list[i]+'.coverage.fits', this_hdr, /silent))[*,*,k]
                 this_hdr = twod_head(this_hdr)
                 if keyword_set(align) then begin     
                    hastrom, this_coverage $
                             , this_hdr, twod_head(target_hdr) $
                             , interp=2 $
                             , cubic=-0.5 $
                             , missing=!values.f_nan
                 endif
              endelse
              coverage_stack[*,i] = reform(this_coverage)
           endfor
        endif else begin
           coverage_stack = finite(data_stack)*1.0
        endelse

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... APPLY NOISE MEASUREMENTS AS WEIGHTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

        if keyword_set(equal_weight) eq 0 then begin
           for i = 0, nf-1 do $
              coverage_stack[*,i] *= 1./(rms_vec[i])^2
        endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... SUM IT UP TO GET THE FINAL COVERAGE
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&        

        coverage_plane = fltarr(sz[1],sz[2])
        coverage_plane[*,*] = reform(total(coverage_stack,2))

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... IDENTIFY WHERE TO COADD
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

        non_zero = where(coverage_plane ne 0, non_zero_ct)

        if non_zero_ct eq 0 then continue

        median_coverage = median(coverage_plane[non_zero])
        coadd_here = coverage_plane / median_coverage ge thresh_coverage
        where_to_coadd = where(coadd_here, coadd_ct)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... COADD
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

        coadd_plane = coadd[*,*,k]
        
        if nf eq 1 then begin
           coadd_plane[where_to_coadd] = data_stack[where_to_coadd]
        endif
        
        if coadd_ct gt 0 and nf gt 1 then begin
           coadd_plane[where_to_coadd] = $
              coadd_data(data_stack[where_to_coadd,*] $
                         , weight = coverage_stack[where_to_coadd,*] $
                         , use_median = use_median)
        endif        

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; ... STORE RESULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

        coadd[*,*,k] = coadd_plane
        coverage[*,*,k] = coverage_plane
     endfor
  endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SET ZERO-COVERAGE DATA TO NAN
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  zero_ind = where(coverage eq 0, zero_ct)
  if zero_ct gt 0 then $
     coadd[zero_ind] = !values.f_nan

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; WRITE THE OUTPUTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if keyword_set(use_disk) then $
     spawn, 'rm -rf '+ scratch_dir     

  writefits, out_root+'.fits', coadd, target_hdr
  sxaddpar, target_hdr, 'BUNIT','COVERAGE'
  writefits, out_root+'.coverage.fits', coverage, target_hdr

end
