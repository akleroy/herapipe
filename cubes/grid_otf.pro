pro grid_otf $
   , data = data $
   , ra = ra $
   , dec = dec $
   , weight = weight $
   , target_hdr = target_hdr $
   , out_root = out_root $
   , show=show $
   , beam_fwhm = beam_fwhm $
   , gauss_kern = gauss_kern $
   , gauss_fwhm = gauss_fwhm $
   , median = use_median $
   , apodize = apodize $
   , apod_fwhm = apod_fwhm

;+
; NAME:
;
; grid_otf
;
; PURPOSE:
;
; Grid individual spectra onto a specified regular grid following the
; recommendations of Mangum et al. (2007). This is written to be of general
; use gridding OTF data (i.e., not telescope specific). To apply this to a
; .fits binary table from IRAM use the IRAM_TABLE_TO_CUBE wrapper.
;
; CATEGORY:
;
; Reduction tool.
;
; CALLING SEQUENCE:
;
; grid_otf, data = data $
;         , ra = ra, dec = dec, weight = weight $
;         , target_hdr = target_hdr $
;         , out_root = 'output' $
;         , beam_fwhm = beam_fwhm $
;         , /gauss_kern $
;         , gauss_fwhm = gauss_fwhm $
;         , /median $         
;         , /show
;
; INPUTS:
;
; The spectra to be gridded and a header defining the target grid. This is
; provided as vectors:
;
; data - an nspec x nchan size two-dimensional array containing the spectra to
;        be gridded.
; ra   - an nspec long array of right ascencions in degrees
; dec  - an nspec long array of declinations in degrees
;
; also needs a target astrometry:
;
; target_hdr - a .fits header (string array) containing the target astromtery
;
; OPTIONAL INPUTS:
;
; weight - an nspec long vector containing weights to be applied to the
;          individual spectra. If not supplied, the porgram assumes equal
;          weights for all data.
;
; beam_fwhm - the fwhm in decimal degrees of the telescope beam. If not
;             supplied, then the program looks for this in the
;             target_header. If it isn't in either place, the program crashes. 
;
; gauss_fwhm - the fwhm in decimal degrees of the gaussian beam to use as a
;              convolution kernel. Requires /GAUSS_KERN to be set.
;
; apod_fwhm - the fwhm in decimal degrees of the gaussian apodization kernel
;             used to fill in missing data and smooth edges.
;
;
; KEYWORD PARAMETERS:
;
; gauss_kern - use a Gaussian convolution kernel. Defaults to FWHM = 1/3
;              BEAM_FWHM, the default convolution kernel used by CLASS.
;
; apodize - as a final step, replace missing data that are near valid
;           measurements with values obtained by smoothing the cube. Default
;           FWHM is 30" but this can be adjusted with the APOD_FWHM parameter.
;
; [ WARNING! Apodization still a work in progress]
;
; median - use a (weighted) median rather than traditional averaging when
;          constructing the cube. Massively slower to compute but potentially
;          reduces sensitivity to outliers at a roughly 20% cost in
;          signal-to-noise. In limited testing, apparent benefits don't
;          outweigh the two OOM increase in computation time. 
;
; OUTPUTS:
;
; out_root+".fits" - a .fits file containing the synthesized cube
;
; out_root+".coverage.fits" - a .fits file containing the coverage of the data
;                             at each pixel
;
; OPTIONAL OUTPUTS:
;
; none
;
; COMMON BLOCKS:
;
; none
;
; SIDE EFFECTS:
;
; none
;
; RESTRICTIONS:
;
; none
;
; PROCEDURE:
;
; Takes the irregularly spaced data and moves them onto the regular .fits grid
; defined by the target header. The procedure involves a convolution. The
; default kernel is a Bessel function tapered by a Gaussian the beam
; width. This is a compromise between practicality and correctness: it
; preserves most spatial scales while being computable in a reasonable amount
; of time. Alternatively, the user can ask the program to use a Gaussian and
; specify a size (default FWHM 1/3 the beam size).
;
; Credit to CLASS's xy_map, E. Rosolowsky, Mangum+ 07.
;
; MODIFICATION HISTORY:
;
; 17dec08 - written following methodology from E. Rosolowsky & Mangum+ 07
;           leroy@mpia.de
;
; 13apr09 - tried my best to speed things up - leroy@mpia.de
;
; 03jun09 - updated documentation and adjusted the final beam size
;           calculations - leroy@mpia.de
;
; 20jul09 - polished, added apodization, updated documentation.
;
; spring14 - some cleanup on pixel size calculation
;
; IF YOU USE THIS AND FIND ERRORS:
;
; email to aleroy@nrao.edu (no promises, but I want to know)
; 
;-

;  on_error, 2

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS, ERROR-CHECKING, AND DEFINITIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
 
; NEED INPUT
  if ((n_elements(data) eq 0) or $
      (n_elements(ra) eq 0) or $
      (n_elements(dec) eq 0)) then $
         message, 'Requires data, ra, and dec vectors.'

; NEED A TARGET HEADER
  if (n_elements(target_hdr) eq 0) then $
     message, 'Requires a target header (TARGET_HDR = ).'
  
; NEED AN OUTPUT FILE 
  if (n_elements(out_root) eq 0) then $
     message, 'Requires an output file root  (OUT_ROOT = ).'

; NEED THE BEAM SIZE
  if (n_elements(beam_fwhm) eq 0) then begin
     beam_fwhm = sxpar(target_hdr,'BMAJ',count=ct)     
     if ct eq 0 then $
        message, 'Requires the beam size via keyword or header.'     
  endif

; THE SIZE OF A PIXEL
  make_axes, target_hdr, ri=ri, di=di  
  pix_scl = sphdist(ri[1,0], di[1,0], ri[0,0], di[0,0], /deg)  

; THE CELL SIZE SHOULD BE 1/3 OF THE FWHM OR SMALLER (WARNING ONLY)
  if abs(sxpar(target_hdr,'CDELT1')) gt beam_fwhm/3.0 then $
     message, 'Warning! Cell size seems to be greater than 1/3 beam.', /info

; APODIZATION KERNEL SIZE
  if n_elements(apod_fwhm) eq 0 then begin
     apod_fwhm = 30./3600.     
  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFINE THE CONVOLVING FUNCTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(gauss_kern) then begin
     if n_elements(gauss_fwhm) eq 0 then $
        gauss_fwhm = beam_fwhm / 3.0
     r_support = 3.0 * gauss_fwhm
     max_conv_fn = 0.5
     cutoff_conv_fn = 0.25 * max_conv_fn
     scale_fwhm = sqrt(beam_fwhm^2 + gauss_fwhm^2)/beam_fwhm
  endif else begin
;    MANGUM, EMERSON, AND GREISEN (2007) ... MODULO TYPOS
     a = 1.55 * (beam_fwhm) / 3.0
     b = 2.52 * (beam_fwhm) / 3.0
     r_support = beam_fwhm
     max_conv_fn = 0.5
     cutoff_conv_fn = 0.25 * max_conv_fn
     scale_fwhm = 1.09          ; APPROXIMATE
  endelse
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE THE OUTPUT AND ASTROMETRY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; GENERATE ASTROMETRY FROM THE TARGET HEADER
  make_axes, target_hdr, ri = ri, di = di, vaxis=vaxis, astrom=astr

; CHECK THAT VAXIS AND SPECTRA HAVE THE SAME LENGTH
  nspec = (size(data))[1]
  nchan = (size(data))[2]
  if (n_elements(vaxis) ne nchan) then begin
     message, 'Velocity info in header and spectra size do not match.', /info
     stop
     return
  endif
  
; FIGURE THE SIZES OF THE OUTPUT
  nx = (size(ri))[1]
  ny = (size(ri))[2]

; IF THERE IS NO WEIGHTING ARRAY, MAKE A DEFAULT (ALL 1s) VERSION
  if n_elements(weight) eq 0 then begin
     weight = fltarr(nspec) + 1.0
  endif

; MAKE TWO EMPTY CUBES: DATA AND COVERAGE
  data_cube = !values.f_nan*fltarr(nx, ny, nchan)
  coverage_cube = fltarr(nx,ny,nchan) ; empty is 0 and not nan
  nspec_map = intarr(nx,ny)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SOME PRECALCULATION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; THE SIZE OF A PIXEL
  pix_scl = abs(sxpar(target_hdr, 'CDELT2'))

; CALCULATE THE X AND Y COORDINATE FOR EACH SPECTRUM
  ad2xy, ra, dec, astr, x_pix, y_pix

; THE SUPPORT RADIUS SQUARED
  r_support_sqrd = r_support^2

; THE SUPPORT RADIUS IN PIXELS
  r_support_pix = r_support / pix_scl
  r_support_pix_sqrd = float(r_support_pix^2)

; THE CONVOLUTION FUNCTION AS A FUNCTION OF DISTANCE SQUARED
  n_pre = 10000
  pre_delta = (1./n_pre*r_support^2)
  pre_delta_pix_sqrd = pre_delta / (pix_scl^2)
  pre_dist_sqrd = findgen(n_pre*1.01)*pre_delta
  pre_dist_sqrd_pix = pre_dist_sqrd / (pix_scl)^2
  pre_dist = sqrt(pre_dist_sqrd)

  if keyword_set(gauss_kern) then begin
     pre_conv_fn = max_conv_fn * $
                   exp(-0.5 * (pre_dist/(gauss_fwhm/2.354))^2)
  endif else begin 
     pre_conv_fn = $
        beselj(!pi*pre_dist/a, 1) / $
        (!pi * pre_dist/a) * exp(-1.0 * (pre_dist/b)^2)
  endelse

; DISTANCE INSIDE WHICH WE CAP THE CONVOLUTION FUNCTION
  cap_dist_sqrd = pre_dist_sqrd[1]
  cap_dist_sqrd_pix = pre_dist_sqrd_pix[1]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER THE SPATIAL AXES OF THE CUBE AND GRID
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; NOTE PRECOMPUTE ON THE R.A. TERM AND CONVOLUTION FUNCTION
; *** ASSUMES THAT SPHERICAL TERMS ARE NEGLIGIBLE OVER THE KERNEL ***        

  for i = 0, nx-1 do begin
     counter, i+1, nx, 'Row '
     for j = 0, ny-1 do begin        

        pix_dist_sqrd = $
           ((x_pix - i)*(x_pix - i) + (y_pix - j)*(y_pix - j))

;       EXTRACT PIXELS WITHIN THE RADIUS OF INTEREST OF THIS PIXEL        
        keep = where(pix_dist_sqrd le r_support_pix_sqrd, keep_ct)
        nspec_map[i,j] = keep_ct

        if keep_ct gt 1 then begin
           
;          INTERPOLATE TO GET THE CONV. FN. FOR EACH DATA POINT
           pre_grid_x = pix_dist_sqrd[keep] / pre_delta_pix_sqrd
           conv_fn = interpolate(pre_conv_fn, pre_grid_x)

;          DATA RIGHT ON TOP OF THE PIXEL TO THE MAX OF THE CONV. FN.
           ind = where(pix_dist_sqrd[keep] lt cap_dist_sqrd_pix, ct)
           if ct gt 0 then conv_fn[ind] = max_conv_fn

;          PLACE A MINIMUM THRESHOLD NEEDED TO CONSIDER A GRID POINT
           coverage = total(conv_fn)
           if coverage gt cutoff_conv_fn then begin

;             GET A VECTOR OF NORMALIZED WEIGHTS
              combined_weight = conv_fn * weight[keep]
              total_weight = total(combined_weight)
              combined_weight /= total_weight

;             IF REQUESTED, PERFORM A "WEIGHTED MEDIAN" ...
              if keyword_set(use_median) then begin
                 spectrum = fltarr(nchan)
                 for k = 0, nchan-1 do begin
;                   GET THE DATA (NOT YET SORTED)
                    sorted_vals = data[keep,k]

;                   SORT IT AND PUT THE WEIGHTS IN THE SAME ORDER
                    order = sort(sorted_vals)
                    sorted_weights = combined_weight[order]
                    sorted_vals = sorted_vals[order]

;                   CALCULATE THE CUMULATIVE DISTRIBUTION FUNCTION
                    cumul = total(sorted_weights, /cumulative) - $
                            sorted_weights*0.5

;                   INTERPOLATE TO THE MIDDLE OF THE CDF
                    spectrum[k] = interpol(sorted_vals,cumul,0.5)
                 endfor
;             ... ELSE JUST GET THE AVERAGE SPECTRUM
              endif else $
                 spectrum = reform(data[keep,*] ## combined_weight)

;             UPDATE THE DATA AND WEIGHTING CUBES
              data_cube[i,j,*] = spectrum
              coverage_cube[i,j,*] = total_weight ; was "coverage"
           endif
        endif

;       HANDLE THE CASE OF ONLY ONE DATA POINT
        if keep_ct eq 1 then begin
           conv_fn = (pix_dist_sqrd[keep] lt cap_dist_sqrd_pix) ? $
                     max_conv_fn : $
                     interpolate(pre_conv_fn $
                                 , pix_dist_sqrd[keep] / pre_delta_pix_sqrd)
           if (conv_fn gt cutoff_conv_fn) then begin
              data_cube[i,j,*] = $
                 data[keep,*]
              coverage_cube[i,j,*] = (conv_fn*weight[keep])[0]
           endif
        endif        

     endfor
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; UPDATE HEADER AND WRITE THE DATA TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  cube_file = out_root + '.fits'
  coverage_file = out_root + '.coverage.fits'
  nspec_file = out_root + '.nspec.fits'

; INCREASE THE BEAM SIZE TO REFLECT THE CONVOLUTION
  beam_fwhm = scale_fwhm*beam_fwhm

; GENERATE THE OUTPUT HEADER
  hdr_out = target_hdr
  sxaddpar, hdr_out, 'HISTORY', 'GRID_OTF (IDL): converted spectra to cube.'
  if keyword_set(gauss_kern) then $
     sxaddpar, hdr_out, 'HISTORY', 'Convolved with Gaussian convolution function.' $
  else $
     sxaddpar, hdr_out, 'HISTORY', 'Convolved with optimized kernel from Mangum+07.'

  if keyword_set(apodize) then $
     sxaddpar, hdr_out, 'HISTORY', 'Missing data and edges apodized.'


  if keyword_set(gauss_kern) then begin
     sxaddpar, hdr_out, 'BMAJ', beam_fwhm
     sxaddpar, hdr_out, 'BMIN', beam_fwhm     
  endif else begin
     sxaddpar, hdr_out, 'BMAJ', beam_fwhm, '*But* not Gaussian.'
     sxaddpar, hdr_out, 'BMIN', beam_fwhm, '*But* not Gaussian.'
  endelse

; WRITE THE DATA TO DISK
  writefits, cube_file, data_cube, hdr_out

; WRITE THE COVERAGE TO DISK
  sxaddpar, hdr_out, 'BUNIT', 'COVERAGE' ; CHANGE FROM K -> COVERAGE
  writefits, coverage_file, coverage_cube, hdr_out

; WRITE THE NUMBER OF SPECTRA COADDED AT EACH POINT TO DISK
  hdr_out = twod_head(hdr_out)
  sxaddpar, hdr_out, 'BUNIT', 'SPECTRA' $
            , 'NUMBER OF SPECTRA' ; CHANGE FROM K -> # SPEC
  writefits, nspec_file, nspec_map, hdr_out

end                             ; of grid_otf
