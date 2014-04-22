pro add_fake_source $
   , infile = infile $
   , outfile = outfile $
   , source_xctr = x_ctr $
   , source_yctr = y_ctr $
   , source_vctr = v_ctr $
   , source_fwhm = fwhm $
   , source_vfwhm = v_fwhm $
   , source_peak = i_peak $
   , source_flux = flux $
   , efficiency = eff

;+
; NAME:
;
; add_fake_source
;
; PURPOSE:
;
; Add a fake gaussian signal to the input data.
;
; CATEGORY:
;
; Reduction tool.
;
; CALLING SEQUENCE:
;
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
;
; MODIFICATION HISTORY:
;
;
; IF YOU USE THIS AND FIND ERRORS:
;
; email to leroy@mpia.de
; 
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS, ERROR-CHECKING, AND DEFINITIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; REQUIRE AN INPUT AND OUTPUT FILE 
  if n_elements(outfile) eq 0 or $
     n_elements(infile) eq 0 then begin
     message, 'Need both INFILE and OUTFILE to work properly.'
  endif
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE VECTORS OF VELOCITY, RA, AND DEC TO GO WITH THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  data = mrdfits(infile,1,hdr)

; MEASURE THE SIZE OF THE DATA
  sz = size(data)
  
; BUILD SKY COORDINATES (COMPLICATED BY CONVENTION CHANGE)
  if sxpar(hdr,'CRVAL3') eq 0.0 then begin
     dec = double(data.crval3 + data.cdelt3)
  endif else begin
     dec = double(sxpar(hdr,'CRVAL3') + data.cdelt3)
  endelse
  
  median_dec = median(dec)
  
  if sxpar(hdr,'CRVAL2') eq 0.0 then begin
     ra = double(data.crval2) + $
          double(data.cdelt2/cos(!dtor*median_dec))
  endif else begin
     ra = double(sxpar(hdr,'CRVAL2')) + $
          double(data.cdelt2/cos(!dtor*median_dec))
  endelse  
; ... UNITS OF RA AND DEC SHOULD NOW BE DECIMAL DEGREES

; PARSE THE HEADER TO CALCULATE THE VELOCITY AXIS
  crval = sxpar(hdr,'VELO-LSR')
  crpix = sxpar(hdr,'CRPIX1')
  cdelt = sxpar(hdr,'DELTAV')
  
  v = findgen(n_elements(data[0].spectrum))
  vdif = v - (crpix-1.0)
  vaxis = vdif * cdelt + crval
  if abs(cdelt) gt 100. then $
     vaxis /= 1e3
; ... UNITS OF VAXIS SHOULD NOW BE V_LSR, KM/S

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DETAILS OF THE FAKE SOURCE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;
; Always a 3D Gaussian with a single FWHM in the spatial direction (i.e., it's
; round). Defaults to 11" FWHM (about the IRAM 30m beam at 230 GHz).
;

; FAKE SOURCE DEFAULTS
  if n_elements(x_ctr) eq 0 then $
     x_ctr = mean(ra)

  if n_elements(y_ctr) eq 0 then $
     y_ctr = mean(dec)

  if n_elements(v_ctr) eq 0 then $
     v_ctr = mean(vaxis)

; DEFAULT TO AN 13 ARCSECOND SOURCE ...
  if n_elements(fwhm) eq 0 then $
     fwhm = 13./3600.

; ... WITH 20 KM/S LINE WIDTH
  if n_elements(v_fwhm) eq 0 then $
     v_fwhm = 20.

; ... AND THE 30M 230GHZ EFFICIENCY
  if n_elements(eff) eq 0 then $
     eff = 0.52/0.91
  
; ... AND THE PEAK INTENSITY TO 1 K UNLESS
  if n_elements(i_peak) eq 0 and n_elements(flux) eq 0 then $
     i_peak = 1.

; WORK OUT THE PEAK INTENSITY IN CASE THE FLUX IS SPECIFIED
; (FLUX IS TAKEN TO BE IN K KM/S ARCSEC^2)
  if n_elements(flux) gt 0 and n_elements(i_peak) eq 0 then begin
     v_area_kms = sqrt(2.*!pi)*v_fwhm / sqrt(8.*alog(2))
     sky_area_as2 = (sqrt(2.*!pi)*(fwhm*3600.) / sqrt(8.*alog(2)))^2
     i_peak = flux / aky_area_as2 / v_area_kms
  endif
                 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ADD THE APPROPRIATE LINE TO THESE SPECTRA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for j = 0, n_elements(x_ctr)-1 do begin
     
;    FIND ALL DATA WITHIN 3 FWHM OF THE FAKE SOURCE CENTER
     dist = sphdist(ra, dec, x_ctr[j], y_ctr[j], /deg)
     ind = where(dist lt 3.0 * fwhm[j], ct)

;    BAIL IF THE SOURCE IS OUTSIDE THE MAP
     if ct eq 0 then begin
        message, 'No data found near specified fake source position.', /info
        continue
     endif  

;    THE SPATIAL SCALING FACTOR FOR EACH SPECTRUM
     spatial_fac = exp(-1.0*dist[ind]^2/2./(fwhm[j]/2.354)^2)

;    THE VELOCITY SCALING FACTOR FOR EACH CHANNEL
     vel_fac = exp(-1.0*(vaxis - v_ctr[j])^2/2./(v_fwhm[j]/2.354)^2)
     
;    ADD THE LINE TO EACH AFFECTED SPECTRUM
     for i = 0L, ct-1 do begin
        spec = data[ind[i]].spectrum
        spec += spatial_fac[i] * vel_fac * i_peak[j] * eff
        data[ind[i]].spectrum = spec
     endfor
     
;    MAKE A NOTE FOR THE HEADER
     sxaddpar, hdr, 'HISTORY', 'ADD_FAKE_SOURCE (IDL) has added a fake source.'
     sxaddpar, hdr, 'HISTORY', '... centered at (ra, dec, vel): '+$
               string(x_ctr[j]) + ', ' + $
               string(y_ctr[j]) + ', ' + $
               string(v_ctr[j])
     sxaddpar, hdr, 'HISTORY', '... with FWHM (sky, vel): ' + $
               string(fwhm[j]) + ', ' + $
               string(v_fwhm[j])
     sxaddpar, hdr, 'HISTORY', '... and peak intensity (main beam): ' + $
               string(i_peak[j])
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO THE OUTPUT FILE (DELETING FIRST TO AVOID APPENDING)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  spawn, 'rm '+outfile
  mwrfits, data, outfile, hdr  

  return

end                             ; of add_fake_source
