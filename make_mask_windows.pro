pro make_mask_windows $
;  INPUT
   , mask_file $
   , gname = gname $
;  OUTPUT FILE STEM
   , out_root = out_root $
;  OPTIONS RELATING TO ASTROMETRY AND UNITS
   , target_hdr = target_hdr $
   , hel_to_lsr = hel_to_lsr $
   , ms_to_kms = ms_to_kms $
;  OPTIONS TO PROCESS THE MASK
   , pad_vel = pad_vel $
   , exp_rad = exp_rad $
   , reject = reject_thresh $
;  TUNING PARAMETERS
   , fit_window = fit_window $
   , hard_min = hard_min $
   , hard_max = hard_max

;+
; NAME:
;
; make_mask_windows
;
; PURPOSE:
;
; Create several .fits files that indicate the location of the line and a
; surrounding "fitting region." This procedure uses a mask as a template.
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
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
; Four fits files: 
;
; out_root_noise_low.fits, out_root_line_low.fits, out_root_line_high.fits,
; out_root_noise_high.fits
;
; These give the locations of the upper and lower boundaries of the fitting
; region (the "noise" files) and the upper and lower boundary of the likely
; location of the line (the "line" files).
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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS, DEFINITIONS, ERROR-CHECKING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; PHYSICAL CONSTANTS BATCH FILE
  @constants.bat

; ROOT OF OUTPUT FILE NAMES
  if n_elements(out_root) eq 0 then $
     out_root = '../windows/'+gname+'_mask_window'
  
; READ THE VELOCITY FIELD
  mask = readfits(mask_file, hdr, /silent)

; GET INFO FROM MY DATABASE ON THIS GALAXY
  s = things_galaxies(gname)

; CONSTRUCT THE VELOCITY AXIS
  make_axes, hdr, vaxis=vaxis, /vonly

; CONVERT M/S TO KM/S
  if keyword_set(ms_to_kms) then begin
     vaxis /= 1000.
  endif

; CONVERT HELIOCENTRIC TO LSR VELOCITY
  if keyword_set(hel_to_lsr) then begin
  
     if n_elements(gname) eq 0 then begin
        message, 'Need galaxy name to convert helio -> lsr', /info
        return
     endif

;    CONVERT FROM HELIOCENTRIC VELOCITY TO LSR
     helio2lsr, 0., vlsr, ra = s.ra_deg, dec = s.dec_deg
     vfield += vlsr
  endif

; MOVE TO TARGET ASTROMETRY 
  if n_elements(target_hdr) gt 0 then begin
     mask = cube_hastrom(cube_in = mask $
                         , hdr_in = hdr $
                         , target_hdr = target_hdr)
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS THE MASK AS REQUESTED
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; REJECT REGIONS BELOW A CERTAIN SIZE
  if n_elements(reject_thresh) gt 0 then $
     reject_regions, mask, reject_thresh

; GROW THE MASK
  if n_elements(exp_rad) gt 0 then $
     mask = exp_mask(mask, rad=exp_rad)
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WORK OUT THE HIGHEST AND LOWEST VELOCITY AS A FUNCTION OF LOCATION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  sz = size(mask)

  mean_vel = mean(vaxis)
  line_high = fltarr(sz[1],sz[2])+mean_vel
  line_low = fltarr(sz[1],sz[2])+mean_vel

  for i = 0, sz[1]-1 do begin
     for j = 0, sz[2]-1 do begin

        spec = mask[i,j,*]        
        ind = where(spec eq 1, ct)
        if ct eq 0 then continue
                
        vel_here = vaxis[ind]
        line_high[i,j] = max(vel_here)
        line_low[i,j]  = min(vel_here)        
     endfor
  endfor

  if n_elements(pad_vel) ne 0 then begin
     line_low -= pad_vel
     line_high += pad_vel
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE CORRESPONDING NOISE WINDOWS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  

  noise_low = line_low - fit_window
  noise_high = line_high + fit_window

; ALLOW A HARD MINIMUM AND MAXIMUM
  if n_elements(hard_max) gt 0 then begin
     noise_high = noise_high < hard_max
     nan_ind = where(finite(noise_high) eq 0, nan_ct)
     if nan_ct gt 0 then noise_high[nan_ind] = hard_max
  endif

  if n_elements(hard_min) gt 0 then $
     noise_low = noise_low > hard_min  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
  
  sxaddpar, hdr, 'BUNIT', 'KM/S', 'LSR VELOCITY'

  writefits, out_root+'_noise_low.fits', noise_low, hdr
  writefits, out_root+'_line_low.fits', line_low, hdr
  writefits, out_root+'_line_high.fits', line_high, hdr
  writefits, out_root+'_noise_high.fits', noise_high, hdr

end                             ; of make_mask_windows
