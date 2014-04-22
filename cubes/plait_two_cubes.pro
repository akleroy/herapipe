pro plait_two_cubes $
   , file_one = file_one $
   , angle_one = angle_one $
   , file_two = file_two $
   , angle_two = angle_two $
   , out_file = out_file $
   , mask_file = mask_file $
   , test_accuracy = test_accuracy $
   , show = show 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ASSIGN ONE CUBE TO BE VERTICAL (ANGLE 0) AND ONE HORIZONTAL (90)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  angle_diff = abs(angle_one - angle_two)
  if angle_diff ne 90 then begin
     message, 'Maps not 90 degrees out of phase.'    
  endif
  
; PHASE SHIFT NEGATIVE ANGLES BY 180 DEGREES
  if (angle_one) lt 0 then angle_one += 180
  if (angle_two) lt 0 then angle_two += 180  
  
; FIND THE ONE THAT'S CLOSER TO VERTICAL (I.E., 0 DEGREES)
  diff_from_zero_one = $
     min([abs(angle_one),abs(angle_one - 180)])
  diff_from_zero_two = $
     min([abs(angle_two),abs(angle_two - 180)])
  
; IDENTIFY A HORIZONTAL AND A VERTICAL CUBE AND WORK OUT HOW MUCH ROTATION, IF
; ANY, NEEDS TO BE APPLIED TO THE DATA
  if diff_from_zero_one le diff_from_zero_two then begin
     vert_file = file_one
     hori_file = file_two
     if (angle_one lt 90) then $
        rot = angle_one $
     else $
        rot = angle_one - 180.
  endif else begin
     vert_file = file_two
     hori_file = file_one
     if (angle_two lt 90) then $        
        rot = angle_two $
     else $
        rot = angle_two - 180.
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  cube_hori = readfits(hori_file, hori_hdr)
  twod_hori_hdr = twod_head(hori_hdr)

  cube_vert = readfits(vert_file, vert_hdr)
  twod_vert_hdr = twod_head(vert_hdr)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THEM ZERO OUTSIDE THE MASK REGION AND WHERE THEY ARE N-A-N
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(mask_file) eq 0 then begin
     have_mask = 0
  endif else begin
     test = file_search(mask_file, count=count)
     have_mask = count gt 0
  endelse
  
  if have_mask then begin
     mask = readfits(mask_file)     
  endif else begin
     mask = (cube_hori eq cube_hori)*0B + 1B
  endelse


  blank_hori = where(finite(cube_hori) eq 0 or (mask eq 0), blank_hori_ct)
  blank_vert = where(finite(cube_vert) eq 0 or (mask eq 0), blank_vert_ct)
  
  if blank_hori_ct gt 0 then cube_hori[blank_hori] = 0.
  if blank_vert_ct gt 0 then cube_vert[blank_vert] = 0.

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ROTATE THE CUBE (IF NEEDED)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if abs(rot) gt 1. then begin
     orig_cube_hori = cube_hori
     orig_cube_vert = cube_vert

     sz = size(cube_hori)
     for i = 0, sz[3]-1 do begin
        hrot, cube_hori[*,*,i], twod_hori_hdr, plane, new_hdr $
              , rot, -1, -1, 2, MISSING = 0. $
              , CUBIC = -0.5
        cube_hori[*,*,i] = plane

        hrot, cube_vert[*,*,i], twod_vert_hdr, plane, new_hdr $
              , rot, -1, -1, 2, MISSING = 0. $
              , CUBIC = -0.5
        cube_vert[*,*,i] = plane

     endfor
  endif 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE THE FFT-PLANE WEIGHTING ARRAYS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; THIS NEEDS IMPROVEMENT. FOR THE MOMENT, WE USE 1e-4 RADIAN, OR ABOUT 21"
  
  sz = size(cube_vert)

  pix_radian = abs(sxpar(vert_hdr, 'CDELT1')*!dtor)
  
; AN IMAGE OF Y-AXIS SPATIAL FREQUENCIES IN UNITS OF radian^-1
  y_img = (fltarr(sz[1])+1.) # findgen(sz[2])
  y_img = (y_img) < abs(sz[2] - y_img)
  y_img /= (sz[2] * pix_radian) 
  
; AN IMAGE OF X-AXIS SPATIAL FREQUENCIES IN UNITS OF radian^-1
  x_img = findgen(sz[1]) # (fltarr(sz[2])+1.)
  x_img = (x_img) < abs(sz[1] - x_img)
  x_img /= (sz[1] * pix_radian)

; SCALE OF FREQUENCIES TO PRUNE
  scale   = 1e-4                ;1./3600.            ; radian
  vert_wt = fltarr(sz[1],sz[1])+1.0
  vert_wt = (1.0 - (cos(y_img*scale*!pi/2.))^2)*(abs(y_img) lt 1./scale) + $
            (abs(y_img) ge 1./scale)*1.0
  
  hori_wt = fltarr(sz[1],sz[1])+1.0
  hori_wt = (1.0 - (cos(x_img*scale*!pi/2.))^2)*(abs(x_img) lt 1./scale) + $
            (abs(x_img) ge 1./scale)*1.0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE AN OUTPUT CUBE AND THEN GO PLANE BY PLANE AND PLAIT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; COPY AND BLANK ONE OF THE CUBES TO BE THE OUTPUT CUBE
  combined_cube = cube_hori*0.
  sz = size(combined_cube)

  for i = 0, sz[3]-1 do begin
     
     plane_hori = cube_hori[*,*,i]
     plane_vert = cube_vert[*,*,i]
     
     fft_hori = fft(plane_hori)
     fft_vert = fft(plane_vert)

     fft_avg = (fft_hori*hori_wt + fft_vert*vert_wt) / (hori_wt + vert_wt)
     fft_avg[0,0] = (fft_vert[0,0] + fft_hori[0,0])*0.5

     combined_cube[*,*,i] = fft(fft_avg,1)
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEROTATE ROTATE THE CUBE (IF NEEDED)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if abs(rot) gt 1. then begin
     sz = size(combined_cube)
     for i = 0, sz[3]-1 do begin
        hrot, combined_cube[*,*,i], twod_hori_hdr, plane, new_hdr $
              , -1.*rot, -1, -1, 2, MISSING = !values.f_nan $
              , CUBIC = -0.5
        combined_cube[*,*,i] = plane        
     endfor
  endif

  if blank_hori_ct gt 0 then combined_cube[blank_hori] = !values.f_nan
  if blank_vert_ct gt 0 then combined_cube[blank_vert] = !values.f_nan

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; OPTIONALLY TEST THE IMPACT OF ROTATION/INTERPOLATION ON THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(test_accuracy) and abs(rot) gt 1. then begin
     sz = size(cube_hori)
     for i = 0, sz[3]-1 do begin
        hrot, cube_hori[*,*,i], twod_hori_hdr, plane, new_hdr $
              , -1.*rot, -1, -1, 2, MISSING = !values.f_nan $
              , CUBIC = -0.5
        cube_hori[*,*,i] = plane
        
        hrot, cube_vert[*,*,i], twod_vert_hdr, plane, new_hdr $
              , -1.*rot, -1, -1, 2, MISSING = !values.f_nan $
              , CUBIC = -0.5
        cube_vert[*,*,i] = plane
     endfor

     if blank_hori_ct gt 0 then cube_hori[blank_hori] = !values.f_nan
     if blank_vert_ct gt 0 then cube_hori[blank_vert] = !values.f_nan
     if blank_hori_ct gt 0 then cube_vert[blank_hori] = !values.f_nan
     if blank_vert_ct gt 0 then cube_vert[blank_vert] = !values.f_nan

     ind = where(finite(cube_hori) and finite(orig_cube_hori))
     message, 'Rotation introduces '+ $
              string(mad((cube_hori/orig_cube_hori)[ind]))+ $
              ' scatter about '+ $
              string(median((cube_hori/orig_cube_hori)[ind])), /info
     message, 'Rotation introduces '+ $
              string(mad((cube_vert/orig_cube_vert)[ind]))+ $
              ' scatter about '+ $
              string(median((cube_vert/orig_cube_vert)[ind])), /info
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  writefits, out_file, combined_cube, hori_hdr

end                             ; of plait_two_cubes
