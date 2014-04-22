pro compare_cube_stack $
   , cube_list $
   , ref_cube = ref_cube_file $
   , bright_map = bright_map_file $
   , window_mask = window_mask_file $
   , psfile = psfile

; Placeholder.
;
;Program needs to take a stack of cubes, a "bright region map," and a window
;cube and compare the cubes to the spectrum of a reference cube or the average
;spectrum. Running some basic statistics makes sense, too. Should be fairly
;straightforward.
;
;

  mask = readfits(window_mask_file, /silent)
  sz = size(mask)

; APPLY THE BRIGHT-REGION MAP TO THE MASK
  if n_elements(bright_map_file) gt 0 then begin
     bright_map = readfits(bright_map_file, /silent)
     for k = 0, sz[3]-1 do $
        mask[*,*,k] = bright_map*mask[*,*,k]
  endif
  weight = total(total(mask,1),1)

  spec_ra = fltarr(sz[3], n_elements(cube_list))

; GRAB A SPECTRUM FROM EACH CUBE
  for i = 0, n_elements(cube_list)-1 do begin
     cube = readfits(cube_list[i], hdr)
     spec = total(total(cube*mask,1,/nan),1) / weight
     spec *= (weight gt 0.25*max(weight))
     spec_ra[*,i] = spec
  endfor

  make_axes, hdr, vaxis=vaxis
  vaxis /= 1e3

; GET THE REFERENCE SPECTRUM
  if n_elements(ref_cube_file) gt 0 then begin
     cube = readfits(ref_cube_file)
     ref_spec = total(total(cube*mask,1,/nan),1) / weight
     ref_spec *= (weight gt 0.25*max(weight))
  endif else begin
     ref_spec = total(spec_ra,2) / total(finite(spec_ra),2)
  endelse

; PLOTTING RANGE
;  plot_max = 1.1*(max(spec_ra,/nan) > max(ref_spec,/nan))
;  plot_min = (min(spec_ra,/nan) < min(ref_spec,/nan)) < 0. 
;  plot_min -= 0.1*abs(plot_min)

; INITIALIZE THE PS PLOTTER
  if n_elements(psfile) gt 0 then begin
     ps, /ps, /def, file=psfile, xsize=7, ysize=10, /color
  endif
  
  !p.multi=[0,1,2]
  
  for k = 0, n_elements(cube_list)-1 do begin
     plot, vaxis, spec_ra[*,k] $
           ;, yrange=[plot_min, plot_max], /ystyle $
           , /xstyle $
           , title = '!6'+cube_list[k] $
           , xtitle = '!6LSR Velocity (km s!u-1!n)' $
           , ytitle = '!6Intensity'
     oplot, vaxis, ref_spec, color=getcolor('red')     
  endfor

; CLOSE THE PS FILE THEN SHOW THE RESULTS
  if n_elements(psfile) gt 0 then begin
     ps, /x
     spawn, 'gv '+psfile+' &'
  endif
  
end                             ; of compare_cube_stack
