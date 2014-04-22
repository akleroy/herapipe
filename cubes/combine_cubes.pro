pro combine_cubes $
   , cube_list = cube_list $
   , noise_cube_list = noise_cube_list $
   , noice_vec = noise_vec $
   , out_file = out_file $
   , wt_out_file = wt_out_file
 
  n_cubes = n_elements(cube_list)

  cube = readfits(cube_list[0], hdr)
  wt_cube = (1./readfits(noise_cube_list[0],wt_hdr))^2

  for i = 1, n_cubes-1 do begin
     to_add = readfits(cube_list[i], add_hdr)
     wt_to_add = (1./readfits(noise_cube_list[i],wt_add_hdr))^2

     to_add = $
        cube_hastrom(cube_in=to_add, hdr_in=add_hdr, target=hdr)

     wt_to_add = $
        cube_hastrom(cube_in=wt_to_add, hdr_in=wt_add_hdr, target=wt_hdr)
     
     only_new = where((wt_to_add gt 0 and finite(to_add)) and $
                      (finite(cube) eq 0 or wt_cube eq 0), only_new_ct)
     
     overlap =  where((wt_to_add gt 0 and finite(to_add)) and $
                      (finite(cube) and wt_cube gt 0), overlap_ct)

;    ONLY AREA COVERED BY THE NEW CUBE
     if only_new_ct gt 0 then begin
        cube[only_new] = to_add[only_new]
        wt_cube[only_new] = wt_to_add[only_new]
     endif

;    OVERLAP BETWEEN NEW AND OLD CUBE
     if overlap_ct gt 0 then begin
        cube[overlap] = (wt_to_add[overlap]*to_add[overlap] + $
                         wt_cube[overlap]*cube[overlap]) / $
                        (wt_to_add[overlap] + wt_cube[overlap])
        wt_cube[overlap] = wt_to_add[overlap] + wt_cube[overlap]
     endif

  endfor

  writefits, out_file, cube, hdr

  if n_elements(wt_out_file) gt 0 then $
     write_fits, wt_out_file, wt_cube, hdr

end
