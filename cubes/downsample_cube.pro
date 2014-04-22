pro downsample_cube $
   , in_file = in_file $
   , out_file = out_file $
   , iterations = iterations

  cube = readfits(in_file, hdr)

  make_axes, hdr, vaxis=vaxis, /vonly

  v0 = vaxis[0]
  
  for k = 0, iterations-1 do begin

     sz = size(cube)
     new_cube = dblarr(sz[1],sz[2],sz[3]/2)

     for i = 0, sz[1]-1 do begin
        
        for j = 0, sz[2]-1 do begin
        
           spec = reform(cube[i,j,*])

           spec = hans(3, spec)

           new_cube[i,j,*] = spec[2*indgen(sz[3]/2)+1]

        endfor

     endfor

     cube = new_cube
     vaxis = vaxis[2*indgen(sz[3]/2)+1]

  endfor

  sxaddpar,hdr, 'CDELT3', vaxis[1]-vaxis[0], 'Updated by downsample_cube.pro'
  sxaddpar,hdr, 'CRPIX3', 1, 'Updated by downsample_cube.pro'
  sxaddpar,hdr, 'CRVAL3', vaxis[0], 'Updated by downsample_cube.pro'
  sxaddpar,hdr, 'NAXIS3', n_elements(vaxis), 'Updated by downsample_cube.pro'
  writefits, out_file, cube, hdr

end
