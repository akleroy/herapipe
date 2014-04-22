function make_iram_cube_header, tab_hdr = tab_hdr $
                                , pix_scale = pix_scale $
                                , xsize=xsize, ysize=ysize $
                                , xctr=xctr, yctr=yctr

  naxis_vec = [xsize, ysize, sxpar(tab_hdr,'MAXIS1')]
  mkhdr, hdr, 3, naxis_vec

; MAKE THE VELOCITY AXIS (ALONG THE THIRD DIMENSION)
  sxaddpar, hdr, 'CTYPE3', 'VELOCITY', 'M/S'
  sxaddpar, hdr, 'CRVAL3', sxpar(tab_hdr,'VELO-LSR')
  sxaddpar, hdr, 'CRPIX3', sxpar(tab_hdr,'CRPIX1')
  sxaddpar, hdr, 'CDELT3', sxpar(tab_hdr,'DELTAV')

; MAKE THE POSITION AXES
  sxaddpar, hdr, 'CTYPE1', 'RA---TAN'
; ... ALLOW A USER-SUPPLIED R.A. CENTER, ELSE USE THE TABLE HEADER
  if n_elements(xctr) gt 0 then $
     sxaddpar, hdr, 'CRVAL1', xctr $
  else $
     sxaddpar, hdr, 'CRVAL1', sxpar(tab_hdr,'CRVAL2')
  sxaddpar, hdr, 'CRPIX1', xsize/2 
  sxaddpar, hdr, 'CDELT1', -1.0*pix_scale

; ... ALLOW A USER-SUPPLIED DECLINATION CENTER, ELSE USE THE TABLE HEADER
  sxaddpar, hdr, 'CTYPE2', 'DEC--TAN'
  if n_elements(xctr) gt 0 then $
     sxaddpar, hdr, 'CRVAL2', yctr $
  else $
     sxaddpar, hdr, 'CRVAL2', sxpar(tab_hdr,'CRVAL3')

; ADD THE J2000 EQUINOX
  sxaddpar, hdr, 'EQUINOX', 2000

  sxaddpar, hdr, 'CRPIX2', ysize/2
  sxaddpar, hdr, 'CDELT2', pix_scale

; PUT THE RESOLUTION OF THE 30M IN THE HEADER
  beam_fwhm = 2460. / (sxpar(tab_hdr,'RESTFREQ')/1d9) / 3600.
  sxaddpar, hdr, 'BMAJ', beam_fwhm
  sxaddpar, hdr, 'BMIN', beam_fwhm

; TRANSFER OTHER ENTRIES OF INTEREST
  sxaddpar, hdr, 'LINE', sxpar(tab_hdr,'LINE')

; ADD UNITS (TA* in K)
  sxaddpar, hdr, 'BUNIT', 'K', 'TA*, apply eff. to get Tmb'

; ADD A STAMP
  sxaddpar,hdr,'HISTORY','Header created by MAKE_IRAM_CUBE_HEADER.pro'
  spawn, 'date', date_str 
  sxaddpar,hdr,'HISTORY',date_str

  return, hdr

end                             ; of make_iram_cube_header
