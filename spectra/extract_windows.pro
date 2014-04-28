function extract_windows, data

  max_nwindows = max(data.nwindows)

  if max_windows eq 0 then $
     return, !values.f_nan

  sz = size(data)

  win_arr = fltarr(max_nwindows*2,sz[1])*!values.f_nan

  for i = 0L, sz[1]-1 do begin
     nwindows = data[i].nwindows
     windows = fltarr(2*nwindows)
     reads, data[i].windows, windows
     win_arr[0:(2*nwindows-1),i] = windows
  endfor

  return, win_arr

end
