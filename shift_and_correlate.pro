function shift_and_correlate $
   , data_in $
   , ind=ind $
   , ngrid=ngrid $
   , cube=cube

  data = data_in
  mask = (data eq data)*0B
  if n_elements(ind) eq 0 then ind = where(finite(data))
  mask[ind] = 1B
  bad = where(mask eq 0, badct)
  if badct gt 0 then data[bad] = !values.f_nan
  if n_elements(ngrid) eq 0 then ngrid = 5

  output = fltarr(2*ngrid+1, 2*ngrid+1)*!values.f_nan

  for i = -ngrid, ngrid do begin
     counter, i+ngrid, 2*ngrid+1, 'Row'
     for j = -ngrid, ngrid do begin        
        x = data[ind]
        y = (keyword_set(cube)) ? (shift(data,i,j,0))[ind] : $
            (shift(data,i,j))[ind]
        keep = where(finite(x) and finite(y))
        x = x[keep]
        y = y[keep]
        output[i+ngrid,j+ngrid] = correlate(x,y)
     endfor
  endfor

  return, output
end
