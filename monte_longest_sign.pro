pro monte_longest_sign, hans=han_size

  niter=100000L
  nspec = 465

  longest_1 = intarr(niter)
  longest_2 = intarr(niter)
  longest_3 = intarr(niter)

  kern = [0.25, 1., 0.25]
  kern /= total(kern)

  for i = 0L, niter-1 do begin
     spec = randomn(seed,nspec)
     spec = convol(spec, kern, /edge_truncate)

     if n_elements(han_size) gt 0 then $
        spec = hans(han_size,spec)

     positive = label_region(spec gt 0)
     length_pos = intarr(max(positive))
     for j = 1, max(positive) do $
        length_pos[j-1] = total(positive eq j)
     
     negative = label_region(spec lt 0)
     length_neg = intarr(max(negative))
     for j = 1, max(negative) do $
        length_neg[j-1] = total(negative eq j)

     length = [length_pos, length_neg]
     length = length[sort(length)]
     length = reverse(length)

     longest_1[i] = length[0]
     longest_2[i] = length[1]
     longest_3[i] = length[2]
  endfor

  !p.multi=[0,2,1]

  x = findgen(50)
  y = fltarr(50)
  for i = 0, 49 do y[i] = total(longest_1 lt x[i])/(niter*1.0)

  fasthist, longest_1, /log
  plot, x, y, xtitle='longest string', ytitle='cumulative distribution'
  
  print, total(longest_1 ge 20.)/(niter*1.0)

  !p.multi=[0,2,1]

  x = findgen(50)
  y = fltarr(50)
  for i = 0, 49 do y[i] = total(longest_2 lt x[i])/(niter*1.0)

  fasthist, longest_2, /log
  plot, x, y, xtitle='longest string', ytitle='cumulative distribution'

  print, total(longest_2 ge 15.)/(niter*1.0)

  !p.multi=[0,2,1]

  x = findgen(50)
  y = fltarr(50)
  for i = 0, 49 do y[i] = total(longest_3 lt x[i])/(niter*1.0)

  fasthist, longest_3, /log
  plot, x, y, xtitle='longest string', ytitle='cumulative distribution'

  print, total(longest_3 ge 13.)/(niter*1.0)

  stop

end
