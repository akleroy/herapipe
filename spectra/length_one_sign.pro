function length_one_sign, spec, longest=longest
  
  tot = total(spec gt 0)
  if tot eq 0 or tot eq n_elements(spec) then $
     return, n_elements(spec)    

  if (total(spec gt 0) gt 0) then begin     
     positive = label_region(spec gt 0)  
     if max(positive) gt 0 then begin
        length_pos = intarr(max(positive))    
        for i = 1, max(positive) do $
           length_pos[i-1] = total(positive eq i)
        length = length_pos
     endif
  endif

  if (total(spec lt 0) gt 0) then begin     
     negative = label_region(spec lt 0)
     if max(negative) gt 0 then begin
        length_neg = intarr(max(negative))
        for i = 1, max(negative) do $
           length_neg[i-1] = total(negative eq i)
        length = n_elements(length) gt 0 ? [length, length_neg] : length_neg
     endif
  endif
    
  length = length[sort(length)]
  length = reverse(length)

  if keyword_set(longest) then begin
     return, length[0]
  endif else begin
     return, length
  endelse

end
