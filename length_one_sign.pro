function length_one_sign, spec, longest=longest
  
  if (total(spec gt 0) gt 0) then begin     
     positive = label_region(spec gt 0)  
     length_pos = intarr(max(positive))
     for i = 1, max(positive) do $
        length_pos[i-1] = total(positive eq i)
  endif

  if (total(spec lt 0) gt 0) then begin     
     negative = label_region(spec lt 0)
     length_neg = intarr(max(negative))
     for i = 1, max(negative) do $
        length_neg[i-1] = total(negative eq i)
  endif
  
  length = [length_pos, length_neg]
  length = length[sort(length)]
  length = reverse(length)

  if keyword_set(longest) then begin
     return, length[0]
  endif else begin
     return, length
  endelse

end
