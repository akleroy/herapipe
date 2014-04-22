pro update_telescop_field, data

  uniq_telescop = (data.telescop)[uniq(data.telescop, sort(data.telescop))]
  needs_replacing = 0B
  for k = 0L, n_elements(uniq_telescop)-1 do begin
     if (strpos(uniq_telescop[k],'WA') ne -1) or $
        (strpos(uniq_telescop[k],'WB') ne -1) then $
           needs_replacing = 1B
  endfor

  if needs_replacing eq 0B then return
  
  message, 'Replacing old TELESCOP field with new values.', /info
  
  for k = 0L, n_elements(data)-1 do begin
     str = data[k].telescop
     str = repstr(str, 'WA', 'W01')
     str = repstr(str, 'WB', 'W02')
     data[k].telescop = strcompress(str,/remove_all)
  endfor

end
