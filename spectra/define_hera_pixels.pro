pixel_list = '30M-'+ [ $
               'W01-1H01','W01-1H02','W01-1H03', $
               'W01-1H04','W01-1H05','W01-1H06', $
               'W01-1H07','W01-1H08','W01-1H09', $
               'W02-2H01','W02-2H02','W02-2H03', $
               'W02-2H04','W02-2H05','W02-2H06', $
               'W02-2H07','W02-2H08','W02-2H09']

if keyword_set(fts) then begin      
   pixel_list = '30M-'+ [ $
                'F02-1H01','F02-1H02','F02-1H03', $
                'F02-1H04','F02-1H05','F02-1H06', $
                'F02-1H07','F02-1H08','F02-1H09', $
                'F02-2H01','F02-2H02','F02-2H03', $
                'F02-2H04','F02-2H05','F02-2H06', $
                'F02-2H07','F02-2H08','F02-2H09']            

endif

; TBD - add the finer FTS option

npix = n_elements(pixel_list)

