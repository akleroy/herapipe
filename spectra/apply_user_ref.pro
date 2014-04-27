pro apply_user_ref $
   , list_file $
   , tag = tag

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SET SOME DEFAULTS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; MINIMUM REFERENCE TIME IN SECONDS (USED WITH THE RULE OF THUMB)
  tref_min = 5.0                

; TOLERANCE FOR A SCAN ANGLE MATCH (IN DEGREES)
  tol_scan_ang = 1.0
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file $
           , working_name, already_ref_sub, ref_defn_file $
           , format='X,A,A,A', /silent $
           , comment="#"
  working_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)
  already_ref_sub = strupcase(strcompress(already_ref_sub, /remove_all))
  ref_defn_file = strcompress(ref_defn_file,/remove_all)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA SETS AND IDENTIFY REFERENCES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin

;    FLAG TELLING US TO WRITE
     need_to_write = 0

;    CHECK IF WE HAVE ANY USER FLAGS TO APPLY
     dummy = file_search(ref_defn_file[i], count=ref_defn_ct)
     if ref_defn_file[i] eq 'none' or ref_defn_ct eq 0 then $
        continue     

;    READ THE DATA
     indir = '../spectra/'
     infile = indir+working_name[i]+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
     
     print, 'ding!'

;    MEASURE THE SIZE OF THE DATA
     sz = size(data)
     
;    READ THE USER-DEFINED FILE
     readcol, ref_defn_file[i] $
              , ref_subscan, ref_scanang, ref_time $
              , format='A,A,A', /silent $
              , comment="#"
     ref_subscan = strcompress(ref_subscan,/remove_all)
     ref_scanang = strcompress(ref_scanang,/remove_all)
     ref_time = strcompress(ref_time,/remove_all)
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    LOOP OVER REFERENCE DEFINITIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     for j = 0, n_elements(ref_subscan)-1 do begin

;       BUILD AN INDEX FOR THIS DEFINITION

;       ... START WITH THE SUBSCAN
        if ref_subscan[j] eq 'all' then $
           is_ref = bytarr(sz[1])+1B $
        else $
           is_ref = data.subscan eq ref_subscan[j]

;       ... NEXT IMPLEMENT THE SCANANGLE
        if ref_scanang[j] ne 'all' then $
           is_ref *= abs(data.scan_ang - ref_scanang[j]) lt tol_scan_ang

;       ... FINALLY THE TIME
        if ref_time[j] ne 'all' then begin
           if ref_time[j] eq 'thumb' then begin
;          ... THE RULE-OF-THUMB CASE
              sub_scan_length = data.subscan_stop_ut - data.subscan_start_ut
              ref_needed = data.obstime * sqrt(sub_scan_length / data.obstime)
              ref_needed = ref_needed > tref_min
              is_ref *= $
                 (((data.ut - data.subscan_start_ut) le ref_needed) or $
                  ((data.subscan_stop_ut - data.ut) le ref_needed))
           endif else begin
;          ... THE USER SPECIFIES A TIME
              is_ref *= $
                 (((data.ut - data.subscan_start_ut) le ref_time[j]) or $
                  ((data.subscan_stop_ut - data.ut) le ref_time[j]))
           endelse
        endif
        
;       EVERYTHING THAT SURVIVES THIS HASH IS A REFERENCE
        ref_ind = where(is_ref, ref_ct)

;       SET THE REFERENCE POSITIONS TO HAVE ON = 0 IN THE DATA STRUCTURE
        if ref_ct gt 0 then begin
           need_to_write = 1
           data[ref_ind].on = 0
        endif

     endfor
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ... WRITE (IF NEEDED)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%    

     if need_to_write then begin

;       WRITE THE DATA (FIRST DELETE TO KEEP FROM APPENDING)
        outfile = infile
        spawn, 'rm '+outfile
        mwrfits, data, outfile, hdr

     endif
  
  endfor

end                             ; of apply_user_ref
