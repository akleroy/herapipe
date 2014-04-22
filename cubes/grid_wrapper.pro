pro grid_wrapper $
   , list_file $
   , tag = tag $
   , outfile = outfile $
   , target_hdr = target_hdr $
   , apply_gain = apply_gain $
   , gain_file = gain_ascii_file $
   , split_by_angle = split_by_angle $
   , scan_angle = scan_angle $
   , grid_by_median=use_median $
   , gauss_kern = gauss_kern $  
   , gauss_fwhm = gauss_fwhm $  
   , show=show                  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN *ALL* DATA (MEMORY INEFFICIENT, BUT SIMPLE AND FAST)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin

;    READ THE DATA
     indir = '../spectra/'
     infile = indir+working_name[i]+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)

;    ISOLATE THE DATA THAT WE ARE INTERESTED IN
     if keyword_set(split_by_angle) then begin
        ind = where(data.flagged eq 0 and $
                    data.on eq 1 and $
                    data.ref_subtracted eq 1 and $
                    round(data.scan_ang) eq scan_angle, ct)
     endif else begin
        ind = where(data.flagged eq 0 and $
                    data.on eq 1 and $
                    data.ref_subtracted eq 1, ct)
     endelse

     if ct eq 0 then continue

     data = data[ind]

;    IF REQUESTED, GAIN CORRECTIONS TO THE DATA
     if keyword_set(apply_gain) then begin
        readcol, gain_file, format='A,A,F,F' $
                 , gain_day, gain_pixel, gain_val, gain_unc $
                 , count = nlines
        
;       HARDCODE SOME PRUNING
        gain_cap = 1.5
        bad_gains = where((gain_val gt gain_cap) or $
                          (gain_val lt 1./gain_cap), bad_ct)
        if bad_ct gt 0 then begin
           message, 'Some gains outside reasonable range. '+$
                    'Capping these at edge of range.', /info
           gain_val[bad_gains] = $
              ((gain_val[bad_gains] > 1./gain_cap) < gain_cap)
        endif

;       APPLY THE GAINS
        for j = 0, nlines-1 do begin
           ind = where((data.telescop eq gain_pixel[j]) and $
                       (working_name[i] eq gain_day[j]), ct)
           if (ct gt 0) then $
              data[ind].spectrum = data[ind].spectrum * gain_val[i]
        endfor
        
     endif

;    IDENTIFY WEIGHTS FOR THE DATA
     if total(finite(data.fit_rms) eq 0) gt 0 or $
        total(data.fit_rms eq 0) gt 0 then begin
        message, 'RMS about the fit not available for '+$
                 'all spectra. Using Tsys-based weighting.', /info
        this_weight = 1./data.tsys^2 
     endif else begin 
        this_weight = 1./data.fit_rms^2
     endelse

;    CHECK THAT THE VELOCITY AXIS IS THE SAME FOR EACH DATA SET
     if total(data.v0 ne data[0].v0) gt 0 or $
        total(data.deltav ne data[0].deltav) gt 0 $
     then begin
        message, 'Velocity axis mismatch. Stopping for inspection.'    
     endif
     
;    SAVE THE COORDINATES, WEIGHT, AND SPECTRA
     if n_elements(ra) eq 0 then begin
        ra = data.ra_deg
        dec = data.dec_deg
        spec = transpose(data.spectrum)
        weight = this_weight
     endif else begin
        ra = [ra,data.ra_deg]
        dec = [dec,data.dec_deg]
        spec = [spec,transpose(data.spectrum)]
        weight = [weight,this_weight]
     endelse

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALL THE GRIDDING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    ADD THE SCAN ANGLE TO THE HEADER IF APPROPRIATE
     if keyword_set(split_by_angle) then $
        sxaddpar, target_hdr, 'SCANANGL', scan_angle, 'ROUNDED SCAN ANGLE'


;    TOSS THESE DATA INTO THE CUBE
     grid_otf, data=spec $
               , ra=ra $
               , dec=dec $
               , weight=weight $
               , target = target_hdr $
               , out_root = outfile $
               , median = use_median $
               , gauss_kern = gauss_kern $
               , gauss_fwhm = gauss_fwhm $
               , show = show
 
end                             ; of grid_wrapper
