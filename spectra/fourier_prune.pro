pro fourier_prune $
   , list_file $
   , tag = tag $
   , bad_channels = bad_channels $
   , show = show $
   , report = report $
   , channel_prof = kern_in $
   , monte = monte

;
; AKL - BE CAREFUL HERE!
;
; IN ADDITION TO NO OVERARCHING DOCUMENTATION YET, I'M REUSING VARIABLES HERE
; TO KEEP THE MEMORY USAGE LIGHT.
;
  
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; SET SOME DEFAULTS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; A TOY CHANNEL PROFILE ... COULD BE REPLACED LATER
  if n_elements(kern_in) eq 0 then $
     kern_in = [0.25, 1.0, 0.25]
  
  kern = kern_in/total(kern_in)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ LIST OF DATA FILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  readcol, list_file, working_name $
           , format='X,A', /silent $
           , comment="#"
  data_name = strcompress(working_name, /remove_all)
  ndata = n_elements(working_name)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER DATA SETS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for i = 0, ndata-1 do begin

;    READ THE DATA
     indir = working_dir+'spectra/'
     infile = indir+working_name[i]+'_'+tag+'.processed.fits'
     dummy = file_search(infile, count=count)
     if count eq 0 then begin
        message, 'File not found '+string(working_name[i])+'. Skipping.', /info
        continue
     endif
     data = mrdfits(infile,1,hdr, /silent)
     
;    MEASURE THE SIZE OF THE DATA
     sz = size(data)

;    WORK ONLY WITH REFERENCE-SUBTRACTED, ON SOURCE, UNFLAGGED DATA
     on_ind = where(data.on eq 1 and $
                    data.ref_subtracted eq 1 and $
                    data.flagged eq 0, on_ct)
     if on_ct eq 0 then begin
        message, 'Only Fourier prune on reference-subtracted data.', /info
        continue
     endif

;    MEASURE THE NUMBER OF CHANNELS AND SPECTRA
     n_chan = n_elements(data[0].spectrum)
     n_spec = on_ct ; (size(data))[1]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS THEN TAKE FFT OF DATA, BUILD AXIS, ETC.
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;    REPLACE EDGES WITH (SLIGHTLY) INNER VALUES. EDGES ARE LATER
;    IGNORED, BUT be CAREFUL!  I'M NOT SURE IF THIS IS DESIRABLE.

     for j = 0L, n_spec-1 do begin
        spec = data[on_ind[j]].spectrum
        spec[0:9] = mean(spec[10:20],/nan)
        spec[n_chan-10:n_chan-1] = mean(spec[n_chan-20:n_chan-11],/nan)
        data[on_ind[j]].spectrum = spec
     endfor

;    SEPARATE OUT AN ARRAY OF SPECTRA
     spec_ra = data[on_ind].spectrum

;    TURN NOT-A-NUMBERS INTO ZEROS FOR THE FFT (DEBATABLE AND CRAPPY)
     nan_ind = where(finite(spec_ra) eq 0, nan_ct)
     if nan_ct gt 0 then spec_ra[nan_ind] = 0.

;    TAKE THE FFT OF EACH SPECTRUM
     fft_ra = fft(spec_ra,-1,dim=1) ;

;    CALCULATE THE X-AXIS FOR THE POWER SPECTRUM
     dfreq_mhz = sxpar(hdr,'CDELT1')
     power_spec_xaxis = findgen(n_chan/2.0+1) / (n_chan*dfreq_mhz)     
     
;    GET THE MEAN POWER SPECTRUM AVERAGED OVER ALL SPECTRA
     mean_power_spec = $
        total((abs(fft_ra)^2)[0:(n_chan/2.0+1),*],2,/nan)/(n_spec*1.0)
;        median(power_spec_ra[0:(n_chan/2+1),*],dim=2)

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; REPLACE BAD FFT CHANNELS WITH INTERPOLATED DATA FROM NEARBY CHANNELS
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     for j = 0L, n_elements(bad_channels) - 1 do begin

        chan_lo = bad_channels[j]
        chan_hi = n_chan - bad_channels[j]

;       INTERPOLATE IN POWER (RATHER THAN MAGNITUDE) ... AN OPEN QUESTION
        amp_vec = $
           sqrt((abs(fft_ra[chan_lo-2,*])^2 + abs(fft_ra[chan_lo+2,*])^2)*0.5)

;       ... AN OPEN QUESTION: RANDOM VS. INTERPOLATED PHASES?
        phase_vec = $
           randomu(seed,n_spec)*2.0*!pi

        fft_ra[chan_lo,*] = $
           complex(amp_vec*cos(phase_vec), amp_vec*sin(phase_vec))

        fft_ra[chan_hi,*] = $
           complex(amp_vec*cos(-1.*phase_vec), amp_vec*sin(-1.*phase_vec))

     endfor

     if n_elements(bad_channels) gt 0 then begin
;       REVERSE THE FFT AND PLACE THE NEW DATA INTO THE OLD TABLE
        data[on_ind].spectrum = real_part(fft(fft_ra,dim=1,/inverse))
        
;       RE-TAKE THE FFT OF EACH SPECTRUM
;       (MAKES SURE THE PHASE DOESN'T SCREW US)
        spec_ra = data[on_ind].spectrum
        fft_ra = fft(spec_ra,dim=1,-1)

;       UPDATE THE POWER SPECTRUM
        pruned_mean_power_spec = $
           total((abs(fft_ra)^2)[0:(n_chan/2+1),*],2,/nan)/(n_spec*1.0)
        pruned_phase_ra = atan(fft_ra, /phase)        
     endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; TEST A CHANNEL PROFILE USING A MONTE CARLO REALIZATION
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

;    IF REQUESTED, MAKE A MONTE CARLO SPECTRUM
     if keyword_set(monte) then begin
        for j = 0L, n_spec-1 do $
           spec_ra[*,j] -= mean(spec_ra[*,j],/nan)
        rms = mad(spec_ra)
        
        sz = size(spec_ra)
        spec_ra = randomn(seed,sz[1],sz[2])
        
        for j = 0L, n_spec-1 do $
           spec_ra[*,j] = convol(spec_ra[*,j], kern, /edge_truncate)
                
        spec_ra *= rms/mad(spec_ra)
        fft_ra = fft(spec_ra, dim=1, -1)
        
        monte_mean_power_spec = $
           total((abs(fft_ra)^2)[0:(n_chan/2+1),*],2,/nan)/(n_spec*1.0)
     endif


; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; PLOT THE VARIOUS POWER SPECTRA IF REQUESTED
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

     if keyword_set(show) then begin

        loadct, 0, /silent
        reversect

;       PLOT THE POWER SPECTRUM, OVERPLOT THE MEAN, THEN RE-PLOT
        plot, power_spec_xaxis, mean_power_spec $
              , xtitle='!6Inverse Frequency [MHz!u-1!n]' $
              , ytitle='!6Power (abs(FFT)!u2!n)' $
              , yrange=[0., max(mean_power_spec[30:*],/nan)] $
              , title='Power Spectrum '+data_name[i]

        if n_elements(pruned_mean_power_spec) gt 0 then begin
           oplot, power_spec_xaxis, pruned_mean_power_spec $
                  , color=getcolor('red') $
                  , thick=5
        endif

        if n_elements(monte_mean_power_spec) gt 0 then begin
           oplot, power_spec_xaxis, monte_mean_power_spec $
                  , color=getcolor('forest') $
                  , thick=5
        endif

        oplot, power_spec_xaxis, mean_power_spec $
               , color=getcolor('black'), thick=3
        legend, lines=[0,0,0], box=0 $
                , /top, /right $
                , ['!6Pow. Spec. Before Pruning', $
                   'Pow. Spec. After Pruning', $
                   'Monte Carlo + Channel Profile'] $
                , color= $
                [getcolor('black'), getcolor('red'), getcolor('forest')] $
                , textcolor=getcolor('black') $
                , thick=[5,5,5]
        
        if keyword_set(report) then begin
           im = tvrd(true=1)
           write_jpeg, working_dir+'reports/fft_'+data_name[i]+'.jpeg', im, true=1
        endif

     endif

;    WRITE THE POWER SPECTRUM TO DISK IN A TEXT FILE
     if keyword_set(report) then begin
        get_lun, lun
        openw, lun, working_dir+'reports/fft_spec'+data_name[i]+'.txt'
        printf,lun,'Column 1: Frequency [MHz^-1]'
        printf,lun,'Power: (abs(FFT)^2 of spectrum in TA*)'
        for jj = 0, n_elements(power_spec_xaxis)-1 do $
           printf,lun,power_spec_xaxis[jj], mean_power_spec[jj]
        close, lun
        free_lun, lun
     endif

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; WRITE PROCESSED DATA TO DISK
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

;    (DELETING THE OLD VERSION TO KEEP MWRFITS FROM APPENDING)
     outfile = infile
     spawn, 'rm '+outfile
     mwrfits, data, outfile, hdr

  endfor

end
