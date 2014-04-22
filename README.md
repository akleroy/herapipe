HERA pipeline

Example directory structure for one of the HERACLES targets:

ngc2903/
|-- code
|   |-- bad_data.txt
|   |-- bad_data.txt~
|   |-- cube_info_beta.idl
|   |-- default_ref.txt
|   |-- default_ref.txt~
|   |-- first_red
|   |-- hi
|   |-- log.txt
|   |-- log.txt~
|   |-- orig_data.txt
|   |-- orig_data.txt~
|   |-- pipe_batch.pro
|   |-- pipe_batch.pro~
|   |-- pixel_gains_beta.txt
|   |-- reports_poly2
|   |-- test_data.txt
|   `-- test_data.txt~
|-- cubes
|   |-- _beta
|   |   |-- mask_blank_beta.fits
|   |   |-- mask_line_beta.fits
|   |   |-- ngc2903_beta.fits
|   |   |-- ngc2903_beta_bright_map.fits
|   |   |-- ngc2903_beta_mask.fits
|   |   |-- ngc2903ang1_beta.fits
|   |   |-- ngc2903ang1_beta_raw.coverage.fits
|   |   |-- ngc2903ang1_beta_raw.fits
|   |   |-- ngc2903ang1_beta_raw.nspec.fits
|   |   |-- ngc2903ang2_beta.fits
|   |   |-- ngc2903ang2_beta_raw.coverage.fits
|   |   |-- ngc2903ang2_beta_raw.fits
|   |   |-- ngc2903ang2_beta_raw.nspec.fits
|   |   |-- noise_beta.fits
|   |   `-- temp.fits
|   |-- reference
|   `-- windows
|-- masks
|   |-- first_red
|   |-- hi
|   |-- reference
|   |   |-- make_on_mask.pro
|   |   |-- make_on_mask.pro~
|   |   `-- ngc2903_on_mask.fits
|   `-- windows
|       |-- make_line_mask.pro
|       |-- make_line_mask.pro~
|       |-- ngc2903_hi_high.fits
|       |-- ngc2903_hi_low.fits
|       |-- ngc2903_line_high.fits
|       |-- ngc2903_line_low.fits
|       |-- ngc2903_loweredge_high.fits
|       |-- ngc2903_loweredge_low.fits
|       |-- ngc2903_mask.fits
|       |-- ngc2903_upperedge_high.fits
|       `-- ngc2903_upperedge_low.fits
|-- other_data
|   |-- first_red
|   |   |-- ngc2903_mask.fits
|   |   |-- ngc2903_prelim.fits
|   |   `-- ngc2903_prelim.heracles.fits
|   |-- hi
|   |   |-- ngc2903_things_mom1.fits
|   |   `-- ngc2903_things_na.fits
|   `-- reports_poly2
|       |-- fft_ngc2903_13feb08.jpeg
|       |-- fft_ngc2903_16feb08.jpeg
|       |-- fft_ngc2903_17feb08.jpeg
|       |-- fft_ngc2903_21feb08.jpeg
|       |-- fft_ngc2903_27nov07.jpeg
|       |-- fft_specngc2903_13feb08.txt
|       |-- fft_specngc2903_16feb08.txt
|       |-- fft_specngc2903_17feb08.txt
|       |-- fft_specngc2903_21feb08.txt
|       |-- fft_specngc2903_27nov07.txt
|       |-- flag_noise_ngc2903_13feb08.jpeg
|       |-- flag_noise_ngc2903_16feb08.jpeg
|       |-- flag_noise_ngc2903_17feb08.jpeg
|       |-- flag_noise_ngc2903_21feb08.jpeg
|       |-- flag_noise_ngc2903_27nov07.jpeg
|       |-- flag_ripple_ngc2903_13feb08.jpeg
|       |-- flag_ripple_ngc2903_16feb08.jpeg
|       |-- flag_ripple_ngc2903_17feb08.jpeg
|       |-- flag_ripple_ngc2903_21feb08.jpeg
|       |-- flag_ripple_ngc2903_27nov07.jpeg
|       |-- flagging_ngc2903_13feb08.txt
|       |-- flagging_ngc2903_16feb08.txt
|       |-- flagging_ngc2903_17feb08.txt
|       |-- flagging_ngc2903_21feb08.txt
|       |-- flagging_ngc2903_27nov07.txt
|       |-- gain_ngc2903_13feb08.ps
|       |-- gain_ngc2903_16feb08.ps
|       |-- gain_ngc2903_17feb08.ps
|       |-- gain_ngc2903_21feb08.ps
|       |-- gain_ngc2903_27nov07.ps
|       |-- noise_report_beta.txt
|       |-- noise_vs_scale_beta.jpeg
|       |-- onoff_ngc2903_13feb08.jpeg
|       |-- onoff_ngc2903_16feb08.jpeg
|       |-- onoff_ngc2903_17feb08.jpeg
|       |-- onoff_ngc2903_21feb08.jpeg
|       `-- onoff_ngc2903_27nov07.jpeg
|-- reports
|   |-- fft_ngc2903_13feb08.jpeg
|   |-- fft_ngc2903_16feb08.jpeg
|   |-- fft_ngc2903_17feb08.jpeg
|   |-- fft_ngc2903_21feb08.jpeg
|   |-- fft_ngc2903_27nov07.jpeg
|   |-- fft_specngc2903_13feb08.txt
|   |-- fft_specngc2903_16feb08.txt
|   |-- fft_specngc2903_17feb08.txt
|   |-- fft_specngc2903_21feb08.txt
|   |-- fft_specngc2903_27nov07.txt
|   |-- flag_noise_ngc2903_13feb08.jpeg
|   |-- flag_noise_ngc2903_16feb08.jpeg
|   |-- flag_noise_ngc2903_17feb08.jpeg
|   |-- flag_noise_ngc2903_21feb08.jpeg
|   |-- flag_noise_ngc2903_27nov07.jpeg
|   |-- flag_ripple_ngc2903_13feb08.jpeg
|   |-- flag_ripple_ngc2903_16feb08.jpeg
|   |-- flag_ripple_ngc2903_17feb08.jpeg
|   |-- flag_ripple_ngc2903_21feb08.jpeg
|   |-- flag_ripple_ngc2903_27nov07.jpeg
|   |-- flagging_ngc2903_13feb08.txt
|   |-- flagging_ngc2903_16feb08.txt
|   |-- flagging_ngc2903_17feb08.txt
|   |-- flagging_ngc2903_21feb08.txt
|   |-- flagging_ngc2903_27nov07.txt
|   |-- gain_ngc2903_13feb08.ps
|   |-- gain_ngc2903_16feb08.ps
|   |-- gain_ngc2903_17feb08.ps
|   |-- gain_ngc2903_21feb08.ps
|   |-- gain_ngc2903_27nov07.ps
|   |-- noise_report_beta.txt
|   |-- noise_report_prev.txt
|   |-- noise_vs_scale_beta.jpeg
|   |-- onoff_ngc2903_13feb08.jpeg
|   |-- onoff_ngc2903_16feb08.jpeg
|   |-- onoff_ngc2903_17feb08.jpeg
|   |-- onoff_ngc2903_21feb08.jpeg
|   `-- onoff_ngc2903_27nov07.jpeg
`-- spectra
    |-- _beta
    |-- _fake
    |-- _off1
    |-- _off2
    |-- ngc2903_13feb08_beta.processed.fits
    |-- ngc2903_16feb08_beta.processed.fits
    |-- ngc2903_17feb08_beta.processed.fits
    |-- ngc2903_21feb08_beta.processed.fits
    `-- ngc2903_27nov07_beta.processed.fits
