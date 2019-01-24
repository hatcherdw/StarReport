;;======================================================================;;
;; Star Report Pipeline, Central Michigan University
;; 
;; spectrum_analyzer.pro 
;; 
;; The main goal of this program is to read in flux and wavelength 
;; data from a text file and store this data in two arrays (flux and 
;; wave) to generate a plot of a single spectrum. All plotting is done
;; elsewhere but the arrays are constructed here. Equivalent width
;; (EW) is also calculated within this procedure.
;;
;;
;; Authors: Christian Hannah & Glennon Fagan 
;;
;;======================================================================;;


PRO spectrum_analyzer, temp_file
RESOLVE_ALL, /QUIET
COMMON spec_block, cc_waves, center, EW, ew_error, flag, flux, fwhm, ha_waves, mask_flux, peak_flux, peak_num, peak_wave, peak1_flux, peak1_wave, peak2_flux, peak2_wave, shell_para, snr, snr_error, symmetry, vr_ratio, wave
COMMON star_block, analysis_logs, root_direc, available_dates, available_frames, analysis_dates, analysis_frames, mult_stars, star_dec, star_hd, star_hr, star_proper_name, star_name, star_ra, star_vmag, start_time, sub_direc, targets_dir_name, trend, ver_num

;; Define a flag for spectra experiencing analytical errors
flag = INTARR(3)

;;Define some color tables
color_table = [60,250,150,220,30,90,110,190,10]
color_text = ['blue','red','green','orange','violet','lt_blue','cyan','yellow','black']

;; Check to see if spectrum_plot has been ran alone
IF ISA(temp_file) THEN spec_plot_alone = 0 ELSE spec_plot_alone = 1

IF spec_plot_alone THEN BEGIN
   SPAWN, 'pwd', root_direc
   root_direc = root_direc+'/' ;;working directory 
   sub_direc = 'Spectra'       ;;directory the spectra files are stored in
ENDIF

;; prompt the user to enter a frame number for spec_plot_alone
IF spec_plot_alone THEN BEGIN
   temp_file = ''
   READ, temp_file, PROMPT='Enter frame number: '
ENDIF

;; get a path to the file location so it can be in a seperate directory 
spec_file = FILEPATH('*'+temp_file+'*', ROOT_DIR=root_direc, SUBDIR=sub_direc)




;;***************************** IMPORT the observed spectrum first ****************************;;

;; print status
text = 'Importing data from an ASCII file - Frame: '+ temp_file
IF spec_plot_alone EQ 1 THEN BEGIN
   PRINT, ' '
   PRINT, text
   PRINT, ' '
ENDIF ELSE BEGIN
   analysis_logs = [analysis_logs,'']
   analysis_logs = [analysis_logs, text]
   analysis_logs = [analysis_logs,'']
ENDELSE

;; determine the number of lines in the ASCII file
check_if_ascii = QUERY_ASCII(spec_file, ascii)
number_of_lines = ascii.lines

;; read in the entire ASCII file into a string array
ascii_strings = STRARR(number_of_lines)

text = 'Number of ASCII lines detected ='+ STRCOMPRESS(STRING(number_of_lines))
IF spec_plot_alone EQ 1 THEN BEGIN
   PRINT, text
ENDIF ELSE BEGIN
   analysis_logs = [analysis_logs,text]
ENDELSE

tmp_string = ''
OPENR, file_unit, spec_file, /GET_LUN
counter = 0L
WHILE NOT EOF(file_unit) DO BEGIN
   ;; read the entire line as a string
   READF, file_unit, tmp_string
   ascii_strings[counter] = tmp_string
   counter += 1
ENDWHILE
CLOSE, file_unit
FREE_LUN, file_unit

;; setup the main variales that will hold the imported values
;; Note: these might have to be trimmed
wave = DBLARR(number_of_lines)
flux = DBLARR(number_of_lines)


;;define some needed variables and formats
tmp_wave = 0D
tmp_flux = 0D
counter = 0L
old_format_1 = '       ???.?????      ?.????????' 
old_format_2 = '       ???.?????       ?.???????'
old_format_3 = '       ???.?????       ??.??????'
new_format = ' ???.??????  ??.??????'

error_lines = INTARR(number_of_lines)
j = 0

;; import the values from each valid string
FOR i = 0, N_ELEMENTS(ascii_strings)-1 DO BEGIN
      ;; determine the format of the ASCII file (old vs new)
      IF STRMATCH(ascii_strings[i], old_format_1) EQ 1 THEN BEGIN
        ;; old format 
        READS, ascii_strings[i], tmp_wave, tmp_flux, FORMAT='(F16.5, F16.8)'
        counter += 1
      ENDIF ELSE IF STRMATCH(ascii_strings[i], old_format_2) EQ 1 THEN BEGIN
        READS, ascii_strings[i], tmp_wave, tmp_flux, FORMAT='(F16.5, F16.7)'
        counter += 1
      ENDIF ELSE IF STRMATCH(ascii_strings[i], new_format) EQ 1 THEN BEGIN
         ;; new format
         READS, ascii_strings[i], tmp_wave, tmp_flux, FORMAT='(F11.6, F11.6)'
         counter += 1
      ENDIF ELSE IF STRMATCH(ascii_strings[i], old_format_3) EQ 1 THEN BEGIN
         READS, ascii_strings[i], tmp_wave, tmp_flux, FORMAT='(F16.5, F16.6)'
         counter += 1
      ENDIF ELSE BEGIN
         ;; count lines that are not of known format
         error_lines[j] = i+1
         j += 1
         tmp_wave = 0D
         tmp_flux = 0D
      ENDELSE
     wave[counter-1] = tmp_wave
     flux[counter-1] = tmp_flux
ENDFOR


;; print how many lines were imported and how many that were not of known format
text  = 'Number of lines imported from the ASCII file ='+ STRCOMPRESS(STRING(counter))
text1 = '**WARNING: File contains'+ STRCOMPRESS(STRING(N_ELEMENTS(WHERE(error_lines NE 0))))+ $
        ' lines of unknown format.**'
error_line_nums_arr = STRCOMPRESS(STRING(error_lines[WHERE(error_lines NE 0)]))
error_line_nums= ''
FOR i=0, N_ELEMENTS(error_line_nums_arr)-1 DO BEGIN
   error_line_nums += error_line_nums_arr[i]
ENDFOR
text2 = 'Line Number(s):'+error_line_nums
IF spec_plot_alone EQ 1 THEN BEGIN
   PRINT, text
   PRINT, ' '
   PRINT, text1
   PRINT, text2
   PRINT, ' '
ENDIF ELSE BEGIN
   analysis_logs = [analysis_logs,text]
   analysis_logs = [analysis_logs,'']
   analysis_logs = [analysis_logs,text1]
   analysis_logs = [analysis_logs,text2]
   analysis_logs = [analysis_logs,'']
ENDELSE


;; remove the elements from the arrays that were not needed
wave = wave[0:counter-1]
flux = flux[0:counter-1]







;;***************************** CALL PROCEDURE Noise_Analyzer *******************************;;

noise_analyzer, color_table, smooth_flux, spec_plot_alone, temp_file, $
                GAUSS_SIGMA=gauss_sigma, INTSIG_WAVE_INT=intsig_wave_int, INTSIG_WAVE_FIN=intsig_wave_fin, $
                INT_SIG=int_sig




;;***************************** CALL PROCEDURE Peak_Analyzer *******************************;;

peak_analyzer, spec_plot_alone, temp_file, color_table







;;*************************** Calculate the equivalent width (EW) ***************************;;

;; start stddev at flux[40] because some spectra have high flux values                                      
;; at the beginning                                                                                         
left_bound = VALUE_LOCATOR(wave, 652.55)
right_bound = VALUE_LOCATOR(wave, 653.09)
std_flux = STDDEV(flux[left_bound :right_bound])

;; determine delta wavelengths and delta flux
delta_wave = DBLARR(N_ELEMENTS(wave)-1)
delta_flux = DBLARR(N_ELEMENTS(flux)-1)
FOR i=0, N_ELEMENTS(wave)-2 DO BEGIN
   delta_wave[i] = wave[i+1] - wave[i]
   delta_flux[i] = (flux[i+1] - flux[i])/2
ENDFOR

;; smooth the spectrum                                                        
smooth_flux = SMOOTH(flux, 30)

;; determine start and stop indecies for integration   
c = 2.0D ;; constant to multiply by std_flux
i = MIN(WHERE(wave GE 654.0)) 
WHILE (smooth_flux[i] LE (1.0 + (c*std_flux)) && smooth_flux[i] GE (1.0 - (c*std_flux))) DO BEGIN
   i += 1
   IF i EQ N_ELEMENTS(flux) - 1 THEN BREAK
ENDWHILE

;; define wave_start (subtract 10 from i to make sure the entire
;; emission is included)
wave_start = i-10


adjusted = 0 ;;this will be used to trigger a print statement informing the user that 
             ;;the start index was set manually

;;make sure the wave_start makes sense
IF wave_start GT value_locator(wave, 657.0) THEN BEGIN 
   wave_start = value_locator(wave, 654.509)
   i = wave_start
   adjusted = 1
ENDIF


i = i+20 ;; so we start looking for the return to continuum from inside the emission/absorption 
WHILE (smooth_flux[i] GE 1.0 + (c*std_flux) || smooth_flux[i] LE 1.0 - (c*std_flux)) || $
      (smooth_flux[i+10] GE 1.0 + (c*std_flux) || smooth_flux[i+10] LE 1.0 - (c*std_flux)) || $
      (smooth_flux[i+15] GE 1.0 + (c*std_flux) || smooth_flux[i+15] LE 1.0 - (c*std_flux)) || $
      (wave[i] LT 656.5) DO BEGIN
   i +=1
   IF i GE value_locator(wave, 658.407) THEN BEGIN ;; avoid illegal subscript error
      wave_end = value_locator(wave, 658.407)
      IF adjusted THEN adjusted = 3 ELSE adjusted = 2
      break
   ENDIF
ENDWHILE


IF ~ISA(wave_end) THEN wave_end = i + 5 ;;add 5 to make sure the whole emission is included

;; warn the user that the start and stop indices were manually set
text  = '** NOTICE: Start index for integration set manually. **'
text1 = '** NOTICE: End index for integration set manually. **'
text2 = '** NOTICE: Start and end index for integration set manually. **'
IF spec_plot_alone EQ 1 THEN BEGIN
   CASE adjusted OF
      1: PRINT, text
      2: PRINT, text1
      3: PRINT, text2
      ELSE:
   ENDCASE
   PRINT, ' '
ENDIF ELSE BEGIN
   CASE adjusted OF
      1: analysis_logs = [analysis_logs,text]
      2: analysis_logs = [analysis_logs,text1]
      3: analysis_logs = [analysis_logs,text2]
      ELSE:
   ENDCASE
analysis_logs = [analysis_logs,'']
ENDELSE


;; do arithmetic for EW
EW = 0.0D
FOR i= wave_start, wave_end DO BEGIN
   EW += ((1-flux[i]) + delta_flux[i]) * delta_wave[i]
ENDFOR
 

;; Check EW value with IDL function
;; subtract the flux values from 1 
flux_test = DBLARR(N_ELEMENTS(flux))
FOR i=wave_start, wave_end DO BEGIN
   flux_test[i] = 1 - flux[i]
ENDFOR
area = INT_TABULATED(wave, flux_test)


;; Propagate error for EW
sum_errors = 0.0D
FOR i=wave_start, wave_end DO BEGIN
   sum_errors += (ABS((delta_flux[i] + 1 - flux[i]) * delta_wave[i]) * (std_flux/flux[i]))^2
   ew_error = SQRT(sum_errors)
ENDFOR

text = '   Equivalent Width:'+ STRCOMPRESS(STRING(EW))+ ' ±'+ STRCOMPRESS(STRING(ew_error))+ ' nm'
IF spec_plot_alone EQ 1 THEN BEGIN
   PRINT, text
   PRINT, ' '
ENDIF ELSE BEGIN
   analysis_logs = [analysis_logs,text]
   analysis_logs = [analysis_logs,'']
ENDELSE


;; double check EW with IDL function and print warnings
text ='**WARNING: Inconsistencies in EW. EW calculation does not match IDL procedure EW.**' 
text1='   Manually calculated EW:'+ STRCOMPRESS(STRING(EW))+ ' nm'
text2='   IDL calculated EW:'+ STRCOMPRESS(STRING(area))+ ' nm'
text3='   Difference between the manual area and the IDL area is:'+ $
             STRCOMPRESS(STRING(ABS(area - EW)))+ ' nm'
IF ABS(area - EW) GT 0.01 THEN BEGIN
   IF spec_plot_alone EQ 1 THEN BEGIN
      PRINT, text
      PRINT, text1
      PRINT, text2
      PRINT, text3
   ENDIF ELSE BEGIN
      analysis_logs = [analysis_logs, text]
      analysis_logs = [analysis_logs,text1]
      analysis_logs = [analysis_logs,text2]
      analysis_logs = [analysis_logs,text3]
   ENDELSE
   flag[0] = 1
ENDIF




;;************************* Calculate EW with fixed wavelength *******************************;;

pm_wave = 2.3

fixed_wave_beg = value_locator(wave, 656.28-pm_wave)
fixed_wave_end = value_locator(wave, 656.28+pm_wave)

;; do arithmetic
ew_fixed = 0.0D
FOR i= fixed_wave_beg, fixed_wave_end DO BEGIN
   EW_fixed += ((1-flux[i]) + delta_flux[i]) * delta_wave[i]
ENDFOR

;; Propagate error for fixed EW
sum_errors = 0.0D
FOR i=fixed_wave_beg, fixed_wave_end DO BEGIN
   sum_errors += (ABS((delta_flux[i] + 1 - flux[i]) * delta_wave[i]) * (std_flux/flux[i]))^2
   fixed_ew_error = SQRT(sum_errors)
ENDFOR



;;************************* Add the HD number to output *********************************;;

IF spec_plot_alone THEN BEGIN ;;only do this for spectrum_analyzer running alone;
   
  ;; read in the obslog with csv_reader                                       
   variables = csv_reader('stellar_obslog.csv')                                
   check_if_ascii = QUERY_ASCII('stellar_obslog.csv', ascii)                   
   number_of_lines = ascii.lines

   ;;define some variables
   frame_num = STRARR(number_of_lines)
   name = ''
   ut_date = ''
   
   FOR i=0,number_of_lines-1 DO BEGIN
      ;;store each part of the line in a respective variable
      frame_num[i] = variables[1,i]
      IF frame_num[i] EQ temp_file THEN BEGIN
         name = variables[2,i]
         name = STRCOMPRESS(name, /REMOVE_ALL)
         ut_date = variables[3,i]
         IF STRMATCH(STRMID(name,0,1),'[!H1234567890]') THEN BEGIN            
            name = 'HR 9999'
         ENDIF ELSE IF STRMATCH(STRMID(name,0,1),'[!H]') THEN BEGIN
            name = 'HD '+name 
         ENDIF ELSE name = 'HR '+STRMID(name,2)
         break
      ENDIF
   ENDFOR

   ;; define star name and call hdhr_locator at end of program to print the hd and
   ;; hr numbers
   star_name = name

   hdhr_locator
  
   ;; report HD/HR/DATE to terminal
   PRINT, '   DATE: '+ut_date
   PRINT, ' '
   
ENDIF



;;************************************* Create Plots ************************************;;

plotter, temp_file, spec_plot_alone, color_table, ew_fixed, fixed_ew_error, fixed_wave_beg, $
         fixed_wave_end, intsig_wave_int, intsig_wave_fin, int_sig, wave_end, wave_start




END
