;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; coordinator.pro
;;
;; This procedure does a number of things including: calling
;; hdhr_locator.pro and match_maker.pro. It will also convert the
;; observation dates to Julian dates. After, it will call
;; spectrum_analyzer.pro in order to conduct an individual analysis
;; on all the frames for the star. Once that analysis is completed it
;; will call the procedure stability_test.pro which look assess all
;; the spectra together to look for trends over time. coordinator.pro
;; also makes calls to halpha_centralizer.pro and
;; cross_correlator.pro which, if the user wanted, will align all the
;; spectra. The rest of the procedure does the following: reports
;; information to the terminal, conducts a consistency check between
;; the frame header information to the information gathered from the
;; stellar observation log, and creates various log files.
;;
;; Author: W. Glennon Fagan
;;======================================================================;;


PRO coordinator, BP_DATE_IND=bp_date_ind, CENTERS=centers, EWS_ERROR=EWs_error, $
                 EWS_TRANS=EWs, FLUXS=fluxs, FRAME_DATES_EXCLUDED = frame_dates_excluded, $
                 FRAMES_EXCLUDED = frames_excluded, FWHMS=fwhms, NUM_ANALYZED=num_analyzed, $
                 N_FILES=n_files, NUM_SING=num_sing, NUM_DOUB=num_doub, PEAK_ONE_WAVES=peak_one_waves,$
                 PEAK_ONE_FLUXS=peak_one_fluxs, PEAK_NUMS=peak_nums, PEAK_TWO_WAVES=peak_two_waves, $
                 PEAK_TWO_FLUXS=peak_two_fluxs, RP_DATE_IND=rp_date_ind, SHELL_PARAS=shell_paras, $
                 SING_PEAK_WAVES=sing_peak_waves, SING_PEAK_FLUXS=sing_peak_fluxs, SNRS=snrs, $
                 SORTED_FILES=sorted_files, SORTED_JDATES=sorted_jdates, SP_DATE_IND=sp_date_ind, $
                 SYMMETRIES=symmetries, VR_RATIOS=vr_ratios, WAVES=waves
COMPILE_OPT DEFINT32, STRICTARR, STRICTARRSUBS
RESOLVE_ALL, /QUIET
COMMON spec_block
COMMON star_block
COMMON abort_block

analysis_logs = ''  ;define the variable that will hold the information from the analysis


;;******************************** CALL HDHR_LOCATOR ************************************;;

full_pipeline = 1 ;;variable used to ensure HDHR_locator will be used as if it were being executed alone
hdhr_locator, full_pipeline



;;************************ SKIP PROGRAM IF NUMBER NOT IN CATALOG *************************;;

IF no_hdhr_num EQ 1 THEN BEGIN
   GOTO, jump
ENDIF




;;********************* CALL MATCH_MAKER.PRO TOGET MATCHING FILES ***********************;;

IF ~mult_stars THEN BEGIN
   PRINT, '=========================== START ==========================='
   PRINT, "Searching for matching frames: "
ENDIF

match_maker, EXPECTED_DATES=expected_dates, EXPECTED_FRAMES=expected_frames, $
             FILE_DATES=file_dates, FRAME_DATES_EXCLUDED = frame_dates_excluded, $
             FRAMES_EXCLUDED = frames_excluded, FRAME_OBSERVER=frame_observer, $
             HEADER_TARGET_NAME=header_target_name, MISSING_DATES=missing_dates,$
             MISSING_FRAMES=missing_frames, MISSING_FRAME_OBS = missing_frame_obs, $
             NAME_STRINGS=name_strings, NOTES_FOR_EXCLUDED = notes_for_excluded, $
             NUM_DIREC_FILES=num_direc_files


num_analyzed = N_ELEMENTS(analysis_frames)

;;************************ SKIP IF NO EXPECTED FRAMES *************************;; 
IF no_exp_frames EQ 1 THEN BEGIN
   GOTO, jump
ENDIF

;;************************ SKIP IF NO MATCHING FRAMES *************************;;
;when no matching frames we wont skip the entire program. This way an
;expected frames list can still be created
IF no_match_frames EQ 1 THEN BEGIN
   GOTO, moveto
ENDIF


;;******************************** EXTRACT OBS DATE **************************************;;

years = INTARR(num_analyzed)
months = INTARR(num_analyzed)
days = INTARR(num_analyzed)
temp_years = 0
temp_months = 0
temp_days = 0
junk1 = ''
junk2 = ''
FOR i=0, num_analyzed-1 DO BEGIN
   READS, analysis_dates[i], temp_years, junk1, temp_months, junk2, temp_days, FORMAT='(F4, A1, F2, A1, F2)'

   years[i] = temp_years
   months[i] = temp_months
   days[i] = temp_days
ENDFOR

;;convert each files date to Julian date for ease when sorting files         
jul_dates = LONARR(num_analyzed)
FOR i=0, num_analyzed-1 DO BEGIN
   jul_dates[i] = JULDAY(months[i], days[i], years[i])
ENDFOR

sorted_jdates = jul_dates[SORT(jul_dates)]
sorted_ind = SORT(jul_dates)

;;use indecies from sorted julian date to places files sorted into a            
;;new variable                                                                               
sorted_files = STRARR(num_analyzed)
j = 0
FOR i=0, num_analyzed-1 DO BEGIN
   sorted_files[i] = analysis_frames[sorted_ind[i]]
ENDFOR




;;******************************** CALL SPECTRUM_ANALYZER.PRO *******************************;;

IF ~mult_stars THEN BEGIN
   PRINT, ""
   PRINT, "Begin analysis:  "
   PRINT, ""
ENDIF

;; define some variables to hold information needed to generate a star report
j = 0
EWs = DBLARR(num_analyzed)
EWs_error = DBLARR(num_analyzed)
waves = DBLARR(512 , num_analyzed)
fluxs = DBLARR(512 , num_analyzed)
sing_peak_waves = DBLARR(num_analyzed)
sing_peak_fluxs = DBLARR(num_analyzed)
peak_one_waves = DBLARR(num_analyzed)
peak_two_waves = DBLARR(num_analyzed)
peak_one_fluxs = DBLARR(num_analyzed)
peak_two_fluxs = DBLARR(num_analyzed)
vr_ratios = DBLARR(num_analyzed)
shell_paras = DBLARR(num_analyzed)
symmetries = DBLARR(num_analyzed)
snrs = DBLARR(num_analyzed)
peak_nums = INTARR(num_analyzed)
centers = DBLARR(num_analyzed)
fwhms = DBLARR(num_analyzed)
flags = INTARR(3,num_analyzed)
sp_date_ind = INTARR(num_analyzed)
bp_date_ind = INTARR(num_analyzed)
rp_date_ind = INTARR(num_analyzed)


;; use counters to know how many double/single peak spectra there are
num_sing = 0
num_doub = 0

FOR i=0, num_analyzed-1 DO BEGIN
   analysis_logs = [analysis_logs, '']
   analysis_logs = [analysis_logs, '=============================================================']
   analysis_logs = [analysis_logs,'File processed: '+sorted_files[i]+' ('+STRCOMPRESS(STRING(i+1))+ $
                    ' of' + STRCOMPRESS(STRING(num_analyzed))+')']

   temp_file = sorted_files[i]

   
   ;; CALL PROCEDURE
   spectrum_analyzer, temp_file
   

   ;; SAVE DATA TO ARRAYS
   EWs[j] = EW
   EWs_error[j] = ew_error
   waves[0,j] = wave
   fluxs[0,j] = flux 
   symmetries[j] = symmetry
   snrs[j] = snr
   peak_nums[j] = peak_num
   centers[j] = center
   flags[0,j] = flag
   sp_date_ind[j] = -1
   bp_date_ind[j] = -1
   rp_date_ind[j] = -1

   IF peak_num EQ 1 THEN BEGIN
      sing_peak_waves[j] = peak_wave
      sing_peak_fluxs[j] = peak_flux
      sp_date_ind[j] = i
      fwhms[j] = fwhm
      num_sing += 1
   ENDIF ELSE IF peak_num EQ 2 THEN BEGIN
      peak_one_waves[j] = peak1_wave
      peak_one_fluxs[j] = peak1_flux
      bp_date_ind[j] = i
      peak_two_waves[j] = peak2_wave
      peak_two_fluxs[j] = peak2_flux
      rp_date_ind[j] = i
      vr_ratios[j] = vr_ratio
      shell_paras[i] = shell_para
      num_doub += 1
   ENDIF ELSE BEGIN
      analysis_logs = [analysis_logs, '**NOTICE: NO PEAK DATA**']
      analysis_logs = [analysis_logs, '']
   ENDELSE

   j += 1


   ;;PRINT OUT STATUS REPORT:
   iteration = '      '+STRING(STRCOMPRESS(i+1, /remove_all))
   print_length = STRCOMPRESS(STRING($
                  STRLEN(" Spectra Analyzed for '"+star_proper_name+"' ("+star_name+')') + $
                  STRLEN(STRING(STRCOMPRESS(STRLEN(iteration+1), /remove_all)))$
                  ),/remove_all)


   PRINT, STRING(STRCOMPRESS(num_analyzed, /remove_all))+" Spectra Analyzed for '"+$
          star_proper_name+"' ("+star_name+')', $                            
          FORMAT = '($, "'+iteration+' / " ,A'+print_length+',"'+STRING(13B)+'")' 

ENDFOR

   PRINT,'';empty print statement to close the open line


;; delete READ_FILES.txt
SPAWN, 'rm -f READ_FILES.txt'




;;************************ CALL PROCEDURE stability_test.pro ************************;;
;stop
stability_test, EWs, EWs_error, sorted_jdates



;stop
;;********************** CALL PROCEDURE halpha_centralizer.pro **********************;; 

halpha_centralizer, waves, centers, SHIFTS=shifts




;;*********************** CALL PROCEDURE cross_correlator.pro ***********************;;

cross_correlator, waves, fluxs, num_analyzed, shifts




;;******************************* REPORTING TO TERMINAL *****************************;;

;;Report to user their numbers here so it comes after the processing
IF ~mult_stars THEN BEGIN
   PRINT, ''
   PRINT, ''
   PRINT, "   Star's Proper Name: "+ star_proper_name
   PRINT, "   Star's HD number: ", STRCOMPRESS(STRING(star_hd))
   PRINT, "   Star's HR number: ", STRCOMPRESS(STRING(star_hr))
   PRINT, ''
ENDIF

num_direc_files = STRCOMPRESS(STRING(num_direc_files), /REMOVE_ALL)
IF frames_excluded[0] EQ -1 THEN BEGIN;if there are no excluded frames
   num_expected = STRCOMPRESS(STRING(N_ELEMENTS(expected_frames)), /REMOVE_ALL)
   num_excluded = '0'
ENDIF ELSE BEGIN;if there are frames to be excluded
   num_expected = STRCOMPRESS(STRING(N_ELEMENTS(expected_frames)), /REMOVE_ALL)
   num_excluded = STRCOMPRESS(STRING(N_ELEMENTS(frames_excluded)), /REMOVE_ALL)
ENDELSE

IF ~mult_stars THEN BEGIN
   PRINT, '   '+num_expected + ' frames expected as per the stellar_obslog.txt.' 
   PRINT, '   '+ STRCOMPRESS(STRING(num_analyzed), /REMOVE_ALL) + ' out of ' + num_direc_files + $
          ' frames were analyzed from the Spectra directory.'
   PRINT, '   '+num_excluded + ' frames excluded as per excluded_frames.txt.'
   PRINT, '' 

ENDIF




;;****************** COMPARING USER INPUTTED STAR TO STAR # IN THE HEADER *****************;;

FOR i=0, N_ELEMENTS(header_target_name)-1 DO BEGIN
   IF header_target_name[i] EQ star_hd || $
      STRLOWCASE(header_target_name[i]) EQ STRCOMPRESS(STRLOWCASE(star_proper_name), /REMOVE_ALL) || $
      STRCOMPRESS(STRLOWCASE(header_target_name[i]), /REMOVE_ALL) EQ 'hr'+star_hr || $
      header_target_name[i] EQ ''THEN BEGIN;nothings wrong so dont do anything
   ENDIF ELSE BEGIN
      analysis_logs = [analysis_logs, '']
      analysis_logs = [analysis_logs, '**WARNING:', analysis_frames[i], $
             ' POSSIBLE INCORRECT FILE FOR REQUESTED STAR. CHECK HEADER IN FILE.**']
      analysis_logs = [analysis_logs, '']

      PRINT, '**WARNING:', analysis_frames[i], $
             ' POSSIBLE INCORRECT FILE FOR REQUESTED STAR. CHECK HEADER IN FILE.**'
   ENDELSE
ENDFOR




;;****************** COMPARING DATE FROM OBSLOG TO DATE FROM HEADER *******************;;

FOR i=0, N_ELEMENTS(analysis_dates)-1 DO BEGIN
   IF STRMID(name_strings[i],0,1) EQ '%' THEN BEGIN;;determines if theres a header
      IF analysis_dates[i] NE file_dates[i] THEN BEGIN;;compares obslog and header
         analysis_logs = [analysis_logs, '']
         analysis_logs = [analysis_logs, '**WARNING:', analysis_frames[i], $
                ' HAS AN INCONSISTENCY IN OBSERVATION DATE BETWEEN FILE AND STELLAR OBS LOG**']
         analysis_logs = [analysis_logs, '']
         
         PRINT, '**WARNING:', analysis_frames[i], $
                ' HAS AN INCONSISTENCY IN OBSERVATION DATE BETWEEN FILE AND STELLAR OBS LOG**'
      ENDIF
   ENDIF
ENDFOR


moveto:;;move here if no matching frames were found. this way expected frames can still be made 
;;*********************************** CREATE LOG FOLDER ************************************;;

IF ~FILE_TEST('Logs', /DIRECTORY) THEN SPAWN, 'mkdir Logs'

;;At this point the logs directory must exist now so we can CD into it
CD, 'Logs'




;;************************ CREATE LOG FILE FOR EXPECTED STAR FRAMES ************************;;

OPENW, new_file_unit, 'expected_frames.log', /GET_LUN

PRINTF, new_file_unit, '% ' +systime(0)+': expected_frames.log created'
PRINTF, new_file_unit, '% '
PRINTF, new_file_unit, '% ' + star_proper_name
PRINTF, new_file_unit, '% ' + star_name
PRINTF, new_file_unit, '% NO. OF FRAMES:', STRCOMPRESS(N_ELEMENTS(expected_frames))
PRINTF, new_file_unit, 'Frame   ',  'UT Date      ', 'Observer'
PRINTF, new_file_unit, '-----   ',  '----------   ', '----------'

FOR i=0, N_ELEMENTS(expected_frames)-1 DO BEGIN
   PRINTF, new_file_unit, expected_frames[i], expected_dates[i], $
           frame_observer[i], FORMAT='(A5, A13, A13)'
ENDFOR

PRINTF, new_file_unit,''
PRINTF, new_file_unit,'----------- List of Input Files for process_sss.pro -----------'
PRINTF, new_file_unit,''

FOR i=0, N_ELEMENTS(expected_frames)-1 DO BEGIN
   PRINTF, new_file_unit, expected_frames[i]+'.nrm'
ENDFOR

PRINTF, new_file_unit,''
PRINTF, new_file_unit,'----------- List of Frames -----------'
PRINTF, new_file_unit,''

FOR i=0, N_ELEMENTS(expected_frames)-1 DO BEGIN
   PRINTF, new_file_unit, expected_frames[i]
ENDFOR

CLOSE, new_file_unit
FREE_LUN, new_file_unit

;copy this log file and place the copy in the stars' directory
SPAWN, '\cp -r expected_frames.log '+root_direc+targets_dir_name+'/HD'+star_hd+'/Logs'

;;The following if statement is used to skip the creation of the rest
;;of the log files. If no matching frames were found but there were
;;expected frames then only the expected frames log file should be made.
IF no_match_frames EQ 1 THEN BEGIN
   CD, '..';move back to working frame
   abort = 1
   GOTO, jump
ENDIF




;;************************ CREATE LOG FILEFOR STAR FRAMES AND DATES ************************;;

OPENW, new_file_unit, 'analyzed_frames.txt', /GET_LUN
PRINTF, new_file_unit, '% ' +systime(0)+': analyzed_frames.txt created'
PRINTF, new_file_unit, '% '
PRINTF, new_file_unit, '% ' + star_proper_name
PRINTF, new_file_unit, '% ' + star_name
PRINTF, new_file_unit, '% NO. OF FRAMES:', STRCOMPRESS(N_ELEMENTS(analysis_frames)) 
PRINTF, new_file_unit, 'Frame   ', 'Obs. Date    ', 'EW       ', 'SNR     ', 'Flags' 
PRINTF, new_file_unit, '-----   ', '----------   ', '------   ', '-----   ', '-----'

FOR i=0, num_analyzed-1 DO BEGIN
   PRINTF, new_file_unit, analysis_frames[i], analysis_dates[i], EWs[i], snrs[i], $
           flags[*, i], FORMAT='(A5, A13, F9.3, F8.1, I4, I2, I2)'
ENDFOR

CLOSE, new_file_unit
FREE_LUN, new_file_unit

;copy this log file and place the copy in the stars' directory
SPAWN, '\cp -r analyzed_frames.txt '+root_direc+targets_dir_name+'/HD'+star_hd+'/Logs'




;;************************ CREATE LOG FILE FOR MISSING FRAMES ************************;;

OPENW, new_file_unit, 'missing_frames.log', /GET_LUN

PRINTF, new_file_unit, '% ' +systime(0)+': missing_frames.log created'
PRINTF, new_file_unit, '% '
PRINTF, new_file_unit, '% ' + star_proper_name
PRINTF, new_file_unit, '% ' + star_name

;;a bunch of if statements to determine what to print according to the
;;existence of missing frames or not 
mf_type = ISA(missing_frames, /NUMBER) ;; if its a number theres no missing frames
IF mf_type THEN BEGIN ;;i.e. no missing frames
   PRINTF, new_file_unit, '% NO. OF FRAMES: 0'
ENDIF ELSE BEGIN
   PRINTF, new_file_unit, '% NO. OF FRAMES:', STRCOMPRESS(N_ELEMENTS(missing_frames))
ENDELSE


PRINTF, new_file_unit, 'Frame   ',  'UT Date      ', 'Observer'
PRINTF, new_file_unit, '-----   ',  '----------   ', '-----------'


IF mf_type THEN BEGIN ;;i.e. no missing frames
   FOR i=0, N_ELEMENTS(missing_frames)-1 DO BEGIN
      PRINTF, new_file_unit, '*** NO MISSING FRAMES ***'
   ENDFOR
ENDIF ELSE BEGIN
   FOR i=0, N_ELEMENTS(missing_frames)-1 DO BEGIN
      PRINTF, new_file_unit, missing_frames[i], missing_dates[i], $
              missing_frame_obs[i], FORMAT='(A5, A13, A14)'
   ENDFOR
ENDELSE

CLOSE, new_file_unit
FREE_LUN, new_file_unit

;copy this log file and place the copy in the stars' directory
SPAWN, '\cp -r missing_frames.log '+root_direc+targets_dir_name+'/HD'+star_hd+'/Logs'




;;************************ CREATE LOG FILE FOR EXCLUDED FRAMES ************************;;

OPENW, new_file_unit, 'frames_excluded.log', /GET_LUN
PRINTF, new_file_unit, '% ' +systime(0)+': frames_excluded.txt created'
PRINTF, new_file_unit, '% '
PRINTF, new_file_unit, '% ' + star_proper_name
PRINTF, new_file_unit, '% ' + star_name

IF frames_excluded[0] EQ -1 THEN BEGIN
   PRINTF, new_file_unit, '% NO. OF FRAMES: 0'
ENDIF ELSE BEGIN
   PRINTF, new_file_unit, '% NO. OF FRAMES:', STRCOMPRESS(N_ELEMENTS(frames_excluded))
ENDELSE

PRINTF, new_file_unit, 'Frame   ', 'Obs. Date    ', 'Notes:'
PRINTF, new_file_unit, '-----   ', '----------   ', '----->'

IF frames_excluded[0] EQ -1 THEN BEGIN
PRINTF, new_file_unit, '*** NO FRAMES EXCLUDED ***'
ENDIF ELSE BEGIN
   FOR i=0, N_ELEMENTS(frames_excluded)-1 DO BEGIN
      str_len = STRLEN(notes_for_excluded[i])+3
      str_len = STRCOMPRESS(STRING(str_len), /REMOVE_ALL)
      PRINTF, new_file_unit, frames_excluded[i], frame_dates_excluded[i],$
              notes_for_excluded[i], FORMAT='(A5, A13, A'+str_len+')'
   ENDFOR
ENDELSE

CLOSE, new_file_unit
FREE_LUN, new_file_unit

;copy this log file and place the copy in the stars' directory
SPAWN, '\cp -r frames_excluded.log '+root_direc+targets_dir_name+'/HD'+star_hd+'/Logs'




;;************************ CREATE LOG FILE FOR SPECTRAL ANALYSIS INFO ************************;;

OPENW, new_file_unit, 'analysis_log.log', /GET_LUN
PRINTF, new_file_unit, '% ' +systime(0)+': analysis_logs.log created'
PRINTF, new_file_unit, '% '
PRINTF, new_file_unit, '% ' + star_proper_name
PRINTF, new_file_unit, '% ' + star_name

FOR i=0, N_ELEMENTS(analysis_logs)-1 DO BEGIN
   PRINTF, new_file_unit, analysis_logs[i]
ENDFOR
CLOSE, new_file_unit
FREE_LUN, new_file_unit

;copy this log file and place the copy in the stars' directory
SPAWN, '\cp -r analysis_log.log '+root_direc+targets_dir_name+'/HD'+star_hd+'/Logs'

CD, '..' ;;move back into the working directory 

jump: ;; jumps here if no matches are found for HR/HD number
END
