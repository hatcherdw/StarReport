;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; targets_loop.pro
;;
;; This is a short and simple procedure whose purpose is to loop
;; star_report.pro for multiple targets. The array titled star_list
;; is the list of targets to be analyzed.  When
;; star_report.pro is being ran in a loop like this, the on screen
;; version of the star report will not be generated. Do not fret
;; though; a copy of each star's updated star report can be
;; found in the star's directory titled HD[star's HD number]_report.ps.
;;
;; Author: W. Glennon Fagan 
;;======================================================================;;


PRO targets_loop
RESOLVE_ALL, /QUIET
COMMON abort_block, no_exp_frames, no_hdhr_num, no_match_frames, user_entry
COMMON star_block, analysis_logs, root_direc, available_dates, available_frames, analysis_dates, analysis_frames, mult_stars, star_dec, star_hd, star_hr, star_proper_name, star_name, star_ra, star_vmag, start_time, sub_direc, targets_dir_name, trend, ver_num

begin_time = systime(1)

;;targets list to run
working_list = 'primary_targets.txt'
;working_list = 'secondary_targets.txt'
;working_list = 'cal_targets.txt'

;; Define directory to store all the stars' data
targets_dir_name = 'Primary_Targets' 
;targets_dir_name = 'Secondary_Targets'
;targets_dir_name = 'Calibration_Collaboration_Targets'

user_answer = ''
WHILE user_answer NE 'y' && user_answer NE 'Y' DO BEGIN
   PRINT, ''
   PRINT, 'The file that contains the list of stars to analyze: ', working_list
   PRINT, 'Your data will be stored in the directory:           ', targets_dir_name
   
   READ, user_answer, PROMPT="If this is incorrect please type 'n' below and 'y' if correct: "
   PRINT, ''
   
   IF user_answer EQ 'N' || user_answer EQ 'n' THEN BEGIN
      temp_file_name = ''
      READ, temp_file_name, PROMPT=$
            "Enter the new name of the file containing a list of stars or 'n' to not change: "
      IF temp_file_name NE 'N' && temp_file_name NE 'n' THEN BEGIN
         working_list = temp_file_name
      ENDIF

      temp_dir_name = ''
      READ, temp_dir_name, PROMPT=$
            "Enter the new title of the directory in which to store star data or 'n' to not change: "
      IF temp_dir_name NE 'N' && temp_dir_name NE 'n' THEN BEGIN
         targets_dir_name = temp_dir_name
      ENDIF
   ENDIF
ENDWHILE




;; ********************************* Read working_list into array **************************************;;

check_if_ascii = QUERY_ASCII(working_list, ascii)
number_of_lines = ascii.lines ;store the number of lines in a variable

;; read in the entire ASCII file into a string array
file_lines = STRARR(number_of_lines);String array that will contain each line of file as an element
tmp_string = '';tmp string array to be used with READF command
OPENR, file_unit, working_list, /GET_LUN;open the working_list file in file_unit
counter = 0L ;define a counter set to a Long Integer

WHILE NOT EOF(file_unit) DO BEGIN
   ;; read the entire line as a string
   READF, file_unit, tmp_string
   file_lines[counter] = tmp_string
   counter += 1
ENDWHILE
CLOSE, file_unit ;cloe the file_unit and then free the LUN
FREE_LUN, file_unit

star_list = STRARR(N_ELEMENTS(file_lines)-1);define array that will hold list of stars for analysis
temp_line = ''
junk = ''
FOR i=0, N_ELEMENTS(file_lines)-1 DO BEGIN;;starts at 1 to ignore first line containing column titles
   IF STRMID(file_lines[i], 0, 1) NE '%' THEN BEGIN   
      READS, file_lines[i], junk, temp_line, FORMAT='(A22, A6)'
      temp_line = STRCOMPRESS(temp_line, /REMOVE_ALL)
      star_list[i-1] = 'HD' + temp_line
      temp_line = ''
   ENDIF
ENDFOR

star_list = star_list[WHERE(star_list NE '')]


;define variables for use in a log file to hold stars with no frames
;or not in the catalog
targets_list_log_pn = STRARR(N_ELEMENTS(star_list))
targets_list_log_hd = STRARR(N_ELEMENTS(star_list))
targets_list_log_hr = STRARR(N_ELEMENTS(star_list))
targets_list_log_issue = STRARR(N_ELEMENTS(star_list))



;*********************************
s_counter = 0 ;stable trend counter
v_counter = 0 ;variable trend counter
dl_counter = 0 ;exhibiting disk loss counter
dg_counter = 0 ;exhibiting disk growth counter
na_counter = 0 ;too few spectra to assess counter
;********************************* 


FOR i=0, N_ELEMENTS(star_list)-1 DO BEGIN;;starts at 1 to ignore first line containing column titles 
num_stars = N_ELEMENTS(star_list)
print, '   Analyzing star '+STRING(STRCOMPRESS(i+1, /remove_all))+' / '+$
       STRING(STRCOMPRESS(num_stars, /remove_all))

   star_report, star_list[i]

   IF no_hdhr_num EQ 1 || no_exp_frames EQ 1 || no_match_frames EQ 1 THEN BEGIN
      targets_list_log_pn[i] = star_proper_name 
      targets_list_log_hd[i] = star_hd
      targets_list_log_hr[i] = star_hr
   ENDIF

   IF no_hdhr_num EQ 1 THEN BEGIN
      targets_list_log_issue[i] = 'No HD/HR number' 
   ENDIF ELSE IF no_exp_frames EQ 1 THEN BEGIN
      targets_list_log_issue[i] = 'No expected frames'
   ENDIF ELSE IF no_match_frames EQ 1 THEN BEGIN
      targets_list_log_issue[i] = 'No matching frames'
   ENDIF

   CASE trend OF
      ;depending on the trend increment the counter accordingly
      'Stable': s_counter += 1
      'Variable': v_counter += 1
      'Disk-loss': dl_counter += 1
      'Disk-growth': dg_counter += 1
      ELSE: na_counter +=1
   ENDCASE

print, ''

ENDFOR



;remove empty spaces in the arrays
targets_list_log_pn = targets_list_log_pn[WHERE(targets_list_log_pn NE '')]
targets_list_log_hd = targets_list_log_hd[WHERE(targets_list_log_hd NE '')]
targets_list_log_hr = targets_list_log_hr[WHERE(targets_list_log_hr NE '')]
targets_list_log_issue = targets_list_log_issue[WHERE(targets_list_log_issue NE '')]

CD, 'Logs'
OPENW, new_file_unit, 'targets_list_info.log', /GET_LUN

PRINTF, new_file_unit, '% ' +systime(0)+': Targets_list_info.log created'
PRINTF, new_file_unit, '% '
PRINTF, new_file_unit, 'Proper Name   HD       HR     Issue'
PRINTF, new_file_unit, '-----------   ------   ----   ------------------'

FOR i=0, N_ELEMENTS(targets_list_log_pn)-1 DO BEGIN
   PRINTF, new_file_unit, targets_list_log_pn[i], targets_list_log_hd[i], targets_list_log_hr[i], $
           targets_list_log_issue[i], FORMAT='(A11, A9, A7, A21)'
ENDFOR
CLOSE, new_file_unit
FREE_LUN, new_file_unit
CD, '..'





PRINT, 'Stable:', s_counter
PRINT, 'Variable:', v_counter
PRINT, 'Disk-loss:', dl_counter
PRINT, 'Disk-growth:', dg_counter
PRINT, 'N/A:', na_counter






run_time = SYSTIME(1) - begin_time
minutes = run_time / 60
runtime_minutes = FLOOR(minutes)
runtime_seconds = FLOOR((minutes - runtime_minutes) * 60)

PRINT, ''
PRINT, '===================='+STRCOMPRESS(STRING(N_ELEMENTS(star_list)))+$
       ' STARS ANALYZED ===================='
PRINT, ''
PRINT, 'RUN TIME:', STRCOMPRESS(STRING(runtime_minutes))+'m '+$
       STRCOMPRESS(STRING(runtime_seconds), /REMOVE_ALL)+'s'
PRINT, ''
PRINT, '==========================================================='


END
