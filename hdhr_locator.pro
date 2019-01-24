;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; hdhr_locator.pro
;;
;; This procedure is meant to prompt the user to enter the HD or HR
;; number for the star that the user wishes to be analyzed. The
;; beginning of this procedure is dedicated to propmting the user and
;; giving them useful feedback while doing it. hdhr_locator.pro will
;; then store all the HD and HR numbers and proper names in respective
;; variables from a file titled master_list.txt. Then depending on whether the user
;; entered an HD or HR number it will find it in the arrays and obtain
;; the stars corresponding classifications. In addition to these
;; various classifications being pulled out, the stars RA and DEC is
;; obtained. The important variables defined here are: 'star_hd',
;; 'star_hr', 'star_proper_name', 'star_ra', and 'star_dec'.
;;
;; In order for this pipeline to work there must be a file titled
;; master_list.txt that will contain the HD and HR numbers, proper
;; name, RA and Dec, and the magnitude. The pipeline will run as long
;; as it at least has the HD. If you have, say, two different lists
;; of targets and you want them all to be in this master_list.txt then
;; you can easily combine them into one with a command in the
;; terminal. Say the two lists are called primary.txt and
;; secondary.txt. To combine them into one list titled
;; master_list.txt do the following: "cat primary.txt secondary.txt >
;; master_list.txt" And as long as there are no empty lines,
;; and all lines with comments begin with a "%" then this will be fine
;; as is, noting that the stars in the master_list.txt do not need to
;; be in any specific order.
;;
;; The last section of this procedure has to do with printing to the
;; terminal if the user executed hdhr_locator.pro by itself.   
;;
;; Author: W. Glennon Fagan
;;======================================================================;;


PRO hdhr_locator, full_pipeline
RESOLVE_ALL, /QUIET
COMMON abort_block, no_exp_frames, no_hdhr_num, no_match_frames, user_entry
COMMON star_block, analysis_logs, root_direc, available_dates, available_frames, analysis_dates, analysis_frames, mult_stars, star_dec, star_hd, star_hr, star_proper_name, star_name, star_ra, star_vmag, start_time, sub_direc, targets_dir_name, trend, ver_num

;;if star_name is defined that means that the pipeline is being ran
;;for multiple targets (loop_the_loop) and therefore
IF ISA(star_name) THEN user_entry = star_name
IF ISA(star_name) THEN prompt_user = 0 ELSE prompt_user = 1

;;if ver_num is not defined then that means that hdhr_locator.pro is
;;being executed alone and these variables will be used to tell the
;;procedure that
IF ~ISA(ver_num) THEN hdhr_locator_alone = 1 ELSE hdhr_locator_alone = 0
IF ~ISA(ver_num) THEN mult_stars = 0




;;******************************** OBTAIN STAR NAME ***********************************;;
  
;;ask for user input to define star name
star_spec = ''
star_num = 0L
run_loop = 1 ;variable that is used to determine whether or not to repeat the following loop
WHILE run_loop EQ 1 DO BEGIN

   ;print statements to notify the user of the directory in
   ;which the star's data will be stored in
   IF ~hdhr_locator_alone THEN BEGIN ;if this is true then hdhr_locator.pro is being executed alone 
      IF prompt_user THEN BEGIN      ;if prompt user is equal to one then targets loop isn't being executed
         PRINT, ''
         PRINT, 'Your data will be stored in the directory: ', targets_dir_name
         PRINT, "If this is incorrect please type 'N' below"
      ENDIF
   ENDIF

   ;;if statement to prompt user to enter star name if star_report is not
   ;;being executed for multiple targets
   IF ~prompt_user THEN star_name = star_name ELSE BEGIN
      star_name = ''
      READ, star_name, PROMPT='Enter star`s HD or HR number (HD/HR ******): '
      user_entry = star_name ;;set the entry to a seperate variable so 'star_name' can be 
                             ;;altered to both HD/HR number in the beginning of coordinator.pro
   ENDELSE

   ;if the user enter 'N' then prompt the user to enter a new 
   ;directory in which to store the star's data
   IF user_entry EQ 'N' || user_entry EQ 'n' THEN BEGIN
      READ, targets_dir_name, PROMPT='Enter the new title of the directory in which to store star data: '
      run_loop = 1
      GOTO, END_FOR
   ENDIF
   
   
   star_name = STRLOWCASE(STRCOMPRESS(star_name, /REMOVE_ALL)) ;;used to get rid of all spaces and fix case
   strlen = STRLEN(star_name)                      ;;find length of string
   
   ;;extract the designator and number and store in individual variables
   star_spec = STRMID(star_name, 0, 2)
   star_num = STRMID(star_name, 2, strlen)
   
   
   ;;************************** CHECK FOR BAD ENTRIES ***************************;;
   
   ;;Begin with making sure the hd/hr designator was entered properly
   run_loop = 0 ;;set equal to 0 so the loop wont be ran again unless there was a bad entry
   IF star_spec EQ 'hr' || star_spec EQ 'hd' THEN BEGIN
   ENDIF ELSE BEGIN
      run_loop = 1 ;;set equal to 1 if the entry had a bad HD/HR designator
      PRINT, '**ERROR: HD/HR DESIGNATION NOT SPECIFIED PROPERLY**'
   ENDELSE

   ;;Check that the length of the number doesn't exceed 6 digits 
   IF STRLEN(star_num) GT 6 THEN BEGIN
      run_loop = 1 ;;set equal to 1 if the entry had a bad HD/HR designator
      PRINT, '**ERROR: HD/HR NUMBER TOO LONG. HD: ≥6 & HR: ≥4**'
   ENDIF
   
   ;case statement that is used to check if entry had included letters to the hd/hr number
   bad_entry = 0
   FOR i=0, STRLEN(star_num)-1 DO BEGIN
      CASE STRMID(star_num, i, 1) OF
         '0':
         '1':
         '2':
         '3':
         '4':
         '5':
         '6':
         '7':
         '8':
         '9':
         ELSE: BEGIN
            run_loop = 1 ;;set equal to 1 if the entry had a bad HD/HR designator
            PRINT, '**ERROR: LETTERS WERE INCLUDED IN THE HD/HR NUMBER**'
            GOTO, end_for;Skip to end of for loop to avoid mult. printing
         END
      ENDCASE
   ENDFOR
   end_for:
   
   IF ~mult_stars THEN BEGIN
      PRINT, ''
   ENDIF
ENDWHILE




;;************************** SET VARIABLE TOMEASURERUN TIME ***************************;;
;;note this is done here so that the time to enter an HD/HR
;;number isn't factored into the run time.
start_time = systime(1)
;;******************************************************* ***************************;;




;;************************** STORE ALL HD/HR NUM IN ARRAY ***************************;;

check_if_ascii = QUERY_ASCII('master_list.txt', ascii)
number_of_lines = ascii.lines ;variable containg number of lines in file

;; read in the entire ASCII file into a string array 
ascii_lines = STRARR(number_of_lines)
tmp_string = ''
OPENR, file_unit, 'master_list.txt', /GET_LUN
counter = 0L

WHILE NOT EOF(file_unit) DO BEGIN
   ;; read the entire line as a string    
   ;; emacs star_finder.pro          
   READF, file_unit, tmp_string
   ascii_lines[counter] = tmp_string
   counter += 1
ENDWHILE
CLOSE, file_unit
FREE_LUN, file_unit

 
;;define temp varibles to be used when reading out values from each line
tmp_hd = ''
tmp_hr = ''
tmp_proper_name = ''

tmp_ra_h = ''
tmp_ra_d = ''
tmp_ra_s = ''
tmp_dec_d = ''
tmp_dec_m = ''
tmp_dec_s = ''

tmp_v_mags = ''

counter = 0;;counter must be long

;;define string arrays to hold hd and hr 
hd_str = STRARR(N_ELEMENTS(ascii_lines)-N_ELEMENTS(WHERE(STRMID(ascii_lines,0,1) EQ '%')))
hr_str = STRARR(N_ELEMENTS(ascii_lines)-N_ELEMENTS(WHERE(STRMID(ascii_lines,0,1) EQ '%')))
proper_name = STRARR(N_ELEMENTS(ascii_lines)-N_ELEMENTS(WHERE(STRMID(ascii_lines,0,1) EQ '%')))

ra_list = STRARR(N_ELEMENTS(ascii_lines)-N_ELEMENTS(WHERE(STRMID(ascii_lines,0,1) EQ '%')))
dec_list = STRARR(N_ELEMENTS(ascii_lines)-N_ELEMENTS(WHERE(STRMID(ascii_lines,0,1) EQ '%')))

v_mags = STRARR(N_ELEMENTS(ascii_lines)-N_ELEMENTS(WHERE(STRMID(ascii_lines,0,1) EQ '%')))

FOR i=0, N_ELEMENTS(ascii_lines)-1 DO BEGIN
   ;;the following if statement is used to ignore lines that begin
   ;;with a comment as they are used to signify a comment
   IF STRMID(ascii_lines[i], 0, 1) NE '%' THEN BEGIN
      READS, ascii_lines[i], tmp_proper_name, tmp_hr, tmp_hd, tmp_ra_h, tmp_ra_d, tmp_ra_s, $
             tmp_dec_d, tmp_dec_m, tmp_dec_s, tmp_v_mags, FORMAT='(A12, A7, A9, A5, A3, A3, A6, A3, A3, A6)'
      counter += 1
      
      
      ;handle RAs and DECs
      tmp_ra_s =  STRCOMPRESS(STRING(ROUND(DOUBLE(tmp_ra_s))), /REMOVE_ALL)
      tmp_dec_s = STRCOMPRESS(STRING(ROUND(DOUBLE(tmp_dec_s))), /REMOVE_ALL)
      
      ;;if the seconds and arcseconds are only one character long then
      ;;put a '0' before it
      IF STRLEN(tmp_ra_s) THEN tmp_ra_s = '0'+tmp_ra_s
      IF STRLEN(tmp_dec_s) THEN tmp_dec_s ='0'+tmp_dec_s
      
      ;;define the degree sign to be used in the RA and DEC Strings
      deg = '!9%!X'
      
      ;;Create a complete string array of the RA and DEC
      ra_list[counter-1] = tmp_ra_h + 'h' + tmp_ra_d + 'm ' + tmp_ra_s + 's'
      dec_list[counter-1] = tmp_dec_d + deg + tmp_dec_m + "' " + tmp_dec_s + '"'
      
      
      ;handle hd/hr numbers
      hd_str[counter-1] = tmp_hd
      hr_str[counter-1] = tmp_hr
      proper_name[counter-1] = tmp_proper_name
      v_mags[counter-1] = tmp_v_mags
      
   ENDIF
   
   ;;set all temp variables back to an empty string
   tmp_proper_name = ''
   tmp_hd = ''
   tmp_hr = ''
   junk = ''
   tmp_v_mags = ''
ENDFOR


;;trim of any white space on either end of the hd and hr number
hd = STRTRIM(hd_str, 2)
hr = STRTRIM(hr_str, 2)


;;search for corresponding  hr or hd number, while 
star_hd = '';will contain the official HD for the entry 
star_hr = '' ;will contain the official HR for the entry
IF star_spec EQ 'hd' THEN BEGIN
   star_hd = star_num
   star_hr = hr[WHERE(hd EQ star_num)]
   star_hr = star_hr[0];;this is to take care of possible duplicates
   star_proper_name = proper_name[WHERE(hd EQ star_num)]
   star_proper_name = star_proper_name[0]

   star_ra = ra_list[WHERE(hd EQ star_num)]
   star_ra = star_ra[0]
   star_dec = dec_list[WHERE(hd EQ star_num)]
   star_dec = star_dec[0]

   star_vmag = v_mags[WHERE(hd EQ star_num)]
   star_vmag = star_vmag[0]
ENDIF ELSE IF star_spec EQ 'hr' THEN BEGIN
   star_hr = star_num
   star_hd = hd[WHERE(hr EQ star_num)]
   star_hd = star_hd[0];;this is to take care of possible duplicates
   star_proper_name = proper_name[WHERE(hr EQ star_num)]
   star_proper_name = star_proper_name[0]

   star_ra = ra_list[WHERE(hr EQ star_num)]
   star_ra = star_ra[0]
   star_dec = dec_list[WHERE(hr EQ star_num)]
   star_dec = star_dec[0]

   star_vmag = v_mags[WHERE(hr EQ star_num)]
   star_vmag = star_vmag[0]
ENDIF



;;******** IF NO HR NUMBER star_hr = '----' ********;;
IF star_hr EQ '' THEN BEGIN
   star_hr = '----'
ENDIF


;;******** IF NO V MAG VALUE star_vmag = 'N/A' ********;;
IF star_vmag EQ '    ' THEN BEGIN
   star_vmag = 'N/A'
ENDIF ELSE BEGIN ;remove white space in beginning
   star_vmag = STRTRIM(star_vmag,1)
ENDELSE


;;********IF NO PROPER NAME star_proper_name = '-----'********;

IF star_proper_name EQ '              ' THEN BEGIN;no proper name in catalog 
   star_proper_name = '-----'
ENDIF ELSE BEGIN
   star_proper_name = STRTRIM(star_proper_name, 2)
ENDELSE




;;************************** SKIP PROGRAM IF NUMBER NOT IN CATALOG ***************************;;

;;if statement that is used to skip program if no hd number is found
;;in the catalog. Because every star should have an HD number,
;;if it's not found it is most likely not in the
;;catalog                                                      

no_hdhr_num = 0
test = 0
test = WHERE(hd EQ star_hd)
IF test[0] EQ -1 || star_hd EQ 0 THEN BEGIN
trend = 'N/A'
no_exp_frames = 0 ;; set other abort variables to 0 so we can skip program and still 
no_match_frames = 0 ;; print abort messages at the end of star_report.pro
no_hdhr_num = 1 ;;set to 1 because the entry was not in the catalog
ENDIF 




;;************************** CODE WHEN EXECUTING HDHR_LOCATOR ALONE ***************************;; 

;;Statement that will print the HD/HR number for a user requested
;;HD/HR number if hdhr_locator or spectrum_analyzer was executed individually
IF ISA(full_pipeline) EQ 0 THEN BEGIN
   IF test[0] EQ -1 || star_hd EQ 0 THEN BEGIN;when the entry isn't in the catalog
      PRINT, ''
      PRINT, '     YOU ENTERED: ' +  user_entry
      PRINT, '     **WARNING: NO HD/HR NUMBER FOUND. STAR MISSING FROM CATALOG**'
      PRINT, ''
      PRINT, ''
   ENDIF ELSE BEGIN;when the entry is in the catalog
      PRINT, ''
      PRINT, 'Corresponding HD/HR numbers:'
      PRINT, ''
      IF star_hr EQ 0 THEN BEGIN
         PRINT, 'Proper Name: ' + star_proper_name
         PRINT, "     HD "+STRCOMPRESS(STRING(star_hd), /REMOVE_ALL)
         PRINT, "     HR ----"
      ENDIF ELSE BEGIN
         PRINT,'Proper Name: '+ star_proper_name
         PRINT, "     HD "+STRCOMPRESS(STRING(star_hd), /REMOVE_ALL)
         PRINT, "     HR "+STRCOMPRESS(STRING(star_hr), /REMOVE_ALL)
      ENDELSE
      PRINT, ''
      PRINT, ''
   ENDELSE
ENDIF


;; define full name string
   star_name = 'HD '+ STRCOMPRESS(STRING(star_hd), /REMOVE_ALL) + ' / HR ' $
               + STRCOMPRESS(STRING(star_hr), /REMOVE_ALL)

END
