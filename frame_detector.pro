;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; frame_detector.pro
;;
;; This procedure is the part of the pipelin that detects all 
;; expected frames for acertain star in a ``.csv" file known as the 
;; stellar observation log. It will look for frames using the HD
;; number, the HR number, and the proper name each
;; individually. Arrays of matching frames will then be concatenated to
;; make one comprehensive list. 
;;
;; Author: W. Glennon Fagan
;;======================================================================;;


PRO frame_detector, FRAME_OBSERVER=frame_observer, EXPECTED_DATES=expected_dates, $
                    EXPECTED_FRAMES=expected_frames
RESOLVE_ALL, /QUIET
COMMON star_block
COMMON abort_block


variables = csv_reader('stellar_obslog_2019-01-08.csv')

check_if_ascii = QUERY_ASCII('stellar_obslog_2019-01-08.csv', ascii)
number_of_lines = ascii.lines


frame_num = STRARR(number_of_lines)
name = STRARR(number_of_lines)
ut_date = STRARR(number_of_lines)
start_ut = STRARR(number_of_lines)
exp_time = STRARR(number_of_lines)
exp_time2 = STRARR(number_of_lines)
obs_year = STRARR(number_of_lines)
t_red = STRARR(number_of_lines)
t_blue = STRARR(number_of_lines)
t_room = STRARR(number_of_lines)
t_coffin = STRARR(number_of_lines)
observer = STRARR(number_of_lines)
notes = STRARR(number_of_lines)

FOR i=0,number_of_lines-1 DO BEGIN  
   ;;store each part of the line in a respective variable
   frame_num[i] = variables[1,i]
   name[i] = variables[2,i]
   name[i] = STRCOMPRESS(name[i], /REMOVE_ALL)
   ut_date[i] = variables[3,i]
   start_ut[i] = variables[4,i]
   exp_time[i] = variables[5,i]
   exp_time2[i] = variables[6,i]
   obs_year[i] = variables[7,i]
   t_red[i] = variables[11,i]
   t_blue[i] = variables[12,i]
   t_room[i] = variables[13,i]
   t_coffin[i] = variables[14,i]
   observer[i] = variables[15,i]
   observer[i] = STRCOMPRESS(observer[i], /REMOVE_ALL)
   notes[i] = variables[16,i]
ENDFOR




;;**************************** SEARCH FOR MATCHING HR FRAMES ****************************;;

;;Search for frames for matching HR number if HR number exists
hr_matches = 0
IF star_hr NE '----' THEN BEGIN
   hr_log_num = name[WHERE(STRMID(name, 0, 2) EQ 'HR')] 
   hr_log_num = STRMID(hr_log_num, 2)               ;;array of stars that used hr in log
   hr_log_ind = WHERE(STRMID(name, 0, 2) EQ 'HR')   ;;array of indices for stars that used hr in log
   match_ind_hr = STRARR(N_ELEMENTS(WHERE(hr_log_num EQ star_hr))) ;;array of indices for stars that matched     user input star
   j = 0
   
   ;;For loop to find all matching frames
   FOR i=0, N_ELEMENTS(hr_log_num)-1 DO BEGIN
      IF hr_log_num[i] EQ star_hr THEN BEGIN
         hr_matches += 1
         match_ind_hr[j] = hr_log_ind[i]
         j++
      ENDIF
   ENDFOR

   
   ;;For loop to use matching indices to fill an array with the frame
   ;;numbers for all matching
   match_frames_hr = STRARR(N_ELEMENTS(match_ind_hr)) ;;array of frame # for matching stars
   IF match_ind_hr[0] EQ '' THEN BEGIN                ;no matches with hr number
      match_frames_hr = ''
   ENDIF ELSE BEGIN
      FOR i=0, N_ELEMENTS(match_ind_hr)-1 DO BEGIN
         match_frames_hr[i] = frame_num[match_ind_hr[i]]
      ENDFOR
   ENDELSE
ENDIF




;;*************************** SEARCH FOR MATCHING HD FRAMES ***************************;;

;;Search for frames for matching HD number
hd_matches = 0 
hd_log_num = name[WHERE(name EQ star_hd)];this variable contains an array of the stars HD number that has a number of elements equal to the amount of HD number entries are present for this star in the obslog 
hd_log_ind = WHERE(name EQ star_hd);this variable contains an array of indices associated with the stellar obslog for HD entries 
match_ind_hd = STRARR(N_ELEMENTS(WHERE(hd_log_num EQ star_hd))) ;array of matching frame's indices
j=0

FOR i=0, N_ELEMENTS(hd_log_num)-1 DO BEGIN
   IF hd_log_num[i] EQ star_hd THEN BEGIN
      hd_matches += 1
      match_ind_hd[j] = hd_log_ind[i]
      j++
   ENDIF
ENDFOR

;;For loop to use matching indices to fill an array with the frame
;;numbers for all matching
match_frames_hd = STRARR(N_ELEMENTS(match_ind_hd))
FOR i=0, N_ELEMENTS(match_ind_hd)-1 DO BEGIN
   match_frames_hd[i] = frame_num[match_ind_hd[i]]
ENDFOR

IF hd_log_ind[0] EQ -1 THEN BEGIN
   match_frames_hd = -1 ;THIS IS IF THERE ARE NO MATCHING FRAMES WITH THE HD NUMBER
ENDIF




;;*************************** SEARCH FOR MATCHING PROPER NAME FRAMES ***************************;; 

;;make two new variables that will set both of these to lower case
;;letters to help ensure all mathces are gathered
lc_targets = STRLOWCASE(name)
lc_pn = STRCOMPRESS(STRLOWCASE(star_proper_name), /REMOVE_ALL)

match_ind_pn = WHERE(lc_targets EQ lc_pn)
IF match_ind_pn[0] EQ -1 THEN BEGIN ;no matches with proper name
   match_frames_pn = ''
   pn_matches = 0
ENDIF ELSE BEGIN
   match_frames_pn = frame_num[WHERE(lc_targets EQ lc_pn)]
   pn_matches = N_ELEMENTS(match_frames_pn)
ENDELSE




;;*************************** MAKE COMPREHENSIVE MATCHING FRAME LIST ***************************;;

;;If statement to determine what and how to print it, determined by
;;the existence of an HR number or if there are no frames at all
no_exp_frames = 0 ;;used to skip pipeline if no expected frames
IF hd_matches EQ 0 && hr_matches EQ 0 && pn_matches EQ 0 THEN BEGIN
   IF ~mult_stars THEN BEGIN
      PRINT, ''
      PRINT, '**NO EXPECTED FRAMES**'
      PRINT, ''
   ENDIF
   no_match_frames = 0;define this as 0 so that it can run through the if statement at the end of star_report
   no_exp_frames = 1 ;;set to 1 (on) if no expected frames
   GOTO, jump                                                       ;;jump to end
ENDIF ELSE IF star_hr EQ '----' || match_frames_hr[0] EQ '' THEN BEGIN ;IF NO HR NUM OR NO FRAMES LISTED UNDER ONE
   IF pn_matches NE 0 THEN BEGIN ;;IF THERE ARE FRAMES LISTED UNDER A PROPER NAME
      IF match_frames_hd[0] EQ -1 THEN BEGIN ;IF THERE ARE NO FRAMES WITH HD ENTRIES 
         expected_frames_tmp = [match_frames_pn]
         ind_dates = [match_ind_pn]
         temp = [match_ind_pn]
      ENDIF ELSE BEGIN;IF THERE ARE FRAMES WITH HD ENTRIES
         expected_frames_tmp = [match_frames_hd, match_frames_pn]
         ind_dates = [match_ind_hd, match_ind_pn]
         temp = [match_ind_hd, match_ind_pn]
      ENDELSE

      int_mflt = LONG(expected_frames_tmp)         ;;convert to int
      expected_frames = int_mflt[sort(int_mflt)]   ;;sort as integers
      expected_frames = STRCOMPRESS(STRING(expected_frames), /REMOVE_ALL)
      
      ind_dates = ind_dates[SORT(ind_dates)]   ;sort the dates to ensure theyre in chronological order
      expected_dates = ut_date[ind_dates]
      
      frame_observer = observer[temp[SORT(temp)]] ;array of corresponding observers of frame
   ENDIF ELSE BEGIN

      expected_frames = match_frames_hd
      expected_dates = ut_date[match_ind_hd]
      frame_observer = observer[match_ind_hd] ;array of corresponding observers of frame 
   ENDELSE
   
   ;;************************ PRINTING ************************;;
   
   IF ~mult_stars THEN BEGIN
      PRINT, 'There are', STRCOMPRESS(N_ELEMENTS(expected_frames)) , ' frames expected for '+star_name
   ENDIF

ENDIF ELSE BEGIN                   ;IF THERE IS HD AND HR FRAMES
   IF pn_matches NE 0 THEN BEGIN   ;HD, HR, frames and PN

      IF match_frames_hd[0] EQ -1 THEN BEGIN ;IF THERE ARE NO FRAMES WITH HD ENTRIES
         expected_frames_tmp = [match_frames_hr, match_frames_pn]
         ind_dates = [match_ind_hr, match_ind_pn];concatinate the date indices for the star
         temp = [match_ind_hr, match_ind_pn]
      ENDIF ELSE BEGIN;IF THERE ARE FRAMES WITH HD ENTRIES
         expected_frames_tmp = [match_frames_hd, match_frames_hr, match_frames_pn]
         ind_dates = [match_ind_hd, match_ind_hr, match_ind_pn];concatinate the date indices for the star
         temp = [match_ind_hd, match_ind_hr, match_ind_pn]
      ENDELSE

      int_mflt = LONG(expected_frames_tmp)         ;;convert to int
      expected_frames = int_mflt[sort(int_mflt)]   ;;sort as integers
      expected_frames = STRCOMPRESS(STRING(expected_frames), /REMOVE_ALL)
      
      ind_dates = ind_dates[SORT(ind_dates)] ;sort the dates to ensure theyre in chronological order
      expected_dates = ut_date[ind_dates]
      
      frame_observer = observer[temp[SORT(temp)]] ;array of corresponding observers of frame
   ENDIF ELSE BEGIN ; HD, HR, frames but no PN
      IF match_frames_hd[0] EQ -1 THEN BEGIN ;IF THERE ARE NO FRAMES WITH HD ENTRIES
         expected_frames_tmp = [match_frames_hr]
         ind_dates = [match_ind_hr]
         temp = [match_ind_hr]
      ENDIF ELSE BEGIN;IF THERE ARE FRAMES WITH HD ENTRIES
         expected_frames_tmp = [match_frames_hd, match_frames_hr]
         ind_dates = [match_ind_hd, match_ind_hr]
         temp = [match_ind_hd, match_ind_hr]
      ENDELSE

      int_mflt = LONG(expected_frames_tmp)         ;;convert to int
      expected_frames = int_mflt[sort(int_mflt)]   ;;sort as integers
      expected_frames = STRCOMPRESS(STRING(expected_frames), /REMOVE_ALL)
      
      ind_dates = ind_dates[SORT(ind_dates)]   ;sort the dates to ensure theyre in chronological order
      expected_dates = ut_date[ind_dates]
      
      frame_observer = observer[temp[SORT(temp)]] ;array of corresponding observers of frame
   ENDELSE
   
   
   ;;************************ PRINTING ************************;;
   IF ~mult_stars THEN BEGIN   
      PRINT, 'There are', STRCOMPRESS(N_ELEMENTS(expected_frames)) , ' frames expected for '+ star_name
   ENDIF   

ENDELSE  

jump:;jump here if there are no expected frames and therefore no matching frames
END
