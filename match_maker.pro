;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; match_maker.pro
;;
;; The purpose of this procedure is to search for frames in the
;; directory containing all spectra and fill an array with the frames
;; that are expected for the requested star according to the stellar
;; observation log. This is done by comparing all the frame numbers in
;; the spectra directory to the list of expected frames gathered from
;; the frame_detector.pro procedure. 
;;
;; Additionally, the frame headers are extracted which will be used in
;; a comparison in the procedure coordinator.pro.
;;
;; Author: W. Glennon Fagan
;;======================================================================;;


PRO match_maker, EXPECTED_DATES=expected_dates, EXPECTED_FRAMES=expected_frames, $
                 FILE_DATES=file_dates, FRAME_DATES_EXCLUDED=frame_dates_excluded, $
                 FRAMES_EXCLUDED=frames_excluded, FRAME_OBSERVER=frame_observer, $
                 HEADER_TARGET_NAME=header_target_name, MISSING_DATES=missing_dates,$
                 MISSING_FRAMES=missing_frames, MISSING_FRAME_OBS=missing_frame_obs, $
                 NAME_STRINGS=name_strings, NOTES_FOR_EXCLUDED=notes_for_excluded, $
                 NUM_DIREC_FILES=num_direc_files
RESOLVE_ALL, /QUIET
COMMON star_block
COMMON abort_block



exclude_poor_frames = 1



;;************************ FIND ALL SPECTRA FOR THIS STAR *************************;;

;;Call frame_detector to obtain a list of all eexpected frame numbers
;;for a given star, named 'expected_frames'
frame_detector, FRAME_OBSERVER=frame_observer, $
                EXPECTED_DATES=expected_dates, EXPECTED_FRAMES=expected_frames

;if no_exp_frames is on (equal to one) then skip to end of program
;this comes from frame_detector if there are no expected frames
IF no_exp_frames THEN GOTO, jump

;; go to the directory containing sprectra and save the names of every
;; spectrum and temporarily move it to the main directory
CD, sub_direc
SPAWN, 'ls *.txt > READ_FILES.txt'
SPAWN, 'mv -f READ_FILES.txt '+root_direc
CD, '..'

fname_temp = ''
n_files = 0
file_unit = 0
OPENR, file_unit, root_direc+'READ_FILES.txt', ERROR = err, /GET_LUN
IF(err NE 0) THEN BEGIN
   PRINT, '**ERROR: Cannot open the READ_FILES.txt file**'
ENDIF
WHILE NOT EOF(file_unit) DO BEGIN
   READF, file_unit, fname_temp
   IF (fname_temp NE '') THEN BEGIN
      n_files=n_files+1
      IF(n_files EQ 1) THEN BEGIN
         file_list=[fname_temp]
      ENDIF ELSE BEGIN
         file_list=[file_list, fname_temp]
      ENDELSE
   ENDIF
ENDWHILE
CLOSE, file_unit
FREE_LUN, file_unit

num_direc_files = N_ELEMENTS(file_list) ;;# of spectra to be tested for matches




;;************************** MAKE ARRAY OF ONLY FRAME #'S ***************************;;

;;Loop to extract only the frame # from ALL file names that are in
;;your spectra directory and store in a new array.
direc_frame_list = STRARR(N_ELEMENTS(file_list)) ;;Array of ALL frames NUMBERS from directory
numbers=['0','1','2','3','4','5','6','7','8','9'];used in following loop to search for non number characters
FOR i=0, N_ELEMENTS(file_list)-1 DO BEGIN
   line_len = STRLEN(file_list[i]) ;number of characters in file name 
   line_string = STRARR(line_len) ;array where each element will hold a character of the file name

   FOR j=4, line_len-1 DO BEGIN
      line_string[j] = STRMID(file_list[i], j,1, /REVERSE_OFFSET)
   ENDFOR
   line_string = line_string[WHERE(line_string NE '')]; remove any empty elements
   
   k=0
   frame_temp =''
   WHILE WHERE(line_string[k] EQ numbers) NE -1 DO BEGIN
      frame_temp = [frame_temp, line_string[k]] ;starting after the '.txt' tack on each number of the frame number until it reaches something other than a number. NOTICE: at this point the frame number is in reverse
      k++ ;increment k
   ENDWHILE
   frame_temp = REVERSE(frame_temp[WHERE(frame_temp NE '')]) ;remove empty elements; reverse order of elements
   direc_frame_list[i] = STRJOIN(frame_temp) ;turn array into a string
ENDFOR
direc_frame_list = direc_frame_list[WHERE(direc_frame_list NE '')] ;remove empty elements




;;*********************** FIND MATCHING FRAMES IN SPECTRA DIRECTORY ************************;;

;;For loop that will compare all expected frames for a given star from
;;the stellar obs log to all the frames in the directory. This will
;;then store all matching file names in an array
available_frames = STRARR(N_ELEMENTS(expected_frames))    ;;has full file title/list of frames for requested star
missing_frames = expected_frames     ;;array of frames that were expected but not found in spectra directory
missing_frame_obs = frame_observer   ;;array of observers corresponding to missing frames 

FOR i=0, N_ELEMENTS(direc_frame_list)-1 DO BEGIN
   ;;nested loop to compare a certain frame the directory
   ;;to all the frames for that known star in the expected_frames
   match = 0
   FOR j=0, N_ELEMENTS(expected_frames)-1 DO BEGIN
      IF direc_frame_list[i] EQ expected_frames[j] THEN BEGIN
         missing_frames[j] = ''                     ;;remove frame if we have it
         available_frames[j] = direc_frame_list[i]   ;;similarly store the frame number alone if matching
         match = 1
      ENDIF
   ENDFOR 
ENDFOR


;;Next is the array that contains just frame numbers that are
;;available in our directory
tempaf = available_frames[WHERE(available_frames NE '')]
available_frames = tempaf

;;now we will fill the array for all of the missing dates
;;corresponding to the missing frames
missing_dates = expected_dates[WHERE(missing_frames NE '')]

;;now we will fill the array for all of the observers
;;corresponding to the missing frames  
missing_frame_obs = frame_observer[WHERE(missing_frames NE '')]

;;next is missing frames, the array containing all missing frames from
;;the expected list
tempmf = missing_frames[WHERE(missing_frames NE '')]
missing_frames = tempmf
IF ~ISA(WHERE(missing_frames NE ''), /ARRAY) THEN BEGIN ;; if no missing frames set them both equal to -1
   missing_frames = -1
   missing_dates = -1
ENDIF 


;;IF statement to skip program if no matches
no_match_frames = 0
IF available_frames[0] EQ '' THEN BEGIN
   PRINT, '**WARNING: NO FILE MATCHES DETECTED**'
   no_match_frames = 1
   SPAWN, 'rm -f READ_FILES.txt';remove here because the program now skips where it would normally be removed
   GOTO, jump
ENDIF



;;********************** EXTRACTING OBS DATE FROM OBS LOG ************************;;

available_frame_dates = STRARR(N_ELEMENTS(available_frames)) ;;array containing obs dates for requested star
FOR i=0, N_ELEMENTS(available_frames)-1 DO BEGIN
   available_frame_dates[i] = expected_dates[WHERE(expected_frames EQ available_frames[i])]
ENDFOR




;;********************** COMPARE FRAMES TO LIST OF EXCLUDED FRAMES **********************;;

;;Import list of frames to be excluded
IF FILE_TEST('exclude_frames.txt') && exclude_poor_frames EQ 1 THEN BEGIN ;test for exclude frames file. If   exists use it:                                                                                                 

   check_if_ascii = QUERY_ASCII('exclude_frames.txt', ascii)
   number_of_lines = ascii.lines
   
   ;; read in the entire ASCII file into a string array
   ascii_lines = STRARR(number_of_lines)
   tmp_string = ''
   OPENR, file_unit, 'exclude_frames.txt', /GET_LUN
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

   
   ;;extract just the file names and store in an array
   exclude_frames = STRARR(N_ELEMENTS(WHERE(STRMID(ascii_lines, 0, 1) NE '%'))) ;list of ALL frames to be excluded                                                                                                    
   exclude_notes = STRARR(N_ELEMENTS(WHERE(STRMID(ascii_lines, 0, 1) NE '%')))  ;list of notes
   j=0
   FOR i=0, N_ELEMENTS(ascii_lines)-1 DO BEGIN
      IF STRMID(ascii_lines[i], 0, 1) NE '%' THEN BEGIN
         str_len = STRLEN(ascii_lines[i])
         exclude_frames[j] = STRMID(ascii_lines[i], 0, 5)
         star = STRMID(ascii_lines[i], 5, 11)
         exclude_notes[j] = STRMID(ascii_lines[i], 19, str_len)
         j++
      ENDIF
      star = ''
   ENDFOR

   analysis_frames = STRARR(N_ELEMENTS(available_frames)) ;array that will contain the updated list of frames for that star that has had the frames to be excluded removed
   analysis_dates = STRARR(N_ELEMENTS(available_frame_dates)) ;array to contain the updated list of corresponding  obs dates for the frames that arent excluded for the star                                           
   frames_excluded = STRARR(N_ELEMENTS(exclude_frames)) ;array to contain frames that WILL be excluded for that star                                                                                                   
   frame_dates_excluded = STRARR(N_ELEMENTS(exclude_frames))
   notes_for_excluded = STRARR(N_ELEMENTS(exclude_frames))
   exclude = 0
   k = 0

   ;nested for loop to compare each frame in the matching frames list to
   ;each one of the frames on the excluded frames list
   FOR i=0, N_ELEMENTS(available_frames)-1 DO BEGIN
      FOR j=0, N_ELEMENTS(exclude_frames)-1 DO BEGIN
         IF available_frames[i] EQ exclude_frames[j] THEN BEGIN
            exclude = 1 ;set this to 1(true) if this frame is in the exclude frames list
         ENDIF
      ENDFOR
      IF exclude NE 1 THEN BEGIN ;if that frame was NOT on the list to be excluded add it to the updated list of frames for the star, along with its corresponding obs date           
         analysis_frames[i] = available_frames[i]
         analysis_dates[i] = available_frame_dates[i]
      ENDIF ELSE BEGIN
         frames_excluded[k] = available_frames[i]
         frame_dates_excluded[k] = available_frame_dates[i]
         k++
      ENDELSE
      exclude = 0;it is vital that this is initialized back to 0
   ENDFOR

   IF k EQ 0 THEN BEGIN
      frames_excluded = -1
   ENDIF ELSE BEGIN
      ;;remove blanks in the array
      analysis_frames = analysis_frames[WHERE(analysis_frames NE '')]
      analysis_dates = analysis_dates[WHERE(analysis_dates NE '')]
      frames_excluded = frames_excluded[WHERE(frames_excluded NE '')]
      frame_dates_excluded = frame_dates_excluded[WHERE(frame_dates_excluded NE '')]

      ;for loop to fill array with corresponding notes as to why frameswere excluded
      FOR i=0, N_ELEMENTS(frames_excluded)-1 DO BEGIN
         notes_for_excluded[i] = exclude_notes[WHERE(frames_excluded[i] EQ exclude_frames)]
      ENDFOR
      notes_for_excluded = notes_for_excluded[WHERE(notes_for_excluded NE '')] ;remove blank indices
   ENDELSE
ENDIF ELSE BEGIN
   ;If the exclude_frames.txt file doesn't exist just act like 
   ;the original list of expected frames is the updated list                       
   frames_excluded = -1
   analysis_frames = available_frames
   analysis_dates = available_frame_dates
ENDELSE



;;*************************** CONSISTENCY CHECKING *******************************;;

;;extract star number and date line from spectrum file
CD, sub_direc
date_strings = STRARR(N_ELEMENTS(analysis_frames))
name_strings = STRARR(N_ELEMENTS(analysis_frames))
FOR i=0, N_ELEMENTS(analysis_frames)-1 DO BEGIN
   star_asp = STRARR(3)
   temp_string = ''
   OPENR, file_unit, '*'+analysis_frames[i]+'*', /GET_LUN
   FOR j=0, 2 DO BEGIN
      READF, file_unit, temp_string
      star_asp[j] = temp_string
   ENDFOR
   name_strings[i] = star_asp[0]
   date_strings[i] = star_asp[1]
   
   CLOSE, file_unit
   FREE_LUN, file_unit
ENDFOR
CD, '..'



;;************** EXTRACTING OBJECT NAME FROM FILE for consistency checking ***************;;

tmp_object = STRARR(N_ELEMENTS(analysis_frames))
header_target_name = STRARR(N_ELEMENTS(analysis_frames))

FOR i=0, N_ELEMENTS(name_strings)-1 DO BEGIN
   ;;extract object name and convert to long integer
   str_length = strlen(name_strings[i])
   junk = ''
   CASE str_length OF
      15: tmp_object[i] = STRMID(name_strings[i], 13, 1)
      16: tmp_object[i] = STRMID(name_strings[i], 13, 2)
      17: tmp_object[i] = STRMID(name_strings[i], 13, 3)
      18: tmp_object[i] = STRMID(name_strings[i], 13, 4)
      19: tmp_object[i] = STRMID(name_strings[i], 13, 5)
      20: tmp_object[i] = STRMID(name_strings[i], 13, 6)
      21: tmp_object[i] = STRMID(name_strings[i], 13, 7)
      ELSE: analysis_logs=[analysis_logs, $
                           'Frame '+ analysis_frames[i]+ ' has unreadable name or is missing one.']
   ENDCASE
   
   header_target_name[i] = STRCOMPRESS(STRLOWCASE(tmp_object[i]), /REMOVE_ALL)
ENDFOR




;;**************** EXTRACTING OBS DATE FROM FILE FOR CONSISTENCY CHECKING *****************;;

;;Extract on date from date string
file_dates = STRARR(N_ELEMENTS(date_strings))
FOR i=0, N_ELEMENTS(date_strings)-1 DO BEGIN
   junk = ''
   temp_date = ''
   READS, date_strings[i], junk, temp_date, FORMAT='(A12, A10)'
   file_dates[i] = temp_date
ENDFOR




;;********************************* CALL FRAME RETRIEVER **********************************;;

frame_retriever

jump: ;;jumps here if no matches are found for HR/HD number or no expected frames
END
