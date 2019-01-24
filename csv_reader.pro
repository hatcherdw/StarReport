;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; csv_reader.pro
;;
;; This function was designed to read in a .csv file and extract all
;; values for each line and return a matrix array with all values for
;; all lines of the .csv file. In this matrix the rows represent the
;; values in each line of the .csv file whereas each line represents
;; the columns of it.
;;
;; Author: W. Glennon Fagan 
;;======================================================================;;


FUNCTION csv_reader, file_name


;;************************ READ LINE FROM LOG INTO ARRAY ************************;;

check_if_ascii = QUERY_ASCII(file_name, ascii)
number_of_lines = ascii.lines

;; read in the entire ASCII file into a string array
file_lines = STRARR(number_of_lines)
tmp_string = ''
OPENR, file_unit, file_name, /GET_LUN
counter = 0L

WHILE NOT EOF(file_unit) DO BEGIN
   ;; read the entire line as a string
   READF, file_unit, tmp_string
   file_lines[counter] = tmp_string
   counter += 1
ENDWHILE
CLOSE, file_unit
FREE_LUN, file_unit




;;************************ ANALYZE AND STORE DATA FROM EACH LINE ************************;;

FOR i=0, N_ELEMENTS(file_lines)-1 DO BEGIN

   line_length = STRLEN(file_lines[i])
   
   ;;Fill an array with each part of the string
   ;;This is done to easily locate commas and extract data
   line_arr = STRARR(line_length)
   FOR j=0, line_length-1 DO BEGIN
      line_arr[j] = STRMID(file_lines[i], j, 1)
   ENDFOR

   num_var = N_ELEMENTS(WHERE(line_arr EQ ',')) + 1 ;;determines how many values are in file 

   IF i EQ 0 THEN BEGIN
      ;;this will be used to compare to subsequent lines
      ;;in case there is inconsistencies in the number of values
      num_var_line1 = N_ELEMENTS(WHERE(line_arr EQ ',')) + 1
   ENDIF

   ;;locate commas and fill an array with the indices of them
   com_loc = INTARR(N_ELEMENTS(WHERE(line_arr EQ ',')) + 2)
   com_loc[0] = 0 ;;start the array with 0
   com_loc[1:N_ELEMENTS(com_loc)-2] = WHERE(line_arr EQ ',')
   com_loc[N_ELEMENTS(com_loc)-1] = line_length - 1 ;;end array with the length of the string   
   ;;For loop to store each csv into an array named 'info'
   info = STRARR(N_ELEMENTS(com_loc)-1)
   FOR j=0, N_ELEMENTS(com_loc)-2 DO BEGIN
      IF j EQ 0 THEN BEGIN
         info[j] = STRMID(file_lines[i], com_loc[j], com_loc[j+1] - com_loc[j])
      ENDIF ELSE IF j EQ N_ELEMENTS(com_loc)-2 THEN BEGIN
         info[j] = STRMID(file_lines[i], com_loc[j]+1, com_loc[j+1] - com_loc[j] + 1)
      END ELSE BEGIN
         info[j] = STRMID(file_lines[i], com_loc[j]+1, com_loc[j+1] - com_loc[j]-1)
      ENDELSE
   ENDFOR
  
   ;;************** STORE ALL VARIABLES IN AN MXN MATRIX **************;;

   IF i EQ 0 THEN BEGIN
      ;;if statement so that the array is only define the first loop
      ;;so it doesn't reset values after each line
      variables = STRARR(num_var, number_of_lines)
   ENDIF


;;***********statement to check if there are any lines w/ the wrong
;;***********amount of values
;   IF num_var NE num_var_line1 THEN BEGIN
;      ;;print statement if wrong different number of values at some
;      ;;line compared to the other lines 
;      PRINT, 'ERROR: Inconsistency in the number of variables at line:', STRCOMPRESS(i+1)
;   ENDIF
;;*****************************************************************


   FOR k=0, num_var_line1-1 DO BEGIN
   variables[k,i] = info[k]
   ENDFOR

ENDFOR


RETURN, variables
END
