;;======================================================================;;
;; Star Report Pipeline, Central Michigan University
;;                            
;; n_extrema.pro
;;                      
;; This function will return an array containg the indices for a user
;; defined number of extrema points in a given array. This function
;; utilizes the IDL procedures MAX and MIN.
;;
;; Parameters: arr = array to be searched
;;             num_ex = number of extrema 
;;             max_min = integer value; 1 for maxima, 0 for minima
;;
;; Calling Sequence: result = n_extrema(arr, num_ex, max_min)
;; 
;;
;; Author: Christian Hannah
;;
;;======================================================================;;


FUNCTION n_extrema, arr, num_ex, max_min


;; define arrays for storing indices and the original array for
;; editing
duplicate_arr = arr
extrema = INTARR(num_ex)

;; if number of elements in input array is less than num_ex, just
;; return every index
IF N_ELEMENTS(duplicate_arr) LT num_ex THEN BEGIN
   FOR i=0,N_ELEMENTS(duplicate_arr)-1 DO BEGIN
      extrema[i] = WHERE(duplicate_arr EQ duplicate_arr[i])
      GOTO, hop ;;skip to return indices
   ENDFOR
ENDIF

;; get extrema indices
FOR i=0,num_ex-1 DO BEGIN
   IF max_min EQ 1 THEN extrema[i] = MIN(WHERE(duplicate_arr EQ MAX(duplicate_arr)))
   IF max_min EQ 0 THEN extrema[i] = MIN(WHERE(duplicate_arr EQ MIN(duplicate_arr)))
   duplicate_arr[extrema[i]] = 0D
ENDFOR

;; put those indices in order
extrema = extrema[SORT(extrema)]

hop: ;;program jumps here if number of elements in input array is less than num_ex

RETURN, extrema

END
