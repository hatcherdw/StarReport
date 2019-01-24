;;======================================================================;;
;; Star Report Pipeline, Central Michigan University 
;;
;; value_locator.pro
;;
;; This function will return an index corresponding to a value in the
;; supplied array that is closest to the specified value.
;;
;; Parameters: arr = array to be searched
;;             value = value for which the closest element in the
;;                     array will be found 
;; 
;; Calling Sequence: result = value_locator(arr, value)
;;
;;
;; Author: Christian Hannah
;;
;;======================================================================;; 


FUNCTION value_locator, arr, value


;; create an array to hold the differences between each element of the
;; array and the value
diffs = DBLARR(N_ELEMENTS(arr))

;; obtain differences
FOR i=0, N_ELEMENTS(arr)-1 DO BEGIN
   diffs[i] = ABS(arr[i]-value)
ENDFOR

;; smallest difference corresponds to the closest value
closest = WHERE(diffs EQ MIN(diffs))
closest = MAX(closest)
RETURN, closest

END
