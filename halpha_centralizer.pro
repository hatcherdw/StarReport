;;======================================================================;;
;; Star Report Pipeline, Central Michigan University 
;;                                      
;; halpha_centralizer.pro
;; 
;; This program takes the wavelength vectors from spectra matching the
;; supplied hd/hr number and shifts the values appropriately to center
;; all spectra on 656.28 nm using the Doppler-shift equation in log space.
;; 
;;
;; Author: Christian Hannah
;;
;;======================================================================;;


PRO halpha_centralizer, waves, centers, SHIFTS=shifts
COMMON spec_block


;; add 650 to the array of central wavelengths because 650 was
;; subtracted in the symmetry calculation in peak_analyzer.pro
center_obs = centers + 650D
center_rest = 656.28D

;; take the log of both values to properly determine the shifts
log_center_obs = ALOG10(center_obs)
log_center_rest = ALOG10(center_rest)

;; subtract the values to get the shifts
shifts = DBLARR(N_ELEMENTS(centers))
FOR i=0, N_ELEMENTS(centers)-1 DO BEGIN
   shifts[i] = log_center_rest - log_center_obs[i]
ENDFOR

;;take the log of the wavelength arrays for each spectrum to 
;;properly apply the shifts
log_waves_obs = ALOG10(waves)

;; apply the shifts
shift_log_waves_obs = DBLARR(512, N_ELEMENTS(centers))
FOR i=0, N_ELEMENTS(centers)-1 DO BEGIN
   FOR j=0,511 DO BEGIN
      shift_log_waves_obs[j,i] = log_waves_obs[j,i] + shifts[i]
   ENDFOR
ENDFOR

;; convert back to regular wavelengths
ha_waves = 10D^(shift_log_waves_obs)


;; Do the same shifting for the individual peak wavelengths
;FOR i=0, N_ELEMENTS(centers)-1 DO BEGIN
;   IF sing_peak_waves[i] NE 0 THEN BEGIN
;      sing_peak_waves[i] = sing_peak_waves[i]*(10D^shifts[i])
;   ENDIF
;   IF peak_one_waves[i] NE 0 THEN BEGIN
;      peak_one_waves[i] = peak_one_waves[i]*(10D^shifts[i])
;   ENDIF
;   IF peak_two_waves[i] NE 0 THEN BEGIN
;      peak_two_waves[i] = peak_two_waves[i]*(10D^shifts[i])
;   ENDIF
;ENDFOR


END
