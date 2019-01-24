;;======================================================================;; 
;; Star Report Pipeline, Central Michigan University 
;;
;; peak_analyzer.pro
;;
;; This procedure is designed to find all (if any) emmission peaks
;; within a spectrum. This program will allow for a maximum of 2 peaks
;; per spectrum. So, if a triple peak spectrum is ran, only 2 peaks will
;; ever be reported. Once the peaks are located, Gaussian fits to the
;; top of each peak provide the peak flux and wavelength information.
;;
;; In addition to finding the peaks, this procedure also performs a
;; level of asymmetry calculation. In this calculation, a gaussian fit
;; is applied to the wings of the emission ONLY. This fit produces a
;; central wavelength value for the emission (assuming symmetric
;; wings). This central wavelength is used to compute a degree of
;; asymmetry for the emission. It is also used in
;; halpha_centralizer.pro to shift the spectra.
;;
;;
;; Author: Christian Hannah
;;
;;======================================================================;;


PRO peak_analyzer, spec_plot_alone, temp_file, color_table
RESOLVE_ALL, /QUIET
COMMON peak_block, gauss_wings, interval, interval_1, interval_2, peak_fit, peak_fit_1, peak_fit_2, w_interval
COMMON spec_block
COMMON star_block

;; create switch for generating diagnostics plots
diagnostics = 1

DEVICE, DECOMPOSED = 0 ; Colors are taken from a color table (not specified by R, G, B)    
LOADCT, 39, /SILENT

SPAWN, 'pwd', og_loc
og_loc = og_loc+'/' ;; original directory location


;; define some printing switches used to control output at end of program
fit_err = 0
print_sing_peak = 0
print_dbl_peak = 0
fwhm_err = 0
improper_wing_err = 0
wing_err = 0


;;define initial smoothing value for smoothing the spectrum based on
;;average flux value in range 655.78 - 656.78 nm
start_i = VALUE_LOCATOR(flux, 655.78)
end_i = VALUE_LOCATOR(flux, 656.78)
avg_f = TOTAL(flux[start_i : end_i])/N_ELEMENTS(flux[start_i : end_i])
IF avg_f LT 1.2 THEN smoothing = 6 ELSE smoothing = 3


;;==========================================================================================;;
;;============================== BEGIN SEARCHING FOR PEAKS =================================;;
;;==========================================================================================;;
;;
;; REPEAT...UNTIL loop to increment smoothing until 2 or less peaks
;; fitting the requirements discussed below are found or a smoothing value of 20
;; is reached (indicates absorption-only spectrum). So, the spectrum
;; gets smoother upon each iteration of this loop. All code related to
;; searching for peaks lies within this loop and is thus repeated each time the
;; spectrum is re-smoothed. 
;;
REPEAT BEGIN

   ;;********************** FINDING ALL PEAKS **************************;;
   ;;
   ;; This section of the procedure will find all local maxima within the 
   ;; spectrum. This is done using the DERIV procedure from IDL that 
   ;; calculates the derivative at each point in the spectrum. The peaks
   ;; are found by storing the indices where the derivative changes from 
   ;; postive to negative. To ensure this point is in fact the largest 
   ;; point in the immediate area, the "peak" flux value is then compared
   ;; to the flux values of the points to the immediate right and
   ;; left. If either point is larger, it will be considered as the peak.
   ;; 
   
   ;; smooth the flux to help ignore some noise
   asmooth_flux = DBLARR(N_ELEMENTS(flux))
   asmooth_flux = SMOOTH(flux[40:*], smoothing)
   
   smoothing+=1 ;;increment smoothing

   absorption = 0 ;; used to print a statement if spectrum is an absorption

   ;;calculate derivatives
   derivs = DERIV(wave[40:*], asmooth_flux)

   ;;flag points with positive derivatives
   pos = WHERE(derivs GE 0)
   pos_der_flag = derivs*0D
   pos_der_flag[pos] = 1D

   ;; flag the peaks (where deriv goes from >0 to <0)
   der_change = derivs*0D
   FOR i=0, N_ELEMENTS(der_change)-2 DO BEGIN
      IF pos_der_flag[i] GT 0 && pos_der_flag[i+1] EQ 0 THEN der_change[i] = 1D
   ENDFOR
   
   peak_ind = WHERE(der_change EQ 1D) ;;array of peak indices
   num_peaks = N_ELEMENTS(peak_ind)   ;;number of peaks found

   ;; check the points to the imediate right and left and uses the
   ;; highest value for the peak index
   FOR i=0, num_peaks-1 DO BEGIN
      k = peak_ind[i]
      max_check = MAX([asmooth_flux[(k-1)>0], asmooth_flux[k], asmooth_flux[(k+1)<(N_ELEMENTS(asmooth_flux)-1)]],$
                      index)
      IF index NE 1 THEN peak_ind[i] = ((peak_ind[i]-1+index)>0)<(N_ELEMENTS(asmooth_flux)-1)
   ENDFOR

   ;;*******************************************************************;;



   ;;********************** COMPUTE PEAK DATA **************************;;
   ;;
   ;; for each peak found compute:
   ;; 1) The number of points with derivarive positive sitting at its
   ;; left (l_npts) plus the number of points with negative derivarive sitting at 
   ;; its right (r_npts). This sum is stored in peak_npts.
   ;; 2) Three weights for each peak, consisting of the addition of absolute value 
   ;; of the derivatives for the points calculated in 1. Stored in
   ;; l_weight, r_weight, and peak_weight (sum of l_weight and
   ;; r_weight)
   ;; 3) the flux difference between peak and left trough (l_flux) and
   ;; the flux difference between peak and right trough (r_flux)
   ;; 4) The distance in wavelength from the peak to 656.28 nm stored
   ;; in diffs.
   ;;
   peak_weight = DBLARR(num_peaks)
   l_weight = DBLARR(num_peaks)
   r_weight = DBLARR(num_peaks)
   peak_npts = INTARR(num_peaks)
   l_npts = INTARR(num_peaks)
   r_npts = INTARR(num_peaks)
   diffs = DBLARR(num_peaks)
   l_flux = DBLARR(num_peaks)
   r_flux = DBLARR(num_peaks)
   l_trough_ind = INTARR(num_peaks)
   r_trough_ind = INTARR(num_peaks)

   FOR i=0, num_peaks-1 DO BEGIN
      peak_npts[i] = 1
      l_npts[i] = 0
      r_npts[i] = 0
      l_weight[i] = 0.0
      r_weight[i] = 0.0
      peak_weight[i] = 0.0
      
      ;;left search
      good = 1
      j = peak_ind[i]
      WHILE good EQ 1 DO BEGIN
         j = (j-1)>0
         IF derivs[j] GE 0 THEN BEGIN
            peak_npts[i] = peak_npts[i]+1
            l_npts[i] = l_npts[i]+1
            l_weight[i] = l_weight[i]+abs(derivs[j])
            peak_weight[i] = peak_weight[i]+abs(derivs[j])
         ENDIF ELSE IF asmooth_flux[j+1]-asmooth_flux[j] GT 0 THEN BEGIN
            peak_npts[i] = peak_npts[i]+1
            r_npts[i] = r_npts[i]+1
         ENDIF ELSE good = 0
         IF j EQ 0 THEN good = 0
      ENDWHILE
      
      ;;right search
      good = 1
      j = peak_ind[i]
      WHILE good EQ 1 DO BEGIN
         j = (j+1)<(N_ELEMENTS(asmooth_flux)-1)
         IF derivs[j] LT 0 THEN BEGIN
            peak_npts[i] = peak_npts[i]+1
            r_npts[i] = r_npts[i]+1
            r_weight[i] = r_weight[i]+abs(derivs[j])
            peak_weight[i]= peak_weight[i]+abs(derivs[j])
         ENDIF ELSE IF asmooth_flux[j-1]-asmooth_flux[j] GT 0 THEN BEGIN
            peak_npts[i] = peak_npts[i]+1
            r_npts[i] = r_npts[i]+1
         ENDIF ELSE good = 0
         IF j EQ N_ELEMENTS(asmooth_flux)-1 THEN good = 0
      ENDWHILE

      ;; calculate diffs
      ref = wave[VALUE_LOCATOR(wave, 656.28)]
      diffs[i] = ref-wave[peak_ind[i]+40]
      
      ;;calculate l_flux and r_flux
      l_trough_ind[i] = N_EXTREMA(asmooth_flux[(peak_ind[i]-l_npts[i])>0:peak_ind[i]],1,0)+(peak_ind[i]-l_npts[i])>0
      r_trough_ind[i] = N_EXTREMA(asmooth_flux[peak_ind[i]:(peak_ind[i]+r_npts[i])<471],1,0)+peak_ind[i]<471
      l_flux[i] = asmooth_flux[peak_ind[i]]-asmooth_flux[l_trough_ind[i]]
      r_flux[i] = asmooth_flux[peak_ind[i]]-asmooth_flux[r_trough_ind[i]]
   ENDFOR
   
   ;;*******************************************************************;;



   ;;********************* COMBINE PEAKS *******************************;;
   ;;
   ;; To avoid having two seperate peaks in one due to noise, combine
   ;; peaks using information regarding the difference in flux
   ;; from peak to each neighboring trough (r_flux and l_flux). 
   ;;
   ;; For any two peaks to be combined the flux at both peaks must 
   ;; be > 1±4*std_flux.
   ;; 
   ;; In addition, the following cases are used to combine peaks:
   ;; 1) left peak has r_flux < percent_f * the flux at the peak,
   ;; right peak has l_flux < percent_f * the flux at the peak, and 
   ;; the peaks are closer than 0.5 nm apart.
   ;; 2) left peak has r_flux < percent_f * the flux at the peak and
   ;; l_flux > percent_f * the flux at the peak and right peak must
   ;; have l_flux > percent_f * the flux at the peak.
   ;; 3) right_peak has l_flux <  percent_f * the flux at the peak and
   ;; r_flux >  percent_f * the flux at the peak and left peak has
   ;; r_flux >  percent_f * the flux at the peak
   ;;
   
   ;; if the first peak_ind is zero, remove this peak as it will cause
   ;; problems when removing zeros from newly constructed temp arrays
   IF peak_ind[0] EQ 0 THEN BEGIN 
      peak_ind = peak_ind[1:*]
      peak_weight = peak_weight[1:*]
      l_weight = l_weight[1:*]
      r_weight = r_weight[1:*]
      peak_npts = peak_npts[1:*]
      l_npts = l_npts[1:*]
      r_npts = r_npts[1:*]
      diffs = diffs[1:*]
      l_flux = l_flux[1:*]
      r_flux = r_flux[1:*]
      l_trough_ind = l_trough_ind[1:*]
      r_trough_ind = r_trough_ind[1:*]
      num_peaks = num_peaks-1
   ENDIF
   
   ;; calculate std_flux
   left_bound = VALUE_LOCATOR(wave, 652.55)
   right_bound = VALUE_LOCATOR(wave, 654.55)
   std_flux = STDDEV(flux[left_bound :right_bound])

   percent_f = .02 ;;percentage of flux
   changed = 1     ;;set to 1 to enter the loop
   
   WHILE changed && num_peaks GT 2 DO BEGIN
      changed = 0 
      
      temp_peak_ind = INTARR(num_peaks)
      temp_peak_weight = DBLARR(num_peaks)
      temp_l_weight = DBLARR(num_peaks)
      temp_r_weight = DBLARR(num_peaks)
      temp_peak_npts = INTARR(num_peaks)
      temp_l_npts = INTARR(num_peaks)
      temp_r_npts = INTARR(num_peaks)
      temp_diffs = DBLARR(num_peaks)
      
      j=0
      FOR i=0, num_peaks-2 DO BEGIN
         ;; Case 1
         ;; move peak_ind between the two peaks and sum characteristics
         ;; from each peak
         IF r_flux[i] LE percent_f*asmooth_flux[peak_ind[i]] && $
            l_flux[i+1] LE percent_f*asmooth_flux[peak_ind[i+1]] && $
            ABS(wave[peak_ind[i]]-wave[peak_ind[i+1]]) LE 0.5 && $
            ABS(1-asmooth_flux[peak_ind[i]]) GT 4*std_flux && $
            ABS(1-asmooth_flux[peak_ind[i+1]]) GT 4*std_flux THEN BEGIN
            ;print, 'case 1'
            temp_peak_ind[j] = peak_ind[i]+FIX((peak_ind[i+1]-peak_ind[i])/2D)
            temp_peak_weight[j] = peak_weight[i]+peak_weight[i+1]
            temp_l_weight[j] = l_weight[i]+l_weight[i+1]
            temp_r_weight[j] = r_weight[i]+r_weight[i+1]
            temp_peak_npts[j] = peak_npts[i]+peak_npts[i+1]
            temp_l_npts[j] = l_npts[i]+l_npts[i+1]
            temp_r_npts[j] = r_npts[i]+r_npts[i+1]
            temp_diffs[j] = ref-wave[temp_peak_ind[j]+40]
            j+=1
            i+=1
            changed = 1

         ;; Case 2
         ;; get rid of this index and combine statistics with peak
         ;; to imediate right
         ENDIF ELSE IF r_flux[i] LE percent_f*asmooth_flux[peak_ind[i]] && $
            l_flux[i] GT percent_f*asmooth_flux[peak_ind[i]] && $
            l_flux[i+1] GT percent_f*asmooth_flux[peak_ind[i+1]] && $
            ABS(1-asmooth_flux[peak_ind[i]]) GT 4*std_flux && $
            ABS(1-asmooth_flux[peak_ind[i+1]]) GT 4*std_flux THEN BEGIN
            ;print, 'case 2'
            temp_peak_ind[j] = peak_ind[i+1]
            temp_peak_weight[j] = peak_weight[i]+peak_weight[i+1]
            temp_l_weight[j] = peak_weight[i]+l_weight[i+1]
            temp_r_weight[j] = r_weight[i]
            temp_peak_npts[j] = peak_npts[i]+peak_npts[i+1]
            temp_l_npts[j] = peak_npts[i]+l_npts[i+1]
            temp_r_npts[j] = r_npts[i+1]
            temp_diffs[j] = ref-wave[temp_peak_ind[j]+40]
            j+=1
            i+=1
            changed = 1

         ;; Case 3   
         ;; get rid of this index and combine statistics with peak
         ;; to imediate left
         ENDIF ELSE IF l_flux[i+1] LE percent_f*asmooth_flux[peak_ind[i+1]] && $
            r_flux[i+1] GT percent_f*asmooth_flux[peak_ind[i+1]] && $
            r_flux[i] GT percent_f*asmooth_flux[peak_ind[i]] && $
            ABS(1-asmooth_flux[peak_ind[i]]) GT 4*std_flux && $
            ABS(1-asmooth_flux[peak_ind[i+1]]) GT 4*std_flux THEN BEGIN
            ;print, 'case 3'
            temp_peak_ind[j] = peak_ind[i]
            temp_peak_weight[j] = peak_weight[i]+peak_weight[i+1]
            temp_l_weight[j] = l_weight[i]
            temp_r_weight[j] = r_weight[i]+peak_weight[i+1]
            temp_peak_npts[j] = peak_npts[i]+peak_npts[i+1]
            temp_l_npts[j] = l_npts[i]
            temp_r_npts[j] = r_npts[i]+peak_npts[i+1]
            temp_diffs[j] = ref-wave[temp_peak_ind[j]+40]
            j+=1
            i+=1
            changed = 1

         ;; no combination needed
         ENDIF ELSE BEGIN
            temp_peak_ind[j] = peak_ind[i]
            temp_peak_weight[j] = peak_weight[i]
            temp_l_weight[j] = l_weight[i]
            temp_r_weight[j] = r_weight[i]
            temp_peak_npts[j] = peak_npts[i]
            temp_l_npts[j]= l_npts[i]
            temp_r_npts[j]= r_npts[i]
            temp_diffs[j] = ref-wave[temp_peak_ind[j]+40]
            j+=1
            ;; ensure that the last peak_ind is also added to the temp arrays  
            IF i EQ num_peaks-2 THEN BEGIN 
               temp_peak_ind[j] = peak_ind[i+1]
               temp_peak_weight[j] = peak_weight[i+1]
               temp_l_weight[j] = l_weight[i+1]
               temp_r_weight[j] = r_weight[i+1]
               temp_peak_npts[j] = peak_npts[i+1]
               temp_l_npts[j]= l_npts[i+1]
               temp_r_npts[j]= r_npts[i+1]
               temp_diffs[j] = ref-wave[temp_peak_ind[j]+40]
               j+=1
            ENDIF
         ENDELSE
      ENDFOR

      ;;remove zeros from temp arrays
      temp_peak_ind = temp_peak_ind[WHERE(temp_peak_ind NE 0)]
      temp_peak_weight = temp_peak_weight[WHERE(temp_peak_ind NE 0)]
      temp_l_weight = temp_l_weight[WHERE(temp_peak_ind NE 0)]
      temp_r_weight = temp_r_weight[WHERE(temp_peak_ind NE 0)]
      temp_peak_npts = temp_peak_npts[WHERE(temp_peak_ind NE 0)]
      temp_l_npts = temp_l_npts[WHERE(temp_peak_ind NE 0)]
      temp_r_npts = temp_r_npts[WHERE(temp_peak_ind NE 0)]
      temp_diffs = temp_diffs[WHERE(temp_peak_ind NE 0)]
      
      ;; define new value for num_peaks
      num_peaks = N_ELEMENTS(temp_peak_ind)

      ;; update l_flux and r_flux arrays
      r_flux = DBLARR(num_peaks)
      l_flux = DBLARR(num_peaks)
      l_trough_ind = INTARR(num_peaks)
      r_trough_ind = INTARR(num_peaks)
      FOR i=0, num_peaks-1 DO BEGIN
         l_trough_ind[i] = N_EXTREMA(asmooth_flux[(temp_peak_ind[i]-temp_l_npts[i])>0:temp_peak_ind[i]],1,0)+$
                        (temp_peak_ind[i]-temp_l_npts[i])>0
         r_trough_ind[i] = N_EXTREMA(asmooth_flux[temp_peak_ind[i]:(temp_peak_ind[i]+temp_r_npts[i])<471],1,0)+$
                        temp_peak_ind[i]
         l_flux[i] = asmooth_flux[temp_peak_ind[i]]-asmooth_flux[l_trough_ind[i]]
         r_flux[i] = asmooth_flux[temp_peak_ind[i]]-asmooth_flux[r_trough_ind[i]]
      ENDFOR

      ;;test plotting
      ;PLOT, wave[40:*], asmooth_flux, PSYM=-3
      ;OPLOT, wave[temp_peak_ind+40], asmooth_flux[temp_peak_ind], PSYM=4, COLOR=150
      ;print, temp_peak_ind
      ;stop

      ;; change name of temp arrays
      peak_ind = temp_peak_ind
      peak_weight = temp_peak_weight
      l_weight = temp_l_weight
      r_weight = temp_r_weight
      peak_npts =temp_peak_npts
      l_npts = temp_l_npts
      r_npts = temp_r_npts
      diffs = temp_diffs

   ENDWHILE
   
   ;;*******************************************************************;;



   ;;******** COMPUTE DATA FOR PEAK SELECTION AND NARROW RANGE *********;;

   ;;calculate a standard deviation for weights, after excluding 2
   ;;largest weights to avoid having a skewed stddev from true peaks
   max_weights = N_EXTREMA(peak_weight,2,1)
   temp_weight = peak_weight
   temp_weight[max_weights] = 0D
   std_weight = STDDEV(temp_weight)
   avg_weight = TOTAL(temp_weight)/(N_ELEMENTS(temp_weight)-2)

   max_l_weights = N_EXTREMA(l_weight,2,1)
   temp_l_weight = l_weight
   temp_l_weight[max_l_weights] = 0D
   std_l_weight =STDDEV(temp_l_weight)
   avg_l_weight =TOTAL(temp_l_weight)/(N_ELEMENTS(temp_l_weight)-2)

   max_r_weights = N_EXTREMA(r_weight,2,1)
   temp_r_weight = r_weight
   temp_r_weight[max_r_weights] = 0D
   std_r_weight =STDDEV(temp_r_weight)
   avg_r_weight =TOTAL(temp_r_weight)/(N_ELEMENTS(temp_r_weight)-2)
   
   ;; Calculate the stddev of l_npts only for peaks left of 656.28 nm (std_ll_npts)
   ;; and stddev of r_npts for peaks right of 656.28 nm (std_rr_npts).
   rr_npts = r_npts[WHERE(diffs LT 0)]
   std_rr_npts = STDDEV(rr_npts)
   avg_rr_npts = TOTAL(rr_npts)/N_ELEMENTS(rr_npts)
   ll_npts = l_npts[WHERE(diffs GT 0)]
   std_ll_npts = STDDEV(ll_npts)
   avg_ll_npts = TOTAL(ll_npts)/N_ELEMENTS(ll_npts)

   ;;narrow peaks by excluding ones outside of the range of 656.28 ± 0.8 nm
   test = WHERE(ABS(diffs) LE 0.8D)
   IF test[0] NE -1 THEN BEGIN
      close_peak_ind = peak_ind[WHERE(ABS(diffs) LE 0.8D)]
      close_peak_weight = peak_weight[WHERE(ABS(diffs) LE 0.8D)]
      close_diffs = diffs[WHERE(ABS(diffs) LE 0.8D)]
      close_l_weight = l_weight[WHERE(ABS(diffs) LE 0.8D)]
      close_r_weight = r_weight[WHERE(ABS(diffs) LE 0.8D)]
      close_peak_npts = peak_npts[WHERE(ABS(diffs) LE 0.8D)]
      close_l_npts = l_npts[WHERE(ABS(diffs) LE 0.8D)]
      close_r_npts = r_npts[WHERE(ABS(diffs) LE 0.8D)]
      close_l_trough_ind = l_trough_ind[WHERE(ABS(diffs) LE 0.8D)]
      close_r_trough_ind = r_trough_ind[WHERE(ABS(diffs) LE 0.8D)]
      close_r_flux = r_flux[WHERE(ABS(diffs) LE 0.8D)]
      close_l_flux = l_flux[WHERE(ABS(diffs) LE 0.8D)]
   ENDIF ELSE BEGIN
      absorption = 1
      peak_num = 0
      GOTO, jump
   ENDELSE
   
   ;;*******************************************************************;;



   ;;******************** BEST PEAK SELECTION **************************;;
   ;;
   ;; Run remaining peak_ind through a set of conditions to determine
   ;; the best peak/peaks (max 2).
   ;; 
   ;; These conditions fall under 2 cases:
   ;; 
   ;; Case A: Less than 4 total peaks found
   ;; -> Conditions: (all apply if only one close_peak_ind)
   ;;      1) l_weight > 4
   ;;      2) r_weight > 4
   ;;      3) l_npts ≥ 3
   ;;      4) r_npts ≥ 3
   ;; -> If there are more than one close_peak_ind, also must satisfy:
   ;;      Left: 2 & 3 & peak_weight > 6
   ;;      Right: 1 & 4 & peak_weight > 6
   ;; 
   ;; Case B: More than 4 total peaks found
   ;; -> Conditions:
   ;;      Left:
   ;;         1) peak_weight > avg_weight+std_weight
   ;;         2a) l_weight > avg_l_weight and r_npts ≥ 2
   ;;               or
   ;;         2b) l_npts ≥ avg_ll_npts+std_ll_npts and r_npts ≥ 3
   ;;         
   ;;      Right:
   ;;         1) peak_weight > avg_weight+std_weight
   ;;         2a) r_weight > avg_r_weight and l_npts ≥ 2
   ;;               or
   ;;         2b) r_npts ≥ avg_rr_npts+std_rr_npts and l_npts ≥ 3
   ;;      
   ;;      *** ALSO, if peak_weight > avg_weight+25*std_weight, peak is
   ;;      counted regardless of previous conditions. ***
   ;;
   
   ;; set up arrays to contain best peak info
   best_peak_ind = close_peak_ind*0D
   best_peak_npts = close_peak_ind*0D
   best_l_npts = close_peak_ind*0D
   best_r_npts = close_peak_ind*0D
   best_l_trough_ind = close_peak_ind*0D
   best_r_trough_ind = close_peak_ind*0D
   k=0

   ;;
   ;; CASE A
   ;;
   ;; average and stddev statistics aren't useful if there aren't 
   ;; enough peaks found. So, if 3 or less peak_ind are found then
   ;; compare values with different conditions.. 
   IF N_ELEMENTS(peak_ind) LT 4 THEN BEGIN
      PRINT, 'NOTICE: LESS THAN 4 TOTAL PEAKS FOUND'
      ;;just one peak
      IF N_ELEMENTS(close_peak_ind) EQ 1 THEN BEGIN
         IF close_l_weight[0] GT 4.0 && close_r_weight[0] GT 4.0 && close_l_npts[0] GE 3 && $
            close_r_npts[0] GE 3 THEN BEGIN
            best_peak_ind[k] = close_peak_ind[0]
            best_peak_npts[k] = close_peak_npts[0]
            best_l_npts[k] = close_l_npts[0]
            best_r_npts[k] = close_r_npts[0]
            best_l_trough_ind[k] = close_l_trough_ind[0]
            best_r_trough_ind[k] = close_r_trough_ind[0]
         ENDIF
      ENDIF ELSE BEGIN
         ;;multiple peaks
         FOR i=0, N_ELEMENTS(close_peak_ind)-1 DO BEGIN
            ;;left peaks
            IF close_diffs[i] GT 0 && close_peak_weight[i] GT 6D THEN BEGIN
               IF close_r_weight[i] GT 4.0 && close_l_npts[i] GE 3 THEN BEGIN
                  best_peak_ind[k] = close_peak_ind[i]
                  best_peak_npts[k] = close_peak_npts[i]
                  best_l_npts[k] = close_l_npts[i]
                  best_r_npts[k] = close_r_npts[i]
                  best_l_trough_ind[k] = close_l_trough_ind[i]
                  best_r_trough_ind[k] = close_r_trough_ind[i]
                  k+=1
               ENDIF
               ;;right peaks
            ENDIF ELSE IF close_diffs[i] LE 0 && close_peak_weight[i] GT 6D THEN BEGIN
               IF close_l_weight[i] GT 4.0 && close_r_npts[i] GE 3 THEN BEGIN
                  best_peak_ind[k] = close_peak_ind[i]
                  best_peak_npts[k] = close_peak_npts[i]
                  best_l_npts[k] = close_l_npts[i]
                  best_r_npts[k] = close_r_npts[i]
                  best_l_trough_ind[k] = close_l_trough_ind[i]
                  best_r_trough_ind[k] = close_r_trough_ind[i]
                  k+=1
               ENDIF
            ENDIF
         ENDFOR
      ENDELSE
   ENDIF

   ;;
   ;; CASE B
   ;;
   ;; use avg and stddev statistics to select at most 2 best peaks
   ;; using the following conditions
   w=0
   FOR i=0, N_ELEMENTS(close_peak_ind)-1 DO BEGIN
      ;;left peaks
      IF close_diffs[i] GT 0 && close_peak_weight[i] GT avg_weight+std_weight && $
         close_r_flux[i] GT percent_f*asmooth_flux[close_peak_ind[i]] THEN BEGIN
         IF (close_l_weight[i] GT avg_l_weight && close_r_npts[i] GE 2) || $
            (close_l_npts[i] GE FLOOR(avg_ll_npts+std_ll_npts) && $
             close_r_npts[i] GE 3) THEN BEGIN
            best_peak_ind[w] = close_peak_ind[i]
            best_peak_npts[w] = close_peak_npts[i]
            best_l_npts[w] = close_l_npts[i]
            best_r_npts[w] = close_r_npts[i]
            best_l_trough_ind[w] = close_l_trough_ind[i]
            best_r_trough_ind[w] = close_r_trough_ind[i]
            w += 1
         ENDIF ELSE IF close_peak_weight[i] GT avg_weight+25*std_weight THEN BEGIN
            best_peak_ind[w] = close_peak_ind[i]
            best_peak_npts[w] = close_peak_npts[i]
            best_l_npts[w] = close_l_npts[i]
            best_r_npts[w] = close_r_npts[i]
            best_l_trough_ind[w] = close_l_trough_ind[i]
            best_r_trough_ind[w] = close_r_trough_ind[i]
            w += 1
         ENDIF
      ENDIF
      ;;right peaks
      IF close_diffs[i] LE 0 && close_peak_weight[i] GT avg_weight+std_weight && $
      close_l_flux[i] GT percent_f*asmooth_flux[close_peak_ind[i]] THEN BEGIN
         IF (close_r_weight[i] GT avg_r_weight && close_l_npts[i] GE 2) || $
            (close_r_npts[i] GE avg_rr_npts+std_rr_npts && $
             close_l_npts[i] GE 3) THEN BEGIN
            best_peak_ind[w] = close_peak_ind[i]
            best_peak_npts[w] = close_peak_npts[i]
            best_l_npts[w] = close_l_npts[i]
            best_r_npts[w] = close_r_npts[i]
            best_l_trough_ind[w] = close_l_trough_ind[i]
            best_r_trough_ind[w] = close_r_trough_ind[i]
            w += 1
         ENDIF ELSE IF close_peak_weight[i] GT avg_weight+25*std_weight THEN BEGIN
            best_peak_ind[w] = close_peak_ind[i]
            best_peak_npts[w] = close_peak_npts[i]
            best_l_npts[w] = close_l_npts[i]
            best_r_npts[w] = close_r_npts[i]
            best_l_trough_ind[w] = close_l_trough_ind[i]
            best_r_trough_ind[w] = close_r_trough_ind[i]
            w += 1
         ENDIF
      ENDIF
   ENDFOR
   
   ;;remove zeros from best_peak_ind and store number of best peaks in n_best
   best_peak_ind = best_peak_ind[WHERE(best_peak_ind NE 0)]
   best_peak_npts = best_peak_npts[WHERE(best_peak_ind NE 0)]
   best_l_npts = best_l_npts[WHERE(best_peak_ind NE 0)]
   best_l_trough_ind = best_l_trough_ind[WHERE(best_peak_ind NE 0)]
   best_r_trough_ind = best_r_trough_ind[WHERE(best_peak_ind NE 0)]
   n_best = N_ELEMENTS(best_peak_ind)
   
   ;;*******************************************************************;;

ENDREP UNTIL n_best LE 2 || smoothing GT 20 ;;maximum 2 peaks to be found and smoothing should not exceed 20

;;==========================================================================================;;
;;=============================== END SEARCHING FOR PEAKS ==================================;;
;;==========================================================================================;; 






;;==========================================================================================;; 
;;======================== MAKE SURE "BEST" PEAKS AREN'T JUST NOISE ========================;;
;;==========================================================================================;; 

;; set absorption to 1 if no peaks were found, 1 peak was found but is
;; farther than 0.5 nm from 656.28 nm, or smoothing exceeded 20.
IF best_peak_ind[-1] EQ 0 || smoothing GT 20 || (n_best EQ 1 && ABS(ref-wave[best_peak_ind+40]) GT 0.5) THEN BEGIN
   absorption = 1
   peak_num = 0
   GOTO, jump
ENDIF ELSE BEGIN ;; define number of peaks for the spectrum if not an absorption-only
   CASE n_best OF 
      1: peak_num = 1
      2: peak_num = 2
   ENDCASE
ENDELSE

;; if 2 peaks are found make sure one of them isn't just noise
;; near the continuum
IF peak_num EQ 2 THEN BEGIN
   l_strength = asmooth_flux[best_peak_ind[0]]
   r_strength = asmooth_flux[best_peak_ind[1]]
   IF (ABS(1-l_strength) LT 0.15*ABS(1-r_strength) XOR ABS(1-r_strength) LT 0.15*ABS(1-l_strength)) || $
   ((l_strength-1 GT 0.25 && r_strength-1 LT 0) XOR (l_strength-1 LT 0 && r_strength-1 GT 0.25)) THEN BEGIN
      peak_num = 1
      temp_bpi = INTARR(1)
      temp_pnpts = INTARR(1)  
      temp_lnpts = INTARR(1)
      temp_rnpts = INTARR(1)
      temp_lti = INTARR(1)
      temp_rti = INTARR(1)
      IF l_strength GT r_strength THEN BEGIN
         temp_bpi[0] = best_peak_ind[0]
         temp_pnpts[0] = best_peak_npts[0]
         temp_lnpts[0] = best_l_npts[0]
         temp_rnpts[0] = best_r_npts[0]
         temp_lti[0] = best_l_trough_ind[0]
         temp_rti[0] = best_r_trough_ind[0]
      ENDIF ELSE BEGIN
         temp_bpi[0] = best_peak_ind[1]
         temp_pnpts[0] = best_peak_npts[1]
         temp_lnpts[0] = best_l_npts[1]
         temp_rnpts[0] = best_r_npts[1]
         temp_lti[0] = best_l_trough_ind[1]
         temp_rti[0] = best_r_trough_ind[1]
      ENDELSE
      best_peak_ind = temp_bpi
      best_peak_npts = temp_pnpts
      best_l_npts = temp_lnpts
      best_r_npts = temp_rnpts
      best_l_trough_ind = temp_lti
      best_r_trough_ind = temp_rti
   ENDIF
ENDIF

;; adjust best_peak_ind, best_l_trough_ind, and best_r_trough_ind back
;; to proper indices for full flux array
best_peak_ind = best_peak_ind+40
best_l_trough_ind = best_l_trough_ind+40 
best_r_trough_ind = best_r_trough_ind+40

;; Test Plotting
;PLOT, wave[40:*], asmooth_flux, PSYM=-3, XRANGE=[654,659]
;IF best_peak_ind[-1] NE 0 THEN BEGIN
;   OPLOT, wave[best_peak_ind], asmooth_flux[best_peak_ind-40], PSYM=4, COLOR=150
;ENDIF ELSE BEGIN
;   XYOUTS, 0.5, 0.5, 'NO PEAKS', CHARSIZE=8, COLOR=250, /NORMAL, ALIGNMENT=0.5
;ENDELSE

;;==========================================================================================;;
;;=================== END MAKE SURE "BEST" PEAKS AREN'T JUST NOISE =========================;;
;;==========================================================================================;;






;;==========================================================================================;;
;;================================= FITTING THE PEAKS ======================================;;
;;==========================================================================================;;
;;
;; This area of the procedure performs the peak fitting necessary to
;; report peak flux values and peak wavelengths.
;;
;; This process also helps to determine if the peak ranges found really
;; have a peak in them. It does this by selecting the num_pts largest
;; flux values in the range of indices given by the trough indices and
;; fitting a gaussian to them. The number of maximum points
;; included is incremented by 1 until the center of the gaussian fit(s) lies within
;; the range(s) or all points in the range are included. If all points in
;; a range are included and the fit's center still does not
;; lie within the range, a warning will be printed. 
;;

IF peak_num EQ 1 THEN BEGIN ;; 1 peak
   ;;define range of indices including all points to the left and
   ;;right of the peak (both ending at the left and right troughs of
   ;;the peak)
   range = INDGEN(best_r_trough_ind[0]-best_l_trough_ind[0]) + best_l_trough_ind[0] 

   ;; determine the proper interval for the top of the peak and fit a
   ;; gaussian to determine the peak flux and wavelength 
   num_pts = 5
   REPEAT BEGIN
      num_pts += 1
      interval = n_extrema(flux[range], num_pts, 1)+range[0]
      peak_fit = GAUSSFIT(wave[interval], flux[interval], CHISQ=chisq, gauss_para, NTERMS=3)
   ENDREP UNTIL gauss_para[1] GE wave[interval[0]] && gauss_para[1] LE wave[interval[-1]] $
      || num_pts GE N_ELEMENTS(flux[range])

   ;;***Test Plotting***
   ;OPLOT, wave[interval], peak_fit, PSYM=-3, COLOR=150
   
   ;; make sure the fits are correct                                                                                
   IF num_pts GE N_ELEMENTS(flux[range]) && (gauss_para[1] LT wave[interval[0]] || $
                                             gauss_para[1] GT wave[interval[-1]]) THEN BEGIN
      fit_err = 1 ;;set switch to print warning message
      flag[1] = 2
      peak_flux = 0D
      peak_wave = 0D
      peak_num = 0
      symmetry = 0D
      center = 6.28D
      GOTO, bypass
   ENDIF

   print_sing_peak = 1 ;; set switch to print single peak information 

   ;; define peak flux and wavelength
   peak_flux = gauss_para[0]
   peak_wave = gauss_para[1]

   ;;calculate FWHM for single peak spectra
   l_rng = INDGEN(best_peak_ind[0]-best_l_trough_ind[0])+best_l_trough_ind[0]
   r_rng = INDGEN(best_r_trough_ind[0]-best_peak_ind[0])+best_peak_ind[0]
   
   try_again: ;;program jumps here if half_max points are incorrect
   IF N_ELEMENTS(WHERE(l_rng EQ 0)) EQ N_ELEMENTS(l_rng) || $
      N_ELEMENTS(WHERE(r_rng EQ 511)) EQ N_ELEMENTS(r_rng) THEN BEGIN
      fwhm_err = 1 ;; set switch to print warning message
      ;PRINT, 'WARNING: FWHM ERROR. Correct half-max points could not be found.'
      fwhm = 0.0
      GOTO, break_out
   ENDIF
   
   IF flux[l_rng[0]] GT 1 && flux[r_rng[-1]] GT 1 THEN BEGIN
      half_max = (peak_flux-1D)/2D + 1D ;;flux value for half peak flux
   ENDIF ELSE BEGIN
      avg_base = (flux[l_rng[0]] + flux[r_rng[-1]])/2D
      half_max = (peak_flux-avg_base)/2D + avg_base
   ENDELSE

   flux_left = VALUE_LOCATOR(flux[l_rng],half_max) + l_rng[0] ;; index of point left of peak at half_max
   flux_right = VALUE_LOCATOR(flux[r_rng],half_max) + r_rng[0] ;; index of point right of peak at half_max

   ;;define a threshold (min_sep) depending on peak_flux value
   IF peak_flux GT 10 THEN min_sep = 2.0 ELSE min_sep = 0.8

   IF ABS(flux[flux_left]-half_max) GT min_sep THEN BEGIN
      l_rng = l_rng-10>0
      GOTO, try_again
   ENDIF
   IF ABS(flux[flux_right]-half_max) GT min_sep THEN BEGIN
      r_rng = r_rng+10<511
      GOTO, try_again
   ENDIF
   
   fwhm = wave[flux_right] - wave[flux_left] ;; FWHM

   break_out: ;; jump here when exiting the above loop due half_max points not being found
   
   ;OPLOT, [wave[flux_left],wave[flux_left]], [flux[flux_left],flux[flux_left]], PSYM=4, COLOR=200
   ;OPLOT, [wave[flux_right],wave[flux_right]], [flux[flux_right],flux[flux_right]],PSYM=4, COLOR=200

ENDIF ELSE BEGIN ;; 2 peaks
   ;; switches used to mark which peak was properly fit in the event
   ;; that one peak does not produce a satisfactory fit
   blue = 0
   red = 0

   ;;define range of indices including all points to the left and 
   ;;right of peak used in weight
   range_1 = INDGEN(best_r_trough_ind[0]-best_l_trough_ind[0]) + best_l_trough_ind[0]
   range_2 = INDGEN(best_r_trough_ind[1]-best_l_trough_ind[1]) + best_l_trough_ind[1]
   
   ;; LEFT PEAK
   ;; determine the proper interval for a fit and fit a gaussian to 
   ;; the top of the peak to determine the true peak flux and
   ;; wavelength                                               
   num_pts = 4
   REPEAT BEGIN
      num_pts += 1
      interval_1 = n_extrema(flux[range_1], num_pts, 1)+range_1[0]
      peak_fit_1 = GAUSSFIT(wave[interval_1], flux[interval_1], CHISQ=chisq, gauss_para_1, NTERMS=3)
   ENDREP UNTIL gauss_para_1[1] GE wave[interval_1[0]] && gauss_para_1[1] LE wave[interval_1[-1]] $
      || num_pts GE N_ELEMENTS(flux[range_1])
   
   ;; make sure the fits are correct                                                                                
   IF num_pts GE N_ELEMENTS(flux[range_1]) && (gauss_para_1[1] LT wave[interval_1[0]] || $
                                               gauss_para_1[1] GT wave[interval_1[-1]]) THEN BEGIN
      fit_err = 1 ;; set switch to print warning message
      flag[1] = 2
      peak1_flux = 0D
      peak1_wave = 0D
      peak_num -= 1
      symmetry = 0D
      center= 6.28D
   ENDIF ELSE BEGIN
      red = 1
      peak1_flux = gauss_para_1[0]
      peak1_wave = gauss_para_1[1]
   ENDELSE

   ;; RIGHT PEAK
   num_pts = 4
   REPEAT BEGIN
      num_pts += 1
      interval_2 = n_extrema(flux[range_2], num_pts, 1)+range_2[0]
      peak_fit_2 = GAUSSFIT(wave[interval_2], flux[interval_2], CHISQ=chisq, gauss_para_2, NTERMS=3)
   ENDREP UNTIL gauss_para_2[1] GE wave[interval_2[0]] && gauss_para_2[1] LE wave[interval_2[-1]] $
      || num_pts GE N_ELEMENTS(flux[range_2])
   
   ;; make sure the fits are correct
   IF num_pts GE N_ELEMENTS(flux[range_2]) && (gauss_para_2[1] LT wave[interval_2[0]] || $
                                               gauss_para_2[1] GT wave[interval_2[-1]]) THEN BEGIN
      fit_err = 1 ;; set switch to print warning message
      flag[1] = 2
      peak2_flux = 0D
      peak2_wave = 0D
      peak_num -= 1
      symmetry = 0D
      center= 6.28D
   ENDIF ELSE BEGIN
      blue = 1
      peak2_flux = gauss_para_2[0]
      peak2_wave = gauss_para_2[1]
   ENDELSE
   
   IF peak_num EQ 0 THEN BEGIN
      GOTO,bypass
   ENDIF ELSE IF peak_num EQ 1 THEN BEGIN
      print_sing_peak = 1 ;; set switch to print single peak data
      IF red EQ 1 THEN BEGIN
         gauss_para = gauss_para_1
         peak_fit = peak_fit_1
         interval = interval_1
         peak_flux = peak1_flux
         peak_wave = peak1_wave

         ;; now that this is a single peak, calculate FWHM
         l_rng = INDGEN(best_peak_ind[0]-best_l_trough_ind[0])+best_l_trough_ind[0]
         r_rng = INDGEN(best_r_trough_ind[0]-best_peak_ind[0])+best_peak_ind[0]

         go_again_1: ;;program jumps here if half_max points are incorrect 
         IF N_ELEMENTS(WHERE(l_rng EQ 0)) EQ N_ELEMENTS(l_rng) || $
            N_ELEMENTS(WHERE(r_rng EQ 511)) EQ N_ELEMENTS(r_rng) THEN BEGIN
            fwhm_err = 1 ;; set switch to print warning message
            fwhm = 0.0
            GOTO, break_out_1
         ENDIF

         IF flux[l_rng[0]] GT 1 && flux[r_rng[-1]] GT 1 THEN BEGIN
            half_max = (peak_flux-1D)/2D + 1D ;;flux value for half peak flux
         ENDIF ELSE BEGIN
            avg_base = (flux[l_rng[0]] + flux[r_rng[-1]])/2D
            half_max = (peak_flux-avg_base)/2D + avg_base
         ENDELSE
         
         flux_left = VALUE_LOCATOR(flux[l_rng],half_max) + l_rng[0] ;; index of point left of peak at half_max
         flux_right = VALUE_LOCATOR(flux[r_rng],half_max) + r_rng[0] ;; index of point right of peak at half_max 
         
         ;;define a threshold (min_sep) depending on peak_flux value
         IF peak_flux GT 10 THEN min_sep = 2.0 ELSE min_sep = 0.8

         IF ABS(flux[flux_left]-half_max) GT min_sep THEN BEGIN
            l_rng = l_rng-10>0
            GOTO, go_again_1
         ENDIF
         IF ABS(flux[flux_right]-half_max) GT min_sep THEN BEGIN
            r_rng = r_rng+10<511
            GOTO, go_again_1
         ENDIF
         
         fwhm = wave[flux_right] - wave[flux_left] ;; FWHM
         
         break_out_1:

         ;OPLOT, [wave[flux_left],wave[flux_left]], [flux[flux_left],flux[flux_left]], PSYM=4, COLOR=200
         ;OPLOT, [wave[flux_right],wave[flux_right]], [flux[flux_right],flux[flux_right]],PSYM=4, COLOR=200

      ENDIF ELSE IF blue EQ 1 THEN BEGIN
         gauss_para = gauss_para_2
         peak_fit = peak_fit_2
         interval = interval_2
         peak_flux = peak2_flux
         peak_wave = peak2_wave
         ;; now that this is a single peak, calculate FWHM
         l_rng = INDGEN(best_peak_ind[0]-best_l_trough_ind[0])+best_l_trough_ind[0]
         r_rng = INDGEN(best_r_trough_ind[0]-best_peak_ind[0])+best_peak_ind[0]

         go_again_2: ;;program jumps here if half_max points are incorrect
         IF N_ELEMENTS(WHERE(l_rng EQ 0)) EQ N_ELEMENTS(l_rng) || $
            N_ELEMENTS(WHERE(r_rng EQ 511)) EQ N_ELEMENTS(r_rng) THEN BEGIN
            fwhm_err = 1 ;; set switch to print warning message
            fwhm = 0.0
            GOTO, break_out_2
         ENDIF

         IF flux[l_rng[0]] GT 1 && flux[r_rng[-1]] GT 1 THEN BEGIN
            half_max = (peak_flux-1D)/2D + 1D ;;flux value for half peak flux 
         ENDIF ELSE BEGIN
            avg_base = (flux[l_rng[0]] + flux[r_rng[-1]])/2D
            half_max = (peak_flux-avg_base)/2D + avg_base
         ENDELSE

         flux_left = VALUE_LOCATOR(flux[l_rng],half_max) + l_rng[0] ;; index of point left of peak at half_max 
         flux_right = VALUE_LOCATOR(flux[r_rng],half_max) + r_rng[0] ;; index of point right of peak at half_max 
         
         ;;define a threshold (min_sep) depending on peak_flux value
         IF peak_flux GT 10 THEN min_sep = 2.0 ELSE min_sep = 0.8

         IF ABS(flux[flux_left]-half_max) GT min_sep THEN BEGIN
            l_rng = l_rng-10>0
            GOTO, go_again_2
         ENDIF
         IF ABS(flux[flux_right]-half_max) GT min_sep THEN BEGIN
            r_rng = r_rng+10<511
            GOTO, go_again_2
         ENDIF

         fwhm = wave[flux_right] - wave[flux_left] ;; FWHM

         break_out_2:

         ;OPLOT, [wave[flux_left],wave[flux_left]], [flux[flux_left],flux[flux_left]], PSYM=4, COLOR=200
         ;OPLOT, [wave[flux_right],wave[flux_right]], [flux[flux_right],flux[flux_right]],PSYM=4, COLOR=200

      ENDIF
   ENDIF ELSE BEGIN
      print_dbl_peak = 1 ;; set switch to print double peak data
      
      ;;***Test Plotting***  
      ;OPLOT, wave[interval_1], peak_fit_1, PSYM=-3, COLOR=150
      ;OPLOT, wave[interval_2], peak_fit_2, PSYM=-3, COLOR=150

      ;;determine central_depth between peaks and average peak value 
      left_ind = VALUE_LOCATOR(wave, peak1_wave)
      right_ind = VALUE_LOCATOR(wave, peak2_wave)>left_ind+1
      central_depth = MIN(flux[left_ind : right_ind])
      avg_peak = (peak1_flux + peak2_flux)/2.0
      
      ;;determine shell parameter and V/R ratio  
      shell_para = avg_peak/central_depth
      vr_ratio = ALOG10(peak1_flux/peak2_flux)

   ENDELSE 

   ;OPLOT, wave[interval_1], peak_fit_1, PSYM=-3, COLOR=150
   ;OPLOT, wave[interval_2], peak_fit_2, PSYM=-3, COLOR=150

ENDELSE

;;==========================================================================================;;
;;=============================== END FITTING THE PEAKS ====================================;;
;;==========================================================================================;;






;;==========================================================================================;;
;;========================= CALCULATING SYMMETRY AND WING FIT ==============================;;
;;==========================================================================================;;
;;
;; This section is used to determine the central wavelength and level
;; of asymmetry(symmetry) of an emission spectrum. This is done by
;; fitting the wings of emission with a gaussian and using its peak
;; data for the central wavelength. Symmetry is determined in the
;; following ways
;;
;; 1. Single peak spectra - The peak wavelength is divided by the
;;                          central wavelength
;;
;; 2. Double peak spectra - The distance from peak 1 to the
;;                          central wavelength is divided by the
;;                          distance from peak 2 to the central
;;                          wavelength    
;;
;; Using this method, symmetry values closer to 1 are considered
;; more symmetric
;;
;; The wing fits do not always work properly and sometimes give the
;; wrong central wavelength. These issues usually come from asymmetric
;; wings or the wings chosen are too small to give a proper fit.
;; 
;; Other failures come from the gauss_fit function not being able to 
;; fit the wings at all. These issues are rare but do occur. Spectra 
;; exhibiting this issue will be flagged 
;; 

;; this smoothed flux will be used in the fitting as to make the plot
;; less noisy to ease the fitting
wsmooth_flux = smooth(flux,20)
start_scan = VALUE_LOCATOR(wave, 654.0)

;;determine if the spectrum contains absorption
smooth_der = DERIV(wave, wsmooth_flux)
sum_der = TOTAL(smooth_der[40:VALUE_LOCATOR(wave, 655.0)])
IF sum_der GE -1.4 THEN abs = 0 ELSE abs = 1 ;;abs is switch for absorption
 
;; define the symmetry variable 
symmetry = 0D

retry:
IF peak_num EQ 2 || peak_num EQ 1 THEN BEGIN
   IF ~abs THEN BEGIN
      ;; no absorption
      
      ;; alter the percentage of the peak used for wing fit based on peak flux
      IF MAX(flux[best_peak_ind]) GE 2 THEN k = 0.3D ELSE k = 0.55D
          
      ;; get indecies for the wings
      wing_ind = WHERE(wsmooth_flux[start_scan :*] LE 1D + (k * (MAX(flux[best_peak_ind]) - $
                                                                     1D))) + start_scan
      ;; check for breaks in continuity
      j = 0
      brk_ind = INTARR(N_ELEMENTS(wing_ind))
      FOR i=0, N_ELEMENTS(wing_ind)-2 DO BEGIN
         IF wing_ind[i]+1 NE wing_ind[i+1] THEN BEGIN
            brk_ind[j] = i
            j += 1
            brk_ind[j] = i+1
            j += 1
         ENDIF
      ENDFOR
      ;; fit the wings
      IF peak_num EQ 2 THEN BEGIN
         IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 2 THEN BEGIN 
            ;; dip doesnt go below 1D+(k*(MAX(flux[best_peak_ind])-1D))
            
            w_interval = wing_ind
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[1,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            
            ;; make sure the center is reasonable 
            IF center GT 8D || center LT 4D THEN BEGIN
               center = 6.28D
               improper_wing_err = 1 ;; set switch to print warning message
               flag[1] = 2
            ENDIF
            
            one_frm_cent = ABS(center - (gauss_para_1[1] - 650D))
            two_frm_cent = ABS(center - (gauss_para_2[1] - 650D))
            symmetry = one_frm_cent/two_frm_cent
            ;PRINT, '   Symmetry of peak(s):', STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
         ENDIF ELSE IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 4 THEN BEGIN 
            ;; dip goes below 1D+(k*(MAX(flux[best_peak_ind])-1D))

            w_interval = [wing_ind[0:brk_ind[0]], wing_ind[brk_ind[3]:*]]
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[2,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            
            ;; make sure the center is reasonable
            IF center GT 8D || center LT 4D THEN BEGIN
               center = 6.28D
               improper_wing_err = 1 ;; set switch to print warning message
               flag[1] = 2
            ENDIF

            one_frm_cent = ABS(center - (gauss_para_1[1] - 650D))
            two_frm_cent = ABS(center - (gauss_para_2[1] - 650D))
            symmetry = one_frm_cent/two_frm_cent
            ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
         ENDIF ELSE BEGIN
            wing_err = 1 ;; set switch to print warning message
            flag[1] = 2
            center = 6.28D
         ENDELSE
      ENDIF ELSE IF peak_num EQ 1 THEN BEGIN
         w_interval = wing_ind
         gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                ESTIMATES=[2,656.28,1,1],NTERMS=4)
         center = wgauss_para[1] - 650D
         
         ;; make sure the center is reasonable                                                       
         IF center GT 8D || center LT 4D THEN BEGIN
            center = 6.28D
            improper_wing_err = 1 ;; set switch to print warning message
            flag[1] = 2
         ENDIF
         
         symmetry = ABS(gauss_para[1] - 650D)/center
         ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
         ;PRINT, ' '
      ENDIF

   ENDIF ELSE BEGIN ;; absorption

      ;; get indecies for the wings
      min_rng = INDGEN(best_r_trough_ind[0]-best_l_trough_ind[0]+20)+(best_l_trough_ind[0]-10)
      min_for_wings = 1D - ABS(0.65D * (1D - MIN(flux[min_rng]))) 
      wing_ind = WHERE(wsmooth_flux[start_scan :*] GE min_for_wings) + start_scan

      ;; check for breaks in continuity                  
      j = 0
      brk_ind = INTARR(N_ELEMENTS(wing_ind))                                                
      FOR i=0, N_ELEMENTS(wing_ind)-2 DO BEGIN
         IF wing_ind[i]+1 NE wing_ind[i+1] THEN BEGIN
            brk_ind[j] = i 
            j += 1
            brk_ind[j] = i+1
            j += 1
         ENDIF
      ENDFOR

      ;; make sure the whole spectrum wasn't included and if it
      ;; was, treat it as an emission with no absorption
      IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 1 && WHERE(brk_ind NE 0) EQ -1 THEN BEGIN
         abs = 0 
         GOTO, retry
      ENDIF
 
      ;; fit the wings
      IF peak_num EQ 2 THEN BEGIN
         IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 2 THEN BEGIN 
            w_interval = wing_ind
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[-2,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            
            ;; make sure the center is reasonable  
            IF center GT 8D || center LT 4D THEN BEGIN
               center = 6.28D
               improper_wing_err = 1 ;; set switch to print warning message
               flag[1] = 2
            ENDIF
            
            one_frm_cent = ABS(center - (gauss_para_1[1] - 650D))
            two_frm_cent = ABS(center - (gauss_para_2[1] - 650D))
            symmetry = one_frm_cent/two_frm_cent
            ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
         ENDIF ELSE IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 4 THEN BEGIN 
            w_interval = [wing_ind[0:brk_ind[0]], wing_ind[brk_ind[3]:*]]
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[-2,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            
            ;; make sure the center is reasonable  
            IF center GT 8D || center LT 4D THEN BEGIN
               center = 6.28D
               improper_wing_err = 1 ;; set switch to print warning message
               flag[1] = 2
            ENDIF

            one_frm_cent = ABS(center - (gauss_para_1[1] - 650D))
            two_frm_cent = ABS(center - (gauss_para_2[1] - 650D))
            symmetry = one_frm_cent/two_frm_cent
            ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
         ENDIF ELSE IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 6 THEN BEGIN 
            w_interval = [wing_ind[0:brk_ind[0]], wing_ind[brk_ind[5]:*]]
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[-2,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            
            ;; make sure the center is reasonable    
            IF center GT 8D || center LT 4D THEN BEGIN
               center = 6.28D
               improper_wing_err = 1 ;; set switch to print warning message
               flag[1] = 2
            ENDIF
            
            one_frm_cent = ABS(center - (gauss_para_1[1] - 650D))
            two_frm_cent = ABS(center - (gauss_para_2[1] - 650D))
            symmetry = one_frm_cent/two_frm_cent
            ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
        
         ENDIF ELSE BEGIN
            wing_err = 1 ;; set switch to print warning message
            center = 6.28D
            flag[1] = 2
         ENDELSE
         
      ENDIF ELSE IF peak_num EQ 1 THEN BEGIN
         IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 4 THEN BEGIN 
            w_interval = [wing_ind[0:brk_ind[0]], wing_ind[brk_ind[3]:*]]
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[-2,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            symmetry = ABS(gauss_para[1] - 650D)/center
            ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
         ENDIF ELSE IF N_ELEMENTS(WHERE(brk_ind NE 0)) EQ 2 THEN BEGIN
            w_interval = [wing_ind[0:brk_ind[0]], wing_ind[brk_ind[1]:*]]
            gauss_wings = GAUSSFIT(wave[w_interval],wsmooth_flux[w_interval],wgauss_para, $
                                   ESTIMATES=[-2,656.28,1,1],NTERMS=4)
            center = wgauss_para[1] - 650D
            symmetry = ABS(gauss_para[1] - 650D)/center
            ;PRINT,'   Symmetry of peak(s):',STRCOMPRESS(STRING(symmetry))
            ;PRINT, ' '
         ENDIF ELSE BEGIN
            ;PRINT, '   **WARNING: Symmetry Calculation Error**'
            center = 6.28D
            symmetry = 0D
         ENDELSE
      ENDIF
   ENDELSE
ENDIF ELSE BEGIN
   ;PRINT, '   **WARNING: Symmetry Calculation Error**'
   flag[1] = 2
   center = 6.28D
   symmetry = 0D
ENDELSE


;;test plotting
;OPLOT, wave[w_interval], flux[w_interval], PSYM=-3, COLOR=150
;k = WHERE(flux EQ MIN(flux[min_rng]))
;OPLOT, [wave[k],wave[k]], [flux[k], flux[k]], PSYM=4
;OPLOT, wave[w_interval], gauss_wings, PSYM=-3, COLOR=250
;stop

;;==========================================================================================;;
;;====================== END CALCULATING SYMMETRY AND WING FIT =============================;;
;;==========================================================================================;;






;;==========================================================================================;;
;;================================= PRINTING OUTPUT ========================================;;
;;==========================================================================================;;
;;
;; see what messages needed to be printed based on switches. Also
;; determine if output should go to log file or terminal 
IF fit_err THEN BEGIN
   text = '**WARNING: Peak fitting malfunction. Incorrect fit.**'
   IF spec_plot_alone THEN BEGIN ;; print to terminal
      PRINT, ' ' 
      PRINT, text
      PRINT, ' '
   ENDIF ELSE BEGIN ;;print to log file
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text]
      analysis_logs = [analysis_logs,'']
   ENDELSE
ENDIF 
IF print_sing_peak THEN BEGIN
   text1 = '   Spectrum has one peak.'
   text2 = '   Max flux of peak:'+STRCOMPRESS(STRING(gauss_para[0]))
   text3 = '   The peak occurs at:'+STRCOMPRESS(STRING(gauss_para[1]))+' nm'
   IF spec_plot_alone THEN BEGIN ;; print to terminal
      PRINT, ' '
      PRINT, text1
      PRINT, ' '
      PRINT, text2
      PRINT, text3
      PRINT, ' '
   ENDIF ELSE BEGIN ;; print to log file
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text1]
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text2]
      analysis_logs = [analysis_logs,text3]
      analysis_logs = [analysis_logs,'']
   ENDELSE
ENDIF
IF print_dbl_peak THEN BEGIN
   text1 = '   Spectrum has two peaks.'
   text2 = '   Max flux of blue peak:'+STRCOMPRESS(STRING(peak1_flux))
   text3 = '   The blue peak occurs at:'+STRCOMPRESS(STRING(peak1_wave))+' nm'
   text4 = '   Max flux of red peak:'+STRCOMPRESS(STRING(peak2_flux))
   text5 = '   The red peak occurs at:'+STRCOMPRESS(STRING(peak2_wave))+' nm'
   IF spec_plot_alone THEN BEGIN ;; print to terminal
      PRINT, ' '
      PRINT, text1
      PRINT, ' '
      PRINT, text2
      PRINT, text3
      PRINT, ' '
      PRINT, text4
      PRINT, text5
      PRINT, ' '
   ENDIF ELSE BEGIN ;; print to log file
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text1]
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text2]
      analysis_logs = [analysis_logs,text3]
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text4]
      analysis_logs = [analysis_logs,text5]
      analysis_logs = [analysis_logs,'']
   ENDELSE
ENDIF
IF fwhm_err THEN BEGIN
   text = '**WARNING: FWHM ERROR. Correct half-max points could not be found.**'
   IF spec_plot_alone THEN BEGIN ;; print to terminal 
      PRINT, ' '
      PRINT, text
      PRINT, ' '
   ENDIF ELSE BEGIN ;; print to log file 
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text]
      analysis_logs = [analysis_logs,'']
   ENDELSE
ENDIF
IF improper_wing_err THEN BEGIN
   text = '** WARNING: IMPROPER WING FIT **'
   IF spec_plot_alone THEN BEGIN ;; print to terminal
      PRINT, ' '
      PRINT, text
      PRINT, ' '
   ENDIF ELSE BEGIN ;; print to log file
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text]
      analysis_logs = [analysis_logs,'']
   ENDELSE
ENDIF
IF wing_err THEN BEGIN
   text = '**WARNING: WING FIT ERROR**'
   IF spec_plot_alone THEN BEGIN ;; print to terminal  
      PRINT, ' '
      PRINT, text
      PRINT, ' '
   ENDIF ELSE BEGIN ;; print to log file
      analysis_logs = [analysis_logs,'']
      analysis_logs = [analysis_logs,text]
      analysis_logs = [analysis_logs,'']
   ENDELSE
ENDIF

;;==========================================================================================;;
;;============================= END PRINTING OUTPUT ========================================;;
;;==========================================================================================;;






;;==========================================================================================;;
;;=========================== HANDLE ABSORPTION JUMP =======================================;;
;;==========================================================================================;;

jump: ;; program 'jumps' here if absorption spectrum is detected 
IF Absorption EQ 1 THEN BEGIN
   IF spec_plot_alone THEN BEGIN ;; print to terminal
      PRINT, '   ** Absorption spectrum. No emissions detected **'
      PRINT, ' '
   ENDIF ELSE BEGIN ;; print to log file
      analysis_logs = [analysis_logs,'   ** Absorption spectrum. No emissions detected **']
      analysis_logs = [analysis_logs,'']
   ENDELSE
      
   symmetry = 0D ;; define symmetry to not have star_report fail   
   center = 6.28D ;; define center to not have star_report fail and have the plots be 
                  ;; unaffected by halpha_centralizer
ENDIF

bypass: ;; program 'jumps' here if the peak fit malfunctions

;;==========================================================================================;;
;;======================== END HANDLE ABSORPTION JUMP ======================================;;
;;==========================================================================================;;






;;==========================================================================================;;
;;============================= GENERATE DIAGNOSTIC REPORT =================================;;
;;==========================================================================================;;

IF ~diagnostics || spec_plot_alone THEN GOTO, skip

;; set device parameters
;!P.FONT = -1
SET_PLOT, 'ps'
DEVICE, /COLOR, BITS_PER_PIXEL=8
DEVICE, FILENAME='peak_analysis_'+temp_file+'.ps', ENCAPSULATED=0,$
        FONT_SIZE=12, XSIZE=8, YSIZE=11, XOFFSET=0.5, YOFFSET=0, /INCHES

;;get index for proper obs date
this_date = WHERE(analysis_frames EQ temp_file)

;; plot the spectrum
PLOT, wave, flux, XTITLE='Wavelength (nm)', YTITLE='F/Fc', XRANGE = [653, 659], XSTYLE=1, $
      CHARSIZE=1.2, CHARTHICK=3, XTHICK=3, YTHICK=3, THICK=3, $
      TITLE='Frame = '+temp_file+' / Obs. Date = '+analysis_dates[this_date], $
      POS=[.1, 0.5, 0.95, 0.9], /DATA, /NOERASE
 
;; overplot a dashed line at the continuum level                                       
OPLOT, [653.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=color_table[1]

;; plot the peak fit for single peak spectra  
IF ISA(interval) && ISA(peak_fit) && peak_num EQ 1 THEN $
   OPLOT, wave[interval], peak_fit, PSYM=-3, COLOR=color_table[2], THICK=4

;; plot the peak fits with multiple peaks 
IF ISA(interval_1) && ISA(interval_2) && ISA(peak_fit_1) && ISA(peak_fit_2) && peak_num EQ 2 THEN BEGIN
   OPLOT, wave[interval_1], peak_fit_1, PSYM=-3, COLOR=color_table[2], THICK=4
   OPLOT, wave[interval_2], peak_fit_2, PSYM=-3, COLOR=color_table[2], THICK=4
ENDIF
   
;; plot wing fit       
IF ISA(w_interval) && ISA(gauss_wings) && peak_num NE 0 THEN $
   OPLOT, wave[w_interval], gauss_wings, PSYM=-3, COLOR=color_table[4],THICK=2.5

;; plot central wavelength
IF ISA(center) && peak_num NE 0 THEN $
OPLOT, [center+650.0, center+650.0], [0,100], PSYM=-3, COLOR=color_table[5], THICK=2.5  


;; add useful peak data to the report
XYOUTS, 0.03, 0.01, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
XYOUTS, 0.5, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.20, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.5, 0.9625, $
        '_________________________________________________________________________________________________',$
        ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL


IF ISA(peak_flux) && peak_num EQ 1 THEN $
   XYOUTS, 0.15, 0.4, 'Peak Flux: '+STRING(peak_flux, FORMAT='(F5.2)'), CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL 

IF ISA(peak_wave) && peak_num EQ 1 THEN $
   XYOUTS, 0.15, 0.38, 'Peak Wavelength: '+STRING(peak_wave, FORMAT='(F7.2)')+' nm', CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL

IF ISA(peak1_flux) && peak_num EQ 2 THEN $
   XYOUTS, 0.15, 0.40, 'Peak Flux: '+STRING(peak1_flux, FORMAT='(F5.2)'), CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL

IF ISA(peak1_wave) && peak_num EQ 2 THEN $
   XYOUTS, 0.15, 0.38, 'Peak Wavelength: '+STRING(peak1_wave, FORMAT='(F7.2)')+' nm', CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL

IF ISA(peak2_flux) && peak_num EQ 2 THEN $
   XYOUTS, 0.15, 0.36, 'Peak Flux: '+STRING(peak2_flux, FORMAT='(F5.2)'), CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL

IF ISA(peak2_wave) && peak_num EQ 2 THEN $
   XYOUTS, 0.15, 0.34, 'Peak Wavelength: '+STRING(peak2_wave, FORMAT='(F7.2)')+' nm', CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL

IF ISA(num_pts) && peak_num NE 0 THEN BEGIN
   XYOUTS, 0.55, 0.38, 'num_pts = '+STRING(num_pts, FORMAT='(I3)'), CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL
ENDIF ELSE IF peak_num NE 0 THEN BEGIN
   XYOUTS, 0.55, 0.38, 'num_pts = N/A', CHARSIZE=1.2, CHARTHICK=2, /NORMAL
ENDIF

IF ISA(center) && peak_num NE 0 THEN $
   XYOUTS, 0.55, 0.4, 'Central Wavelength: '+STRING((center+650.0), FORMAT='(F7.2)')+' nm', CHARSIZE=1.2, $
           CHARTHICK=2, /NORMAL 

XYOUTS, 0.1, 0.22, 'Green line shows Gaussian fit applied to peak(s)', $
        CHARSIZE=1.2, CHARTHICK=2, /NORMAL
XYOUTS, 0.1, 0.2, 'Violet line shows Gaussian fit applied to wings.', $ 
        CHARSIZE=1.2, CHARTHICK=2, /NORMAL

IF ISA(center) && peak_num NE 0 THEN $
   XYOUTS, 0.1, 0.18, 'Vertical blue line shows the central wavelength given by the wing fit.',$
           CHARSIZE=1.2, CHARTHICK=2, /NORMAL


DEVICE, /CLOSE_FILE
SET_PLOT, 'X'

;;move the file to the diagonostics directory for that star 
SPAWN, 'mv -f '+'peak_analysis_'+temp_file+'.ps ' + og_loc + targets_dir_name + '/HD' + $
       star_hd + '/Diagnostics/'

;;==========================================================================================;;
;;========================= END GENERATE DIAGNOSTIC REPORT =================================;;
;;==========================================================================================;;



skip: ;; jump here if skipping diagnostic report

END
