;;======================================================================;;
;; Star Report Pipeline, Central Michigan University
;;
;; cross_correlator.pro
;;
;; The purpose of this procedure is to shift the wavelength array for
;; each spectrum to align it with the first recorded spectrum when all
;; spectra are plotted together. This is done using cross correlation.
;;
;;
;; Author: Christian Hannah
;;
;;======================================================================;;


PRO cross_correlator, waves, fluxs, num_analyzed, shifts, SORTED_JDATES=sorted_jdates
COMMON spec_block
COMMON star_block


;;define the number of data points per spectrum
data_num = 512

;;define the starting and ending values for the model wavelength array used in resampling
start_wave = 655.0
end_wave = 658.0

;; switch for diagnostic plots 
diagnostics = 1




;;===================================================================================;;
;;===================== RESAMPLE FLUXS W/ CONSTANT SPACING ==========================;;
;;===================================================================================;;
;;
;; Because cross correlation relies on equally spaced data points, we
;; must first construct a model wavelength array and resample each
;; flux array with the model wavelength array. The resampling is done
;; using the INTERPOL IDL procedure which utilizes linear
;; interpolation.
;;
;;*** ATTENTION: *****************************************************
;;* It is crucial to resample all flux arrays with the same model    *
;;* wavelength array as cross correlation is blind to the wavelength *
;;* scale and only cares about index                                 *
;;********************************************************************
;;

;; specify the number of desired points in sample wavelength
i_num = 512 

;; take the log of the wavelength arrays.
l_waves = ALOG10(waves[[0:data_num-1],*])

;;define equally-spaced model wavelength array
l_start = ALOG10(start_wave)
l_end = ALOG10(end_wave)
les_waves = (DINDGEN(i_num)/DOUBLE(i_num-1))*(l_end-l_start)+l_start

;;get wavelength spacing
spacing = (les_waves[i_num-1]-les_waves[0])/(i_num-1)

;;resample (using interpolation) fluxs with es_waves
es_fluxs = DBLARR(i_num, num_analyzed) 
FOR i=0, num_analyzed-1 DO BEGIN
   es_fluxs[*,i] = INTERPOL(fluxs[*,i], l_waves[*,i], les_waves)
ENDFOR     

;;===================================================================================;;
;;=================== END RESAMPLE FLUXS W/ CONSTANT SPACING ========================;;
;;===================================================================================;;






;;===================================================================================;;
;;=========================== CONDUCT CROSS CORRELATION =============================;;
;;===================================================================================;;

;; define some arrays needed
cc_shifts = DBLARR(num_analyzed) ;;array to contain the shift values (in lags) for each spectrum

sl_waves = DBLARR(data_num, num_analyzed) ;;matrix to contain shifted l_waves

ccs = DBLARR(i_num, num_analyzed) ;;matrix to contain the cross correlation values 

lags = FIX(DINDGEN(i_num)-(i_num/2)) ;;lags to be used in cross correlation

;;plot each edited spectrum one by one if desired
;FOR i=0, num_analyzed-1 DO BEGIN
;   PLOT, les_waves, es_fluxs[*,i], PSYM=-3
;   stop
;ENDFOR

;; cross correlate the resampled fluxs and store results in ccs. **Each
;; row of ccs represents a seperate cross correlation**.
FOR i=0, num_analyzed-1 DO BEGIN
   ccs[*,i] = C_CORRELATE(es_fluxs[*,i],es_fluxs[*,0], lags)
ENDFOR 

;; plot the ccs vs. lags one by one if desired
;FOR i=0, num_analyzed-1 DO BEGIN
;   PLOT,  lags, ccs[*,i], PSYM=-3
;   stop
;ENDFOR

;;===================================================================================;;
;;======================== END CONDUCT CROSS CORRELATION ============================;;
;;===================================================================================;;






;;===================================================================================;;
;;==================== OBTAIN PEAK AND TROUGH INDICES FOR CCS =======================;;
;;===================================================================================;;
;;
;; This section of the procedure will obtain indices for every local
;; max and min in each cross correlation.
;;

;;define lag range to search for cross correlation peaks
start_ind = MIN(WHERE(lags GE -250))
end_ind = MAX(WHERE(lags LE 250))

;;define a matrix to hold the indices of elements in each cross
;;correlation corresponding to a negative slope. (Each row represents
;;a cross correlation) 
neg_ind = INTARR(N_ELEMENTS(lags[start_ind : end_ind]), num_analyzed)

;; populate neg_ind
FOR j=0, num_analyzed-1 DO BEGIN
   f=0
   FOR i=start_ind, end_ind-1 DO BEGIN
      IF ccs[i+1,j] LT ccs[i,j] || ccs[i+2,j] LT ccs[i+1,j] THEN BEGIN
         neg_ind[f,j] = i+1
         f++
      ENDIF
   ENDFOR 
ENDFOR

;; check for breaks in continuity in each row and store indices associated with
;; discontinuity as temp_break_ind 
temp_break_ind = INTARR(N_ELEMENTS(lags[start_ind : end_ind]), num_analyzed) 
FOR j=0, num_analyzed-1 DO BEGIN
   ;;the first neg_ind for each cross correlation corresponds to a
   ;;"peak" so store it as a temp_break_ind
   temp_break_ind[0,j] = neg_ind[0,j]
   f=1
   FOR i=1, end_ind-start_ind-1 DO BEGIN
      IF neg_ind[i+1,j]-1 NE neg_ind[i,j] && neg_ind[i,j] NE 0 THEN BEGIN
         ;;add both indecies to temp_break_ind
         temp_break_ind[f,j] = neg_ind[i,j] 
         temp_break_ind[f+1,j] = neg_ind[i+1,j]
         f += 2
      ENDIF
   ENDFOR
ENDFOR

;;remove any duplicates from each row of temp_break_ind (the method used to fill the
;;matrix may reslt in duplicates) by creating a new matrix called break_ind
break_ind = INTARR(N_ELEMENTS(lags[start_ind : end_ind]), num_analyzed)
FOR j=0, num_analyzed-1 DO BEGIN
   f=0
   FOR i=0, end_ind-start_ind-1 DO BEGIN
      IF temp_break_ind[i,j] NE temp_break_ind[i+1,j] THEN BEGIN
         break_ind[f,j] = temp_break_ind[i,j]
         f += 1
      ENDIF
   ENDFOR
ENDFOR

;;seperate break_ind into peak_ind and trough_ind as each index
;;corresponds to a peak or trough
trough_ind = INTARR(N_ELEMENTS(lags[start_ind : end_ind]),num_analyzed)
peak_ind = INTARR(N_ELEMENTS(lags[start_ind : end_ind]), num_analyzed) 
FOR j=0, num_analyzed-1 DO BEGIN
   n=0
   FOR i=0, end_ind-start_ind-1, 2 DO BEGIN
      peak_ind[n,j] = break_ind[i,j]
      trough_ind[n,j] = break_ind[i+1,j]
      n += 1
   ENDFOR
ENDFOR

;;===================================================================================;;
;;================== END OBTAIN PEAK AND TROUGH INDICES FOR CCS =====================;;
;;===================================================================================;;






;;===================================================================================;;
;;=============================== DETERMINE SHIFTS ==================================;;
;;===================================================================================;;
;;
;; To determine the shift needed for each spectrum, an initial gauss fit is
;; first applied to each cross correlation in the lag range of
;; [-80:80]. In the event that there are no peak_ind within the lag
;; range secified above, the shift from this initial gauss fit will be
;; used if the peak of the fit lies within the lag range. Otherwise, a
;; shift of 0 is assigned.
;;
;; To select the peak for the shift, the list of peak indices
;; (peak_ind) is narrowed to conatain only peaks within the lag range.
;; The "best" peak is then determined by comparing ccs values for each
;; peak. The largest ccs value will be considered as the peak. A
;; Gaussian fit is then applied to the top of this peak to obtain the
;; correct shift. If the fit is unsuccessful, the original peak index
;; is used to extract a shift.
;;

;;define some needed values and arrays
no_shift_ind = WHERE(lags EQ 0) ;;index for a shift of 0
num_pts = 6 ;;number of points to include in the gauss fit of the CC peak
gauss_fits = DBLARR(num_pts,num_analyzed) ;;matrix to store gauss fits of the CC peaks

;;redefine start_ind and end_ind for a more narrow range to apply
;;initial gaussfit
start_ind = MIN(WHERE(lags GE -80))
end_ind = MAX(WHERE(lags LE 80))

;;define a matrix to store initial gauss fits of the entire CC plot 
init_gauss_fits = DBLARR(N_ELEMENTS(lags[start_ind :end_ind]),num_analyzed) 

;;array to hold shift indices for plotting
shift_inds = INTARR(num_analyzed)

;;determine shifts 
FOR i=0, num_analyzed-1 DO BEGIN
   
   ;;apply an intial gaussfit of the entire CC plot to obtain an
   ;;initial ("guess") shift 
   init_gauss_fits[*,i] = GAUSSFIT(lags[start_ind :end_ind], ccs[[start_ind :end_ind],i], init_gauss_para, $
                                   CHISQ=chisq, ESTIMATES=[100.0,0,5,0], NTERMS=4)
   
   ;;store index of the shift
   init_shift_index = VALUE_LOCATOR(lags, init_gauss_para[1])
   
   ;; narrow list of peak_ind to include only those within the lag
   ;; range specified by start_ind and end_ind
   j=0
   possible_peaks_ind = INTARR(200)
   FOR q=0, N_ELEMENTS(peak_ind[*,i])-1 DO BEGIN
      IF peak_ind[q,i] GT start_ind && peak_ind[q,i] LT end_ind THEN BEGIN
         possible_peaks_ind[j] = peak_ind[q,i]
         j+=1
      ENDIF
   ENDFOR

   ;;remove zeros from possible_peaks_ind
   possible_peaks_ind = possible_peaks_ind[WHERE(possible_peaks_ind NE 0)]
   
   ;; use shift from initial gaussfit if no peak_ind within the lag
   ;; range specified by start_ind and end_ind and it lies within the range
   IF possible_peaks_ind[0] EQ 0 THEN BEGIN
      IF lags[init_shift_index] LT 80 && lags[init_shift_index] GT -80 THEN BEGIN
         cc_shifts[i] = lags[init_shift_index]
         shift_index = init_shift_index
         shift_inds[i] = shift_index
      ENDIF ELSE BEGIN
         cc_shifts[i] = 0
         shift_index = WHERE(lags EQ 0)
         shift_inds[i] = shift_index
      ENDELSE
      GOTO, leap ;;skip past the rest of this code that requires at least one possible_peaks_ind
   ENDIF

   ;; define peak_shift_ind using the max ccs value of the possible_peaks_ind
   possible_peaks_ccs = DBLARR(N_ELEMENTS(possible_peaks_ind)) ;; array to store ccs values for peaks
   FOR q=0, N_ELEMENTS(possible_peaks_ind)-1 DO BEGIN
      possible_peaks_ccs[q] = ccs[possible_peaks_ind[q],i]
   ENDFOR
   peak_shift_ind = possible_peaks_ind[WHERE(possible_peaks_ccs EQ MAX(possible_peaks_ccs))]
   peak_shift_ind = peak_shift_ind[0]


   ;; determine the trough_ind just to the left of peak_shift_ind
   IF MIN(trough_ind[WHERE(trough_ind[*,i] NE 0),i]) LT peak_shift_ind THEN BEGIN
      ;;use MAX of trough_ind that are less than peak_shift_ind
      trough_left = MAX(trough_ind[WHERE(trough_ind[WHERE(trough_ind[*,i] NE 0),i] LT peak_shift_ind),i])
   ENDIF ELSE BEGIN
      ;; use the start_ind for the range used to find the peaks and troughs
      trough_left = start_ind
   ENDELSE
   
   ;; determine the trough_ind just to the right of peak_shift_ind
   ;;use MIN of trough_ind that are greater than peak_shift_ind
   trough_ind_adj = trough_ind[WHERE(trough_ind[*,i] NE 0),i]
   trough_right = MIN(trough_ind[WHERE(trough_ind_adj GT peak_shift_ind),i])
   
   ;;define the range of indices to be fit (indices of the num_pts
   ;;maximum points between the two troughs)
   range = N_EXTREMA(ccs[[trough_left : trough_right],i], num_pts, 1)+trough_left
   
   ;;apply gaussfit to the top of the peak
   gauss_fits[*,i] = GAUSSFIT(lags[range], ccs[range,i], $
                              gauss_para, NTERMS=3)

   ;;store index of the shift given by the gaussfit if it lies within
   ;;the range of points the fit was applied to
   IF gauss_para[1] LE lags[range[num_pts-1]] && gauss_para[1] GE lags[range[0]] THEN BEGIN
      shift_index = VALUE_LOCATOR(lags, gauss_para[1]) 
      shift_inds[i] = shift_index
   ENDIF ELSE BEGIN
      shift_index = peak_shift_ind
      shift_inds[i] = shift_index
   ENDELSE

   ;;store final shift
   cc_shifts[i] = lags[shift_index] 

   leap: ;; jump here if no possible_peaks_ind

   ;;some test plotting
   ;PLOT, lags, ccs[*,i], PSYM=-3, YRANGE=[-.5,1.2]
   ;OPLOT, lags[range], gauss_fits[*,i], PSYM=-3, COLOR=200 
   ;OPLOT, [lags[0],lags[i_num-1]], [mean_ccs, mean_ccs], PSYM=-3, COLOR=200
   ;OPLOT, lags[range], ccs[range,i], COLOR=250, PSYM=-3
   ;OPLOT, [lags[shift_inds[i]],lags[shift_inds[i]]], [ccs[shift_inds[i],i],ccs[shift_inds[i],i]], $
   ;       PSYM=4, COLOR=100
   ;OPLOT, [lags[trough_left],lags[trough_left]], [ccs[trough_left,i],ccs[trough_left,i]], $
   ;       PSYM=1, COLOR=100
   ;OPLOT, [lags[trough_right],lags[trough_right]], [ccs[trough_right,i],ccs[trough_right,i]], $
   ;       PSYM=1, COLOR=100
   ;OPLOT, lags[neg_ind[*,i]], ccs[neg_ind[*,i],i], PSYM=1, COLOR=150
   ;OPLOT, lags[trough_ind[*,i]], ccs[trough_ind[*,i],i], PSYM=1, COLOR=220 
   ;OPLOT, lags[peak_ind[*,i]], ccs[peak_ind[*,i],i], PSYM=1, COLOR=150
   ;OPLOT, lags[start_ind :end_ind], init_gauss_fits[*,i], PSYM=-3, COLOR=150
   ;OPLOT, [lags[init_shift_index], lags[init_shift_index]], $
   ;       [ccs[init_shift_index,i], ccs[init_shift_index,i]], PSYM=4, COLOR=200
   ;stop
ENDFOR


;; plot the ccs vs lags with gauss fits one by one if desired
;; FOR i=0, num_analyzed-1 DO BEGIN
;;    PLOT,  lags, ccs[*,i], PSYM=-3
;;    OPLOT, lags, gauss_fits[*,i], PSYM=-3, COLOR=150
;;    OPLOT, lags, init_gauss_fits[*,i], PSYM=-3, COLOR=250
;;    stop
;; ENDFOR         

;;===================================================================================;;
;;============================= END DETERMINE SHIFTS ================================;;
;;===================================================================================;;






;;===================================================================================;;
;;================================ APPLY SHIFTS =====================================;;
;;===================================================================================;;

;; apply shifts
FOR i=0, num_analyzed-1 DO BEGIN
      sl_waves[*,i] = l_waves[*,i] + cc_shifts[i]*spacing
ENDFOR

;;anti-log the shifted wavelengths
cc_waves = 10D^(sl_waves)

;;===================================================================================;;
;;============================== END APPLY SHIFTS ===================================;;
;;===================================================================================;;





;;********************************* VISUAL CHECK PLOTTING *****************************;;
;; plot the unshifted spectra
;PLOT, waves[*,0], fluxs[*,0], XRANGE=[653.0,659.0], YRANGE=[0,15], PSYM=-3
;FOR i=1, num_analyzed-1 DO BEGIN
;   OPLOT, waves[*,i], fluxs[*,i], PSYM=-3
;ENDFOR

;stop  

;; plot the shifted spectra
;PLOT, cc_waves[*,0], fluxs[*,0], XRANGE=[653.0,659.0], YRANGE=[0,15], PSYM=-3
;FOR i=1, num_analyzed-1 DO BEGIN
;   OPLOT, cc_waves[*,i], fluxs[*,i], PSYM=-3
;ENDFOR

;stop

;; ;;quick plot to preview the post-shift result
;; ;;plot all spectra shifted
;; PLOT, [654, 659],[0, MAX(adj_fluxs)], XTITLE='Wavelength (nm)', YTITLE='F/Fc', XSTYLE=1, $
;;       CHARSIZE=1.3, TITLE='All Spectra', /NODATA

;; ;; overplot a dashed line at the continuum level                                                        
;; OPLOT, [653.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=60

;; ;; plot the spectra
;; FOR i=0, num_analyzed-1 DO BEGIN
;;    OPLOT, cc_waves[*,i], adj_fluxs[*,i], PSYM=-3
;; ENDFOR
;; stop
;;************************************************************************************;;





;;===================================================================================;;
;;========================= GENERATE DIAGNOSTIC REPORTS =============================;;
;;===================================================================================;;

;;define lambda symbol
lambda = '!7k!X'

IF ~diagnostics THEN GOTO, skip

DEVICE, DECOMPOSED = 0 ;; Colors are taken from a color table (not specified by R, G, B)
LOADCT, 39, /SILENT

SPAWN, 'pwd', og_loc
og_loc = og_loc+'/' ;; original directory location

;;Define some color tables
color_table = [60,250,150,220,30,90,110,190,10]
color_text = ['blue','red','green','orange','violet','lt_blue','cyan','yellow','black']


FOR i=0, num_analyzed-1 DO BEGIN
   ;; set device parameters
   SET_PLOT, 'ps'
   DEVICE, /COLOR, BITS_PER_PIXEL=8
   DEVICE, FILENAME='CC_plot_'+analysis_frames[i]+'.ps', ENCAPSULATED=0
   DEVICE, XSIZE=8, YSIZE=11, XOFFSET=0.5, YOFFSET=0, /INCHES

  ;; plot cross correlation plot
   PLOT, lags, ccs[*,i], XTITLE='Lags', YTITLE='Cross Correlation', XRANGE = [MIN(lags), MAX(lags)], $
         YRANGE = [-0.5, 1.5], XSTYLE=1, CHARSIZE=1.2, CHARTHICK=3, XTHICK=3, YTHICK=3, THICK=3, $
         TITLE='Cross Correlation Between '+analysis_frames[0]+' and '+analysis_frames[i], $
         POS=[.13, 0.50, 0.95, 0.90], /DATA, /NOERASE

   ;; overplot the gaussian fit used to narrow the search for the fit
   OPLOT, lags[start_ind :end_ind], init_gauss_fits[*,i], PSYM=-3, COLOR=color_table[1], THICK=3
   
   ;;overplot the index used for the shift
   OPLOT, [lags[shift_inds[i]],lags[shift_inds[i]]], [ccs[shift_inds[i],i],ccs[shift_inds[i],i]], $
          PSYM=4, COLOR=100

   ;; add useful peak data to the report
   XYOUTS, 0.03, 0.01, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
   XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
   XYOUTS, 0.20, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
   XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
   XYOUTS, 0.5, 0.9625, $
           '_________________________________________________________________________________________________',$
           ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL


   XYOUTS, 0.15, 0.4, 'Cross Correlation:', CHARSIZE=1.2, CHARTHICK=1.5, /NORMAL
   XYOUTS, 0.15, 0.38, 'Shift = '+STRING(cc_shifts[i], FORMAT='(F6.1)')+' steps', CHARSIZE=1.2, $
           CHARTHICK=1.5, /NORMAL
   XYOUTS, 0.15, 0.36, 'Shift [in log('+lambda+')] = '+STRING((cc_shifts[i]*spacing)/10D^(-5), $
                                         FORMAT='(F8.4)')+'e-5', CHARSIZE=1.2, CHARTHICK=1.5, /NORMAL

   ;; report differences in halpha_centralizer shifts (shifts) and cc_shifts
   XYOUTS, 0.15, 0.3, 'Halpha Centralizer:', CHARSIZE=1.2, CHARTHICK=1.5, /NORMAL
   XYOUTS, 0.15, 0.28, 'Shift [in log('+lambda+')] = '+ STRING((shifts[i]/10D^(-5)), FORMAT='(F8.4)')+'e-5', $
           CHARSIZE=1.2, CHARTHICK=1.5, /NORMAL

   DEVICE, /CLOSE_FILE
   SET_PLOT, 'X'
   
   ;;move the file to the diagonostics directory for that star
   SPAWN, 'mv -f '+'CC_plot_'+analysis_frames[i]+'.ps ' + og_loc + targets_dir_name +'/HD' + $
          star_hd + '/Diagnostics/'

ENDFOR



;; set device parameters                                                                                
SET_PLOT, 'ps'
DEVICE, /COLOR, BITS_PER_PIXEL=8
DEVICE, FILENAME='Shifted_HD'+star_hd+'.ps', ENCAPSULATED=0,$
        FONT_SIZE=12, XSIZE=8, YSIZE=11, XOFFSET=0.5, YOFFSET=0, /INCHES

;; add useful peak data to the report                                                                               
XYOUTS, 0.03, 0.01, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.20, 0.97, star_name, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.5, 0.9625, $
        '_________________________________________________________________________________________________',$
        ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL

;;plot all spectra unshifted
PLOT, waves[*,0], fluxs[*,0], XRANGE=[653.0,659.0], YRANGE=[0,MAX(fluxs[[50:500],*])], $
      XTITLE=lambda+' (nm)', YTITLE='F/Fc', XSTYLE=1, CHARSIZE=0.8, CHARTHICK=3, XTHICK=3, $
      YTHICK=3, THICK=3, TITLE='Unshifted Spectra for HD'+star_hd, POS=[.1, 0.69, 0.95, 0.94], $
      /DATA, /NOERASE
FOR i=1, num_analyzed-1 DO BEGIN
   OPLOT, waves[*,i], fluxs[*,i], PSYM=-3
ENDFOR

;; plot the CC shifted spectra 
PLOT, cc_waves[*,0], fluxs[*,0], XRANGE=[653.0,659.0], YRANGE=[0,MAX(fluxs[[50:500],*])], $
      XTITLE=lambda+' (nm)', YTITLE='F/Fc', XSTYLE=1, CHARSIZE=0.8, CHARTHICK=3, XTHICK=3, $
      YTHICK=3, THICK=3, TITLE='CC Shifted Spectra', POS=[.1, 0.37, 0.95, 0.63], $
      /DATA, /NOERASE
FOR i=1, num_analyzed-1 DO BEGIN
   OPLOT, cc_waves[*,i], fluxs[*,i], PSYM=-3
ENDFOR

;; plot the H-alpha shifted spectra
PLOT, ha_waves[*,0], fluxs[*,0], XRANGE=[653.0,659.0], YRANGE=[0,MAX(fluxs[[50:500],*])], $
      XTITLE=lambda+' (nm)', YTITLE='F/Fc', XSTYLE=1, CHARSIZE=0.8, CHARTHICK=3, XTHICK=3, $
      YTHICK=3, THICK=3, TITLE='Halpha_centralizer Shifted Spectra', POS=[.1, 0.05, 0.95, 0.31], $
      /DATA, /NOERASE
FOR i=1, num_analyzed-1 DO BEGIN
   OPLOT, ha_waves[*,i], fluxs[*,i], PSYM=-3
ENDFOR

DEVICE, /CLOSE_FILE
SET_PLOT, 'X'

;;move the file to the diagonostics directory for that star                                             
SPAWN, 'mv -f '+'Shifted_HD'+star_hd+'.ps ' + og_loc + targets_dir_name +'/HD' + $
       star_hd + '/Diagnostics/'

;;===================================================================================;;
;;======================= END GENERATE DIAGNOSTIC REPORTS ===========================;;
;;===================================================================================;;




skip: ;; program jumps here if no diagnostics are desired 

END

