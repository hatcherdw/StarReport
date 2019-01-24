;; Star Report Pipeline. Central Michigan University
;;
;; stability_test.pro
;;
;; This procedure is meant to search for trends in the EW
;; vs. time plot for each star. The various outcomes for trends are
;; disk stable and variable and those exhibiting either disk growth or
;; disk loss. 
;; 
;; Author: W. Glennon Fagan

PRO stability_test, EWs, EWs_error, sorted_jdates
RESOLVE_ALL, /QUIET
COMMON star_block


ps_plot = 0 ;set this to 1 if you want a post scipt plot generated
ps_pure_ew_plot = 0; OR set this to 1 if you want a post scipt plot generated for jsut an ew plot 


IF N_ELEMENTS(EWs) LE 5 THEN BEGIN
   trend = 'Too few spectra to assess'
   GOTO, ope_sneakrightpastya
ENDIF


;************USER DEFINED THRESHOLDS************;
user_max_trend_value = 0.66;this value approaches 1 as a trend becomes LESS clear
user_max_stable_stddev = 0.272264;this value approaches 1 as it becomes more stable;;;
user_max_stable_slope_value = 6.8e-5
;***********************************************;


;start by converting the jdates to double
jdates = DOUBLE(sorted_jdates)

jdates = jdates - 2450000.0

;import a color table to be used in plotting
color_table = [60,250,150,220,30,90,110,190,10]
color_text = ['blue','red','green','orange','violet','lt_blue','cyan','yellow','black']


stddev_ew = STDDEV(EWs) ;calculate the standard deviation of all the EWs
mean_ew = MEAN(EWs) ;calculate the average EW
;pos_dev = mean_ew + stddev_ew ;add the STDDEV to the average EW for plotting a line
;neg_dev = mean_ew - stddev_ew ;subtract the STDDEV from the average EW for plotting a line


linfit_para = LINFIT(jdates, EWs) ;fit a line the distibution of EWs vs jdates
slope = linfit_para[1] ;define variable for slope from LINFIT parameters
y_int = linfit_para[0] ;define variable for y intercept from LINFIT parameters
theta_rot = ATAN(slope) ;calculate the angle between the x-axis and the fitted line


;create x and y values for the fitted line
linfit_x = (DINDGEN(50)/(50-1)*(jdates[N_ELEMENTS(jdates)-1]-jdates[0]))+jdates[0]
linfit_y = linfit_x*linfit_para[1] + linfit_para[0]




;************************* Distance from trendline to point manually calc  *************************;

ews_rotated = DBLARR(N_ELEMENTS(ews))
FOR i=0, N_ELEMENTS(ews)-1 DO BEGIN
   h = ((slope*jdates[i]) + y_int) - ews[i]
   ews_rotated[i] = h*COS(theta_rot)
ENDFOR

stddev_ew_rotated = STDDEV(ews_rotated)

;the closer this value is to 1 the less of a clear increasing or
;decreasing trend is evident in the EW vs. Time plot
trend_value = stddev_ew_rotated/stddev_ew

;the closer this value is to 1 the less variability is evident in the
;EW vs Time plot
;stability_value = 1/ (1+stddev_ew)




;****************************** FORK 1: LINEAR TREND OR NOT ******************************;

Ew_trend = 0;will be set to an integer from 1-4 in the following if statement, depending on its trend

IF trend_value GT user_max_trend_value THEN BEGIN ;There IS NO clear trend. Check for degree of stability
   IF stddev_ew LT user_max_stable_stddev THEN BEGIN ;THEN IS STABLE
      ew_trend = 1
      trend = 'Stable'
   ENDIF ELSE BEGIN
      ew_trend = 2
      trend = 'Variable'
   ENDELSE 
ENDIF ELSE BEGIN ;There IS a clear trend. Look for + or - trend 
   IF ABS(slope) GT user_max_stable_slope_value THEN BEGIN;Slope is large enough to be considered unstable
      IF linfit_para[1] GT 0 THEN BEGIN
         ew_trend = 3
         trend = 'Disk-loss'
      ENDIF ELSE BEGIN
         ew_trend = 4
         trend = 'Disk-growth'
      ENDELSE
   ENDIF ELSE BEGIN             ;the slope is not significant enough over the timescale to be 
      ew_trend = 1
      trend = 'Stable'
   ENDELSE
ENDELSE




;************************************* PLOTTING *************************************;

;;****************Calculations for plotting purposes****************
pos_dev = mean_ew + stddev_ew ;add the STDDEV to the average EW for plotting a line
neg_dev = mean_ew - stddev_ew ;subtract the STDDEV from the average EW for plotting a line

;fork 1 calculations
d = stddev_ew_rotated/(cos(theta_rot)) ;displacement in y direc for the lin fit's stddev
beg_upper = linfit_y[0] + d
fin_upper = linfit_y[N_ELEMENTS(linfit_y)-1] + d
beg_lower = linfit_y[0] - d
fin_lower = linfit_y[N_ELEMENTS(linfit_y)-1] - d

;fork 2 calculations
pos_thres_dev = mean_ew + user_max_stable_stddev
neg_thres_dev = mean_ew - user_max_stable_stddev



; Make a vector of 16 points, A[i] = 2pi/16:
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points,
; and set the filled flag:
USERSYM, COS(A), SIN(A), /FILL
;;******************************************************************

;WINDOW, 2
;;****************GENERATE PLOT****************
;PLOT, jdates, EWs, XTITLE='JD (-2540000)', YTITLE='EW (nm)', XSTYLE=1, $
;      CHARSIZE=1.1, CHARTHICK=1.5, YRANGE=[MAX(ews)+0.3, MIN(ews)-0.3], $
;      XTHICK=1.5, YTHICK=1.5, THICK=1.5,SYMSIZE=1.2, PSYM=8, TITLE="'"+star_proper_name+"'  "+star_name
;OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [mean_ew, mean_ew], COLOR=30
;OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [pos_dev, pos_dev], LINESTYLE=5, COLOR=30
;OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [neg_dev, neg_dev], LINESTYLE=5, COLOR=30
;;;*********************************************
;OPLOT, linfit_x, linfit_y, PSYM=-3, COLOR=250
;OPLOT, [linfit_x[0], linfit_x[N_ELEMENTS(linfit_x)-1]], [beg_upper, fin_upper], LINESTYLE=4, COLOR=250
;OPLOT, [linfit_x[0], linfit_x[N_ELEMENTS(linfit_x)-1]], [beg_lower, fin_lower], LINESTYLE=4, COLOR=250
;;***************************************************




;************************************* POST SCRIPT PLOTTING *************************************;

IF ps_plot THEN BEGIN

   alpha = '!7a!X'
   CD, 'poster_plots'
   linethic = 5


   ;;***************************FORK 1 PLOT***************************
   !P.FONT= -1            ;;Set to -1 for Hershey fonts
   SET_PLOT, 'ps'         ;select post script 
   DEVICE, DECOMPOSED = 0 ;Colors are taken from a color table (not specified by R, G, B)
   LOADCT, 39, /SILENT    ;load color table
   ;;***************************************************
   DEVICE, FILENAME='HD'+star_hd+'_stability.ps', ENCAPSULATED=0 , /LANDSCAPE;,$
   ;XSIZE=10, YSIZE=6, XOFFSET=0, YOFFSET=0, /INCHES, /LANDSCAPE

   PLOT, jdates, EWs, XTITLE='JD (-2540000)', YTITLE='H'+alpha+' EW (nm)', XSTYLE=1, $
         CHARSIZE=1.5, CHARTHICK=4.5, XRANGE=[jdates[0]-100, jdates[N_ELEMENTS(jdates)-1]+100],$
         YRANGE=[MAX(ews)+0.3, MIN(ews)-0.3], $
         XTHICK=4.5, YTHICK=4.5, THICK=1.5,SYMSIZE=1, PSYM=8, TITLE="'"+star_proper_name+"'  "+star_name
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [mean_ew, mean_ew], COLOR=30, THICK=linethic
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [pos_dev, pos_dev], $
          LINESTYLE=5, COLOR=30, THICK=linethic
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [neg_dev, neg_dev], $
          LINESTYLE=5, COLOR=30, THICK=linethic
   ;;***************************************************
   OPLOT, linfit_x, linfit_y, PSYM=-3, COLOR=250, THICK=linethic
   OPLOT, [linfit_x[0], linfit_x[N_ELEMENTS(linfit_x)-1]], [beg_upper, fin_upper], $
          LINESTYLE=4, COLOR=250, THICK=linethic
   OPLOT, [linfit_x[0], linfit_x[N_ELEMENTS(linfit_x)-1]], [beg_lower, fin_lower], $
          LINESTYLE=4, COLOR=250, THICK=linethic
   ;;*************************************************** 
   DEVICE, /CLOSE_FILE
   SET_PLOT, 'X'


   ;;***************************FORK 2 PLOT***************************
   !P.FONT= -1            ;;Set to -1 for Hershey fonts
   SET_PLOT, 'ps'         ;select post script
   DEVICE, DECOMPOSED = 0 ;Colors are taken from a color table (not specified by R, G, B)
   LOADCT, 39, /SILENT    ;load color table
   DEVICE, FILENAME='HD'+star_hd+'_fork2.ps', ENCAPSULATED=0 , /LANDSCAPE;,$
   ;XSIZE=10, YSIZE=6, XOFFSET=0, YOFFSET=0, /INCHES,
   ;/LANDSCAPE          

   PLOT, jdates, EWs, XTITLE='JD (-2540000)', YTITLE='H'+alpha+' EW (nm)', XSTYLE=1, $
         CHARSIZE=1.5, CHARTHICK=4.5, XRANGE=[jdates[0]-100, jdates[N_ELEMENTS(jdates)-1]+100],$
         YRANGE=[MAX(ews)+0.3, MIN(ews)-0.3], $
         XTHICK=4.5, YTHICK=4.5, THICK=1.5,SYMSIZE=1, PSYM=8, TITLE="'"+star_proper_name+"'  "+star_name
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [mean_ew, mean_ew], THICK=linethic
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [pos_dev, pos_dev], $
          LINESTYLE=5, THICK=linethic
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [neg_dev, neg_dev], $
          LINESTYLE=5, THICK=linethic
   ;;now add the sigma that is defined by the user in the form of the user_min_stable_value
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [pos_thres_dev, pos_thres_dev], $
          LINESTYLE=5, COLOR=250, THICK=linethic
   OPLOT, [jdates[0], jdates[N_ELEMENTS(jdates)-1]], [neg_thres_dev, neg_thres_dev], $
          LINESTYLE=5, COLOR=250, THICK=linethic

   
   DEVICE, /CLOSE_FILE
   SET_PLOT, 'X'


   CD, '..'
ENDIF ELSE IF ps_pure_ew_plot THEN BEGIN
   alpha = '!7a!X'
   CD, 'poster_plots'
   linethic = 5

   
   ;;***************************FORK 1 PLOT***************************
   !P.FONT= -1            ;;Set to -1 for Hershey fonts
   SET_PLOT, 'ps'         ;select post script
   DEVICE, DECOMPOSED = 0 ;Colors are taken from a color table (not specified by R, G, B)
   LOADCT, 39, /SILENT   ;load color table                                                                   
   ;;*****************************************************************
   DEVICE, FILENAME='HD'+star_hd+'_jdate_ew.ps', ENCAPSULATED=0 , /LANDSCAPE;,$
   ;XSIZE=10, YSIZE=6, XOFFSET=0, YOFFSET=0, /INCHES, /LANDSCAPE

   PLOT, jdates, EWs, XTITLE='JD (-2540000)', YTITLE='H'+alpha+' EW (nm)', XSTYLE=1, $
         CHARSIZE=1.5, CHARTHICK=4.5, XRANGE=[jdates[0]-100, jdates[N_ELEMENTS(jdates)-1]+100],$
         YRANGE=[MAX(ews)+0.3, MIN(ews)-0.3], $
         XTHICK=4.5, YTHICK=4.5, THICK=1.5,SYMSIZE=1, PSYM=8, TITLE="'"+star_proper_name+"'  "+star_name

   DEVICE, /CLOSE_FILE
   SET_PLOT, 'X'

   CD, '..'
ENDIF 


ope_sneakrightpastya:;sneak on past the curve fitting if the star has been stable

END
