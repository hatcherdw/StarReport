;;======================================================================;;
;; Star Report Pipeline, Central Michigan University
;;
;; plotter.pro
;;
;; This procedure was designed to conduct all plotting resulting from the
;; spectroscopic analysis. It can be divided into two main parts. (1)
;; Generate an on-screen plot of the spectrum when
;; spectrum_analyzer.pro is ran alone. (2) Generate Halpha and
;; Continnum postscript plots when star_report.pro is being
;; executed. Two directories are created during this procedure in part
;; (2). These are titled 'Halpha_plots' and 'Continuum_plots' and
;; can be found within the star's directory. 
;;
;; Authors: Christian Hannah & Glennon Fagan
;;
;;======================================================================;;

PRO plotter,temp_file, spec_plot_alone, color_table, ew_fixed, fixed_ew_error, fixed_wave_beg, $
            fixed_wave_end, intsig_wave_int, intsig_wave_fin, int_sig, wave_end, wave_start
RESOLVE_ALL, /QUIET
COMMON spec_block
COMMON star_block
COMMON peak_block

;; define thickness variables 
plothick = 3
lineThick = 2



;;===============================================================================;;
;;=========================== SPECTRUM_ANALYZER ALONE ===========================;;
;;===============================================================================;;

IF spec_plot_alone THEN BEGIN

   ;;**************** Plot Spectrum ******************;;
   
   dummyTitle='!3'
   WINDOW, 0, XSIZE=100, YSIZE=100
   PLOT, [0, 1], [0, 1], TITLE=dummyTitle
   WDELETE, 0

   PLOT, wave, flux, XTITLE='Wavelength (nm)', YTITLE='F/Fc', XRANGE = [653, 659], XSTYLE=1, $
         CHARSIZE=1.1, CHARTHICK=plotThick, XTHICK=plotThick, YTHICK=plotThick, THICK=lineThick, $
         TITLE='Frame = '+temp_file
   
   ;; overplot a dashed line at the continuum level
   OPLOT, [653.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=70, THICK=1.5



   ;;*********** Excluded Region Plotting ************;;
   ;;
   ;;Loop used to plot all regions removed
   ;;This method of plotting must be used in order to plot a continuous
   ;;line of the points excluded and no have each of those regions be
   ;;connected by a line in the graph. Don't know what I mean?
   ;;un-comment the plot below and stop the program right after it, then
   ;;run spectrum_plot alone.
   ;PLOT, wave, flux, XRANGE = [650, 659], YRANGE = [0, 2], THICK=2.5
   ;OPLOT, wave[WHERE(mask_flux NE 0)], flux[WHERE(mask_flux NE 0)], $
   ;PSYM = -3,THICK=2.5, COLOR=color_table[7]   

   temp_flux = DBLARR(N_ELEMENTS(flux))
   new_start = 0
   counter = 0
   j=0
   WHILE(j LT N_ELEMENTS(flux)) DO BEGIN
      FOR i=j,N_ELEMENTS(flux) - 2 DO BEGIN
         
         IF mask_flux[i] EQ 0 THEN BEGIN
            temp_flux[i] = flux[i]
            new_start++
            IF (mask_flux[i] EQ 0 && mask_flux[i+1] NE 0) THEN BREAK
         ENDIF
      ENDFOR
      
      IF N_ELEMENTS(WHERE(temp_flux eq 0D )) NE N_ELEMENTS(flux) THEN BEGIN
         IF spec_plot_alone THEN OPLOT, wave[WHERE(temp_flux NE 0)], flux[WHERE(temp_flux NE 0)], $
                                        PSYM = -3,THICK=linethick, COLOR=250
      ENDIF ELSE BEGIN
         BREAK
      ENDELSE
      
      temp_flux = DBLARR(N_ELEMENTS(flux))
      j = new_start+1
      counter ++
   ENDWHILE



   
   ;;**************** PEAK INFO PLOTTING ******************;;
   ;; uncomment this code to show peak fits in on-screen plot when
   ;; spectrum_analyzer.pro is executed alone. 

   ;; plot the peak fit for single peak spectra
   ;IF ISA(interval) && ISA(peak_fit) && peak_num EQ 1 THEN $
   ;   OPLOT, wave[interval], peak_fit, PSYM=-3, COLOR=color_table[2], THICK=linethick + 1
   
   ;; plot the peak fits with multiple peaks
   ;IF ISA(interval_1) && ISA(interval_2) && ISA(peak_fit_1) && ISA(peak_fit_2) && peak_num EQ 2 THEN BEGIN
   ;   OPLOT, wave[interval_1], peak_fit_1, PSYM=-3, COLOR=color_table[2], THICK=linethick + 1
   ;   OPLOT, wave[interval_2], peak_fit_2, PSYM=-3, COLOR=color_table[2], THICK=linethick + 1
   ;ENDIF
   GOTO, jump
ENDIF






;;===============================================================================;;
;;============================== STAR REPORT ====================================;;
;;===============================================================================;;

;;********** Define general variables ************;;

obs_date = analysis_dates[WHERE(analysis_frames EQ temp_file)]

;;define some symbols
plus_minus = '!M+!X'
del = '!7D!X'
lam = '!7k!X'

;; set the lines to be thick for ps output
plothick = 3
lineThick = 4

;;calculate delta lambda for the dynamic range
del_lambda = STRING(wave[wave_end] - wave[wave_start], FORMAT='(F6.3)')

;;calculate delta lambda for the fixed range
fixed_del_lambda = STRING(wave[fixed_wave_end] - wave[fixed_wave_beg], FORMAT='(F6.3)')

;;change them from indices to corresponding wavelengths
fixed_wave_beg = wave[fixed_wave_beg]
fixed_wave_end = wave[fixed_wave_end]



;;****** Create directories for Halpha and Continuum ****;;

;;CD into the stars directory
CD, targets_dir_name
CD, 'HD'+star_hd

;;Create directories for the two types of plots
IF ~FILE_TEST('Halpha_plots', /DIRECTORY) THEN SPAWN, 'mkdir Halpha_plots'
IF ~FILE_TEST('Continuum_plots', /DIRECTORY) THEN SPAWN, 'mkdir Continuum_plots'




;;************** Continuum Plot *****************;;

!P.FONT= -1 ;;Set to -1 for Hershey fonts
SET_PLOT, 'ps' ;select post script
DEVICE, DECOMPOSED = 0 ;Colors are taken from a color table (not specified by R, G, B)
LOADCT, 39, /SILENT ;load color table

DEVICE, FILENAME='Continuum'+temp_file+'_plot.ps', ENCAPSULATED=0,$
        XSIZE=8.5, YSIZE=10, XOFFSET=0, YOFFSET=0.5, /INCHES

PLOT, wave, flux, XTITLE='Wavelength (nm)', YTITLE='F/Fc', XRANGE = [649, 659], $
      YRANGE = [0.85, 1.15], XSTYLE=1, CHARSIZE=1.1, CHARTHICK=plotThick,$
      XTHICK=plotThick, YTHICK=plotthick, THICK=linethick, $
      TITLE='Continuum, Frame = '+temp_file, POS=[.10, 0.4, 0.90, 0.9]

;;while loop to only plot excluded regions
temp_flux = DBLARR(N_ELEMENTS(flux))
new_start = 0
counter = 0
j=0
WHILE(j LT N_ELEMENTS(flux)) DO BEGIN
   FOR i=j,N_ELEMENTS(flux) - 2 DO BEGIN
      IF mask_flux[i] EQ 0 THEN BEGIN
         temp_flux[i] = flux[i]
         new_start++
         IF (mask_flux[i] EQ 0 && mask_flux[i+1] NE 0) THEN BREAK
      ENDIF
   ENDFOR
   IF N_ELEMENTS(WHERE(temp_flux eq 0D )) NE N_ELEMENTS(flux) THEN BEGIN
      IF ~spec_plot_alone THEN OPLOT, wave[WHERE(temp_flux NE 0)], flux[WHERE(temp_flux NE 0)], $
                                      PSYM = -3,THICK=linethick, COLOR=color_table[3]
   ENDIF ELSE BEGIN
      BREAK
   ENDELSE
   temp_flux = DBLARR(N_ELEMENTS(flux))
   j = new_start+1
   counter ++
ENDWHILE

;; overplot a dashed line at the continuum level                                                 
OPLOT, [649.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=color_table[1]

;;plot vertical lines and wavelengths corresponding to initial sigma 
OPLOT,  [wave[intsig_wave_int],wave[intsig_wave_int]], [0, 10]
OPLOT,  [wave[intsig_wave_fin],wave[intsig_wave_fin]], [0, 10]
delta_lambda = STRING(wave[intsig_wave_fin] - wave[intsig_wave_int], FORMAT='(F5.3)')

;;print some outputs                                                                             
XYOUTS, 0.03, 0.01, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.20, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.5, 0.9625, $
        '_________________________________________________________________________________________________',$
        ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.3, 'Obs. Date: '+STRCOMPRESS(obs_date, /REMOVE_ALL),$
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.25, 'EW (dynamic): '+ STRING(ew, FORMAT='(F8.5)')+ ' '+ plus_minus +' ' + $
        STRING(ew_error, FORMAT='(F7.5)')+ ' nm ['+STRING(wave[wave_start], FORMAT='(F7.3)')$
        +': '+ STRING(wave[wave_end],FORMAT='(F7.3)')+', '+del+lam+' ='+del_lambda+']', $
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.2, 'Cont. SNR (gauss): '+ STRING(snr, FORMAT='(F6.2)')+ ' '+ plus_minus+' ' + $
        STRING(snr_error, FORMAT='(F4.2)'),$
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.15, 'Initial Sigma = '+STRING(int_sig, FORMAT='(F6.4)')+ ' ' + $
        '['+STRING(wave[intsig_wave_int],FORMAT='(F7.3)')+ ': '+$
        STRING(wave[intsig_wave_fin],FORMAT='(F7.3)')+', '+del+lam+' = '+ delta_lambda+']', $
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
;XYOUTS, 0.11, 0.10, 'Symmetry:  '+STRING(ew, FORMAT='(F5.2)'), $
;        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL


DEVICE, /CLOSE_FILE
SET_PLOT, 'X'

;;move the newly created postscript
SPAWN, 'mv -f Cont*.ps ' +root_direc+ targets_dir_name+'/HD' + star_hd + '/Continuum_plots'





;;******************* Halpha Plot ********************;;

!P.FONT= -1 ;;Set to -1 for Hershey fonts
SET_PLOT, 'ps' ;select post script
DEVICE, DECOMPOSED = 0 ;Colors are taken from a color table (not specified by R, G, B)
LOADCT, 39, /SILENT ;load color table

DEVICE, FILENAME='Halpha'+temp_file+'_plot.ps',$
        ENCAPSULATED=0, XSIZE=8.5, YSIZE=10, XOFFSET=0, YOFFSET=0.5, /INCHES

PLOT, wave, flux, XTITLE='Wavelength (nm)', YTITLE='F/Fc', XRANGE = [653, 659], XSTYLE=1, $
      CHARSIZE=1.1, CHARTHICK=plotThick, XTHICK=plotThick, YTHICK=plotThick, THICK=lineThick, $
      TITLE='Halpha Emission, Frame = '+temp_file, POS=[0.10, 0.4, 0.9, 0.9]

;; overplot a dashed line at the continuum level                                                 
OPLOT, [649.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=color_table[1]

;; Overplot lines to represent the dynamic starting and stopping waves
;; for EW calc
OPLOT, [wave[wave_start], wave[wave_start]], [0,30]
OPLOT, [wave[wave_end], wave[wave_end]], [0,30]

;; Overplot lines to represent the dynamic starting and stopping waves
;; for EW calc
OPLOT, [fixed_wave_beg, fixed_wave_beg], [0,30], LINESTYLE=2
OPLOT, [fixed_wave_end, fixed_wave_end], [0,30], LINESTYLE=2

OPLOT, [656.28, 656.28], [0,30], LINESTYLE=1

;;print some outputs
XYOUTS, 0.03, 0.01, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.20, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.5, 0.9625, $
        '_________________________________________________________________________________________________',$
        ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.3, 'Obs. Date: '+STRCOMPRESS(obs_date, /REMOVE_ALL),$
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.25, 'EW (dynamic): '+ STRING(ew, FORMAT='(F8.5)')+ ' '+ plus_minus +' ' + $
        STRING(ew_error, FORMAT='(F7.5)')+ ' nm ['+STRING(wave[wave_start], FORMAT='(F7.3)')$
        +': '+ STRING(wave[wave_end],FORMAT='(F7.3)')+', '+del+lam+' ='+del_lambda+']', $
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.235, '(SOLID)', ALIGNMENT=0, CHARSIZE=0.75, CHARTHICK=1.5, /NORMAL
XYOUTS, 0.11, 0.20, 'EW (fixed):   ' + STRING(ew_fixed, FORMAT='(F9.5)') +' '+ plus_minus +' '+$
        STRING(fixed_ew_error, FORMAT='(F7.5)') + ' nm ' + $
        '['+STRING(fixed_wave_beg, FORMAT='(F7.3)') + ': '+ STRING(fixed_wave_end ,FORMAT='(F7.3)')+$
        ', '+del+lam+' ='+fixed_del_lambda+']',$
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS,0.11, 0.185, '(LONG DASHED)', ALIGNMENT=0, CHARSIZE=0.75, CHARTHICK=1.5, /NORMAL
XYOUTS, 0.11, 0.15, 'Cont. SNR (gauss): '+ STRING(snr, FORMAT='(F6.2)')+ ' '+ plus_minus+' ' + $
        STRING(snr_error, FORMAT='(F4.2)'),$
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
;XYOUTS, 0.11, 0.10, 'Symmetry:  '+STRCOMPRESS(ew, /REMOVE_ALL),$
;        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.11, 0.1, 'Central '+lam+':  '+STRING(center+650, FORMAT='(F6.2)')+' nm', $
        ALIGNMENT=0, CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS,0.11, 0.085, '(SHORT DASHED)', ALIGNMENT=0, CHARSIZE=0.75, $
       CHARTHICK=1.5, /NORMAL

DEVICE, /CLOSE_FILE
SET_PLOT, 'X'

;;move the newly created postscript
SPAWN, 'mv -f Halpha*.ps ' +root_direc+ targets_dir_name+'/HD' + star_hd + '/Halpha_plots'




CD, root_direc ;;cd back to the original location

jump:;jump here (skipping star_report section) when spec analyzer is ran alone

END
