;;======================================================================;;
;; Star Report Pipeline. Central Michigan University
;;
;; noise_analyzer.pro
;;
;; The goal of this procedure is to calculate a noise value for the
;; continuum of the spectrum and use it to find the signal to noise
;; ratio (SNR). The procedure begins by finding the initial sigma
;; initial sigma using a histogram and Gaussian fit
;; technique for the continuum of a user defined space. The default is
;; wave indices 475 - 511. In order to get an accurate SNR value,
;; everything that isn't photon noise must be
;; excluded from the continuum, e.g. absorptions from various
;; sources. With the value for initial sigma, a boxing method
;; is used to exclude groups of points rather than individual
;; points. The box size has a default size of 17 indices but can change
;; depending on the value of initial sigma. The program will
;; then analyze the points in that box and count how many points were
;; not contained within 1*(initial sigma). If the number of
;; points not contained is greater than a default value of 3 (this
;; number also varies between certain values for
;; initial sigma, all of the points within that box are
;; excluded from the continuum. A continuum mask array called
;; 'mask_flux' is created containing only 1's and 0's with the
;; 1's representing indices to be excluded. The array of flux
;; values, flux, is then multiplied by this mask creating an array
;; where all flux values to be included in the calculation are now
;; 0. This new array is used to over-plot the excluded regions on the
;; spectrum. Finally to find the SNR, all the points that make up the
;; continuum are plotted in a histogram. The value for sigma given by
;; the histogram (gauss_sigma) is then divided into the mean value of
;; all the points to provided a value for SNR.
;;
;; The end of this program is only used when noise_analyzer.pro is
;; being used as a part of star_report.pro. It creates a postscript
;; file for general diagnostics for the SNR calculation. It includes
;; things such as, the histogram plot of the continuum, SNR values
;; calculated with both standard deviations and the Gaussian
;; fit. Remember that the SNR calculated using the gauss_sigma is the
;; one that is used for further purposes. The SNR calculated using the
;; initial sigma is shown here purely for comparative purposes. These
;; files get stored in a sub-directory called 'Diagnostics' in
;; the star's directory.
;;
;;Author: W. Glennon Fagan
;;======================================================================;;


PRO noise_analyzer, color_table, smooth_flux, spec_plot_alone, temp_file, $
                    GAUSS_SIGMA=gauss_sigma, INTSIG_WAVE_INT=wave_int, $
                    INTSIG_WAVE_FIN=wave_fin, INT_SIG=int_sig
RESOLVE_ALL, /QUIET
COMMON spec_block
COMMON star_block

;;set this variable to 1 for diagnostics and 0 for no diagnostics.
diagnostics = 1
 

;;********************************* CREATE INITIAL SIGMA *********************************;;

;;hard define region will be used to determine an initial sigma
;;define wave indecies for initial standard deviation
wave_int = value_locator(wave, 652.55)
wave_fin = value_locator(wave, 653.09)
delta_lambda = STRING(wave[wave_fin] - wave[wave_int], FORMAT='(F5.3)')

;;define box size when scanning using the initial standard deviation
box_size = 0
;;define the maximum amount of uncontained points before the group
;;will be excluded
max_uncon = 0


;;Determine a standard deviation for the specified region
;int_sig = STDDEV(flux[wave_int : wave_fin])

;;int sig found using a histogram for the specified region
min = .945D
binsize = 0.0085D
nbins = 13

hist_x = (DINDGEN(nbins)/DOUBLE(nbins-1)) * (DOUBLE(nbins-1)*binsize) + $
         (min + 0.5*binsize)
int_hist = DOUBLE(HISTOGRAM(flux[wave_int : wave_fin], MIN = min, BINSIZE=binsize, $
                            NBINS = nbins))
int_gauss = GAUSSFIT(hist_x, int_hist, int_para, SIGMA=int_para_error, NTERMS=3)
int_sig = int_para[2]

;;define box size when scanning using the initial standard deviation
;;and define the maximum amount of uncontained points before the group
;;will be excluded 
IF int_sig GT 0.0125 THEN BEGIN
   box_size = 17
   max_uncon = 4
ENDIF ELSE IF int_sig GT 0.01 THEN BEGIN
   box_size = 15
   max_uncon = 4
ENDIF ELSE IF int_sig GT 0.007 THEN BEGIN
   box_size = 14
   max_uncon = 5
ENDIF ELSE IF int_sig GT 0.0051 THEN BEGIN
   box_size = 13
   max_uncon = 5
ENDIF ELSE BEGIN
   box_size = 12
   max_uncon = 9
ENDELSE




;;****************************** CREATE THE CONTINUUM MASK ********************************;;

deg_sig = 1 ;;To what degree sigma will points be excluded
cont_mask = DBLARR(N_ELEMENTS(wave))

con = 0
uncon = 0
exclude_box = 0
;;for loop to execute this box exclusion method. It is note a running
;;exclusion. For this method to work it needs to move at an increment
;;of the box size
FOR i=0, N_ELEMENTS(flux)-box_size-1,box_size DO BEGIN
   FOR j=i, i+box_size-1 DO BEGIN
      IF flux[j] GT 1+(deg_sig*int_sig) || flux[j] LT 1-(deg_sig*int_sig) THEN BEGIN
         uncon++
      ENDIF ELSE BEGIN
         con++
      ENDELSE
      
      ;;if there is a cosmic ray exclude the entire box.
      IF flux[j] GT 1.1 THEN BEGIN
         exclude_box = 1 
      ENDIF
   ENDFOR
   
   ;;Create a mask array where 0 are points that will NOT be
   ;;included in continuum calculation
   IF uncon GT max_uncon || exclude_box EQ 1 THEN BEGIN
      cont_mask[i : i + box_size - 1] = 0
   ENDIF ELSE BEGIN
      cont_mask[i : i + box_size - 1] = 1
   ENDELSE
;stop
con = 0
uncon = 0
exclude_box = 0
ENDFOR

;;for loop to do the same thing as above but for the very last box of
;;points. The for loop above can't analyze the last box number
;;worth of points so thats what is happening below.
FOR i=N_ELEMENTS(flux)-box_size-1, N_ELEMENTS(flux)-1 DO BEGIN
   IF flux[j] GT 1+(deg_sig*int_sig) || flux[j] LT 1-(deg_sig*int_sig) THEN BEGIN
         uncon++
      ENDIF ELSE BEGIN
         con++
      ENDELSE
ENDFOR
   IF uncon GT max_uncon THEN BEGIN
      cont_mask[N_ELEMENTS(flux)-box_size-1 : N_ELEMENTS(flux)-1] = 0
   ENDIF ELSE BEGIN
      cont_mask[N_ELEMENTS(flux)-box_size-1 : N_ELEMENTS(flux)-1] = 1
   ENDELSE

;;multiply our cont_mask with the flux. Mask flux now holds flux
;;values for points that WILL be included in the continuum noise
mask_flux = cont_mask * flux

;;fill variable with the # of points included in continuum. This is
;;used below in the IF statement. If the number of points in the mask
;;flux is less than the number of points in one box then there was
;;most likely in the calculation. print N/A.
num_points = N_ELEMENTS(WHERE(mask_flux NE 0))
IF num_points LT box_size THEN BEGIN
   snr=-1
   snr_error=-1
   
   IF spec_plot_alone EQ 1 THEN BEGIN
      PRINT, 'SNR IS N/A'
   ENDIF ELSE BEGIN
      analysis_logs = [analysis_logs,'SNR IS N/A']
   ENDELSE

   flag[2] = 3
   GOTO, skip
ENDIF




;;**************************** DETERMINE THE ADJUSTED SIGMA *******************************;;

;; histogram parameters
min = .945D
binsize = 0.0085D
nbins = 13

;; create histogram x values for histogram, centered on unity (one).
hist_x = (DINDGEN(nbins)/DOUBLE(nbins-1)) * (DOUBLE(nbins-1)*binsize) + $
         (min + 0.5*binsize)

;;create a histogram for the flux values included in the
;;continuum. and then fit a gauss to them.
hist_flux = HISTOGRAM(flux[WHERE((mask_flux) NE 0)], MIN = min, BINSIZE=binsize, NBINS = nbins)
gauss_flux = GAUSSFIT(hist_x, hist_flux, para, SIGMA=para_error, NTERMS=3)

;;store values from guass fit in various variables.
height = para[0]/N_ELEMENTS(flux)
hist_center = para[1]
gauss_sigma = para[2]




;;******************************** SNR CALCULATION ********************************;; 

;;**********************************
;;Plot showing the regions that are used for the photon noise
;WINDOW, 1
;plot, wave,(flux*cont_mask),PSYM=-4
;;**********************************

;;calculation using sigma from gauss fit
snr = MEAN(mask_flux[WHERE(mask_flux NE 0)]) / gauss_sigma
snr_error = (para_error[2] / gauss_sigma) * snr

;;calculation using std dev
snr2 = MEAN(mask_flux[WHERE(mask_flux NE 0)]) / STDDEV(mask_flux[WHERE(mask_flux NE 0)])

point_inc = N_ELEMENTS(where(mask_flux NE 0))

;;calculation of the number of regions that were used for continuum
reg_inc = 1
FOR i=0,N_ELEMENTS(flux)-2 DO BEGIN
   IF (mask_flux[i] NE 0 && mask_flux[i+1] EQ 0) THEN BEGIN
      reg_inc++
   ENDIF
ENDFOR




;;******************************* PRINTING ************************************;;

text = '   The number of points included in Photon noise:'+ STRCOMPRESS(STRING(point_inc))
IF spec_plot_alone EQ 1 THEN BEGIN
   PRINT, text
ENDIF ELSE BEGIN
   analysis_logs = [analysis_logs,text]
ENDELSE


warn_num = FLOOR(N_ELEMENTS(FLUX)* 0.1)
IF point_inc LT warn_num THEN BEGIN
   
   text = '   **CAUTION - Number of points included in noise calculation under 10%***'
   IF spec_plot_alone EQ 1 THEN BEGIN
      PRINT, text
   ENDIF ELSE BEGIN
      analysis_logs = [analysis_logs,text]
   ENDELSE
   
   flag[2] = 3
ENDIF

text = '   Number of regions used:'+ STRCOMPRESS(STRING(reg_inc))
IF spec_plot_alone EQ 1 THEN BEGIN
   PRINT, text
   PRINT, ' '
ENDIF ELSE BEGIN
   analysis_logs = [analysis_logs,text]
   analysis_logs = [analysis_logs,'']
ENDELSE

skip:




;;******************************* CREATE DIAGNOSTIC REPORT ************************************;; 

;;if statement to skip diagonostics if necessary
IF ~diagnostics || spec_plot_alone THEN GOTO, move

IF snr NE -1 THEN BEGIN
;;Create a model gauss fit using the gauss parameters found earlier
   increment = (1.06 - 0.96)/100
   model_x = DINDGEN(100, INCREMENT=increment, START=0.96)
   model_gauss = DBLARR(100)
   FOR i=0, 99 DO BEGIN
      z = (model_x[i] - para[1]) / (para[2])
      model_gauss[i] = (para[0]) * EXP((-(z)^2)/2) 
   ENDFOR 
ENDIF

plothick = 3
lineThick = 2


!P.FONT= -1 ;;Set to -1 for Hershey fonts
SET_PLOT, 'ps' ;select post script
DEVICE, DECOMPOSED = 0 ;Colors are taken from a color table (not specified by R, G, B)
LOADCT, 39, /SILENT ;load color table

DEVICE, FILENAME='SNR_calc_'+temp_file+'.ps', ENCAPSULATED=0, XSIZE=8.5, YSIZE=10, $
        XOFFSET=0, YOFFSET=0, /INCHES

;;get index for proper obs date
this_date = WHERE(analysis_frames EQ temp_file)

;;begin by plotting the entire contiuum for the frame
PLOT, wave, flux, XTITLE='Wavelength (nm)', YTITLE='F/Fc', XRANGE = [649, 659], yrange = [0.85,1.15],$
      XSTYLE=1, CHARSIZE=1.2, CHARTHICK=plotthick, XTHICK=plotthick, YTHICK=plotthick, THICK=linethick, $
      TITLE='Frame = '+temp_file+' / Obs. Date = '+analysis_dates[this_date], $
      POS=[.1, 0.53, 0.95, 0.92], /DATA, /NOERASE
OPLOT, [649.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=color_table[2]
OPLOT,  [wave[wave_int],wave[wave_int]], [0, 10]
OPLOT,  [wave[wave_fin],wave[wave_fin]], [0, 10]

;;use the same method as above to OPLOT points not in continuum with red
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
      OPLOT, wave[WHERE(temp_flux NE 0)], flux[WHERE(temp_flux NE 0)], $
                                     PSYM = -3,THICK=linethick, COLOR=color_table[1]
   ENDIF ELSE BEGIN
      BREAK
   ENDELSE
   temp_flux = DBLARR(N_ELEMENTS(flux))
   j = new_start+1
   counter ++
ENDWHILE

IF snr NE -1 THEN BEGIN
;;Next we will plot the histogram with gauss fit and parameters
plot, hist_x, hist_flux, TITLE = 'Noise Histogram', XTITLE = 'FLUX',$                           
      YTITLE = 'COUNTS', PSYM = 10,  XSTYLE=1, CHARSIZE=1.2, CHARTHICK=plotThick, XTHICK=plotthick, $
      YTHICK=plotthick, THICK=linethick, POS=[.1, 0.1, 0.5, 0.45], /DATA, /NOERASE   
oplot, hist_x, gauss_flux, PSYM = 4, THICK=linethick, COLOR=color_table[6] 
OPLOT, model_x, model_gauss, PSYM=-3, THICK=linethick, COLOR=color_table[1]
OPLOT, [para[1], para[1]], [0, para[0]+ 100]
ENDIF

;;define some symbols
plus_minus = '!M+!X'
del = '!7D!X'
lam = '!7k!X'

;;Put other information on the report

XYOUTS, 0.03, 0.01, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.20, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.5, 0.9625, $
        '_________________________________________________________________________________________________',$
        ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL

XYOUTS, 0.57, 0.4, 'Initial Sigma = '+STRING(int_sig, FORMAT='(F6.4)'), $
        ALIGNMENT=0, CHARSIZE=1.2, CHARTHICK=1.8, /NORMAL

XYOUTS, 0.57, 0.37, '['+STRING(wave[wave_int],FORMAT='(F7.3)')+ ': '+$
        STRING(wave[wave_fin],FORMAT='(F7.3)')+', '+del+lam+' = '+delta_lambda+']', $
        ALIGNMENT=0, CHARSIZE=1.15, CHARTHICK=1.8, /NORMAL

IF snr NE -1 THEN XYOUTS, 0.57, 0.28, 'GAUSS SIGMA = ' + STRING(gauss_sigma, FORMAT='(F7.5)') $
                          + ' '+ plus_minus+' ' + STRING(para_error[2], FORMAT='(F7.5)'), $
                          ALIGNMENT=0, CHARSIZE=1.2, CHARTHICK=1.8, /NORMAL

XYOUTS, 0.57, 0.25, 'SNR = '+STRING(snr, FORMAT='(F5.1)')+ ' '+ plus_minus+' ' + $
        STRING(snr_error, FORMAT='(F3.1)'), ALIGNMENT=0,$
        CHARSIZE=1.2, CHARTHICK=1.8, /NORMAL

XYOUTS, 0.57, 0.16, 'STDDEV = ' + STRCOMPRESS(STDDEV(mask_flux[WHERE(mask_flux NE 0)]), /REMOVE_ALL), $ 
        ALIGNMENT=0, CHARSIZE=1.2, CHARTHICK=1.8, /NORMAL
IF snr NE -1 THEN BEGIN
XYOUTS, 0.57, 0.13, 'SNR = '+STRCOMPRESS(snr2, /REMOVE_ALL), ALIGNMENT=0,$
        CHARSIZE=1.2, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.1, 0.03, 'Peak Height: '+ STRING(para[0], FORMAT='(F6.2)'), ALIGNMENT=0, $
        CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.5, 0.03, 'Peak Center: '+ STRING(para[1], FORMAT='(F4.2)'), ALIGNMENT=1, $
        CHARSIZE=1.1, CHARTHICK=1.8, /NORMAL
ENDIF

DEVICE, /CLOSE_FILE
SET_PLOT, 'X'

;;move the file to the diagonostics directory for that star
SPAWN, 'mv -f SNR_calc_'+temp_file+'.ps ' + root_direc + targets_dir_name + '/HD' + $
       star_hd + '/Diagnostics/' 

move:
END
