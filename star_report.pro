;;======================================================================;;
;; Star Report Pipeline, Central Michigan University
;;
;; star_report.pro
;;
;; The purpose of this procedure is to collect all neccessary data
;; from the rest of the pipeline,and use it to generate the Star
;; Report. The procedure will generate an on-screen report and a
;; postscript version that will be copied to both the star's
;; directory and the 'Star_Reports' directory.  
;;
;; Authors: Christian Hannah & Glennon Fagan
;;
;;======================================================================;;


PRO star_report,star_list
COMMON abort_block, no_exp_frames, no_hdhr_num, no_match_frames, user_entry
COMMON spec_block, cc_waves, center, EW, ew_error, flag, flux, fwhm, ha_waves, mask_flux, peak_flux, peak_num, peak_wave, peak1_flux, peak1_wave, peak2_flux, peak2_wave, shell_para, snr, snr_error, symmetry, vr_ratio, wave
COMMON star_block, analysis_logs, root_direc, available_dates, available_frames, analysis_dates, analysis_frames, mult_stars, star_dec, star_hd, star_hr, star_proper_name, star_name, star_ra,  star_vmag, start_time, sub_direc, targets_dir_name, trend, ver_num


;;======================= VERSION NUMBER =======================;;
                          ver_num = '1.0'
;;==============================================================;;



;;*******************************************
;; Define directory to store your star's data
IF ~ISA(targets_dir_name) THEN BEGIN
   ;targets_dir_name = 'Targets';default
   targets_dir_name = 'Primary_Targets'
   ;targets_dir_name = 'Secondary_Targets'
   ;targets_dir_name = 'Misc_Targets'
ENDIF
;;*******************************************

;;if statement to define working star name from star list if the
;;program is being ran for multiple targets in a loop
mult_stars = 0
IF ISA(star_list) THEN mult_stars = 1
IF mult_stars THEN star_name = star_list


SPAWN, 'pwd', root_direc
root_direc = root_direc+'/' ;;working directory
sub_direc = 'Spectra' ;;directory the spectra files are stored in 


;; switch for using the shifted spectra from halpha_centralizer.pro
;; or cross_correlator.pro *** MUST BE DIFFERENT or both are ignored ***  
h_align = 0
cc_align = 1




;; ********************************* CALL PROCEDURE Coordinator.pro **************************************;;

coordinator, BP_DATE_IND=bp_date_ind, CENTERS=centers, EWS_ERROR=EWs_error, $
             EWS_TRANS=EWs, FLUXS=fluxs, FRAME_DATES_EXCLUDED = frame_dates_excluded, $
             FRAMES_EXCLUDED = frames_excluded, FWHMS=fwhms, NUM_ANALYZED=num_analyzed, $
             N_FILES=n_files, NUM_SING=num_sing, NUM_DOUB=num_doub, PEAK_ONE_WAVES=peak_one_waves,$
             PEAK_ONE_FLUXS=peak_one_fluxs, PEAK_NUMS=peak_nums, PEAK_TWO_WAVES=peak_two_waves, $
             PEAK_TWO_FLUXS=peak_two_fluxs, RP_DATE_IND=rp_date_ind, SHELL_PARAS=shell_paras, $
             SING_PEAK_WAVES=sing_peak_waves, SING_PEAK_FLUXS=sing_peak_fluxs, SNRS=snrs, $
             SORTED_FILES=sorted_files, SORTED_JDATES=sorted_jdates, SP_DATE_IND=sp_date_ind, $
             SYMMETRIES=symmetries, VR_RATIOS=vr_ratios, WAVES=waves



;; If there were no file matches skip EVERYTHING  
IF no_exp_frames EQ 1 ||  no_hdhr_num EQ 1 ||  no_match_frames EQ 1 THEN BEGIN
   GOTO, jump
ENDIF





;;********************************** Set up the window *************************************;;
 
SET_PLOT, 'X' ; Select X-windows as a default display 
DEVICE, DECOMPOSED = 0                                                  
LOADCT, 39, /SILENT
color_table = [60,250,150,220,30,90,110,190,10]
color_text = ['blue','red','green','orange','violet','lt_blue','cyan','yellow','black']
!P.FONT = -1

IF ~ISA(star_list) THEN BEGIN
   dummyTitle='!3'
   WINDOW, 0, XSIZE=100, YSIZE=100
   PLOT, [0, 1], [0, 1], TITLE=dummyTitle
   WDELETE, 0
   WINDOW, 0, XSIZE=1300, YSIZE=1300
ENDIF



;;=========================================================================================;;
;;========================= BEGIN GENERATING ON-SCREEN REPORT =============================;; 
;;=========================================================================================;;

;; define some symbols to use in plots
plus_minus = '!M+!X'
lambda = '!7k!X'
alpha = '!7a!X'
delta = '!7D!X'




;;************ AVERAGE SPECTRA CHARACTERISTICS AND DEVELOP STRINGS FOR PRINTING ***********;;

;; average the EWs
avg_ew = TOTAL(EWs)/N_ELEMENTS(EWs)

;; average peak positions
avg_sing_peak_wave = TOTAL(sing_peak_waves[WHERE(sing_peak_waves NE 0D)])/N_ELEMENTS(WHERE(sing_peak_waves NE 0D))
avg_peak_one_wave = TOTAL(peak_one_waves[WHERE(peak_one_waves NE 0D)])/N_ELEMENTS(WHERE(peak_one_waves NE 0D))
avg_peak_two_wave = TOTAL(peak_two_waves[WHERE(peak_two_waves NE 0D)])/N_ELEMENTS(WHERE(peak_two_waves NE 0D))

;; average peak fluxs
avg_sing_peak_flux = TOTAL(sing_peak_fluxs[WHERE(sing_peak_fluxs NE 0D)])/N_ELEMENTS(WHERE(sing_peak_fluxs NE 0D))
avg_peak_one_flux = TOTAL(peak_one_fluxs[WHERE(peak_one_fluxs NE 0D)])/N_ELEMENTS(WHERE(peak_one_fluxs NE 0D))
avg_peak_two_flux = TOTAL(peak_two_fluxs[WHERE(peak_two_fluxs NE 0D)])/N_ELEMENTS(WHERE(peak_two_fluxs NE 0D))

;; average symmetries
avg_symmetry = TOTAL(symmetries)/num_analyzed

;; average snr
avg_snr = TOTAL(snrs[WHERE(snrs NE 0)])/ N_ELEMENTS(WHERE(snrs NE 0))

;; average fwhm
avg_fwhm = TOTAL(fwhms[WHERE(fwhms NE 0D)])/N_ELEMENTS(WHERE(fwhms NE 0D))

;; average shell_para
avg_shell_para = TOTAL(shell_paras[WHERE(shell_paras NE 0D)])/N_ELEMENTS(WHERE(shell_paras NE 0D))

;; average vr_ratio
avg_vr_ratio = TOTAL(vr_ratios[WHERE(vr_ratios NE 0D)])/N_ELEMENTS(WHERE(vr_ratios NE 0D))

;; convert the averages to strings 
aspw_str = STRING(avg_sing_peak_wave, FORMAT='(-F7.3)')+' nm'
apow_str = STRING(avg_peak_one_wave, FORMAT='(-F7.3)')+' nm'
aptw_str = STRING(avg_peak_two_wave, FORMAT='(-F7.3)')+' nm'

aspf_str = STRING(avg_sing_peak_flux, FORMAT='(F-5.2)')
apof_str = STRING(avg_peak_one_flux, FORMAT='(F-5.2)')
aptf_str = STRING(avg_peak_two_flux, FORMAT='(F-5.2)')

aew_str = STRING(avg_ew, FORMAT='(F-7.3)')+'nm'
asym_str = STRING(avg_symmetry, FORMAT='(F-8.2)')
asnr_str = STRING(avg_snr, FORMAT='(F-10.1)')
afwhm_str = STRING(avg_fwhm, FORMAT='(F-5.2)')+'nm'
ashell_para_str = STRING(avg_shell_para, FORMAT='(F-7.3)')
avr_ratio_str = STRING(avg_vr_ratio, FORMAT='(F-7.3)')


;; get rid of zeros in peak_num dependent arrays            
sing_peak_waves = sing_peak_waves[WHERE(sing_peak_waves NE 0D)]
sing_peak_fluxs = sing_peak_fluxs[WHERE(sing_peak_fluxs NE 0D)]
peak_one_waves = peak_one_waves[WHERE(peak_one_waves NE 0D)]
peak_one_fluxs = peak_one_fluxs[WHERE(peak_one_fluxs NE 0D)]
peak_two_waves = peak_two_waves[WHERE(peak_two_waves NE 0D)]
peak_two_fluxs = peak_two_fluxs[WHERE(peak_two_fluxs NE 0D)]
fwhms = fwhms[WHERE(fwhms NE 0D)]

;; Put max and min values for each characteristics into strings for ranges
max_ew_str = STRING(MAX(EWs), FORMAT='(F-10.3)')
min_ew_str = STRING(MIN(EWs), FORMAT='(F-10.3)')

max_spw_str = STRING(MAX(sing_peak_waves), FORMAT='(F-10.3)')
min_spw_str = STRING(MIN(sing_peak_waves), FORMAT='(F-10.3)')
max_pow_str = STRING(MAX(peak_one_waves), FORMAT='(F-10.3)')
min_pow_str = STRING(MIN(peak_one_waves), FORMAT='(F-10.3)')
max_ptw_str = STRING(MAX(peak_two_waves), FORMAT='(F-10.3)')
min_ptw_str = STRING(MIN(peak_two_waves), FORMAT='(F-10.3)')

max_spf_str = STRING(MAX(sing_peak_fluxs), FORMAT='(F-10.3)')
min_spf_str = STRING(MIN(sing_peak_fluxs), FORMAT='(F-10.3)')
max_pof_str = STRING(MAX(peak_one_fluxs), FORMAT='(F-10.3)')
min_pof_str = STRING(MIN(peak_one_fluxs), FORMAT='(F-10.3)')
max_ptf_str = STRING(MAX(peak_two_fluxs), FORMAT='(F-10.3)')
min_ptf_str = STRING(MIN(peak_two_fluxs), FORMAT='(F-10.3)')

max_sym_str = STRING(MAX(symmetries), FORMAT='(F-10.2)')
min_sym_str = STRING(MIN(symmetries), FORMAT='(F-10.2)')

max_snr_str = STRING(MAX(snrs), FORMAT='(F-10.1)')
min_snr_str = STRING(MIN(snrs), FORMAT='(F-10.1)')

max_fwhm_str = STRING(MAX(fwhms), FORMAT='(F-10.2)')
min_fwhm_str = STRING(MIN(fwhms), FORMAT='(F-10.2)')

max_shell_str = STRING(MAX(shell_paras), FORMAT='(F-10.3)')
min_shell_str = STRING(MIN(shell_paras), FORMAT='(F-10.3)')

max_vr_str = STRING(MAX(vr_ratios), FORMAT='(F-10.3)')
min_vr_str = STRING(MIN(vr_ratios), FORMAT='(F-10.3)')


;; make titles for displaying the data
spw_title = 'Single Peak '+lambda+': '
pow_title = 'Blue Peak '+lambda+': '
ptw_title = 'Red Peak '+lambda+': '

spf_title = 'Single Peak Flux: '
pof_title = 'Blue Peak Flux: '
ptf_title = 'Red Peak Flux: '

ew_title = 'H'+alpha+' EW: '
sym_title = 'Symmetry: '
snr_title = 'Cont. SNR: '
fwhm_title = 'FWHM: '
shell_para_title = 'Shell-Parameter: '
vr_ratio_title = 'log(V/R): '


;; average peak_nums to tell if the number of peaks changed over time.                                      
;; first omit any absorptions
apeak_nums = peak_nums[WHERE(peak_nums NE 0)]
avg_peak_num = TOTAL(apeak_nums)/N_ELEMENTS(WHERE(peak_nums NE 0))

;; use avg_peak_num to set up the string arrays for printing to the window
IF avg_peak_num EQ 0D THEN BEGIN ;;only absoprtion spectra 
   prnt_str_1 = STRARR(3)
   prnt_str_1 = [aew_str, '', asnr_str]
   prnt_str_title_1 = STRARR(3)
   prnt_str_title_1 = [ew_title, 'Absorption Spectra', snr_title]
   prnt_str_1_range = STRARR(3)
   prnt_str_1_range = ['['+max_ew_str+' : '+min_ew_str+']', '', '['+min_snr_str+' : '+max_snr_str+']']

   prnt_str_1_open = ['[','','[']
   prnt_str_1_close = [']','',']']
   prnt_str_1_colons = [':','',':']
   prnt_str_1_mins = [max_ew_str,'',min_snr_str]
   prnt_str_1_maxs = [min_ew_str,'',max_snr_str]

ENDIF ELSE IF avg_peak_num EQ 1 THEN BEGIN
   prnt_str_1 = STRARR(5)
   prnt_str_1 = [aew_str, aspw_str, aspf_str, asnr_str, afwhm_str]
   prnt_str_title_1 = STRARR(5)
   prnt_str_title_1 = [ew_title, spw_title, spf_title, snr_title, fwhm_title]
   prnt_str_1_range = STRARR(5)
   prnt_str_1_range = ['['+max_ew_str+' : '+min_ew_str+']', '['+min_spw_str+' : '+max_spw_str+']', $ 
                       '['+min_spf_str+' : '+max_spf_str+']', '['+min_snr_str+' : '+max_snr_str+']', $
                       '['+min_fwhm_str+' : '+max_fwhm_str+']']

   prnt_str_1_open = ['[','[','[','[','[']
   prnt_str_1_close = [']',']',']',']',']']
   prnt_str_1_colons = [':',':',':',':',':']
   prnt_str_1_mins = [max_ew_str,min_spw_str,min_spf_str,min_snr_str,min_fwhm_str]
   prnt_str_1_maxs = [min_ew_str,max_spw_str,max_spf_str,max_snr_str,max_fwhm_str]



ENDIF ELSE IF avg_peak_num EQ 2 THEN BEGIN
   prnt_str_1 = STRARR(8)
   prnt_str_1 = [aew_str, apow_str, apof_str, aptw_str, aptf_str, ashell_para_str, avr_ratio_str, asnr_str]
   prnt_str_title_1 = STRARR(8)
   prnt_str_title_1 = [ew_title, pow_title, pof_title, ptw_title, ptf_title, shell_para_title, vr_ratio_title, $
                       snr_title]
   prnt_str_1_range = STRARR(8)
   prnt_str_1_range = ['['+max_ew_str+' : '+min_ew_str+']', '['+min_pow_str+' : '+max_pow_str+']', $
                       '['+min_pof_str+' : '+max_pof_str+']', '['+min_ptw_str+' : '+max_ptw_str+']', $
                       '['+min_ptf_str+' : '+max_ptf_str+']', '['+min_shell_str+' : '+max_shell_str+']', $
                       '['+min_vr_str+' : '+max_vr_str+']', '['+min_snr_str+' : '+max_snr_str+']']

   prnt_str_1_open = ['[','[','[','[','[','[','[','[']
   prnt_str_1_close = [']',']',']',']',']',']',']',']']
   prnt_str_1_colons = [':',':',':',':',':',':',':',':']
   prnt_str_1_mins = [max_ew_str,min_pow_str,min_pof_str,min_ptw_str,min_ptf_str,min_shell_str,$
                      min_vr_str,min_snr_str]
   prnt_str_1_maxs = [min_ew_str,max_pow_str,max_pof_str,max_ptw_str,max_ptf_str,max_shell_str,$
                      max_vr_str,max_snr_str]


ENDIF ELSE BEGIN
   prnt_str_1 = STRARR(11)
   prnt_str_1 = [aew_str, aspw_str, aspf_str, apow_str, apof_str, aptw_str, aptf_str, ashell_para_str, avr_ratio_str, asnr_str, afwhm_str]
   prnt_str_title_1 = STRARR(11)
   prnt_str_title_1 = [ew_title, spw_title, spf_title, pow_title, pof_title, ptw_title, ptf_title, $
                       shell_para_title, vr_ratio_title, snr_title, fwhm_title]
   prnt_str_1_range = STRARR(11)
   prnt_str_1_range = ['['+max_ew_str+' : '+min_ew_str+']', '['+min_spw_str+' : '+max_spw_str+']', $
                       '['+min_spf_str+' : '+max_spf_str+']', '['+min_pow_str+' : '+max_pow_str+']', $
                       '['+min_pof_str+' : '+max_pof_str+']', '['+min_ptw_str+' : '+max_ptw_str+']', $
                       '['+min_ptf_str+' : '+max_ptf_str+']', '['+min_shell_str+' : '+max_shell_str+']', $
                       '['+min_vr_str+' : '+max_vr_str+']','['+min_snr_str+' : '+max_snr_str+']', $
                       '['+min_fwhm_str+' : '+max_fwhm_str+']']

   prnt_str_1_open = ['[','[','[','[','[','[','[','[','[','[','[']
   prnt_str_1_close = [']',']',']',']',']',']',']',']',']',']',']']
   prnt_str_1_colons = [':',':',':',':',':',':',':',':',':',':',':']
   prnt_str_1_mins = [max_ew_str,min_spw_str,min_spf_str,min_pow_str,min_pof_str,min_ptw_str,$
                      min_ptf_str,min_shell_str,min_vr_str,min_snr_str,min_fwhm_str]
   prnt_str_1_maxs = [min_ew_str,max_spw_str,max_spf_str,max_pow_str,max_pof_str,max_ptw_str,$
                      max_ptf_str,max_shell_str,max_vr_str,max_snr_str,max_fwhm_str]

ENDELSE




;;******************************** OBTAIN MOST RECENT SPECTRUM CHARACTERISTICS ***************************;;

;; save most recent spectrum info into an array, convert the values to
;; strings
spw = sing_peak_waves[N_ELEMENTS(sing_peak_waves)-1]
spw_str = STRING(spw, FORMAT='(F-7.3)')+' nm'
pow = peak_one_waves[N_ELEMENTS(peak_one_waves)-1]
pow_str = STRING(pow, FORMAT='(F-7.3)')+' nm'
ptw = peak_two_waves[N_ELEMENTS(peak_two_waves)-1]
ptw_str = STRING(ptw, FORMAT='(F-7.3)')+' nm'

spf = sing_peak_fluxs[N_ELEMENTS(sing_peak_fluxss)-1]
spf_str = STRING(spf, FORMAT='(F-5.2)')
pof = peak_one_fluxs[N_ELEMENTS(peak_one_fluxs)-1]
pof_str = STRING(pof, FORMAT='(F-5.2)')
ptf = peak_two_fluxs[N_ELEMENTS(peak_two_fluxs)-1]
ptf_str = STRING(ptf, FORMAT='(F-5.2)')

ew = EWs[num_analyzed-1]
ew_err = EWs_error[num_analyzed-1]
ew_str = STRING(ew, FORMAT='(F-7.3)')+plus_minus+' '+STRING(ew_err, FORMAT='(F-6.3)')+'nm'
sym = symmetries[num_analyzed-1]
sym_str = STRING(sym, FORMAT='(F-8.2)')
snr = snrs[num_analyzed-1]
snr_str = STRING(snr, FORMAT='(F-10.1)')
fwhm = fwhms[N_ELEMENTS(fwhms)-1]
fwhm_str = STRING(fwhm, FORMAT='(F-5.2)')+'nm'

;; store above values into arrays based on the peak_num
IF peak_nums[num_analyzed-1] EQ 1 THEN BEGIN
   prnt_str_2 = STRARR(6)
   prnt_str_2 = [ew_str, spw_str, spf_str, sym_str, snr_str, fwhm_str]
   prnt_str_title_2 = STRARR(6)
   prnt_str_title_2 = [ew_title, spw_title, spf_title, sym_title, snr_title, fwhm_title]
ENDIF ELSE IF peak_nums[num_analyzed-1] EQ 2 THEN BEGIN
   prnt_str_2 = STRARR(7)
   prnt_str_2 = [ew_str, pow_str, pof_str, ptw_str, ptf_str, sym_str, snr_str]
   prnt_str_title_2 = STRARR(7)
   prnt_str_title_2 = [ew_title, pow_title, pof_title, ptw_title, ptf_title, sym_title, snr_title]
ENDIF ELSE BEGIN
   prnt_str_2 = STRARR(3)
   prnt_str_2 = [ew_str,' ', snr_str]
   prnt_str_title_2 = STRARR(3)
   prnt_str_title_2 = [ew_title, 'Absorption Spectrum', snr_title]
ENDELSE


;; define the y locations used for adding the characteristics
prnt_str_yloc = DBLARR(11)
prnt_str_yloc = [0.25, 0.23, 0.21, 0.19, 0.17, 0.15, 0.13, 0.11, 0.09, 0.07, 0.05]





;;********************************* ADJUST PEAK DATA FOR EASE OF PLOTTING *****************************;;

;; adjust the julien dates for easier reading
jdates = sorted_jdates - 2450000

sp_jdates = jdates[sp_date_ind[WHERE(sp_date_ind NE -1)]]
bp_jdates = jdates[bp_date_ind[WHERE(bp_date_ind NE -1)]]
rp_jdates = jdates[rp_date_ind[WHERE(rp_date_ind NE -1)]]


;; single peaks
IF N_ELEMENTS(WHERE(sing_peak_fluxs NE 0)) GT 1 || WHERE(sing_peak_fluxs NE 0) NE -1 THEN BEGIN
   sp_fluxs = sing_peak_fluxs[WHERE(sing_peak_fluxs NE 0)]
   max_sp = MAX(sp_fluxs)
   min_sp = MIN(sp_fluxs)
ENDIF ELSE BEGIN
   max_sp = 0D
   min_sp = 0D
ENDELSE

;;double peaks
IF N_ELEMENTS(WHERE(peak_one_fluxs NE 0)) GT 1 || WHERE(peak_one_fluxs NE 0) NE -1 THEN BEGIN
   bp_fluxs = peak_one_fluxs[WHERE(peak_one_fluxs NE 0)]
   max_bp = MAX(bp_fluxs)
   min_bp = MIN(bp_fluxs)
ENDIF ELSE BEGIN
   max_bp = 0D
   min_bp = 0D
ENDELSE

IF N_ELEMENTS(WHERE(peak_two_fluxs NE 0)) GT 1 || WHERE(peak_two_fluxs NE 0) NE -1 THEN BEGIN
   rp_fluxs = peak_two_fluxs[WHERE(peak_two_fluxs NE 0)]
   max_rp = MAX(rp_fluxs)
   min_rp = MIN(rp_fluxs)
ENDIF ELSE BEGIN
   max_rp = 0D
   min_rp = 0D
ENDELSE


;; obtain some max and min values for specifying yrange in the plot
maxes = [max_sp, max_bp, max_rp]
max_peak = MAX(maxes)

mins = [min_sp, min_bp, min_rp]
min_peak = MIN(mins[WHERE(mins NE 0)])





;;******************************* OBTAIN REST OF VALUES TO ADD TO THE REPORT *******************************;;

;; make num_analyzed a string
num_spec = STRTRIM(num_analyzed, 2)


;; convert most recent julian date back to calendar
mr_date = LONG(sorted_jdates[num_analyzed-1])
CALDAT, mr_date, Ml, Dl, Yl
;; convert Ml, Dl, Yl to strings
Ml_str = STRTRIM(Ml, 2)
Dl_str =STRTRIM(Dl, 2)
Yl_str =STRTRIM(Yl, 2)


;; convert oldest julian date back to calendar
first_date = LONG(sorted_jdates[0])
CALDAT, first_date, Mf, Df, Yf
;; convert Mf, Df, Yf to strings 
Mf_str = STRTRIM(Mf, 2)
Df_str =STRTRIM(Df, 2)
Yf_str =STRTRIM(Yf, 2)


;;convert sorted_jdates to fractional years
cal_dates = DBLARR(N_ELEMENTS(sorted_jdates)) ;;array to store fractional calendar year dates

FOR i=0, N_ELEMENTS(sorted_jdates)-1 DO BEGIN
   CALDAT, sorted_jdates[i], M, D, Y ;;convert JD to calendar 
   year_jdate = JULDAY(1,0,Y) ;;convert the year back to JD
   difference = sorted_jdates[i] - year_jdate ;;subtract the JD's to get number of days
   fraction = difference/365.25 ;;convert days to a fractional year
   frac_year = Y + fraction 
   cal_dates[i] = frac_year
ENDFOR


;;determine/add number of nights observed to the report
num_nights = 1
FOR i=0, num_analyzed-2 DO BEGIN
   IF jdates[i] NE jdates[i+1] THEN num_nights+=1
ENDFOR


;;specify number of days used to pad the time axis in plots dealing
;;with time.
axis_padding = 365
;;convert padding to its equivalent in fractional calendar years
cal_pad = 365/365.25




;;****************************************** PLOT ALL SPECTRA ****************************************;;

;;If statement to skip plotting to window if star_report is being looped
IF ~ISA(star_list) THEN BEGIN

   ;; get max and min spectra using EWs 
   max_ind = MIN(WHERE(EWs EQ MIN(EWs)))
   min_ind = MIN(WHERE(EWs EQ MAX(EWs)))
   
   ;; get the proper scale for y axis
   heights = DBLARR(3)
   
   IF ISA(sing_peak_fluxs) THEN heights[0] = MAX(sing_peak_fluxs)
   IF ISA(peak_one_fluxs) THEN heights[1] = MAX(peak_one_fluxs)
   IF ISA(peak_two_fluxs) THEN heights[2] = MAX(peak_two_fluxs)
   
   ;; if all spectra are absorptions the max height should be 1
   IF N_ELEMENTS(WHERE(heights EQ 0)) EQ N_ELEMENTS(heights) THEN heights[0] = 1D
   
   height = MAX(heights)
   
   
   ;; Setup plot for all spectra  
   PLOT, [654, 659],[0, height+0.5], XTITLE=lambda+' (nm)', YTITLE='F/Fc', XSTYLE=1, $
         CHARSIZE=1.3, TITLE='All Spectra', POS=[0.06, 0.69, 0.477, 0.935], /NODATA
   
   ;; overplot a dashed line at the continuum level  
   OPLOT, [653.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=60
   
   ;; plot the spectra 
   IF ~h_align && cc_align THEN BEGIN
      FOR i=0, num_analyzed-1 DO BEGIN
         IF i EQ max_ind || i EQ min_ind THEN GOTO, hop_1 ;; so the max and min can be overplotted last
         OPLOT, cc_waves[*,i], fluxs[*,i], PSYM=-3
         hop_1:
      ENDFOR
   
      ;; overplot max and min spectra                                         
      OPLOT, cc_waves[*,min_ind], fluxs[*,min_ind], PSYM=-3, COLOR=250
      OPLOT, cc_waves[*,max_ind], fluxs[*,max_ind], PSYM=-3, COLOR=150

   ENDIF ELSE IF h_align && ~cc_align THEN BEGIN
      FOR i=0, num_analyzed-1 DO BEGIN
         IF i EQ max_ind || i EQ min_ind THEN GOTO, hop_2 ;; so the max and min can be overplotted last
         OPLOT, ha_waves[*,i], fluxs[*,i], PSYM=-3
         hop_2:
      ENDFOR
      
      ;; overplot max and min spectra
      OPLOT, ha_waves[*,min_ind], fluxs[*,min_ind], PSYM=-3, COLOR=250
      OPLOT, ha_waves[*,max_ind], fluxs[*,max_ind], PSYM=-3, COLOR=150
   ENDIF ELSE BEGIN
      FOR i=0, num_analyzed-1 DO BEGIN
         IF i EQ max_ind || i EQ min_ind THEN GOTO, hop_3 ;; so the max and min can be overplotted last
         OPLOT, waves[*,i], fluxs[*,i], PSYM=-3
         hop_3:
      ENDFOR
      
      ;; overplot max and min spectra
      OPLOT, waves[*,min_ind], fluxs[*,min_ind], PSYM=-3, COLOR=250
      OPLOT, waves[*,max_ind], fluxs[*,max_ind], PSYM=-3, COLOR=150
   ENDELSE
   
   
   ;;************************************** ADD PLOT OF EW vs. TIME ***************************************;;
   
   PLOT, jdates, EWs, XTITLE='JD (-2450000)', YTITLE='H'+alpha+' EW (nm)', $ 
         XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $ 
         YRANGE=[MAX(EWs)+0.2, MIN(EWs)-0.5], XSTYLE=9, CHARSIZE=1.3, $ 
         PSYM=3, POS=[0.541, 0.69, 0.959, 0.935], /NOERASE
   ;; add error bars
   down_err = EWs - EWs_error
   up_err = EWs + EWs_error
   ERRPLOT, jdates, down_err, up_err, /DATA
  
   ;;add an axis to the plot displaying fractional years
   AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=1.3
   XYOUTS, 0.267, 0.64, 'Years', ALIGNMENT=0.5, CHARSIZE=1.3, /NORMAL



   
   ;;********************************** ADD MOST RECENT SPECTRUM PLOT ************************************;;
   
   PLOT, waves[*,num_analyzed-1], fluxs[*,num_analyzed-1], XTITLE=lambda+' (nm)', YTITLE='F/Fc', $ 
         XRANGE = [654, 659], XSTYLE=1, CHARSIZE=1.3, $
         TITLE='Most Recent Spectrum ('+Ml_str+'/'+Dl_str+'/'+Yl_str+')', $ 
         PSYM=-3, POS=[0.06, 0.37, 0.477, 0.62], /NOERASE
   
   ;; oplot continuum
   OPLOT, [649.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=60
   
   

   
   ;;********************************** PEAK STRENGTH VS TIME PLOT ************************************;;   
   
   IF ~avg_peak_num EQ 0D THEN BEGIN ;;do not add plot if all absorptions
      ;; adjust the maximum value of the yrange used in the plot to avoid
      ;; the overlapping of text and points
      sep = max_peak-min_peak
      CASE 1 OF 
         (sep LE 1.2): const = 0.3
         (sep GT 1.2) AND (sep LE 2.8): const = 0.65
         (sep GT 2.8) AND (sep LE 5.0): const = 1.0
         (sep GT 5.0): const = 1.2
         ELSE: PRINT, '**WARNING: max_peak/min_peak not properly defined**'
      ENDCASE
      
      
      ;; set up the plot
      PLOT, [0,0], [0,0], XTITLE='JD (-2450000)', YTITLE='H'+alpha+' Peak Strength (F/Fc)', $
            XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
            YRANGE=[min_peak-0.2 , max_peak+const], XSTYLE=9, CHARSIZE=1.3, $
            PSYM=3, POS=[0.541, 0.37, 0.959, 0.62], /NOERASE
      
      ;; add single peak strengths
      IF ISA(sp_fluxs) THEN $
         OPLOT, sp_jdates, sp_fluxs, PSYM=5
      
      ;; add blue peak strengths
      IF ISA(bp_fluxs) THEN $
         OPLOT, bp_jdates, bp_fluxs, PSYM=4, COLOR=color_table[0]
      
      ;; add red peak strengths
      IF ISA(rp_fluxs) THEN $
         OPLOT, rp_jdates, rp_fluxs, PSYM=1, COLOR=color_table[1]
      
      ;;add an axis to the plot displaying fractional years
      AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=1.3
      XYOUTS, 0.749, 0.64, 'Years', ALIGNMENT=0.5, CHARSIZE=1.3, /NORMAL
   ENDIF





   ;;************************ PLOT SHELL PARAMETER AND V/R VS TIME OR FWHM VS TIME *************************;;

   ;;define percent 
   percent_dbl = N_ELEMENTS(rp_fluxs)/DOUBLE(num_spec)
   percent_sng = N_ELEMENTS(sp_fluxs)/DOUBLE(num_spec)

   IF ISA(rp_fluxs) && percent_sng LE percent_dbl THEN BEGIN 
      
      ;;remove zeros in shell_paras and vr_ratios
      shell_paras = shell_paras[WHERE(shell_paras NE 0D)]
      vr_ratios = vr_ratios[WHERE(vr_ratios NE 0D)]


      ;;set up plot for shell parameters
      PLOT, [0,0], [0,0], XTITLE='JD (-2450000)', $
            XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
            YRANGE=[0.0,MAX(shell_paras)+0.25*MAX(shell_paras)],$
            XSTYLE=9, YSTYLE=5, CHARSIZE=1.3, PSYM=3, POS=[0.541, 0.05, 0.959, 0.3], /NOERASE
      
      ;;plot a dashed line at minimum shell phase value (1.5)
      OPLOT, [jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], [1.5, 1.5], LINESTYLE=2

      ;;plot shell parameters
      OPLOT, bp_jdates, shell_paras, PSYM=4

      
      ;;add a yaxis to the plot for shell parameter 
      AXIS, YAxis=0, YRANGE=[0.0,MAX(shell_paras)+0.25*MAX(shell_paras)], YSTYLE=1, CHARSIZE=1.3, $
            YTITLE='Shell-Parameter'

      ;;set up pot for V/R
      PLOT, [0,0], [0,0], XTITLE='JD (-2450000)', $
            XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
            YRANGE=[MIN(vr_ratios)-.1,MAX(vr_ratios)+.3],$
            XSTYLE=9, YSTYLE=5, CHARSIZE=1.3, PSYM=3, POS=[0.541, 0.05, 0.959, 0.3], /NOERASE
      
      ;;oplot a dashed line at 0 for even peak strengths
      OPLOT, [jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], [0.0, 0.0], LINESTYLE=2, COLOR=150

      ;;plot V/R
      OPLOT, bp_jdates, vr_ratios, PSYM=1, COLOR=150

      ;;add a yaxis to the plot for shell parameter
      AXIS, YAxis=1, YRANGE=[MIN(vr_ratios)-.1,MAX(vr_ratios)+.3], YSTYLE=1, CHARSIZE=1.3, COLOR=150
      XYOUTS, 0.995, 0.175, 'log(V/R)', ALIGNMENT=0.5, ORIENTATION=90, CHARSIZE=1.3, COLOR=150, /NORMAL


      ;;add an axis to the plot displaying fractional years
      AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=1.3
      XYOUTS, 0.749, 0.32, 'Years', ALIGNMENT=0.5, CHARSIZE=1.3, /NORMAL

   ENDIF ELSE IF ~avg_peak_num EQ 0D THEN BEGIN ;;do not add plot if all absorptions 
      max_fwhm = MAX(fwhms)
      min_fwhm = MIN(fwhms)
      
      ;;make fwhms an array if there is only one
      IF N_ELEMENTS(fwhms) EQ 1 THEN fwhms = [fwhms]

      ;; set up plot for fwhms
      PLOT, sp_jdates, fwhms, XTITLE='JD (-2450000)', YTITLE='FWHM (nm)', $
            XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
            YRANGE=[MIN(min_fwhm)-.1,MAX(max_fwhm)+.1], $
            XSTYLE=9, CHARSIZE=1.3, PSYM=4, POS=[0.541, 0.05, 0.959, 0.3], /NOERASE

      ;;add an axis to the plot displaying fractional years
      AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=1.3
      XYOUTS, 0.749, 0.32, 'Years', ALIGNMENT=0.5, CHARSIZE=1.3, /NORMAL
   ENDIF
   

   ;;************************************ ADD TEXT TO THE WINDOW ****************************************;; 
   
   XYOUTS, 0.03, 0.004, 'Ver. '+ver_num+', '+systime(0), CHARSIZE=1.3, CHARTHICK=1.7, /NORMAL
   XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=3, CHARTHICK=1.8, /NORMAL
   XYOUTS, 0.18, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=2, CHARTHICK=1.8, /NORMAL
   XYOUTS, 0.82, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=2, CHARTHICK=1.8, /NORMAL

   XYOUTS, 0.12, 0.92, '# of Spectra: '+num_spec, ALIGNMENT=0.5, CHARSIZE=1.5, /NORMAL
   XYOUTS, 0.03, 0.28, 'Average Characteristics:', CHARSIZE=2.4, /NORMAL
   XYOUTS, 0.03, prnt_str_yloc, prnt_str_title_1, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.17, prnt_str_yloc, prnt_str_1, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.31, 0.28, '[Range]', CHARSIZE=2.0, /NORMAL
   XYOUTS, 0.29, prnt_str_yloc, prnt_str_1_open, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.30, prnt_str_yloc, prnt_str_1_mins, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.37, prnt_str_yloc, prnt_str_1_colons, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.39, prnt_str_yloc, prnt_str_1_maxs, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.45, prnt_str_yloc, prnt_str_1_close, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.555, 0.92, '# of Obs. Nights:'+STRCOMPRESS(num_nights), COLOR=220, CHARSIZE=1.5, /NORMAL
   XYOUTS, 0.785, 0.92, 'First Observed: '+Mf_str+'/'+Df_str+'/'+Yf_str, COLOR=220, CHARSIZE=1.5, /NORMAL
   XYOUTS, 0.785, 0.905, 'Last Observed: '+Ml_str+'/'+Dl_str+'/'+Yl_str, COLOR=220, CHARSIZE=1.5, /NORMAL

   XYOUTS, 0.03, 0.315, 'Obs. Disk Trend: '+ trend, CHARSIZE=1.8, /NORMAL
   XYOUTS, 0.03, 0.3125, '______________', CHARSIZE=1.8, /NORMAL

   ;;do not add these labels if all absorptions because there will be no plot 
   IF ~avg_peak_num EQ 0D THEN BEGIN  
      XYOUTS, 0.88, 0.6, 'Single Peak', CHARSIZE =1.5, /NORMAL
      XYOUTS, 0.88, 0.585, 'Blue Peak', COLOR=color_table[0], CHARSIZE =1.5, /NORMAL
      XYOUTS, 0.88, 0.57, 'Red Peak', COLOR=color_table[1], CHARSIZE =1.5, /NORMAL
   ENDIF


   XYOUTS, 0.0, 0.9625, $
      '________________________________________________________________________________________________', $
      CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL
   XYOUTS, 0.03, 0.275, '_____________________', CHARSIZE=2.4, /NORMAL
   
   ;;Printing the excluded frames takes place in this if statement                                             
   XYOUTS, 0.03, 0.016, 'Frames excluded:', CHARSIZE =1.5, /NORMAL
   IF frames_excluded[0] EQ -1 THEN BEGIN
      XYOUTS, 0.15, 0.016, 'None', CHARSIZE =1.5, /NORMAL
   ENDIF ELSE BEGIN
      j=0.15
      FOR i=0, N_ELEMENTS(frames_excluded)-1 DO BEGIN
         IF i EQ N_ELEMENTS(frames_excluded)-1 THEN BEGIN
            XYOUTS, j, 0.016, frames_excluded[i], CHARSIZE =1.5, /NORMAL
         ENDIF ELSE BEGIN
            XYOUTS, j, 0.016, frames_excluded[i]+',', CHARSIZE =1.5, /NORMAL
         ENDELSE
         j = j + 0.048
      ENDFOR
   ENDELSE

   ;; Signify to the user that the spectra are being shifted  
   IF h_align && ~cc_align THEN XYOUTS,  0.12, 0.905, 'H'+alpha+' Shifted', ALIGNMENT=0.5, $
                                         CHARSIZE=1.5, /NORMAL
   IF cc_align && ~h_align THEN XYOUTS,  0.12, 0.905, 'CC Shifted', ALIGNMENT=0.5, CHARSIZE=1.5, /NORMAL
   

   ;; If some peaks are single and some are double, indicate how many of each    
   IF num_sing NE 0 && num_doub NE 0 THEN BEGIN
      XYOUTS, 0.35, 0.92, '# of SP Spectra:'+STRCOMPRESS(num_sing), CHARSIZE=1.5, /NORMAL
      XYOUTS, 0.348, 0.905, '# of DP Spectra:'+STRCOMPRESS(num_doub), CHARSIZE=1.5, /NORMAL
      XYOUTS, 0.37, 0.89, '__ Max Emission', CHARSIZE=1.5, COLOR=150, /NORMAL
      XYOUTS, 0.37, 0.875, '__ Min Emission', CHARSIZE=1.5, COLOR=250, /NORMAL
   ENDIF ELSE BEGIN 
      XYOUTS, 0.37, 0.92, '__ Max Emission', CHARSIZE=1.5, COLOR=150, /NORMAL
      XYOUTS, 0.37, 0.905, '__ Min Emission', CHARSIZE=1.5, COLOR=250, /NORMAL
   ENDELSE
   

ENDIF;;if statement to stop printing to window if star_report is looped

;;=========================================================================================;;
;;================================ END ON-SCREEN REPORT ===================================;;
;;=========================================================================================;;






;;=========================================================================================;;
;;=============================== GENERATE POSTSCRIPT =====================================;;
;;=========================================================================================;;


plotThick = 3
lineThick = 4

SET_PLOT, 'ps'
DEVICE, /COLOR, BITS_PER_PIXEL=8
DEVICE, FILENAME='HD'+star_hd+'_report.ps', ENCAPSULATED=0
DEVICE, XSIZE=8, YSIZE=10.5, XOFFSET=0.25, YOFFSET=0.25, /INCHES




;;****************************************** PLOT ALL SPECTRA ****************************************;;
;; get max and min spectra using EWs
max_ind = MIN(WHERE(EWs EQ MIN(EWs)))
min_ind = MIN(WHERE(EWs EQ MAX(EWs)))

;; get the proper scale for y axis
heights = DBLARR(3)

IF ISA(sing_peak_fluxs) THEN heights[0] = MAX(sing_peak_fluxs)
IF ISA(peak_one_fluxs) THEN heights[1] = MAX(peak_one_fluxs)
IF ISA(peak_two_fluxs) THEN heights[2] = MAX(peak_two_fluxs)

;; if all spectra are absorptions the max height should be 1
IF N_ELEMENTS(WHERE(heights EQ 0)) EQ N_ELEMENTS(heights) THEN heights[0] = 1D

height = MAX(heights)


;; Setup plot for all spectra
PLOT, [654, 659],[0, height+0.5], XTITLE=lambda+' (nm)', YTITLE='F/Fc', XSTYLE=1, $
      CHARSIZE=0.75, TITLE='All Spectra', POS=[0.056, 0.69, 0.472, 0.94], /NODATA
;; overplot a dashed line at the continuum level
OPLOT, [653.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=60



;; plot the spectra                                                                                                   
IF ~h_align && cc_align THEN BEGIN
   FOR i=0, num_analyzed-1 DO BEGIN
      IF i EQ max_ind || i EQ min_ind THEN GOTO, hop_11 ;; so the max and min can be overplotted last
      OPLOT, cc_waves[*,i], fluxs[*,i], PSYM=-3
      hop_11:
   ENDFOR
   
   ;; overplot max and min spectra 
   OPLOT, cc_waves[*,min_ind], fluxs[*,min_ind], PSYM=-3, COLOR=250
   OPLOT, cc_waves[*,max_ind], fluxs[*,max_ind], PSYM=-3, COLOR=150
   
ENDIF ELSE IF h_align && ~cc_align THEN BEGIN
   FOR i=0, num_analyzed-1 DO BEGIN
      IF i EQ max_ind || i EQ min_ind THEN GOTO, hop_22 ;; so the max and min can be overplotted last
      OPLOT, ha_waves[*,i], fluxs[*,i], PSYM=-3
      hop_22:
   ENDFOR
   
   ;; overplot max and min spectra
   OPLOT, ha_waves[*,min_ind], fluxs[*,min_ind], PSYM=-3, COLOR=250
   OPLOT, ha_waves[*,max_ind], fluxs[*,max_ind], PSYM=-3, COLOR=150

ENDIF ELSE BEGIN
   FOR i=0, num_analyzed-1 DO BEGIN
      IF i EQ max_ind || i EQ min_ind THEN GOTO, hop_33 ;; so the max and min can be overplotted last
      OPLOT, waves[*,i], fluxs[*,i], PSYM=-3
      hop_33:
   ENDFOR
   
   ;; overplot max and min spectra
   OPLOT, waves[*,min_ind], fluxs[*,min_ind], PSYM=-3, COLOR=250
   OPLOT, waves[*,max_ind], fluxs[*,max_ind], PSYM=-3, COLOR=150

ENDELSE





;;************************************** ADD PLOT OF EW vs. TIME ***************************************;;

;; adjust the julian dates for easier reading
jdates = sorted_jdates - 2450000

PLOT, jdates, EWs, XTITLE='JD (-2450000)', YTITLE='H'+alpha+' EW (nm)', $
      XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
      YRANGE=[MAX(EWs)+0.3, MIN(EWs)-0.5], XSTYLE=9, CHARSIZE=0.75, $
      PSYM=3, POS=[0.535, 0.69, 0.954, 0.94], /NOERASE

;;add error bars
down_err = EWs - EWs_error
up_err = EWs + EWs_error
ERRPLOT, jdates, down_err, up_err, /DATA

;;add an axis to the plot displaying fractional years
AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=0.75
XYOUTS, 0.264, 0.64, 'Years', ALIGNMENT=0.5, CHARSIZE=0.75, /NORMAL







;;********************************** ADD MOST RECENT SPECTRUM PLOT ************************************;; 

PLOT, waves[*,num_analyzed-1], fluxs[*,num_analyzed-1], XTITLE=lambda+' (nm)', YTITLE='F/Fc', $
      XRANGE = [654, 659], XSTYLE=1, CHARSIZE=0.75, $
      TITLE='Most Recent Spectrum ('+Ml_str+'/'+Dl_str+'/'+Yl_str+')', $
      PSYM=-3, POS=[0.056, 0.37, 0.472, 0.62], /NOERASE


;; oplot continuum
OPLOT, [649.0, 659.0], [1.0, 1.0], LINESTYLE=2, COLOR=60




;;********************************** PEAK STRENGTH VS TIME PLOT ************************************;; 

IF ~avg_peak_num EQ 0D THEN BEGIN ;;do not add plot if all absorptions
   ;; adjust the maximum value of the yrange use din the plot to avoid
   ;; the overlapping of text and points
   IF max_peak-min_peak LE 1.2 THEN BEGIN
      const = 0.3
   ENDIF ELSE IF max_peak-min_peak GT 1.2 && max_peak-min_peak LE 2.8 THEN BEGIN
      const = 0.7
   ENDIF ELSE IF max_peak-min_peak GT 2.8 && max_peak-min_peak LE 5.0 THEN BEGIN
      const = 1.0
   ENDIF ELSE BEGIN
      const = 1.2
   ENDELSE
   
   
   ;; set up the plot
   PLOT, [0,0], [0,0], XTITLE='JD (-2450000)', YTITLE='H'+alpha+' Peak Strength (F/Fc)', $
         XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
         YRANGE=[min_peak-0.3, max_peak+const], XSTYLE=9, CHARSIZE=0.75, $
         PSYM=3, POS=[0.535, 0.37, 0.954, 0.62], /NOERASE
   
   ;; add single peak strengths
   IF ISA(sp_fluxs) THEN $
      OPLOT, sp_jdates, sp_fluxs, PSYM=5, SYMSIZE=0.5
   
   ;; add blue peak strengths
   IF ISA(bp_fluxs) THEN $
      OPLOT, bp_jdates, bp_fluxs, PSYM=4, SYMSIZE=0.5, COLOR=color_table[0]
   
   ;; add red peak strengths
   IF ISA(rp_fluxs) THEN $
      OPLOT, rp_jdates, rp_fluxs, PSYM=1, SYMSIZE=0.5, COLOR=color_table[1]
   
   ;;add an axis to the plot displaying fractional years 
   AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=0.75
   XYOUTS, 0.746, 0.64, 'Years', ALIGNMENT=0.5, CHARSIZE=0.75, /NORMAL
ENDIF



;;************************ PLOT SHELL PARAMETER AND V/R VS TIME OR FWHM VS TIME *************************;; 
IF ISA(rp_fluxs) THEN BEGIN
   
   ;;remove zeros in shell_paras and vr_ratios
   shell_paras = shell_paras[WHERE(shell_paras NE 0D)]
   vr_ratios = vr_ratios[WHERE(vr_ratios NE 0D)]
      
   
   ;;set up plot for shell parameters
   PLOT, [0,0], [0,0], XTITLE='JD (-2450000)', $
            XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
            YRANGE=[0.0,MAX(shell_paras)+0.25*MAX(shell_paras)],$
            XSTYLE=9, YSTYLE=5, CHARSIZE=0.75, PSYM=3, POS=[0.535, 0.05, 0.954, 0.3], /NOERASE
      
   ;;plot a dashed line at minimum shell phase value (1.5)
   OPLOT, [jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], [1.5, 1.5], LINESTYLE=2
      
   ;;plot shell parameters
   OPLOT, bp_jdates, shell_paras, PSYM=4, SYMSIZE=0.5

      
   ;;add a yaxis to the plot for shell parameter
   AXIS, YAxis=0, YRANGE=[0.0,MAX(shell_paras)+0.25*MAX(shell_paras)], YSTYLE=1, CHARSIZE=0.75, $
         YTITLE='Shell-Parameter'
      
   ;;set up pot for V/R
   PLOT, [0,0], [0,0], XTITLE='JD (-2450000)', $
         XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
         YRANGE=[MIN(vr_ratios)-.1,MAX(vr_ratios)+.3],$
         XSTYLE=9, YSTYLE=5, CHARSIZE=0.75, PSYM=3, POS=[0.535, 0.05, 0.954, 0.3], /NOERASE
   
   ;;oplot a dashed line at 0 for even peak strengths
   OPLOT, [jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], [0.0, 0.0], LINESTYLE=2, COLOR=150
   
   ;;plot V/R
   OPLOT, bp_jdates, vr_ratios, PSYM=1, COLOR=150, SYMSIZE=0.5
   
   ;;add a yaxis to the plot for shell parameter
   AXIS, YAxis=1, YRANGE=[MIN(vr_ratios)-.1,MAX(vr_ratios)+.3], YSTYLE=1, CHARSIZE=0.7, COLOR=150
   XYOUTS, 0.998, 0.175, 'log(V/R)', ALIGNMENT=0.5, ORIENTATION=90, CHARSIZE=0.75, COLOR=150, /NORMAL

   ;;add an axis to the plot displaying fractional years 
   AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=0.75
   XYOUTS, 0.746, 0.32, 'Years', ALIGNMENT=0.5, CHARSIZE=0.75, /NORMAL

ENDIF ELSE IF ~avg_peak_num EQ 0D THEN BEGIN ;;do not add plot if all absorptions
   max_fwhm = MAX(fwhms)
   min_fwhm = MIN(fwhms)
   
   ;;make fwhms an array if there is only one                                                                   
   IF N_ELEMENTS(fwhms) EQ 1 THEN fwhms = [fwhms]
   
   ;; set up plot for fwhms
   PLOT, sp_jdates, fwhms, XTITLE='JD (-2450000)', YTITLE='FWHM (nm)', $
         XRANGE=[jdates[0]-axis_padding, jdates[num_analyzed-1]+axis_padding], $
         YRANGE=[MIN(min_fwhm)-.1,MAX(max_fwhm)+.1], $
         XSTYLE=9, CHARSIZE=0.75, PSYM=4, POS=[0.535, 0.05, 0.954, 0.3], SYMSIZE=0.5, /NOERASE
   
   ;;add an axis to the plot displaying fractional years
   AXIS, XAxis=1, XRANGE=[cal_dates[0]-cal_pad,cal_dates[num_analyzed-1]+cal_pad], XSTYLE=1, CHARSIZE=0.75
   XYOUTS, 0.746, 0.32, 'Years', ALIGNMENT=0.5, CHARSIZE=0.75, /NORMAL
   
ENDIF





;;************************************ ADD TEXT TO THE WINDOW ****************************************;;

XYOUTS, 0.03, 0.004, 'Ver. '+ver_num+', ' + systime(0), CHARSIZE=0.65, /NORMAL
XYOUTS, 0.50, 0.97, star_proper_name, ALIGNMENT=0.5, CHARSIZE=1.8, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.20, 0.97, star_name+'  (V='+star_vmag+')', ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.80, 0.97, star_ra + '  ' + star_dec, ALIGNMENT=0.5, CHARSIZE=1.05, CHARTHICK=1.8, /NORMAL
XYOUTS, 0.12, 0.925, '# of Spectra: '+num_spec, ALIGNMENT=0.5, CHARSIZE=0.7, /NORMAL
XYOUTS, 0.03, 0.28, 'Average Characteristics:', CHARSIZE=1.1, /NORMAL
XYOUTS, 0.03, 0.275, '_____________________', CHARSIZE=1.1, /NORMAL
XYOUTS, 0.03, prnt_str_yloc, prnt_str_title_1, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.17, prnt_str_yloc, prnt_str_1, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.31, 0.28, '[Range]', CHARSIZE=1.1, /NORMAL
XYOUTS, 0.29, prnt_str_yloc, prnt_str_1_open, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.30, prnt_str_yloc, prnt_str_1_mins, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.37, prnt_str_yloc, prnt_str_1_colons, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.39, prnt_str_yloc, prnt_str_1_maxs, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.45, prnt_str_yloc, prnt_str_1_close, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.545, 0.925, '# of Obs. Nights:'+STRCOMPRESS(num_nights), COLOR=220, CHARSIZE=0.7, /NORMAL
XYOUTS, 0.762, 0.925, 'First Observed: '+Mf_str+'/'+Df_str+'/'+Yf_str, COLOR=220, CHARSIZE=0.7, /NORMAL
XYOUTS, 0.762, 0.91, 'Last Observed: '+Ml_str+'/'+Dl_str+'/'+Yl_str, COLOR=220, CHARSIZE=0.7, /NORMAL

XYOUTS, 0.03, 0.315, 'Obs. Disk Trend: '+ trend, CHARSIZE=0.8, /NORMAL
XYOUTS, 0.03, 0.3125, '______________', CHARSIZE=0.8, /NORMAL

;;do not add these labels if all absorptions because there will be no
;;plot
IF ~avg_peak_num EQ 0D THEN BEGIN
   XYOUTS, 0.875, 0.605, 'Single Peak', CHARSIZE =0.7, /NORMAL
   XYOUTS, 0.875, 0.59, 'Blue Peak', COLOR=color_table[0], CHARSIZE =0.7, /NORMAL
   XYOUTS, 0.875, 0.575, 'Red Peak', COLOR=color_table[1], CHARSIZE =0.7, /NORMAL
ENDIF


XYOUTS, 0.5, 0.9625, $
        '_________________________________________________________________________________________________',$
        ALIGNMENT=0.5, CHARSIZE=2.4, CHARTHICK=1.8, /NORMAL

;;Printing the excluded frames takes place in this if statement
XYOUTS, 0.03, 0.015, 'Frames excluded:', CHARSIZE =0.65, /NORMAL
IF frames_excluded[0] EQ -1 THEN BEGIN
   XYOUTS, 0.145, 0.015, 'None', CHARSIZE =0.65, /NORMAL
ENDIF ELSE BEGIN
   j=0.145
   FOR i=0, N_ELEMENTS(frames_excluded)-1 DO BEGIN
      IF i EQ N_ELEMENTS(frames_excluded)-1 THEN BEGIN
         XYOUTS, j, 0.014, frames_excluded[i], CHARSIZE =0.75, /NORMAL
      ENDIF ELSE BEGIN
         XYOUTS, j, 0.014, frames_excluded[i]+',', CHARSIZE =0.75, /NORMAL
      ENDELSE
      j = j + 0.055
   ENDFOR
ENDELSE

;; Signify to the user that the spectra are being shifted                                                   
IF h_align && ~cc_align THEN XYOUTS,  0.12, 0.91, 'H'+alpha+' Shifted', ALIGNMENT=0.5, $
                                      CHARSIZE=0.7, /NORMAL
IF cc_align && ~h_align THEN XYOUTS,  0.12, 0.91, 'CC Shifted', ALIGNMENT=0.5, CHARSIZE=0.7, /NORMAL

;; If some peaks are single and some are double, indicate how many of
;; each
IF num_sing NE 0 && num_doub NE 0 THEN BEGIN
   XYOUTS, 0.33, 0.925, '# of SP Spectra:'+STRCOMPRESS(num_sing), CHARSIZE=0.7, /NORMAL
   XYOUTS, 0.328, 0.91, '# of DP Spectra:'+STRCOMPRESS(num_doub), CHARSIZE=0.7, /NORMAL
   XYOUTS, 0.358, 0.895, '__ Max Emission', CHARSIZE=0.7, COLOR=150, /NORMAL
   XYOUTS, 0.358, 0.88, '__ Min Emission', CHARSIZE=0.7, COLOR=250, /NORMAL
ENDIF ELSE BEGIN
   XYOUTS, 0.358, 0.925, '__ Max Emission', CHARSIZE=0.7, COLOR=150, /NORMAL
   XYOUTS, 0.358, 0.91, '__ Min Emission', CHARSIZE=0.7, COLOR=250, /NORMAL
ENDELSE


;; close star_report.ps/eps file
DEVICE, /CLOSE_FILE
SET_PLOT, 'X'

;;=========================================================================================;;
;;================================ END POSTSCRIPT =========================================;;
;;=========================================================================================;;


;; get working directory
SPAWN, 'pwd', og_loc
og_loc = og_loc+'/'
;;create a directory within the your targets directory to hold all
;;the stars' star reports

CD, targets_dir_name
IF ~FILE_TEST('Star_Reports', /DIRECTORY) THEN SPAWN, 'mkdir '+'Star_Reports'
CD, og_loc

;;Now we will move a copy of that postscript file over to the
;;star's directory and a directory that'll hold all reports
SPAWN, '\cp HD'+star_hd+'_report.ps '+og_loc+targets_dir_name+'/HD'+star_hd ;move to stars directory
SPAWN, '\mv HD'+star_hd+'_report.ps '+og_loc+targets_dir_name+'/Star_Reports' ;move to star_reports direc


jump: ;; jump here if no matches were detected


;; if statement to determine which errors (if any) should be printed
IF no_exp_frames EQ 1 && ISA(star_list) EQ 0 THEN BEGIN
   WINDOW, 0, XSIZE=1300, YSIZE=1300
   XYOUTS, 0.5, 0.55, 'YOU ENTERED: '+ user_entry, ALIGNMENT=0.5, COLOR=250, $
           CHARSIZE=5.0, CHARTHICK=5, /NORMAL
   XYOUTS, 0.5, 0.5, 'NO FRAMES FOUND IN OBSERVATION LOG :(', ALIGNMENT=0.5, COLOR=250, $
           CHARSIZE=5.0, CHARTHICK=5, /NORMAL
ENDIF ELSE IF no_match_frames EQ 1 && ISA(star_list) EQ 0 THEN BEGIN
   WINDOW, 0, XSIZE=1300, YSIZE=1300
   XYOUTS, 0.5, 0.55, 'YOU ENTERED: '+ user_entry, ALIGNMENT=0.5, COLOR=250, $
           CHARSIZE=5.0, CHARTHICK=5, /NORMAL
   XYOUTS, 0.5, 0.5, 'NO SPECTRA FOUND IN YOUR DIRECTORTY :(', ALIGNMENT=0.5, COLOR=250, $
           CHARSIZE=5.0, CHARTHICK=5, /NORMAL
ENDIF ELSE IF no_hdhr_num EQ 1 && ISA(star_list) EQ 0 THEN BEGIN
   WINDOW, 0, XSIZE=1300, YSIZE=1300
   XYOUTS, 0.5, 0.55, 'YOU ENTERED: '+ user_entry, ALIGNMENT=0.5, COLOR=250, $
           CHARSIZE=5.0, CHARTHICK=5, /NORMAL
   XYOUTS, 0.5, 0.5, 'NO HD/HR NUMBER FOUND IN CATALOG :(', ALIGNMENT=0.5, COLOR=250, $
           CHARSIZE=5.0, CHARTHICK=5, /NORMAL
ENDIF


;; compute and report run time
run_time = SYSTIME(1) - start_time
minutes = run_time / 60
runtime_minutes = FLOOR(minutes)
runtime_seconds = FLOOR((minutes - runtime_minutes) * 60)

;;start time is defined near the beginning of hdhr_locator.pro
IF ~mult_stars THEN BEGIN
   PRINT, ''
   PRINT, 'RUN TIME:', STRCOMPRESS(STRING(runtime_minutes))+'m '+$
          STRCOMPRESS(STRING(runtime_seconds), /REMOVE_ALL)+'s'
   PRINT, '============================ END ============================'
ENDIF


END
