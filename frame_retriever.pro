;;======================================================================;;
;; Star Report Pipeline. Central Michigan University  
;;
;; frame_retriever.pro
;;
;; This program's main function is to retrieve all the frames
;; for the requested star and store them in a sub-directory of the
;; directory for that star. It is best describes as the completion of
;; a series of tasks that are described as follows, create, if
;; it doesn't already exist, a directory titled "Targets" then
;; move to it. Create, if it doesn't already exist, a
;; sub-directory titled "HD[star's HD number]" then move to
;; it. Then create two other sub-sub-directories titled
;; "Diagnostics" and "Ascii_files". The Diagnostics directory will
;; be used to store various files that contain various information on
;; how the analysis was done and values associated with them. The
;; Ascii_files directory is the directory that holds are the frames
;; for the star the report is being generated on. 
;;
;; Author: W. Glennon Fagan
;;======================================================================;;


PRO frame_retriever
RESOLVE_ALL, /QUIET
COMMON star_block

IF ~mult_stars THEN BEGIN
   PRINT, "Transfering " +STRCOMPRESS(STRING(N_ELEMENTS(available_frames)), /REMOVE_ALL)$
          +" files to Star's directory: "
ENDIF

;;create a directory for targets if it doesnt already exist. then CD
;;into it
IF ~FILE_TEST(targets_dir_name, /DIRECTORY) THEN SPAWN, 'mkdir '+targets_dir_name
CD, targets_dir_name

;;create a directory for the requested star if it doesn't
;;already exist then CD into it
IF ~FILE_TEST('HD'+star_hd, /DIRECTORY) THEN SPAWN, 'mkdir HD'+star_hd
CD, 'HD'+star_hd


;;make directory to stor various diagnostic reports for analysis that
;;will be used later for spectral analysis
IF ~FILE_TEST('Diagnostics', /DIRECTORY) THEN SPAWN, 'mkdir Diagnostics'


;;create a directory that will contain all the Log files for the star
IF ~FILE_TEST('Logs', /DIRECTORY) THEN SPAWN, 'mkdir Logs'

;;Create a directory that will contain all the ascii files for the
;;star then CD into it
IF ~FILE_TEST('Ascii_files', /DIRECTORY) THEN SPAWN, 'mkdir Ascii_files'
CD, 'Ascii_files'


;;Define some variables while in the ascii files directory
SPAWN, 'pwd', file_loc
file_loc = file_loc + '/'
SPAWN, 'ls -1 | wc -l', dir_num_files ;;number of files currently in ascii files directory for star

jump: ;;jump here if overwrite_check is equal to 0
;;CD over to the spectra directory to copy over ascii files                                                 
CD, root_direc+sub_direc
FOR i=0, N_ELEMENTS(available_frames)-1 DO BEGIN
   SPAWN, '\cp *' + available_frames[i] + '* ' + file_loc
ENDFOR


CD, root_direc ;;CD back to the working directory

END
