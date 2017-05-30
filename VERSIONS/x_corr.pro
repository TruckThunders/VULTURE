pro x_corr

;READCOL,"J005758-264314.data",WAV,DUM,FLUX,SIG,ERR,CONT
RDFLOAT,'J005758-264314.data',WAV,DUM,FLUX,SIG,ERR,CONT,/DOUBLE
READCOL,"zem.dat",ZEM

; Transition Wavelengths
LYA = 1215.6
MGII2796 = 2796.353786 
MGII2803 = 2803.530982
CIV1548 = 1548.195
CIV1551 = 1550.770
; +/- 0.5\AA = 8 pix

; Search redshift for positive detection
REDSHIFT = 1.5337

; Search redshift for negative detection 
;REDSHIFT = 1.6

; Create flux arrays for data and then for filter
NORM_FLUX = FLUX/CONT
FILT_FLUX = MAKE_ARRAY(n_elements(FLUX),1,/DOUBLE,VALUE=1)
FILT_FLUX[where(WAV gt (MGII2796*(1.+REDSHIFT))-0.75 and WAV lt (MGII2796*(1.+REDSHIFT))+0.75)] = 0.

FILT_FLUX[where(WAV gt (MGII2803*(1.+REDSHIFT))-0.75 and WAV lt (MGII2803*(1.+REDSHIFT))+0.75)] = 0.

FILT_FLUX[where(WAV gt (CIV1548*(1.+REDSHIFT))-0.75 and WAV lt (CIV1548*(1.+REDSHIFT))+0.75)] = 0.

FILT_FLUX[where(WAV gt (CIV1551*(1.+REDSHIFT))-0.75 and WAV lt (CIV1551*(1.+REDSHIFT))+0.75)] = 0.

; Search blueward of the Lya emission
BLUE_LIMIT = LYA * (1. + ZEM[0])
SEARCH_FLUX = NORM_FLUX[where(WAV gt BLUE_LIMIT)]
SEARCH_WAV = WAV[WHERE(WAV GT BLUE_LIMIT)]
SEARCH_FILT = FILT_FLUX[where(WAV gt BLUE_LIMIT)]

; CREATE LAG ARRAY
LAG = DINDGEN(n_elements(SEARCH_FLUX)) - 0.5*n_elements(SEARCH_FLUX)
;LAG = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
RESULT = C_CORRELATE(SEARCH_FLUX,SEARCH_FILT,LAG,/DOUBLE)




stop
end
