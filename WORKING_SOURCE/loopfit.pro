
DATAFILE ='/Users/coletamburri/Desktop/VBIframe1.fits' ;input image file
LOOPFILE ='/Volumes/VBI_External/VBI_LOOPTRACE2.dat' ;filename for output data
IMAGE1 =READFITS(datafile,header)
NSM1 =10 ;lowpass filter
;NSM2 = NSM1+2 ;highpass filter
QMED=0
RMIN =40 ;minimum curvature radius of loop (pixels)
LMIN =50 ;minimum loop length (in pixels)
NSTRUC =2000 ;maximum limit of traced structures used in array dimension
NLOOPMAX =2000 ;maximum number of detected loops
NGAP =0 ;number of pixels in loop below flux threshold (0,...3)
THRESH1 =0.0 ;ratio of image base flux to median flux
THRESH2 =2 ;threshold in detected structure (number of significance ;levels with respect to the median flux in loop profile
TEST =1001 ;option for display of traced structures if TEST < NSTRUC
WID=1
PARA =[NSM1,RMIN,LMIN,NSTRUC,NLOOPMAX,NGAP,THRESH1,THRESH2,WID,QMED]
LOOPTRACING_AUTO4,IMAGE1,IMAGE2,LOOPFILE,PARA,OUTPUT,TEST
READCOL,LOOPFILE,ILOOP,X,Y,Z,S ;read output data
N=MAX(ILOOP)+1 ;maximum loop number
WINDOW,0,XSIZE=800,YSIZE=800
PLOT,[0,0],[0,0],XRANGE=[0,4096],YRANGE=[0,4096],XSTYLE=1,YSTYLE=1
LMIN=0
FOR I=0,N-1 DO BEGIN &IND=WHERE(I eq ILOOP,NS) &IF (MAX(S(IND)) ge LMIN) THEN BEGIN &SPLINE_P,X(IND),Y(IND),XX,YY &OPLOT,XX,YY,THICK=1 &ENDIF &ENDFOR