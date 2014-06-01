; Lyra_flaresandgoes3.pro
; ------------------------------------------------------------------------------
; IED  29 Mar 2012, 23 Oct, 17 Dec 2013
; ------------------------------------------------------------------------------
; IDL program to download GOES FITS file from NASA website, LYRA FITS file
;  from P2SC website, event list from NOAA website, and create LYRA flare
;  list for a given day.
; update 26Mar2012: new data sources for GOES and Events
; added 28Mar2012: former independent program LyraGoescurve.pro:
;  IDL program for a daily plot including the GOES curve with LYRA channel 3 
;  and channel 4 proxies: (LYRA - LYRA min)*factor + GOES min
; change 29Mar2012: program LyraFlarelistGoescurve.pro turned into a procedure
;  that can be used everywhere on azrael. External files are downloaded and
;  removed later with updates; *.rst and *.html and *.png files are stored in
;  their respective directories
; change 23Oct2013: remove occultation and LAR gaps with LYTAF
; change 17Dec2013: add calib_lev4a to automatically produce image without gaps
;  in directory /homepage/Level4calibrated/
; ------------------------------------------------------------------------------

pro lyra_flaresandgoes3,DateString

;-------------------------------------------------------------------------------
; Get day in the format YYYYMMDD
;-------------------------------------------------------------------------------

DateYear =strmid(DateString,0,4)
DateMonth=strmid(DateString,4,2)
DateDay  =strmid(DateString,6,2)
if (DateMonth eq '01') then mon='Jan'
if (DateMonth eq '02') then mon='Feb'
if (DateMonth eq '03') then mon='Mar'
if (DateMonth eq '04') then mon='Apr'
if (DateMonth eq '05') then mon='May'
if (DateMonth eq '06') then mon='Jun'
if (DateMonth eq '07') then mon='Jul'
if (DateMonth eq '08') then mon='Aug'
if (DateMonth eq '09') then mon='Sep'
if (DateMonth eq '10') then mon='Oct'
if (DateMonth eq '11') then mon='Nov'
if (DateMonth eq '12') then mon='Dec'

;-------------------------------------------------------------------------------
; Get GOES, LYRA, and EventList files from their various websites 
; after removing older versions (but only for current year). 
;-------------------------------------------------------------------------------

GoesFitsFileName=DateString+'_Gp_xr_1m.txt'
LyraFitsFileName='lyra_'+DateString+'-000000_lev3_std.fits'
EventListFileName=DateString+'events.txt'

if (DateYear eq '2014') then begin
spawn,'rm *_Gp_xr_1m.txt'
spawn,'wget http://www.swpc.noaa.gov/ftpdir/lists/xray/'+GoesFitsFileName
spawn,'rm lyra_*-000000_lev3_std.fits'
spawn,'wget http://proba2.sidc.be/lyra/data/bsd/'+DateYear+'/'+DateMonth+'/'+DateDay+'/'+LyraFitsFileName
spawn,'rm *events.txt'
spawn,'wget http://www.swpc.noaa.gov/ftpdir/indices/events/'+EventListFileName
endif

;----------------------------------------------------------------------------------
; Save data content from GOES, LYRA, and EventsList files to restore or text files.
;----------------------------------------------------------------------------------

ingoesfile='/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/goes'+DateString+'.rst'
inlyrafile='/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/over'+DateString+'.rst'
inevntfile='/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/'+DateString+'events.txt'

res_fs_GoesFits=file_search(GoesFitsFileName)
if (res_fs_GoesFits ne '') then begin
close,1
openr,1,GoesFitsFileName
gtim=dblarr(1440)
gch1=dblarr(1440)
in1=0 & in2=0 & in3=0 & in4=0 & in5=0. & in6=0. & in7=0. & in8=0.
inline=''
for i=1,19 do begin
 readf,1,inline
 endfor
i=0
while(not(eof(1))) do begin
 readf,1,in1,in2,in3,in4,in5,in6,in7,in8
 if (in1 eq fix(DateYear)) then begin
  gtim(i)=float(in6)/60.
  gch1(i)=in8
  i=i+1
  endif
 endwhile 
save,gch1,gtim,filename=ingoesfile
close,1
endif

res_fs_LyraFits=file_search(LyraFitsFileName)
if (res_fs_LyraFits ne '') then begin
fxread,LyraFitsFileName,data,headerp
fxbopen,unit,LyraFitsFileName,'IRRAD LEVEL 3',headerb
length=fxpar(headerb,'NAXIS2')
if (length eq 0) then begin
 tim=[0.,1439./60.]
 ch1=[0.,0.]
 ch2=[0.,0.]
 ch3=[0.,0.]
 ch4=[0.,0.]
 endif else begin
 fxbread,unit,IRRAD_TIME,'TIME'
 fxbread,unit,IRRAD_CHANNEL1,'CHANNEL1'
 fxbread,unit,IRRAD_CHANNEL2,'CHANNEL2'
 fxbread,unit,IRRAD_CHANNEL3,'CHANNEL3'
 fxbread,unit,IRRAD_CHANNEL4,'CHANNEL4'
 fxbread,unit,IRRAD_WARNING,'WARNING'
 fxbclose,unit
 ch1=IRRAD_CHANNEL1
 ch2=IRRAD_CHANNEL2
 ch3=IRRAD_CHANNEL3
 ch4=IRRAD_CHANNEL4
 tim=IRRAD_TIME
 endelse
save,ch1,ch2,ch3,ch4,tim,filename=inlyrafile
endif

if (DateYear eq '2014') then begin
res_fs_EventList=file_search(EventListFileName)
if (res_fs_EventList eq '') then begin
 close,3
 openw,3,EventListFileName
 printf,3,' '
 close,3
 endif
spawn,'cp '+EventListFileName+' '+inevntfile
endif

;-------------------------------------------------------------------------------
; Remove LAR etc, and correct occultation gaps 
; according to former program "gaps_newmask.pro"
;-------------------------------------------------------------------------------

restore,inlyrafile

signvek=intarr(1440)
ctim=fltarr(1440)
cch1=fltarr(1440)
cch2=fltarr(1440)
cch3=fltarr(1440)
cch4=fltarr(1440)
for i=0,n_elements(tim)-1 do begin
ctim(tim(i))=tim(i)
cch1(tim(i))=ch1(i)
cch2(tim(i))=ch2(i)
cch3(tim(i))=ch3(i)
cch4(tim(i))=ch4(i)
signvek(tim(i))=1
endfor

restore,'/home/dammasch/lyra2013over/Jan2013/LYTAF-Template.rst'
spawn,'wget --output-document=/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/'+Datestring+'lytaf.txt "http://proba2.oma.be/lyra/data/lytaf/requestAnnotation.php?submit=download&begin_time='+DateYear+'-'+DateMonth+'-'+DateDay+'T00:00:00Z&end_time='+DateYear+'-'+DateMonth+'-'+DateDay+'T23:59:59Z"'
result=read_ascii('/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/'+Datestring+'lytaf.txt',count=a,template=sTemplate)
list_b=[]
list_e=[]
for i=0,a-1 do begin
if result.field1(i) eq 'UV occ.' then begin
 hou1=float(strmid(result.field2(i),11,2))
 min1=float(strmid(result.field2(i),14,2))
 itv1=long(hou1*60.+min1)
 if (strmid(result.field2(i),8,2) ne DateDay) then itv1=long(0) else list_b=[list_b,itv1-1]
 hou2=float(strmid(result.field4(i),11,2))
 min2=float(strmid(result.field4(i),14,2))
 itv2=long(hou2*60.+min2)
 if (strmid(result.field4(i),8,2) ne DateDay) then itv2=long(1439) else list_e=[list_e,itv2+1]
 signvek(itv1:itv2)=0
 cch1(itv1:itv2)=0.
 cch2(itv1:itv2)=0.
 cch3(itv1:itv2)=0.
 cch4(itv1:itv2)=0.
 endif
if result.field1(i) eq 'Offpoint' then begin
 hou1=float(strmid(result.field2(i),11,2))
 min1=float(strmid(result.field2(i),14,2))
 itv1=long(hou1*60.+min1)
 if (strmid(result.field2(i),8,2) ne DateDay) then itv1=long(0)
 hou2=float(strmid(result.field4(i),11,2))
 min2=float(strmid(result.field4(i),14,2))
 itv2=long(hou2*60.+min2)
 if (strmid(result.field4(i),8,2) ne DateDay) then itv2=long(1439)
 signvek(itv1:itv2)=0
 cch1(itv1:itv2)=0.
 cch2(itv1:itv2)=0.
 cch3(itv1:itv2)=0.
 cch4(itv1:itv2)=0.
 endif
; if result.field1(i) eq 'Recovery' then begin
;  hou1=float(strmid(result.field2(i),11,2))
;  min1=float(strmid(result.field2(i),14,2))
;  itv1=long(hou1*60.+min1)
;  if (strmid(result.field2(i),8,2) ne DateDay) then itv1=long(0)
;  hou2=float(strmid(result.field4(i),11,2))
;  min2=float(strmid(result.field4(i),14,2))
;  itv2=long(hou2*60.+min2)
;  if (strmid(result.field4(i),8,2) ne DateDay) then itv2=long(1439)
;  signvek(itv1:itv2)=0 
;  cch1(itv1:itv2)=0.
;  cch2(itv1:itv2)=0.
;  cch3(itv1:itv2)=0.
;  cch4(itv1:itv2)=0. 
;  endif
if result.field1(i) eq 'Calibration' then begin
 hou1=float(strmid(result.field2(i),11,2))
 min1=float(strmid(result.field2(i),14,2))
 itv1=long(hou1*60.+min1)
 if (strmid(result.field2(i),8,2) ne DateDay) then itv1=long(0)
 hou2=float(strmid(result.field4(i),11,2))
 min2=float(strmid(result.field4(i),14,2))
 itv2=long(hou2*60.+min2)
 if (strmid(result.field4(i),8,2) ne DateDay) then itv2=long(1439)
 signvek(itv1:itv2)=0 
 cch1(itv1:itv2)=0.
 cch2(itv1:itv2)=0.
 cch3(itv1:itv2)=0.
 cch4(itv1:itv2)=0. 
 endif
if result.field1(i) eq 'LAR' then begin
 hou1=float(strmid(result.field2(i),11,2))
 min1=float(strmid(result.field2(i),14,2))
 itv1=long(hou1*60.+min1)
 if (strmid(result.field2(i),8,2) ne DateDay) then itv1=long(0)
 hou2=float(strmid(result.field4(i),11,2))
 min2=float(strmid(result.field4(i),14,2))
 itv2=long(hou2*60.+min2)
 if (strmid(result.field4(i),8,2) ne DateDay) then itv2=long(1439)
 signvek(itv1:itv2)=0 
 cch1(itv1:itv2)=0.
 cch2(itv1:itv2)=0.
 cch3(itv1:itv2)=0.
 cch4(itv1:itv2)=0. 
 endif
endfor
nongap=where(signvek eq 1)

cch10=cch1
cch20=cch2
cch30=cch3
cch40=cch4
if (n_elements(list_e) eq 0) then goto,without_occ
restore,'/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/mask_'+mon+DateYear+'.rst'
mask1=mask(0,*,fix(DateDay-1))
mask2=mask(1,*,fix(DateDay-1))
mask3=mask(2,*,fix(DateDay-1))
mask4=mask(3,*,fix(DateDay-1))
mask1=mask1-mask1(50)
mask2=mask2-mask2(50)
mask3=mask3-mask3(50)
mask4=mask4-mask4(50)
for i=0,n_elements(list_e)-2 do begin
itv1=list_e(i)
itv2=list_e(i+1)-1
mask10=mask1(0:itv2-itv1)
cch10(itv1:itv2)=cch1(itv1:itv2)-mask10
mask20=mask2(0:itv2-itv1)
cch20(itv1:itv2)=cch2(itv1:itv2)-mask20
mask30=mask3(0:itv2-itv1)
cch30(itv1:itv2)=cch3(itv1:itv2)-mask30
mask40=mask4(0:itv2-itv1)
cch40(itv1:itv2)=cch4(itv1:itv2)-mask40
endfor
i=n_elements(list_e)-1
itv1=list_e(i)
itv2=(itv1+98)<1439
if (itv1 le 1439) then begin
mask10=mask1(0:itv2-itv1)
cch10(itv1:itv2)=cch1(itv1:itv2)-mask10
mask20=mask2(0:itv2-itv1)
cch20(itv1:itv2)=cch2(itv1:itv2)-mask20
mask30=mask3(0:itv2-itv1)
cch30(itv1:itv2)=cch3(itv1:itv2)-mask30
mask40=mask4(0:itv2-itv1)
cch40(itv1:itv2)=cch4(itv1:itv2)-mask40
endif
i=0
itv1=(list_e(i)-99)>0
itv2=(list_e(i)-1)>0
if (itv2 ge 0) then begin
mask10=mask1(99-itv2+itv1:99)
cch10(itv1:itv2)=cch1(itv1:itv2)-mask10
mask20=mask2(99-itv2+itv1:99)
cch20(itv1:itv2)=cch2(itv1:itv2)-mask20
mask30=mask3(99-itv2+itv1:99)
cch30(itv1:itv2)=cch3(itv1:itv2)-mask30
mask40=mask4(99-itv2+itv1:99)
cch40(itv1:itv2)=cch4(itv1:itv2)-mask40
endif

without_occ:
tim=ctim(nongap)
ch1=cch10(nongap)
ch2=cch20(nongap)
ch3=cch30(nongap)
ch4=cch40(nongap)
if (n_elements(tim) le 1) then begin
tim=[0.,1439./60.]
ch1=[0.,0.]
ch2=[0.,0.]
ch3=[0.,0.]
ch4=[0.,0.]
endif
save,tim,ch1,ch2,ch3,ch4,filename='/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/gaps'+DateString+'.rst'

;-------------------------------------------------------------------------------
; Continue with old program "lyra_flaresandgoes.pro"
; First plot images for Flare List and produce daily HTML page.
;-------------------------------------------------------------------------------

inlyrafile='/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/gaps'+DateString+'.rst'
res_fs=file_search(inlyrafile)
if (res_fs eq '') then begin
oldtim=[0.,1439./60.]
oldch1=[0.,0.]
oldch2=[0.,0.]
oldch3=[0.,0.]
oldch4=[0.,0.]
endif else begin
restore,inlyrafile
oldtim=tim/60.
oldch1=(ch1>0.)*1000.
oldch2=(ch2>0.)*1000.
oldch3=(ch3>0.)*1000.
oldch4=(ch4>0.)*1000. 
endelse

datestr2=strtrim(string(long(DateString)+long(1)),1)
if (datestr2 eq '20100132') then datestr2='20100201'
if (datestr2 eq '20100229') then datestr2='20100301'
if (datestr2 eq '20100332') then datestr2='20100401'
if (datestr2 eq '20100431') then datestr2='20100501'
if (datestr2 eq '20100532') then datestr2='20100601'
if (datestr2 eq '20100631') then datestr2='20100701'
if (datestr2 eq '20100732') then datestr2='20100801'
if (datestr2 eq '20100832') then datestr2='20100901'
if (datestr2 eq '20100931') then datestr2='20101001'
if (datestr2 eq '20101032') then datestr2='20101101'
if (datestr2 eq '20101131') then datestr2='20101201'
if (datestr2 eq '20101232') then datestr2='20110101'
if (datestr2 eq '20110132') then datestr2='20110201'
if (datestr2 eq '20110229') then datestr2='20110301'
if (datestr2 eq '20110332') then datestr2='20110401'
if (datestr2 eq '20110431') then datestr2='20110501'
if (datestr2 eq '20110532') then datestr2='20110601'
if (datestr2 eq '20110631') then datestr2='20110701'
if (datestr2 eq '20110732') then datestr2='20110801'
if (datestr2 eq '20110832') then datestr2='20110901'
if (datestr2 eq '20110931') then datestr2='20111001'
if (datestr2 eq '20111032') then datestr2='20111101'
if (datestr2 eq '20111131') then datestr2='20111201'
if (datestr2 eq '20111232') then datestr2='20120101'
if (datestr2 eq '20120132') then datestr2='20120201'
if (datestr2 eq '20120230') then datestr2='20120301'
if (datestr2 eq '20120332') then datestr2='20120401'
if (datestr2 eq '20120431') then datestr2='20120501'
if (datestr2 eq '20120532') then datestr2='20120601'
if (datestr2 eq '20120631') then datestr2='20120701'
if (datestr2 eq '20120732') then datestr2='20120801'
if (datestr2 eq '20120832') then datestr2='20120901'
if (datestr2 eq '20120931') then datestr2='20121001'
if (datestr2 eq '20121032') then datestr2='20121101'
if (datestr2 eq '20121131') then datestr2='20121201'
if (datestr2 eq '20121232') then datestr2='20130101'
if (datestr2 eq '20130132') then datestr2='20130201'
if (datestr2 eq '20130229') then datestr2='20130301'
if (datestr2 eq '20130332') then datestr2='20130401'
if (datestr2 eq '20130431') then datestr2='20130501'
if (datestr2 eq '20130532') then datestr2='20130601'
if (datestr2 eq '20130631') then datestr2='20130701'
if (datestr2 eq '20130732') then datestr2='20130801'
if (datestr2 eq '20130832') then datestr2='20130901'
if (datestr2 eq '20130931') then datestr2='20131001'
if (datestr2 eq '20131032') then datestr2='20131101'
if (datestr2 eq '20131131') then datestr2='20131201'
if (datestr2 eq '20131232') then datestr2='20140101'
yea2=strmid(datestr2,0,4) & mn2=strmid(datestr2,4,2) & day2=strmid(datestr2,6,2)
if (mn2 eq '01') then mon2='Jan'
if (mn2 eq '02') then mon2='Feb'
if (mn2 eq '03') then mon2='Mar'
if (mn2 eq '04') then mon2='Apr'
if (mn2 eq '05') then mon2='May'
if (mn2 eq '06') then mon2='Jun'
if (mn2 eq '07') then mon2='Jul'
if (mn2 eq '08') then mon2='Aug'
if (mn2 eq '09') then mon2='Sep'
if (mn2 eq '10') then mon2='Oct'
if (mn2 eq '11') then mon2='Nov'
if (mn2 eq '12') then mon2='Dec'
inrstfile2='/home/dammasch/lyra'+yea2+'over/'+mon2+yea2+'/gaps'+datestr2+'.rst'
res_fs2=file_search(inrstfile2)
if (res_fs2 eq '') then begin
tim=[oldtim,24.,24.+1439./60.]
ch1=[oldch1,0.,0.]
ch2=[oldch2,0.,0.]
ch3=[oldch3,0.,0.]
ch4=[oldch4,0.,0.]
endif else begin   
restore,inrstfile2
tim=[oldtim,tim/60.+24.]
ch1=[oldch1,(ch1>0.)*1000.]
ch2=[oldch2,(ch2>0.)*1000.]
ch3=[oldch3,(ch3>0.)*1000.]
ch4=[oldch4,(ch4>0.)*1000.]
endelse

res_fs_goes=file_search(ingoesfile)
if (res_fs_goes eq '') then begin
oldgtim=[0.,1439./60.]
oldgch1=[0.,0.]
endif else begin
restore,ingoesfile
oldgtim=gtim/60.
oldgch1=gch1*1000.
endelse
ingoesfile2='/home/dammasch/lyra'+yea2+'over/'+mon2+yea2+'/goes'+datestr2+'.rst'
res_fs_goes2=file_search(ingoesfile2)
if (res_fs_goes2 eq '') then begin
gtim=[oldgtim,24.,24.+1439./60.]
gch1=[oldgch1,0.,0.]
endif else begin
restore,ingoesfile2
gtim=[oldgtim,gtim/60.+24.]
gch1=[oldgch1,gch1*1000.]
endelse

close,1 & close,2
openr,1,inevntfile
openw,2,'/home/dammasch/homepage/flares/flare'+DateString+'.html'
printf,2,'<HTML>'
printf,2,'<HEAD>'
printf,2,'<TITLE>'+DateDay+' '+mon+' '+DateYear+' Flare List</TITLE>'
printf,2,'</HEAD>'
printf,2,'<BODY BGCOLOR="#ffffff" TEXT="#000000">'
printf,2,'<H3>'+DateDay+' '+mon+' '+DateYear+' Flare List</H3>'
printf,2,'<IMG SRC="flare'+DateString+'.png">'
printf,2,'<HR>'
printf,2,'<PRE>'
printf,2,'event  begin   max     end     class  region'
printf,2,'-----  -----   -----   -----   -----  ------'
set_plot,'z'
tvlct,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
erase
device,set_resolution=[650,400],set_pixel_depth=8,decomposed=1
device,z_buffering=0
!p.background=5
!p.color=4
!xmin=0 & !xmax=24
!p.charsize=0.75

;if (min(ch4) lt 0.8) then lo04=median(ch4)-0.06 else lo04=min(ch4)
;hi04=max(ch4)<10.

if (min(oldch4) lt 0.8) then lo04=median(oldch4)-0.06 else lo04=min(oldch4)
hi04=max(oldch4)<10.

if (hi04 le lo04) then hi04=lo04+1.
df04=hi04-lo04
!xmin=0 & !xmax=24 & !ymin=lo04 & !ymax=hi04+0.2*df04
ymid1=hi04+0.06*df04
; ------------------------------------------------------------------------------
; FIRST PLOT ROUTINE (per day)
; ------------------------------------------------------------------------------
plot,indgen(10)
plot,tim,ch4,psym=3,xstyle=1,ystyle=1,/nodata, $
title='Channel 2-4   Zirconium-Filter   6-20nm + Xray   (1 minute averages)', $
ytitle='Calibrated Signal / (mW m!U-2!N)', $
xtitle='Time / h UTC, '+DateDay+' '+mon+' '+DateYear
oplot,tim,ch4,psym=2,symsize=0.2,color=3
if (total(ch4) lt 2) then xyouts,2,0.5*lo04+0.5*hi04,'(NO DATA AVAILABLE)',size=1.5
inline=''
while(not(eof(1))) do begin 
readf,1,inline
if (strmid(inline,43,3) eq 'XRA') then begin
 evtfl=strmid(inline,0,4)
 evtfl=strtrim(evtfl,1)
 if (strlen(evtfl) eq 3) then evtfl='0'+evtfl
 if (strlen(evtfl) eq 2) then evtfl='00'+evtfl
 if (strlen(evtfl) eq 1) then evtfl='000'+evtfl
 begfl=strmid(inline,11,2)+':'+strmid(inline,13,2)
 begflt=float(strmid(inline,11,2))+float(strmid(inline,13,2))/60.
 maxfl=strmid(inline,18,2)+':'+strmid(inline,20,2)
 maxflt=float(strmid(inline,18,2))+float(strmid(inline,20,2))/60.
 endfl=strmid(inline,28,2)+':'+strmid(inline,30,2)
 endflt=float(strmid(inline,28,2))+float(strmid(inline,30,2))/60.
 if (endflt lt begflt) then begin
  if (maxflt lt begflt) then maxflt=maxflt+24.
  endflt=endflt+24.
  endif
 clafl=strmid(inline,58,4)
 regfl=strmid(inline,76,4)
 lotm=begflt & hitm=endflt
 ind0=where((tim ge lotm) and (tim le hitm))
 if (ymid1 eq hi04+0.06*df04) then begin
  ymid1=hi04+0.13*df04
  ymid2=hi04+0.10*df04
  endif else begin
  ymid1=hi04+0.06*df04
  ymid2=hi04+0.03*df04
  endelse  
 xyouts,begflt,ymid1,evtfl
 xyouts,begflt,ymid2,clafl
 if (n_elements(ind0) eq 1) then begin
  printf,2,evtfl+'   '+begfl+'   '+maxfl+'   '+endfl+'   '+clafl+'   '+regfl
  endif else begin
printf,2,'<A HREF="flare'+DateString+'_'+evtfl+'.png">'+evtfl+'</A>   '+begfl+'   '+maxfl+'   '+endfl+'   '+clafl+'   '+regfl
  endelse
 endif
endwhile
printf,2,'</PRE>'
printf,2,'</BODY>'
printf,2,'</HTML>'
close,1 & close,2
array=tvrd()
write_png,'/home/dammasch/homepage/flares/flare'+DateString+'.png',array, $
[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]

openr,1,inevntfile
set_plot,'z'
erase
device,set_resolution=[550,800],set_pixel_depth=8,decomposed=1
device,z_buffering=0
!p.background=5
!p.color=4
inline=''
while(not(eof(1))) do begin 
readf,1,inline
if (strmid(inline,43,3) eq 'XRA') then begin
 evtfl=strmid(inline,0,4)
 evtfl=strtrim(evtfl,1)
 if (strlen(evtfl) eq 3) then evtfl='0'+evtfl
 if (strlen(evtfl) eq 2) then evtfl='00'+evtfl
 if (strlen(evtfl) eq 1) then evtfl='000'+evtfl
 begfl=strmid(inline,11,2)+':'+strmid(inline,13,2)
 begflt=float(strmid(inline,11,2))+float(strmid(inline,13,2))/60.
 maxfl=strmid(inline,18,2)+':'+strmid(inline,20,2)
 maxflt=float(strmid(inline,18,2))+float(strmid(inline,20,2))/60.
 endfl=strmid(inline,28,2)+':'+strmid(inline,30,2)
 endflt=float(strmid(inline,28,2))+float(strmid(inline,30,2))/60.
 if (endflt lt begflt) then begin
  if (maxflt lt begflt) then maxflt=maxflt+24.
  endflt=endflt+24.
  endif 
 clafl=strmid(inline,58,4)
 regfl=strmid(inline,76,4)
 !xmin=maxflt-1. & !xmax=maxflt+2.
 ind0=where((tim ge !xmin) and (tim le !xmax))
 if (n_elements(ind0) eq 1) then goto,weiter
 lo1=median(ch1(ind0))-0.04 & hi1=median(ch1(ind0))+0.06 & df1=hi1-lo1
 lo2=median(ch2(ind0))-4. & hi2=median(ch2(ind0))+6. & df2=hi2-lo2
 if (min(ch3(ind0)) lt 2.2) then lo3=median(oldch3)-0.18 else lo3=min(ch3(ind0))
 hi3=max(ch3(ind0))<15.  
 if (hi3 le lo3) then hi3=lo3+1. & df3=hi3-lo3
 if (min(ch4(ind0)) lt 0.66) then lo4=median(oldch4)-0.06 else lo4=min(ch4(ind0))
 hi4=max(ch4(ind0))<10.  
 if (hi4 le lo4) then hi4=lo4+1. & df4=hi4-lo4
 ind5=where((gtim ge !xmin) and (gtim le !xmax))
 lo5=min(gch1(ind5))>1.e-5     & hi5=max(gch1(ind5))<1.     & df5=hi5-lo5
 gentickv3=[ 2., 3., 4., 5.,6.,7.,8.,9.,10.,15.]
 gentickv4=[0.6,0.7,0.8,0.9,1.,2.,3.,4., 5., 6., 7., 8., 9.,10.]
 ind3=where((gentickv3 ge lo3) and (gentickv3 le hi3))
 ind4=where((gentickv4 ge lo4) and (gentickv4 le hi4))
 if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
  tickv3=gentickv3(ind3)
  tickv4=gentickv4(ind4)
  ticks3=n_elements(ind3)
  ticks4=n_elements(ind4)
  endif
; ------------------------------------------------------------------------------
; SECOND PLOT ROUTINE (per flare)
;-------------------------------------------------------------------------------
 plot,indgen(10)

 !ymin=lo1 & !ymax=hi1
 plot,tim,ch1,psym=3,xstyle=1,ystyle=1, $
  /nodata,/normal,position=[0.08,0.80,0.98,0.98],xcharsize=0.000001
 oplot,tim,ch1,psym=2,symsize=0.4,color=1
 plots,[maxflt,maxflt],[lo1+0.9*df1,hi1],thick=2
 xyouts,begflt>!xmin,lo1+0.85*df1,evtfl
 xyouts,begflt>!xmin,lo1+0.80*df1,clafl
 xyouts,0.14,0.96,'Channel 2-1',/normal,color=1
 xyouts,0.44,0.96,'Lyman-alpha Filter',/normal,color=1
 xyouts,0.74,0.96,'120-123nm',/normal,color=1

 !ymin=lo2 & !ymax=hi2
 plot,tim,ch2,psym=3,xstyle=1,ystyle=1, $
  /nodata,/normal,/noerase,position=[0.08,0.62,0.98,0.80],xcharsize=0.000001
 oplot,tim,ch2,psym=2,symsize=0.4,color=0
 xyouts,0.14,0.78,'Channel 2-2',/normal,color=0
 xyouts,0.44,0.78,'Herzberg Filter',/normal,color=0
 xyouts,0.74,0.78,'190-222nm',/normal,color=0

 !ymin=lo3 & !ymax=hi3+0.1*df3
 if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
  plot_io,tim,ch3,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
  /nodata,/normal,/noerase,position=[0.08,0.44,0.98,0.62],xcharsize=0.000001,yticks=ticks3-1,ytickv=tickv3
  endif else begin
  plot,tim,ch3,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
  /nodata,/normal,/noerase,position=[0.08,0.44,0.98,0.62],xcharsize=0.000001
  endelse
 oplot,tim,ch3,psym=2,symsize=0.4,color=2
 xyouts,0.14,0.60,'Channel 2-3',/normal,color=2
 xyouts,0.44,0.60,'Aluminium Filter',/normal,color=2
 xyouts,0.74,0.60,'17-80nm + Xray',/normal,color=2  

 !ymin=lo4 & !ymax=hi4+0.1*df4
 if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
  plot_io,tim,ch4,psym=3,xstyle=1,ystyle=1, $
  /nodata,/normal,/noerase,position=[0.08,0.26,0.98,0.44],xcharsize=0.000001,yticks=ticks4-1,ytickv=tickv4
  endif else begin
  plot,tim,ch4,psym=3,xstyle=1,ystyle=1, $
  /nodata,/normal,/noerase,position=[0.08,0.26,0.98,0.44],xcharsize=0.000001
  endelse
 oplot,tim,ch4,psym=2,symsize=0.4,color=3
 xyouts,0.14,0.42,'Channel 2-4',/normal,color=3
 xyouts,0.44,0.42,'Zirconium Filter',/normal,color=3
 xyouts,0.74,0.42,'6-20nm + Xray',/normal,color=3
 if (total(ch4) lt 2) then xyouts,6,0.5*lo4+0.5*hi4,'(NO DATA AVAILABLE)',size=1

 !ymin=lo5 & !ymax=hi5+0.1*df5
 if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
  plot_io,gtim,gch1,psym=3,xstyle=1,ystyle=1, $
  xtitle='Time / h UTC, '+DateDay+' '+mon+' '+DateYear, $
  /nodata,/normal,/noerase,position=[0.08,0.08,0.98,0.26]
  endif else begin
  plot,gtim,gch1,psym=3,xstyle=1,ystyle=1, $
  xtitle='Time / h UTC, '+DateDay+' '+mon+' '+DateYear, $
  /nodata,/normal,/noerase,position=[0.08,0.08,0.98,0.26]
  endelse
 oplot,gtim,gch1,psym=2,symsize=0.4,color=1
 xyouts,0.14,0.24,'GOES',/normal,color=1
 xyouts,0.44,0.24,'Xray Sensor',/normal,color=1
 xyouts,0.74,0.24,'0.1-0.8nm',/normal,color=1 
 xyouts,0.75,0.01,'(1 minute averages)',/normal,color=4   

 !xmin=0 & !xmax=0 & !ymax=0 & !ymin=0  
 array=tvrd()
 write_png,'/home/dammasch/homepage/flares/flare'+DateString+'_'+evtfl+'.png',array, $
  [0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
 weiter:
 endif
endwhile
close,1

; -------------------------------------------------------------------
; Execute former program "LyraGoescurve.pro"
; Plot images for GOES vs. LYRA Proxies.
; -------------------------------------------------------------------

xsize=650
ysize=500
set_plot,'Z'
tvlct,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
erase
device,set_resolution=[xsize,ysize],set_pixel_depth=8,decomposed=1
device,z_buffering=0
!p.background=5
!p.color=4
!xmin=0 & !xmax=24 

res_fs=file_search(inlyrafile)
if (res_fs eq '') then begin
tim=[0.,1439./60.]
ch1=[0.,0.]
ch2=[0.,0.]
ch3=[0.,0.]
ch4=[0.,0.]
endif else begin
restore,inlyrafile
tim=tim/60.
endelse

res_fs_goes=file_search(ingoesfile)
if (res_fs_goes eq '') then begin
gtim=[0.,1439./60.]
gch1=[0.,0.]
endif else begin
restore,ingoesfile
gtim=gtim/60.
endelse

close,1
openr,1,'/home/dammasch/lyra'+DateYear+'over/'+mon+DateYear+'/MinGoesLyra.txt'
in_day=0 & in_a=0.0 & in_b=0.0 & in_c=0.0
v_day=intarr(31) & v_a=fltarr(31) & v_b=fltarr(31) & v_c=fltarr(31)
for i=0,30 do begin
readf,1,in_day,in_a,in_b,in_c
v_day(i)=in_day & v_a(i)=in_a & v_b(i)=in_b & v_c(i)=in_c
if (fix(DateDay) eq in_day) then begin a=in_a & b=in_b & c=in_c & endif
endfor
close,1

; ------------------------------------------------------------------
; try with a,b,c = simply (almost) daily minimum of significant data
; (this is only possible outside occultation season)
; ------------------------------------------------------------------

if (DateYear ne '2010') then begin
!ymin=0.0000001 & !ymax=0.001
ch30=ch3(where(ch3 gt 0.00215))
ch40=ch4(where(ch4 gt 0.00065))

; limit changed from 0.0022 to 0.00215 after one month of very low values. IED
; limit changed from 0.0008 to 0.00065 after one month of very low values. IED

gch10=gch1(where(gch1 gt 1.0e-07))
ch3s=ch30(sort(ch30))
ch4s=ch40(sort(ch40))
gch1s=gch10(sort(gch10))
if (n_elements(ch3s) ge 11) then a=ch3s(10) else a=0.00215
if (n_elements(ch4s) ge 11) then b=ch4s(10) else b=0.00065

;-------------------------------------------------------------------------
; a changed from 0.0022 to 0.00215 after one month of very low values. IED
; b changed from 0.0008 to 0.00065 after one month of very low values. IED
;-------------------------------------------------------------------------

if (n_elements(gch1s) ge 11) then c=gch1s(10) else c=1.0e-07

;--------------------------------------------------------------------
; some exceptions:
;--------------------------------------------------------------------

if (DateString eq '20131120') then begin a=0.0024987255 & b=0.0013790854 & end
if (DateString eq '20131204') then begin a=0.0024833879 & b=0.0013432010 & end
if (DateString eq '20131218') then begin a=0.0026069226 & b=0.0015260200 & end
if (DateString eq '20140115') then begin a=0.0023985442 & b=0.0011855319 & end
if (DateString eq '20140129') then begin a=0.0024677196 & b=0.0011790076 & end

plot,indgen(10)
plot_io,tim,0.015*(ch3-a)+c,psym=3,xstyle=1,ystyle=1,ytitle='GOES Irradiance / (W m!U-2!N)', $
 title='GOES 0.1-0.8nm (red), LYRA Al (blue) & Zr (grey) proxy', $
 xtitle='Time / h UTC, '+DateDay+' '+mon+' '+DateYear, $
 /nodata,/normal,position=[0.08,0.11,0.96,0.95],ycharsize=0.8
oplot,tim,0.015*(ch3-a)+c,psym=3,color=2
oplot,tim,0.018*(ch4-b)+c,psym=3,color=3
oplot,gtim,gch1,psym=3,color=1
plots,[0,24],[0.000001,0.000001]
plots,[0,24],[0.00001,0.00001]
plots,[0,24],[0.0001,0.0001]
flarelabels=['Y','X','M','C','B','A']
for i=1,4 do xyouts,24.3,3.*10.^(-3-i),flarelabels(i),/data
xyouts,0.75,0.01,'ROB/SIDC, Brussels, Belgium',color=4,/normal,charsize=0.8
endif

if (DateYear eq '2010') then begin
!ymin=0.00000001 & !ymax=0.0001
ch30=ch3(where(ch3 gt 0.0020))
ch40=ch4(where(ch4 gt 0.0006))
gch10=gch1(where(gch1 gt 0.5e-07))
ch3s=ch30(sort(ch30))
ch4s=ch40(sort(ch40))
gch1s=gch10(sort(gch10))
if (n_elements(ch3s) ge 11) then a=ch3s(10) else a=0.0020
if (n_elements(ch4s) ge 11) then b=ch4s(10) else b=0.0006
if (n_elements(gch1s) ge 11) then c=gch1s(10) else c=0.5e-07

;--------------------------------------------------------------------
; some exceptions:
;--------------------------------------------------------------------

if (DateString eq '20101015') then begin a=0.0022592962 & b=0.00075498296 & end
if (DateString eq '20101020') then begin a=0.0023230347 & b=0.00084713358 & end
if (DateString eq '20101021') then begin a=0.0022862532 & b=0.00081330907 & end

plot,indgen(10)
plot_io,tim,0.015*(ch3-a)+c,psym=3,xstyle=1,ystyle=1,ytitle='GOES Irradiance / (W m!U-2!N)', $
 title='GOES 0.1-0.8nm (red), LYRA Al (blue) & Zr (grey) proxy', $
 xtitle='Time / h UTC, '+DateDay+' '+mon+' '+DateYear, $
 /nodata,/normal,position=[0.08,0.11,0.96,0.95],ycharsize=0.8
oplot,tim,0.015*(ch3-a)+c,psym=3,color=2
oplot,tim,0.018*(ch4-b)+c,psym=3,color=3
oplot,gtim,gch1,psym=3,color=1
plots,[0,24],[0.0000001,0.0000001]
plots,[0,24],[0.000001,0.000001]
plots,[0,24],[0.00001,0.00001]
flarelabels=['Y','X','M','C','B','A']
for i=2,5 do xyouts,24.3,3.*10.^(-3-i),flarelabels(i),/data
xyouts,0.75,0.01,'ROB/SIDC, Brussels, Belgium',color=4,/normal,charsize=0.8
endif

!xmin=0 & !xmax=0 & !ymax=0 & !ymin=0
array=tvrd()
outpathfile='/home/dammasch/homepage/GoesLyra'+DateString+'.png'
write_png,outpathfile,array,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]

; -------------------------------------------------------------------------------------
; former program calib_lev4a.pro
; -------------------------------------------------------------------------------------

date=DateString
yea=strmid(date,0,4)
mon=strmid(date,4,2)
day=strmid(date,6,2)
if (mon eq '10') then mon1='Oct'
if (mon eq '11') then mon1='Nov'
if (mon eq '12') then mon1='Dec'
if (mon eq '01') then mon1='Jan'
if (mon eq '02') then mon1='Feb'

outpathfile='/home/dammasch/homepage/Level4calibrated/LyraL4C'+date+'.png'

; ------------------------------------------------------------------------------
; prepare plotting
; ------------------------------------------------------------------------------

xsize=650
ysize=800
set_plot,'Z'
tvlct,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
erase
device,set_resolution=[xsize,ysize], set_pixel_depth=24, decomposed=1
device, z_buffering=0
!p.background=5
!p.color=4
!xmin=0 & !xmax=24

; ------------------------------------------------------------------------------
; get data from formerly saved "gaps" file
; ------------------------------------------------------------------------------

restore,'/home/dammasch/lyra'+yea+'over/'+mon1+yea+'/gaps'+date+'.rst'
tim=tim/60.
ch1=ch1*1000.
ch2=ch2*1000.
ch3=ch3*1000.
ch4=ch4*1000.

; ------------------------------------------------------------------------------
; define upper and lower plotting limits for all four channels
; and decide about linear or logarithmic scaling
; ------------------------------------------------------------------------------

lo1=6.2  & hi1=6.4
lo2=688.   & hi2=708.
lo3=(0.95*min(ch3))>2.0
hi3=(1.1*max(ch3))<15.
if (hi3 le lo3) then hi3=lo3+1.
lo4=(0.95*min(ch4))>0.6
hi4=(1.1*max(ch4))<10.
if (hi4 le lo4) then hi4=lo4+1.
gentickv3=[ 2., 3., 4., 5.,6.,7.,8.,9.,10.,15.]
gentickv4=[0.6,0.7,0.8,0.9,1.,2.,3.,4., 5., 6., 7., 8., 9.,10.]
ind3=where((gentickv3 ge lo3) and (gentickv3 le hi3))
ind4=where((gentickv4 ge lo4) and (gentickv4 le hi4))
if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
tickv3=gentickv3(ind3)
tickv4=gentickv4(ind4)
ticks3=n_elements(ind3)
ticks4=n_elements(ind4)
endif

; ------------------------------------------------------------------------------
; plot PNG image
; ------------------------------------------------------------------------------

plot,indgen(10)
!ymin=lo1 & !ymax=hi1
plot,tim,ch1,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
/nodata,/normal,position=[0.08,0.76,0.98,0.99],xcharsize=0.000001,ycharsize=0.8
oplot,tim,ch1,psym=3,color=1
xyouts,0.14,0.97,'Channel 2-1',/normal,color=1
xyouts,0.44,0.97,'Lyman-alpha Filter',/normal,color=1
xyouts,0.74,0.97,'120-123nm',/normal,color=1

!ymin=lo2 & !ymax=hi2
plot,tim,ch2,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
/nodata,/normal,/noerase,position=[0.08,0.53,0.98,0.76],xcharsize=0.000001,ycharsize=0.8
oplot,tim,ch2,psym=3,color=0
xyouts,0.14,0.74,'Channel 2-2',/normal,color=0
xyouts,0.44,0.74,'Herzberg Filter',/normal,color=0
xyouts,0.74,0.74,'190-222nm',/normal,color=0

!ymin=lo3 & !ymax=hi3
if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim,ch3,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
/nodata,/normal,/noerase,position=[0.08,0.30,0.98,0.53],xcharsize=0.000001,ycharsize=0.8,yticks=ticks3-1,ytickv=tickv3
endif else begin
plot,tim,ch3,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
/nodata,/normal,/noerase,position=[0.08,0.30,0.98,0.53],xcharsize=0.000001,ycharsize=0.8
endelse
oplot,tim,ch3,psym=3,color=2
xyouts,0.14,0.51,'Channel 2-3',/normal,color=2
xyouts,0.44,0.51,'Aluminium Filter',/normal,color=2
xyouts,0.74,0.51,'17-80nm + Xray',/normal,color=2

!ymin=lo4 & !ymax=hi4
if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim,ch4,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
xtitle='Time / h UTC, '+day+' '+mon1+' '+yea,/nodata, $
/normal,/noerase,position=[0.08,0.07,0.98,0.30],yticks=ticks4-1,ytickv=tickv4,ycharsize=0.8
endif else begin
plot,tim,ch4,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
xtitle='Time / h UTC, '+day+' '+mon1+' '+yea,/nodata, $
/normal,/noerase,position=[0.08,0.07,0.98,0.30],ycharsize=0.8
endelse
oplot,tim,ch4,psym=3,color=3
xyouts,0.14,0.28,'Channel 2-4',/normal,color=3
xyouts,0.44,0.28,'Zirconium Filter',/normal,color=3
xyouts,0.74,0.28,'6-20nm + Xray',/normal,color=3
xyouts,0.75,0.01,'ROB/SIDC, Brussels, Belgium',/normal,color=4,charsize=0.8
if (total(ch4) lt 0.1) then xyouts,6,0.5*lo4+0.5*hi4,'(NO DATA AVAILABLE)',size=2.5

!xmin=0 & !xmax=0 & !ymax=0 & !ymin=0

; ------------------------------------------------------------------------------
; create PNG image
; ------------------------------------------------------------------------------

array=tvrd()
write_png,outpathfile,array,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
outpathfile2='/home/dammasch/homepage/Level4calibrated/LyraL4Clatest.png'
write_png,outpathfile2,array,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]

; -------------------------------------------------------------------------------------
; former program calib_lev4b.pro
; -------------------------------------------------------------------------------------

xsize=650
ysize=500
set_plot,'Z'
tvlct,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
erase
device,set_resolution=[xsize,ysize], set_pixel_depth=24, decomposed=1
device, z_buffering=0
!p.background=5
!p.color=4
!xmin=0 & !xmax=24

; ------------------------------------------------------------------------------
; read first lev3 FITS header, get data from formerly saved "gaps" file
; ------------------------------------------------------------------------------

date=DateString

date1=long(date)-2
if (date1 eq 20100200) then date1=20100131
if (date1 eq 20100199) then date1=20100130
if (date1 eq 20101000) then date1=20100930
if (date1 eq 20100999) then date1=20100929
if (date1 eq 20101100) then date1=20101031
if (date1 eq 20101099) then date1=20101030
if (date1 eq 20101200) then date1=20101130
if (date1 eq 20101199) then date1=20101129
if (date1 eq 20110100) then date1=20101231
if (date1 eq 20110099) then date1=20101230
if (date1 eq 20110200) then date1=20110131
if (date1 eq 20110199) then date1=20110130
if (date1 eq 20131000) then date1=20130930
if (date1 eq 20130999) then date1=20130929
if (date1 eq 20131100) then date1=20131031
if (date1 eq 20131099) then date1=20131030
if (date1 eq 20131200) then date1=20131130
if (date1 eq 20131199) then date1=20131129
if (date1 eq 20140100) then date1=20131231
if (date1 eq 20140099) then date1=20131230
if (date1 eq 20140200) then date1=20140131
if (date1 eq 20140199) then date1=20140130
date1=strmid(string(date1),4,8)
yea1=strmid(date1,0,4)
mon1=strmid(date1,4,2)
day1=strmid(date1,6,2)
if (mon1 eq '09') then mn1='Sep'
if (mon1 eq '10') then mn1='Oct'
if (mon1 eq '11') then mn1='Nov'
if (mon1 eq '12') then mn1='Dec'
if (mon1 eq '01') then mn1='Jan'
if (mon1 eq '02') then mn1='Feb'
restore,'/home/dammasch/lyra'+yea1+'over/'+mn1+yea1+'/gaps'+date1+'.rst'
tim1=tim/60.
ch11=ch1*1000.
ch21=ch2*1000.
ch31=ch3*1000.
ch41=ch4*1000.

; ------------------------------------------------------------------------------
; read second lev3 FITS header, get data from formerly saved "gaps" file
; ------------------------------------------------------------------------------

date2=long(date)-1
if (date2 eq 20100200) then date2=20100131
if (date2 eq 20101000) then date2=20100930
if (date2 eq 20101100) then date2=20101031
if (date2 eq 20101200) then date2=20101130
if (date2 eq 20110100) then date2=20101231
if (date2 eq 20110200) then date2=20110131
if (date2 eq 20131000) then date2=20130930
if (date2 eq 20131100) then date2=20131031
if (date2 eq 20131200) then date2=20131130
if (date2 eq 20140100) then date2=20131231
if (date2 eq 20140200) then date2=20140131
date2=strmid(string(date2),4,8)
yea2=strmid(date2,0,4)
mon2=strmid(date2,4,2)
day2=strmid(date2,6,2)
if (mon2 eq '09') then mn2='Sep'
if (mon2 eq '10') then mn2='Oct'
if (mon2 eq '11') then mn2='Nov'
if (mon2 eq '12') then mn2='Dec'
if (mon2 eq '01') then mn2='Jan'
if (mon2 eq '02') then mn2='Feb'
restore,'/home/dammasch/lyra'+yea2+'over/'+mn2+yea2+'/gaps'+date2+'.rst'
tim2=tim/60.
ch12=ch1*1000.
ch22=ch2*1000.
ch32=ch3*1000.
ch42=ch4*1000.

; ------------------------------------------------------------------------------
; read third lev3 FITS header, get data from formerly saved "gaps" file
; ------------------------------------------------------------------------------

date3=date
yea3=strmid(date3,0,4)
mon3=strmid(date3,4,2)
day3=strmid(date3,6,2)
if (mon3 eq '10') then mn3='Oct'
if (mon3 eq '11') then mn3='Nov'
if (mon3 eq '12') then mn3='Dec'
if (mon3 eq '01') then mn3='Jan'
if (mon3 eq '02') then mn3='Feb'
restore,'/home/dammasch/lyra'+yea3+'over/'+mn3+yea3+'/gaps'+date3+'.rst'
tim3=tim/60.
ch13=ch1*1000.
ch23=ch2*1000.
ch33=ch3*1000.
ch43=ch4*1000.

outpathfile='/home/dammasch/homepage/3DayQuicklook/LyraCalSWC'+date3+'.png'

; ------------------------------------------------------------------------------
; define upper and lower plotting limits and decide about lin or log scaling
; ------------------------------------------------------------------------------

lo3=(0.95*min([ch31,ch32,ch33]))>2.0
hi3=(1.1*max([ch31,ch32,ch33]))<15.
if (hi3 le lo3) then hi3=lo3+1.
lo4=(0.95*min([ch41,ch42,ch43]))>0.6
hi4=(1.1*max([ch41,ch42,ch43]))<10.
if (hi4 le lo4) then hi4=lo4+1.
gentickv3=[ 2., 3., 4., 5.,6.,7.,8.,9.,10.,15.]
gentickv4=[0.6,0.7,0.8,0.9,1.,2.,3.,4., 5., 6., 7., 8., 9.,10.]
ind3=where((gentickv3 ge lo3) and (gentickv3 le hi3))
ind4=where((gentickv4 ge lo4) and (gentickv4 le hi4))
if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
tickv3=gentickv3(ind3)
tickv4=gentickv4(ind4)
ticks3=n_elements(ind3)
ticks4=n_elements(ind4)
endif

; ------------------------------------------------------------------------------
; plot PNG image
; ------------------------------------------------------------------------------

plot,indgen(10)
!ymin=lo3 & !ymax=hi3
if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim1,ch31,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
/nodata,/normal,position=[0.08,0.53,0.38,0.95],xcharsize=0.000001,ycharsize=0.8,yticks=ticks3-1,ytickv=tickv3
endif else begin
plot,tim1,ch31,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
/nodata,/normal,position=[0.08,0.53,0.38,0.95],xcharsize=0.000001,ycharsize=0.8
endelse
oplot,tim1,ch31,psym=3,color=2
xyouts,0.10,0.92,'Channel 2-3',color=2,/normal

if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim2,ch32,psym=3,xstyle=1,ystyle=1, $
title='LYRA short-wavelength channels (1 minute averages)', $
/nodata,/normal,/noerase,position=[0.38,0.53,0.68,0.95],xcharsize=0.000001,ycharsize=0.000001
endif else begin
plot,tim2,ch32,psym=3,xstyle=1,ystyle=1, $
title='LYRA short-wavelength channels (1 minute averages)', $
/nodata,/normal,/noerase,position=[0.38,0.53,0.68,0.95],xcharsize=0.000001,ycharsize=0.000001
endelse
oplot,tim2,ch32,psym=3,color=2
xyouts,0.40,0.92,'Aluminium Filter',color=2,/normal

if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim3,ch33,psym=3,xstyle=1,ystyle=1, $
/nodata,/normal,/noerase,position=[0.68,0.53,0.98,0.95],xcharsize=0.000001,ycharsize=0.000001
endif else begin
plot,tim3,ch33,psym=3,xstyle=1,ystyle=1, $
/nodata,/normal,/noerase,position=[0.68,0.53,0.98,0.95],xcharsize=0.000001,ycharsize=0.000001
endelse
oplot,tim3,ch33,psym=3,color=2
xyouts,0.70,0.92,'17-80nm + Xray',color=2,/normal

!ymin=lo4 & !ymax=hi4
if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim1,ch41,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
xtitle='Time / h UTC, '+day1+' '+mn1+' '+yea1,/nodata, $
/normal,/noerase,position=[0.08,0.11,0.38,0.53],yticks=ticks4-1,ytickv=tickv4,ycharsize=0.8
endif else begin
plot,tim1,ch41,psym=3,xstyle=1,ystyle=1,ytitle='Calibrated Signal / (mW m!U-2!N)', $
xtitle='Time / h UTC, '+day1+' '+mn1+' '+yea1,/nodata, $
/normal,/noerase,position=[0.08,0.11,0.38,0.53],ycharsize=0.8
endelse
oplot,tim1,ch41,psym=3,color=3
xyouts,0.10,0.50,'Channel 2-4',color=3,/normal
if (total(ch41) lt 0.1) then xyouts,2,0.5*lo4+0.5*hi4,'(NO DATA AVAILABLE)',size=1.2

if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim2,ch42,psym=3,xstyle=1,ystyle=1, $
xtitle=day2+' '+mn2+' '+yea2,/nodata, $
/normal,/noerase,position=[0.38,0.11,0.68,0.53],ycharsize=0.000001
endif else begin
plot,tim2,ch42,psym=3,xstyle=1,ystyle=1, $
xtitle=day2+' '+mn2+' '+yea2,/nodata, $
/normal,/noerase,position=[0.38,0.11,0.68,0.53],ycharsize=0.000001
endelse
oplot,tim2,ch42,psym=3,color=3
xyouts,0.40,0.50,'Zirconium Filter',color=3,/normal
if (total(ch42) lt 0.1) then xyouts,2,0.5*lo4+0.5*hi4,'(NO DATA AVAILABLE)',size=1.2

if ((n_elements(ind3) gt 1) and (n_elements(ind4) gt 1)) then begin
plot_io,tim3,ch43,psym=3,xstyle=1,ystyle=1, $
xtitle=day3+' '+mn3+' '+yea3,/nodata, $
/normal,/noerase,position=[0.68,0.11,0.98,0.53],ycharsize=0.000001
endif else begin
plot,tim3,ch43,psym=3,xstyle=1,ystyle=1, $
xtitle=day3+' '+mn3+' '+yea3,/nodata, $
/normal,/noerase,position=[0.68,0.11,0.98,0.53],ycharsize=0.000001
endelse
oplot,tim3,ch43,psym=3,color=3
xyouts,0.70,0.50,'6-20nm + Xray',color=3,/normal
if (total(ch43) lt 0.1) then xyouts,2,0.5*lo4+0.5*hi4,'(NO DATA AVAILABLE)',size=1.2
xyouts,0.75,0.01,'ROB/SIDC, Brussels, Belgium',color=4,/normal,charsize=0.8

!xmin=0 & !xmax=0 & !ymax=0 & !ymin=0

; ------------------------------------------------------------------------------
; create PNG image
; ------------------------------------------------------------------------------

array=tvrd()
write_png,outpathfile,array,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]
outpathfile2='/home/dammasch/homepage/3DayQuicklook/LyraCalSWClatest.png'
write_png,outpathfile2,array,[0,255,0,128,0,255],[255,0,0,128,0,255],[0,0,255,128,0,255]

end
