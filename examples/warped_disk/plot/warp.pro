pro warp

  oldn=!D.name & set_plot,'ps'
  
 ; !p.font=0
 
lx=801 

 Openr,lunAo10,'warpdiskg.txt',/Get_Lun
  Point_lun,lunAo10,0
  diskr=fltarr(lx, lx)
  ReadF,lunAo10,diskr

 ;Openr,lunAo10,'warpdisk1_0045.txt',/Get_Lun
 ; Point_lun,lunAo10,0
;  diskr2=fltarr(801,801)
;  ReadF,lunAo10,diskr2

; Openr,lunAo10,'warpdisk1_0090.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr3=fltarr(801,801)
;  ReadF,lunAo10,diskr3

; Openr,lunAo10,'warpdisk1_0135.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr4=fltarr(801,801)
;  ReadF,lunAo10,diskr4
  
; Openr,lunAo10,'warpdisk1_0180.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr5=fltarr(801,801)
;  ReadF,lunAo10,diskr5
  
; Openr,lunAo10,'warpdisk1_0225.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr6=fltarr(801,801)
;  ReadF,lunAo10,diskr6

; Openr,lunAo10,'warpdisk1_0270.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr7=fltarr(801,801)
;  ReadF,lunAo10,diskr7

; Openr,lunAo10,'warpdisk1_0315.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr8=fltarr(801,801)
;  ReadF,lunAo10,diskr8

; Openr,lunAo10,'warpdisk1_0360.txt',/Get_Lun
;  Point_lun,lunAo10,0
;  diskr9=fltarr(801,801)
;  ReadF,lunAo10,diskr9

  nnn=300
  ;diskr[0,0]=2.541	
  ;maxr=max(diskr,min=minr)
  ;print,maxr,minr

 ;======================================================

  ;===============================================
  NN=100
  MM=450
  A=1.2
  B=1.05
    Lc=5.
  ratio=1.;3.529
   l=22. & xxss=l*(ratio) & yyss=l
   ;xoff=(LL-xxss)/2.,yoff=(3*LL/2.-yyss)/2.,
   device,filename='warp1_9.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,/color,xoff=(21.5-xxss)/2.0,yoff=(30-yyss)/2.,$
   set_font='Helvetica Italic', /tt_font
  ; set_font='isolatin1'
   XX=Findgen(101)
   YY=Findgen(101)
   XXX=findgen(11)*(0.1)
  rr=fltarr(100,100)*0+1
  cc=findgen(100)*(7.1416/99)
  xc=cos(findgen(100))
  yc=sin(findgen(100))
 
; loadct,33
  RRR=bytscl(findgen(256))
  GGG=bytscl(findgen(256))
  BBB=bytscl(findgen(256))
  RRR[248:255]=[0,  0,255,255,0, 0,255,255]
  GGG[248:255]=[255,0,0,  255,0,255,0,255]
  BBB[248:255]=[0,  0,0,  0  ,255,255,255,255]
  ;248=green,249=black,250=red,251=yellow,252=blue,253=cyan,254=magenta,255=white
  TVLCT,RRR,GGG,BBB
  
  ;CCol=[248,249,250,251,252,253,254]
  CCol=252
  
 ;plot,[0,800],[0,800],pos=[(xxss/(2*Lc))*1000,(yyss/(2*Lc))*1000,$
  ;(xxss*(2*Lc-1)/(2*Lc))*1000,(yyss*(2*Lc-1)/(2*Lc))*1000],/noerase,/nodata,/device,/ynozero,charsize=2,yticklen=0.001;,yminor=5  

		;index=where(diskr4 ne 0. and diskr4 ne 0.1,nzero)
		;if nzero ne 0 then begin
  		;   maxr=max(diskr4(index),min=minr)
		;endif
		;print,maxr,minr

 		;n_one=floor((maxr-1.)/(maxr-minr)*255)
;print,n_one
 ;loadct,33
 ;tvlct,originalr,originalg,originalb,/get
 ;originalr[n_one-1:n_one+1]=255
 ;originalg[n_one-1:n_one+1]=255
 ;originalb[n_one-1:n_one+1]=255

 ;originalr[0]=255 ;wihte color for blank region
 ;originalg[0]=255
 ;originalb[0]=255

 ;originalr[255]=0 ;black color for black hole shadow
 ;originalg[255]=0
 ;originalb[255]=0
 ;TVLCT,originalr,originalg,originalb
;diskr=rotate(diskr,2)
 ;D=bytscl(diskr)
 xl=xxss*1000.*2./3.
 yl=yyss*1000.*2./3.
	maxg=1.4
	colors=255
	chsize=0.7
	chth=2
	tickth=2
	alp1=-55
	alp2=55
	beta1=-55
	beta2=55
	xlen=0.26*3+0.06*2
	ylen=0.26*3+0.08*2
	xslitlen=0.06
	yslitlen=0.08
	mx=1 ;figure numbers on horizon direction
	my=1 ;figure numbers on vertical direction.
	figlenx=(xlen-xslitlen*(mx-1))/mx
	figleny=(ylen-yslitlen*(my-1))/my
	delx=figlenx+xslitlen
	dely=figleny+yslitlen
	deltax=[delx,0,delx,0]
	deltay=[0,dely,0,dely]
	deltxy=[delx,dely,delx,dely]
	pos_llft=[(1.-xlen)/2.,0.01,(1.-xlen)/2.+figlenx,0.01+figleny]
   for ny=0,my-1 do begin	
     	for nx=0,mx-1 do begin	
		posxy=[pos_llft+deltax*nx+deltay*ny]
		plot,[alp1,alp2],[beta1,beta2],pos=posxy,xrange=[alp1,alp2],yrange=[beta1,beta2],$
		/noerase,/device,/ynozero,/normal,xstyle=4+1,ystyle=4+1,charsize=0.5,xtickv=[0,0.5,1,1.5],$
		xticks=10,xminor=10,xtickname=replicate(' ',10),/nodata

                posit = [(1.-xlen)/2.,0.01+figleny+yslitlen/4.,(1.-xlen)/2.+figlenx,$
		          0.01+figleny+yslitlen/4.*2.]+nx*deltax+ny*deltay

	   If ny*mx+nx eq 0 then begin
		disk=diskr   		
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		;wh1=where(disk gt 1.2)
		;disk(wh1)=1
		deltc=0.05			
		texts='a'
	   endif
	   If ny*mx+nx eq 7 then begin
		disk=diskr2
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05			
		texts='b'
	   endif
	   If ny*mx+nx eq 8 then begin
		disk=diskr3
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='c'
	   endif
	   If ny*mx+nx eq 3 then begin
		disk=diskr4
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='d'
	   endif
	   If ny*mx+nx eq 4 then begin
		disk=diskr5
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='e'
	   endif
	   If ny*mx+nx eq 5 then begin
		disk=diskr6
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='f'
	   endif
	   If ny*mx+nx eq 1 then begin
		disk=diskr7
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='d'
	   endif
	   If ny*mx+nx eq 1 then begin
		disk=diskr8
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='e'
	   endif
	   If ny*mx+nx eq 2 then begin
		disk=diskr9
                index=where(disk ne 0. and disk ne 0.1,nzero)
  		If (nzero ne 0) then begin
		    	maxr1 = max(disk(index),min = minr1)
  		endif  
	  	disk[0,0]=maxg
                ncolor1 = floor(maxr1/disk[0,0]*256)
		deltc=0.05	
		texts='f'
	   endif
                print,maxr1,minr1,ncolor1	
		wh1=where(abs(abs(disk)-1.) lt 0.0031,nzero)
		If nzero ne 0 then begin
			;disk(wh1)=0.05
		endif
		print,max(disk),min(disk)
		index=where(disk ne 0. and disk ne 0.1 and disk ne 0.05,nzero)
		if nzero ne 0 then begin
  		   maxr=max(disk(index),min=minr)
		endif
		;print,maxr,minr,texts

  		index1=where(disk eq 0.,nzero)
  		If (nzero ne 0) then begin
		    	disk(index1)=minr-deltc
  		endif 
  		index2=where(disk eq 0.1,nhori)
  		If (nhori ne 0)then begin
  			disk(index2)=minr-deltc;-0.001;0.9090909
  		endif
  		index3=where(disk eq 0.05,nhori)
  		If (nhori ne 0)then begin
  			;disk(index3)=maxr;deltc;-0.001;0.9090909
  		endif
 	;==========************************
	;set the color table.

 		;n_one=floor((maxr-1.)/(maxr-minr)*255)
		;print,n_one
 		loadct,33
 		tvlct,originalr,originalg,originalb,/get
 		;originalr[n_one-1:n_one+1]=255
 		;originalg[n_one-1:n_one+1]=0
 		;originalb[n_one-1:n_one+1]=0

 		originalr[0]=255 ;wihte color for blank region
 		originalg[0]=255
 		originalb[0]=255

 		originalr[255]=0 ;black color for black hole shadow
 		originalg[255]=0
 		originalb[255]=0
 		TVLCT,originalr,originalg,originalb
 	;*********************************
		D=bytscl(-disk)
		tvscl,rotate(D,0),posxy[0],posxy[1],/NaN,/normal,xsize=figlenx,ysize=figleny

		axis,xaxis=0,xticks=6,xminor=4,xrange=[alp1,alp2],xstyle=1,charsize=chsize,charthick=chth,font=-1,$
		xthick=tickth,xtitle=textoidl('\alpha'),color=colors;xtickname=replicate(' ',2),xtickname=['-4','-2','0','2','4','6','8'],

		axis,xaxis=1,xticks=6,xminor=4,xrange=[alp1,alp2],xstyle=1,xtickname=replicate(' ',7),font=-1,$
		charsize=chsize,charthick=chth,xthick=tickth,color=colors

		axis,yaxis=0,yticks=6,yminor=4,yrange=[beta1,beta2],ystyle=1,charsize=chsize,charthick=chth,font=-1,$
		ythick=tickth,ytitle=textoidl('\beta'),color=colors;ytickname=replicate(' ',2),ytickname=['-6','-4','-2','0','2','4','6'],
		;,ytickv=[0,0.5,1]

	 	axis,yaxis=1,yticks=6,yminor=4,yrange=[beta1,beta2],ystyle=1,ytickname=replicate(' ',7),font=-1,$
		charsize=chsize,charthick=chth,ythick=tickth,color=colors;

		xyouts,-55,-55,texts,charthick=3,charsize=1,color=colors,font=1

		 ;************************************************************
		  colorbar,color=colors, bottom=2,DIVISIONS=6,NCOLORS=ncolor1,CHARSIZE=0.6, FORMAT='(F9.2)',$
                     POSITION=posit,$;VERTICAL=544,$
		  PSCOLOR=pscolor,TITLE=' ',RIGHT=32,Max=maxr1,Min=minr1
		  print,maxr,minr
		 ;************************************************************
	endfor
   endfor 
   ;tvscl,rotate(D,0),(xxss)/3.0/2.*1000,(yyss)/3./2.*1000,/NaN,xsize=xl,ysize=yl
 ;tv,diskr

  ;axis,xaxis=0,xtitle='',xthick=5,charthick=3,font=1;,charsize=2
  ;axis,yaxis=0,ytitle='',xticklen=-0.2,ythick=3,charthick=3,font=1;,font='Helvetica Bold italic'   
  
 ;xyouts,3.3,2.6,'',charthick=3,charsize=2
   device,/close 
   
   free_lun,lunAo10
   set_plot,oldn
   
end
