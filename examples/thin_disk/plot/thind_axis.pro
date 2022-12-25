pro thind_axis

  oldn=!D.name & set_plot,'ps'
  
 ; !p.font=0
 
lx=801
ly=801 

 Openr,lunAo10,'diskg.txt',/Get_Lun
  Point_lun,lunAo10,0
  diskr=fltarr(lx,ly)
  ReadF,lunAo10,diskr
   
  

  nnn=300
  maxr=max(diskr,min=minr)
  print,maxr,minr

  index=where(abs(diskr) ne 0.,nzero)
  If (nzero ne 0) then begin
	maxr1 =max(diskr(index),min=minr1)
  endif 

  diskr=abs(diskr)	

  ;maxr=max(diskr,min=minr)
  print,maxr1,minr1

  index=where(diskr eq 0.,nzero)
  If (nzero ne 0) then begin
	diskr(index)=maxr;-0.05
  endif 

  index=where(diskr eq 0.1,nhori)
  If (nhori ne 0)then begin
  	diskr(index)=maxr;minr+0.05;0.9090909
  endif


  ;index=where(diskr eq 1.,nzero)
  ;If (nzero ne 0) then begin
;	diskr(index)=minr-0.05
  ;endif 
  ;index=where(diskr eq 0.9090909,nhori)
  ;If (nhori ne 0)then begin
  ;	diskr(index)=minr
  ;endif
	
  LLevell=Findgen(nnn+1.0)*((maxr-minr)/nnn)+minr
;print,LLevell
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
   device,filename='weaklens12_11.ps',xsize=xxss,ysize=yyss,bits_per_pixel=8,/color,xoff=(21.5-xxss)/2.0,yoff=(30-yyss)/2.,$
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

 n_one=floor((maxr-1.)/(maxr-minr)*255)
print,n_one
 loadct,0
 tvlct,originalr,originalg,originalb,/get
 ;originalr[n_one]=255
 ;originalg[n_one]=255
; originalb[n_one]=255

 originalr[255]=255
 originalg[255]=255
 originalb[255]=255

 TVLCT,originalr,originalg,originalb
;diskr=rotate(diskr,2)
 ;diskr[0,0] = -5.
 D=bytscl(diskr)
 xl=xxss*1000.*2./3.
 yl=yyss*1000.*2./3.
	ratio= 30./60.
	xlen=0.8
	ylen=xlen*ratio
	mx=1 ;figure numbers on horizon direction
	my=1 ;figure numbers on vertical direction.
	deltax=[xlen/mx,0,xlen/mx,0]
	deltay=[0,ylen/my,0,ylen/my]
	deltxy=[xlen/mx,ylen/my,xlen/mx,ylen/my]
	pos_llft=[(1.-xlen)/2.,(1.-ylen)/2.,(1.-xlen)/2.+xlen/mx,(1.-ylen)/2.+ylen/my]
   for ny=0,my-1 do begin	
     	for nx=0,mx-1 do begin	
		plot,[-10,10],[-10,10],pos=[pos_llft+deltax*nx+deltay*ny],xrange=[-30,30],yrange=[-15,15],$
		/noerase,/device,/ynozero,/normal,xstyle=4+1,ystyle=4+1,charsize=0.5,xtickv=[0,0.5,1,1.5],$
		xticks=10,xminor=10,xtickname=replicate(' ',10),/nodata

		tvscl,rotate(D,0),(1.-xlen)/2.,(1.-ylen)/2.,/NaN,/normal,xsize=xlen,ysize=ylen

		axis,xaxis=0,xticks=6,xminor=5,xrange=[-30,30],xstyle=1,charsize=1.3,charthick=3.3,font=-1,$
		xthick=3,xtitle=textoidl('\alpha');xtickname=replicate(' ',2)xtickname=['-4','-2','0','2','4','6','8'],

		axis,xaxis=1,xticks=6,xminor=5,xrange=[-30,30],xstyle=1,xtickname=replicate(' ',7),font=-1,$
		charsize=2.5,charthick=3,xthick=3

		axis,yaxis=0,yticks=4,yminor=6,yrange=[-13,15],ystyle=1,charsize=1.3,charthick=3.3,font=-1,$
		ythick=3,ytitle=textoidl('\beta');ytickname=replicate(' ',2)ytickname=['-6','-4','-2','0','2','4','6'],
		;,ytickv=[0,0.5,1]

	 	axis,yaxis=1,yticks=4,yminor=6,yrange=[-13,15],ystyle=1,ytickname=replicate(' ',7),font=-1,$
		charsize=2.5,charthick=3,ythick=3;
	endfor
   endfor 
 ;************************************************************
  colorbar,DIVISIONS=6,CHARSIZE=0.8, FORMAT='(F9.3)',POSITION=[0.10,0.71,0.9,0.75],$;VERTICAL=544,$
  PSCOLOR=pscolor,TITLE=' ',TOP=32,RIGHT=32,Max=maxr1,Min=minr1
  print,maxr,minr
 ;************************************************************
   ;tvscl,rotate(D,0),(xxss)/3.0/2.*1000,(yyss)/3./2.*1000,/NaN,xsize=xl,ysize=yl
 ;tv,diskr

  ;axis,xaxis=0,xtitle='',xthick=5,charthick=3,font=1;,charsize=2
  ;axis,yaxis=0,ytitle='',xticklen=-0.2,ythick=3,charthick=3,font=1;,font='Helvetica Bold italic'   
  
 ;xyouts,3.3,2.6,'',charthick=3,charsize=2
   device,/close 
   
   free_lun,lunAo10 
   set_plot,oldn
   
end
