; Make a vector of 16 points, A[i] = 2pi/16:
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, COS(A), SIN(A), /FILL

  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;  window,10,xsize=600,ysize=600
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
   for l=0,TauRanges-1 do begin
     for m=1,Block_prev[l] do begin
        Structures[l,m]=m
     endfor
  endfor 
  SemiejeMaxArray=dblarr(Struct_max)
  SemiejeMinArray=dblarr(Struct_max)
;Ciclos para calculo de densidad media
  densMed=dblarr(nrad)
  tot=dblarr(nrad)
  gapBorder=dblarr(2)
  for r=0,nrad-1 do begin
     tot[r]=0.0
     for l=0,nsec-1 do begin
        tot[r]=dens[l,r]+tot[r]
     endfor
  endfor
  for r=0,nrad-1 do begin
     densMed[r]=tot[r]/nsec
  endfor
  densCGSmed = densMed*2./2.25*1.e7
  tauMed=9.6*densCGSmed
  for r=nrad-1,0,-1 do begin
     if rad[r] lt Rplanet[0] then begin
        if tauMed[r] ge 2. then begin
           gapBorder[0]=rad[r]
           break
        endif
     endif
  endfor
  for r=0,nrad-1 do begin
     if rad[r] gt Rplanet[nPlanets-1] then begin
        if tauMed[r] ge 2. then begin
           gapBorder[1]=rad[r]
           break
        endif
     endif
  endfor
print,gapBorder
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  SemiejeMaxArray(*)=SemiejeMax+(2*HillPlanet[nPlanets-1])
;  SemiejeMinArray(*)=SemiejeMin-(2*HillPlanet[0])
  SemiejeMaxArray(*)=gapBorder[1]
  SemiejeMinArray(*)=gapBorder[0]
  Struct_maxArray=findgen(Struct_max)
  plot,Structures(3,*),BlockRad_ini(3,*),background=255,$                            
       color=0,ytitle='inner and outer edge by structure ( AU )',yrange=[1.,500],xrange=[0,8],psym=6,xstyle=1,$
       xtitle='Structure #',charsize=1.5,/ylog,ystyle=1,symsize=1.1,thick=2
;  f=strtrim(f,2)
  nPlanets=strtrim(nPlanets,2)
;  xyouts,150,530,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,530,nPlanets+' planets',charsize=1.5,color=0,/device
  xyouts,140,500,strcompress(string(ulong(output*1000.))+' years'),charsize=1.5,color=0,/device
  oplot,Structures(3,*),BlockRad_fin(3,*),color=70,psym=6,symsize=1.1,thick=2
  for l=0,TauRanges-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if l eq 3 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      simb=6
      col=40
      oplot,Structures(l,*),BlockRad_ini(l,*),color=col,psym=simb,symsize=1.1,thick=2
      oplot,Structures(l,*),BlockRad_fin(l,*),color=col,psym=simb,symsize=1.1,thick=2
    endif
    if l eq 4 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      simb=6
      col=70
      oplot,Structures(l,*),BlockRad_ini(l,*),color=col,psym=simb,symsize=1.1,thick=2
      oplot,Structures(l,*),BlockRad_fin(l,*),color=col,psym=simb,symsize=1.1,thick=2  
    endif
    if l eq 5 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      simb=7
      col=100
      oplot,Structures(l,*),BlockRad_ini(l,*),color=col,psym=simb,symsize=1.1,thick=2  
      oplot,Structures(l,*),BlockRad_fin(l,*),color=col,psym=simb,symsize=1.1,thick=2  
    endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     if l eq 0 then begin
;        simb=6
;        col=70
;     endif
;     if l eq 1 then begin 
;        simb=6
;        col=160
;     endif
;     if l eq 2 then begin
;        simb=6
;        col=250
;     endif
;    oplot,Structures(l,*),BlockRad_ini(l,*),color=col,psym=8   
;    oplot,Structures(l,*),BlockRad_fin(l,*),color=col,psym=simb  
  endfor
  oplot,Struct_maxArray(*),SemiejeMinArray(*),color=0,psym=0
  oplot,Struct_maxArray(*),SemiejeMaxArray(*),color=0,psym=0
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;  xyouts,390,530,'      !4s!X < 0.5',charsize=1.5,color=70,/device
  xyouts,390,530,'0.5 < !4s!X < 1.0',charsize=1.5,color=40,/device
  xyouts,390,500,'1.0 < !4s!X < 2.0',charsize=1.5,color=70,/device
  xyouts,390,470,'2.0 < !4s!X      ',charsize=1.5,color=100,/device
  write_png,strcompress(dischargeImagesFolder+$
                        '/Structures-Sizes.png',/remove_all),tvrd(0,0,600,600,0,true=1)
