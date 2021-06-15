     Set_plot,'z'
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  MassPerBlock=dblarr(TauRanges,Struct_max)
  DustMassPerBlock=dblarr(TauRanges,Struct_max)
  for l=0,TauRanges-1 do begin
     for m=1,Block_prev[l] do begin
        MassPerBlock[l,m]=total(MassRadPerBlock(*,l,m))
        DustMassPerBlock[l,m]=total(MassRadPerBlock(*,l,m))*0.005168   
     endfor
  endfor
  plot,Structures(3,*),DustMassPerBlock(3,*)/0.001*317.,background=255,$                            
       color=0,ytitle='Dust mass!L0.005-0.25 [!4l!Xm]!N ( M!L!20S!3!N )',yrange=[0.0000002,99],$;M!L!9n!3!N
       xrange=[0,8],psym=6,xstyle=1,$
       xtitle='Structure #',charsize=1.5,/ylog,ystyle=1,symsize=1.1,thick=2
;         f=strtrim(f,2)
  nPlanets=strtrim(nPlanets,2)
;  xyouts,150,530,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,530,nPlanets+' planets',charsize=1.5,color=0,/device
  xyouts,140,500,strcompress(string(ulong(output*1000.))+' years'),charsize=1.5,color=0,/device
  smallMass=dblarr(Struct_max)
  simpleArray=findgen(Struct_max)
;  smallMass(*)=0.0055772788
  smallMass(*)=0.001
;  oplot,simpleArray(*),smallMass(*),color=0,psym=0,linestyle=1
  for l=0,TauRanges-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if l eq 3 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      simb=6
      col=40
      oplot,Structures(l,*),DustMassPerBlock(l,*)/0.001*317.,color=col,psym=simb,symsize=1.1,thick=2
    endif
    if l eq 4 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      simb=6
      col=70
     oplot,Structures(l,*),DustMassPerBlock(l,*)/0.001*317.,color=col,psym=simb,symsize=1.1,thick=2
    endif
    if l eq 5 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      simb=7
      col=100
     oplot,Structures(l,*),DustMassPerBlock(l,*)/0.001*317.,color=col,psym=simb,symsize=1.1,thick=2
    endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;     if l eq 0 then begin
;        simb=1
;        col=70
;     endif
;     if l eq 1 then begin 
;        simb=1
;        col=160
;     endif
;     if l eq 2 then begin
;        simb=1
;        col=250
;     endif
;     oplot,Structures(l,*),DustMassPerBlock(l,*)/0.001*317.,color=col,psym=simb
  endfor
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
;  xyouts,390,530,'      !4s!X < 0.5',charsize=1.5,color=70,/device
  xyouts,390,530,'0.5 < !4s!X < 1.0',charsize=1.5,color=40,/device
  xyouts,390,500,'1.0 < !4s!X < 2.0',charsize=1.5,color=70,/device
  xyouts,390,470,'2.0 < !4s!X      ',charsize=1.5,color=100,/device
  write_png,strcompress(dischargeImagesFolder+$
                        '/Structures-Mass.png',/remove_all),tvrd(0,0,600,600,0,true=1)
