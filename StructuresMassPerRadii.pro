  device,get_decompose=old_decomposed,decomposed=0
  loadct,12
;window,12,xsize=600,ysize=600
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  plot,Structures(0,*),MassTauRadPerBlock(0,*,0),background=255,$                    
       color=0,ytitle='MasaPerRadiiPerStructure ( M!L!9n!3!N )',yrange=[0.0000000001,0.0001],$
       xrange=[0,10],psym=1,xstyle=1,$
       xtitle='Estructura #',charsize=1.5,/ylog,ystyle=1
;  f=strtrim(f,2)
  nPlanets=strtrim(nPlanets,2)
;    xyouts,200,520,'f='+f,charsize=1.5,color=0,/device
  xyouts,200,490,nPlanets+' planets',charsize=1.5,color=0,/device
  for l=0,TauRanges-1 do begin
     for n=0,nPlanets do begin
        if n eq 0 then simb=1
        if n eq 1 then simb=6
        if n eq 2 then simb=7
        if l eq 0 then col=0
        if l eq 1 then col=100
        if l eq 2 then col=200
        oplot,Structures(l,*),MassTauRadPerBlock(l,*,n),color=col,psym=simb
     endfor
  endfor
  xyouts,390,530,'      tau < 0.5',charsize=1.5,color=0,/device
  xyouts,390,500,'0.5 < tau < 3.0',charsize=1.5,color=100,/device
  xyouts,390,470,'3.0 < tau      ',charsize=1.5,color=200,/device
  xyouts,390,440,'+ = 10 AU',charsize=1.5,color=0,/device
  xyouts,390,410,'square =15 AU',charsize=1.5,color=0,/device
  xyouts,390,380,'X =30 AU',charsize=1.5,color=0,/device
  write_png,strcompress(dischargeImagesFolder+$
                        '/Structures-MassPerRadiiPerStructure.png',/remove_all),tvrd(0,0,600,600,0,true=1)

