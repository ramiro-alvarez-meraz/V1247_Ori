;Identifica distribucion de area vs tamaño radial de cada estructura
;  window,6,xsize=600,ysize=600
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  plot,BlockRad_fin(0,*)-BlockRad_ini(0,*),AreaTotalPerBlock(0,*),background=255,$                            
       color=1,xtitle='Tama!Z(00f1)o de estructura',psym=1,xstyle=1,yrange=[1.e-6,1.e5],$
       ytitle='Area de estructura ( UA!E2!N )',charsize=1.5,/xlog,/ylog,xrange=[1.e-4,300]
;         f=strtrim(f,2)
  nPlanets=strtrim(nPlanets,2)
;    xyouts,200,520,'f='+f,charsize=1.5,color=0,/device
  xyouts,200,490,nPlanets+' planets',charsize=1.5,color=0,/device
  for l=0,TauRanges-1 do begin
     oplot,BlockRad_fin(l,*)-BlockRad_ini(l,*),AreaTotalPerBlock(l,*),psym=1,color=(l+1)*60
     if l eq 0 then xyouts,130,530,'       tau < 0.5',charsize=1.5,color=(l+1)*60,/device
     if l eq 1 then xyouts,130,500,'0.5 < tau < 3.0',charsize=1.5,color=(l+1)*60,/device
     if l eq 2 then xyouts,130,470,'3.0 < tau',charsize=1.5,color=(l+1)*60,/device
  endfor
  write_png,strcompress(dischargeImagesFolder+'/AreaStructures-Tamanho'+string(l)+$
                        '.png',/remove_all),tvrd(0,0,600,600,0,true=1)
