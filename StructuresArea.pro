  AccPlusTotalPerBlock=dblarr(TauRanges,Struct_max)
  AccMinusTotalPerBlock=dblarr(TauRanges,Struct_max)
  for l=0,TauRanges-1 do begin
     for m=1,Block_prev[l] do begin
        if (AccTotalPerBlock[l,m] gt 0) then begin
           AccPLusTotalPerBlock[l,m] = AccTotalPerBlock[l,m]
        endif
        if (AccTotalPerBlock[l,m] lt 0) then begin
           AccMinusTotalPerBlock[l,m] = AccTotalPerBlock[l,m]
        endif
     endfor
  endfor

;  window,7,xsize=600,ysize=600
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,2,0,0]
  plot,BlockRad_fin(0,*)-BlockRad_ini(0,*),abs(AccMinusTotalPerBlock(0,*)),$
       background=255,color=1,xtitle='Tama!Z(00f1)o de estructura ( UA )',$
       psym=1,xstyle=1,yrange=[1.e-13,1.e-6],$
       ytitle='Acrecion masa por estructura (M!L!9n!3!N/a!Z(00f1)o)',$
       charsize=1.0,/xlog,/ylog,xrange=[1.e-4,300]
;                f=strtrim(f,2)
  nPlanets=strtrim(nPlanets,2)
;    xyouts,200,520,'f='+f,charsize=1.5,color=0,/device
  xyouts,200,490,nPlanets+' planets',charsize=1.5,color=0,/device
  for l=0,TauRanges-1 do begin
     oplot,BlockRad_fin(l,*)-BlockRad_ini(l,*),abs(AccMinusTotalPerBlock(l,*)),$
           psym=1,color=(l+1)*60
     if l eq 0 then xyouts,130,530,'       tau < 0.5',charsize=1.5,color=(l+1)*60,/device
     if l eq 1 then xyouts,130,500,'0.5 < tau < 3.0',charsize=1.5,color=(l+1)*60,/device
     if l eq 2 then xyouts,130,470,'3.0 < tau',charsize=1.5,color=(l+1)*60,/device
  endfor
  plot,BlockRad_fin(0,*)-BlockRad_ini(0,*),abs(AccPlusTotalPerBlock(0,*)),$
       background=255,color=1,xtitle='Tama!Z(00f1)o de estructura ( UA )',$
       psym=1,xstyle=1,yrange=[1.e-13,1.e-6],$
       ytitle='Salida de masa por estructura (M!L!9n!3!N/a!Z(00f1)o)',$
       charsize=1.0,/xlog,xrange=[1.e-4,300],/ylog
  for l=0,TauRanges-1 do begin
     oplot,BlockRad_fin(l,*)-BlockRad_ini(l,*),abs(AccPlusTotalPerBlock(l,*)),$
           psym=1,color=(l+1)*60
     if l eq 0 then xyouts,130,230,'       tau < 0.5',charsize=1.5,$
                           color=(l+1)*60,/device
     if l eq 1 then xyouts,130,200,'0.5 < tau < 3.0',charsize=1.5,$
                           color=(l+1)*60,/device
     if l eq 2 then xyouts,130,170,'3.0 < tau',charsize=1.5,$
                           color=(l+1)*60,/device
  endfor
  write_png,strcompress(dischargeImagesFolder+'/AreaStructures-Tamanho'+string(l)+'.png',$
                        /remove_all),tvrd(0,0,600,600,0,true=1)  

