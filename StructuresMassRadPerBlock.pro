;Identifica contribucion de masa de cada estructura
  for l=0,TauRanges-1 do begin
;     window,l,xsize=600,ysize=600
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
     plot,rad,MassRadPerBlock(*,l,0),background=255,$                            
          color=1,xtitle='Semimajor axis ( AU )',yrange=[1.e-11,1.e-1],psym=0,xstyle=1,$
          ytitle='Radial mass per structure ( M!L!9n!3!N )',charsize=1.5,/xlog,/ylog
;  f=strtrim(f,2)
  nPlanetsStr=strtrim(nPlanets,2)
;    xyouts,200,520,'f='+f,charsize=1.5,color=0,/device
  xyouts,200,490,nPlanetsStr+' planets',charsize=1.5,color=0,/device
     for m=1,Block_prev[l] do begin
        oplot,rad,MassRadPerBlock(*,l,m),psym=0,color=m*35
     endfor
     if l eq 0 then xyouts,420,500,'       tau < 0.5',charsize=1.5,color=0,/device
     if l eq 1 then xyouts,420,500,'0.5 < tau < 3.0',charsize=1.5,color=0,/device
     if l eq 2 then xyouts,420,500,'3.0 < tau',charsize=1.5,color=0,/device
     write_png,strcompress(dischargeImagesFolder+'/MassRadPerBlock-SemiejeMayor'+$
                           string(l)+'.png',/remove_all),tvrd(0,0,600,600,0,true=1)
  endfor
