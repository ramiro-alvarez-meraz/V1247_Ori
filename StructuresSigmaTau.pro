;Define new palette of colors
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;Radial density, radial acretion and radial space ratio:
;  window,8,xsize=600,ysize=600
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  plot,rad,2./2.25*1.e7*SigmaTau(*,0)/nsec,background=255,$                            
       color=0,xtitle='Radius ( AU )',yrange=[0.001,200.],psym=0,xstyle=1,ystyle=1,$
       ytitle='<  !N!4R!3!N (g cm!E-2!N)  >',charsize=1.5,/xlog,/ylog
  nPlanets=strtrim(nPlanets,2)
;  xyouts,150,530,'f='+f,charsize=1.5,color=0,/device
;  xyouts,150,500,nPlanets+' planets',charsize=1.5,color=0,/device
;  xyouts,140,470,strcompress(string(ulong(output*1000.))+' years'),charsize=1.5,color=0,/device
  oplot,rad,2./2.25*1.e7*SigmaTau(*,0)/nsec,color=70,psym=0
  xyouts,390,530,'       !4s!X < 0.5',charsize=1.5,color=70,/device
  oplot,rad,2./2.25*1.e7*SigmaTau(*,1)/nsec,color=160,psym=0
  xyouts,390,500,'0.5 < !4s!X < 2.0',charsize=1.5,color=160,/device
  oplot,rad,2./2.25*1.e7*SigmaTau(*,2)/nsec,color=250,psym=0
  xyouts,390,470,'2.0 < !4s!X',charsize=1.5,color=250,/device
  write_png,strcompress(dischargeImagesFolder+$
                        '/SigmaTau-SemiejeMayor.png',/remove_all),tvrd(0,0,600,600,0,true=1)
