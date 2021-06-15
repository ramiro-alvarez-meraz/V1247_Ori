        alpha=dblarr(nsec)
        x=dblarr(nsec,nrad)
        y=dblarr(nsec,nrad)
     rmed = (rad(1:nrad)-rad(0:nrad-1))/2+rad(0:nrad-1)
        for ii=0,nrad-1 do begin
           for jj=0,nsec-1 do begin
              alpha[jj]=((jj+0.5)/nsec)*2*!PI
        if (jj eq 0) then  alpha[jj]=0.001
        if (jj eq nsec-1) then  alpha[jj]=2*!PI
              x[jj,ii]=rmed[ii]*cos(alpha[jj])
              y[jj,ii]=rmed[ii]*sin(alpha[jj])              
           endfor
        endfor
  planet=dblarr(10,output+1)
  xplanet=dblarr(nPlanets)
  yplanet=dblarr(nPlanets)
  vxplanet=dblarr(nPlanets)
  vyplanet=dblarr(nPlanets)
  mplanet=dblarr(nPlanets)
  Rplanet=dblarr(nPlanets)
  Hillplanet=dblarr(nPlanets)
  ThetaPlanet=dblarr(nPlanets)
  for n=0,nPlanets-1 do begin
     planet_=strcompress(dataFolder+'/planet'+string(n)+'.dat',/remove_all)
     openr,1,planet_
     s=dblarr(10)
     oldS=0
     for l=0,1100 do begin
        readf,1,s 
        if ((s[0] eq oldS) and (l gt 0)) then continue 
        nn=abs(s[0])
        planet[*,nn]=s
        oldS=s[0]
        if (s[0] eq output) then break
     endfor
     close,1
     xplanet[n]=planet(1,output)
     yplanet[n]=planet(2,output)
     vxplanet[n]=planet(4,output)
     vyplanet[n]=planet(5,output)
     h=xPlanet[n]*vyPlanet[n]-yPlanet[n]*vxPlanet[n]
     d=sqrt(xPlanet[n]*xPlanet[n]+yPlanet[n]*yPlanet[n])
     Ax=xPlanet[n]*vyPlanet[n]*vyPlanet[n]-yPlanet[n]*vxPlanet[n]*vyPlanet[n]-xPlanet[n]/d
     Ay=yPlanet[n]*vxPlanet[n]*vxPlanet[n]-xPlanet[n]*vxPlanet[n]*vyPlanet[n]-yPlanet[n]/d
     ee=sqrt(Ax*Ax+Ay*Ay)
     a=h*h/(1-ee*ee)
     mplanet[n]=planet(7,output)*1.86;/Msolarcgs;/MstarMsolar
     Rplanet[n]=sqrt(xplanet[n]^2+yplanet[n]^2)
     Hillplanet[n]=a*((1./3)*mplanet[n]/1.86)^(1./3.)
     ThetaPlanet[n]=atan(xplanet[n],yplanet[n])+(1.0/nsec*2*!PI)
     if (ThetaPlanet[n] < 0.) then ThetaPlanet[n]=(2.*!PI)+ThetaPlanet[n];+(1.0/nsec*2*!PI)
;print,n,'mplanet=',mplanet[n],ThetaPlanet[n]

   endfor

;The plots
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;Calculo de profundidad optica con Xi=9.6 cm^2 g^⁻1
  tau=9.6*dens*2./2.25*1.e7
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  contour,tau,x,y,charsize=1.5,levels=[0.001,0.01,0.1,0.5,1.0,2.0],$
          c_color=[50,120,130,150,170,240],$
;          c_color=[50,100],$
          xrange=[-50,50],yrange=[-50,50],$
          /xstyle,/fill,/ystyle,xtitle='AU',ytitle='AU'
;loadct,3
;  contour,alog10(tau),x,y,levels=[alog10(0.75),alog10(1.0),alog10(1.25)],$
;          c_color=[50,100,150],/overplot
  nPlanets=strtrim(nPlanets,2)
;  xyouts,150,530,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,520,nPlanets+' planets',charsize=1.5,color=0,/device
;  xyouts,150,100,strcompress(string(abs(output*1000.))+' years'),charsize=1.5,color=0,/device
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;  xyouts,390,540,strcompress('       !4s!X < '+tau_flagStr[0]),charsize=1.0,color=70,/device ;ñ=!Z(00f1)os
;  xyouts,390,520,strcompress(tau_flagStr[0]+' < !4s!X < '+tau_flagStr[1]),charsize=1.0,color=160,/device
;  xyouts,390,500,strcompress(tau_flagStr[1]+' < !4s!X < '+tau_flagStr[2]),charsize=1.0,color=160,/device
;  xyouts,390,480,strcompress(tau_flagStr[2]+' < !4s!X < '+tau_flagStr[3]),charsize=1.0,color=160,/device
;  xyouts,390,460,strcompress(tau_flagStr[3]+' < !4s!X'),charsize=1.0,color=240,/device
;  device,get_decompose=old_decomposed,decomposed=0
;  loadct,0
;  oplot,xcero,ycero,psym=2
  oplot,-xplanet,-yplanet,psym=1,color=250
  write_png,strcompress(dischargeEmisionFolder+'/TauTransition<50.png',/remove_all),tvrd(0,0,600,600,0,true=1)

  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  contour,alog10(tau),x,y,charsize=1.5,levels=[alog10(0.1),alog10(0.5),alog10(2.0)],$;-1.0,-0.30103,0.0,0.30103],$
          c_color=[50,100,150],$
          xrange=[-50,50],yrange=[-50,50],$
          /xstyle,/fill,/ystyle,xtitle='AU',ytitle='AU' ;,$
  nPlanets=strtrim(nPlanets,2)
;  xyouts,150,530,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,520,nPlanets+' planets',charsize=1.5,color=0,/device
;  xyouts,150,100,strcompress(string(abs(output*1000.))+' years'),charsize=1.5,color=0,/device
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;  xyouts,390,540,strcompress('       !4s!X < '+tau_flagStr[0]),charsize=1.0,color=70,/device ;ñ=!Z(00f1)os
;  xyouts,390,520,strcompress(tau_flagStr[0]+' < !4s!X < '+tau_flagStr[1]),charsize=1.0,color=160,/device
;  xyouts,390,500,strcompress(tau_flagStr[1]+' < !4s!X < '+tau_flagStr[2]),charsize=1.0,color=160,/device
;  xyouts,390,480,strcompress(tau_flagStr[2]+' < !4s!X < '+tau_flagStr[3]),charsize=1.0,color=160,/device
;  xyouts,390,460,strcompress(tau_flagStr[3]+' < !4s!X'),charsize=1.0,color=240,/device
;  device,get_decompose=old_decomposed,decomposed=0
;  loadct,0
;  oplot,xcero,ycero,psym=2
;  oplot,-xplanet,-yplanet,psym=1
  write_png,strcompress(dischargeEmisionFolder+'/TauTransition<1.png',/remove_all),tvrd(0,0,600,600,0,true=1)








