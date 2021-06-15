Pro Currents

  nPlanets=7                    ;Number of planets
  output=30                    ;The output
;  TauRanges=3

  if nPlanets eq 5 then npStr='five'
  if nPlanets eq 6 then npStr='six'
  if nPlanets eq 7 then npStr='seven'
  diskcase='alpha_visc-'+npStr+'_planets'		;The case
  outputinterval=1000.		;In years
  AUcgs=1.49597871e13
  Ycgs=31536000.0
  MstarMsolar=1.87
  Msolarcgs=1.9891e33

  dataFolder='~/fargo3d_outputs/V1247_Ori-'+diskcase
  previus_dischargeImagesFolder='~/Dropbox/visualization_tools/images/V1247_Ori-'+diskcase
  file_mkdir,strcompress(previus_dischargeImagesFolder,/remove_all)
  file_mkdir,strcompress(previus_dischargeImagesFolder+'/output'+string(output),/remove_all)
  dischargeImagesFolder=previus_dischargeImagesFolder+'/output'+string(output)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Cells numbers
  dims=dblarr(8)
  dims_=strcompress(dataFolder+'/dims.dat',/remove_all)
  openr,1,dims_
  readf,1,dims
  close,1
  nsec = uint(dims(7))
  nrad = uint(dims(6))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Number of planets
  ;; numberAsString = STRTRIM(k, 2)
  ;; spawn,'ls -l out'+numberAsString+'/planet*.dat | wc -l',planets
  ;; for i=2,4 do begin
  ;;    a=STRTRIM(i, 2)
  ;;    if (a eq planets) then planets=i
  ;; endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Physical quantities arrays 
  acret=dblarr(nsec,nrad)
  DensMaxima=dblarr(nrad)
  tot=dblarr(nrad)
  radii=dblarr(nsec,nrad)
  CurrentsArray=dblarr(nsec,nrad)
  CurrentsWakel=dblarr(nsec)
  CurrentsWaker=dblarr(nsec,nrad)
  acretion=dblarr(nrad)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Coloca los espacios radiales en forma de arreglos
  rad=dblarr(nrad+1)
  rad_=strcompress(dataFolder+'/used_rad.dat',/remove_all)
  openr,1,rad_
  readf,1,rad
  close,1
  rad=rad;/AUcgs
  rmin=rad(0)
  rmax=rad(nrad)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Coloca los espacios radiales en forma de arreglos
  planet=dblarr(10,output+1)
  xcero=dblarr(1)
  ycero=dblarr(1)
  xplanet=dblarr(nPlanets)
  yplanet=dblarr(nPlanets)
  vxplanet=dblarr(nPlanets)
  vyplanet=dblarr(nPlanets)
  mplanet=dblarr(nPlanets)
  Rplanet=dblarr(nPlanets)
  Hillplanet=dblarr(nPlanets)
  ThetaPlanet=dblarr(nPlanets)
  alpha=dblarr(nsec)
  xHill=dblarr(nPlanets,nsec)
  yHill=dblarr(nPlanets,nsec)
  HillplanetTheta=dblarr(nPlanets,nsec)
  HillplanetR=dblarr(nPlanets,nsec)
  gapBorder=dblarr(2)
;  planets=1
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

;     print,xplanet[n],yplanet[n],rad

     
     ;; for i=0,nrad-1 do begin
     ;;    if (xplanet[n] gt rad[i]) and (xplanet[n] lt rad[i+1]) then begin
     ;;       print,i
     ;;       ;; for j=0,nsec-1 do begin
     ;;       ;;    if (yplanet[n] gt rad[i]) and (yplanet[n] lt rad[i+1]) then print,j
     ;;       ;; endfor
     ;;    endif
     ;; endfor
        x=dblarr(nsec,nrad)
        y=dblarr(nsec,nrad)
     rmed = (rad(1:nrad)-rad(0:nrad-1))/2+rad(0:nrad-1)
        for ii=0,nrad-1 do begin
           for jj=0,nsec-1 do begin
              x[jj,ii]=rmed[ii]*cos(alpha[jj])
              y[jj,ii]=rmed[ii]*sin(alpha[jj])              
           endfor
        endfor

     mplanet[n]=planet(7,output)*1.86;/Msolarcgs;/MstarMsolar
print,'mplanet=',mplanet[n],n
     Rplanet[n]=sqrt(xplanet[n]^2+yplanet[n]^2)
;     Rplanet[n]=a
     Hillplanet[n]=a*((1./3)*mplanet[n])^(1./3.);Rplanet[n]*((1./3)*mplanet[n])^(1./3.)
     ThetaPlanet[n]=atan(yplanet[n],xplanet[n])+(1.0/nsec*(2.*!PI))+!PI
     if (ThetaPlanet[n] < 0.) then ThetaPlanet[n]=!PI+ThetaPlanet[n];+(1.0/nsec*2*!PI)
;     xplanet[n]=(Rplanet[n]*cos(ThetaPlanet[n]))
;     yplanet[n]=(Rplanet[n]*sin(ThetaPlanet[n]))

	distPlanetCell=dblarr(nsec,nrad)
;        for i=0,nrad-1 do begin
;           for j=0,nsec-1 do begin
;              xdiff=rmed[i]-xPlanet[n]
;              ydiff=rmed[i]-yPlanet[n] 
;              distPlanetCell[j,i]=sqrt((xdiff*xdiff)+(ydiff*ydiff))
;           endfor
;        endfor

;	dens_max=0
;        for i=0,nrad-1 do begin
;           for j=0,nsec-1 do begin
;              if (distPlanetCell[j,i] le Hillplanet[n]) and (n ge 50) then begin;
;		if dens_max lt dens[j,i] then begin;
;		  dens_max=dens[j,i]	
;		  xPlanet[n]=rmed[i]*cos(alpha[j])
;		  yPlanet[n]=rmed[i]*sin(alpha[j])
;		endif
;	      endif
;           endfor
;        endfor
;     Rplanet[n]=sqrt(xplanet[n]^2+yplanet[n]^2)
;     ThetaPlanet[n]=atan(yplanet[n],xplanet[n])
;     if (ThetaPlanet[n] < 0.) then ThetaPlanet[n]=2.*!PI+ThetaPlanet[n]
;     xplanet[n]=planet(1,output)
;     yplanet[n]=planet(2,output)
;     vxplanet[n]=planet(3,output)
;     vyplanet[n]=planet(4,output)
     h=xPlanet[n]*vyPlanet[n]-yPlanet[n]*vxPlanet[n]
     d=sqrt(xPlanet[n]*xPlanet[n]+yPlanet[n]*yPlanet[n])
     Ax=xPlanet[n]*vyPlanet[n]*vyPlanet[n]-yPlanet[n]*vxPlanet[n]*vyPlanet[n]-xPlanet[n]/d
     Ay=yPlanet[n]*vxPlanet[n]*vxPlanet[n]-xPlanet[n]*vxPlanet[n]*vyPlanet[n]-yPlanet[n]/d
     ee=sqrt(Ax*Ax+Ay*Ay)
     a=h*h/(1-ee*ee)
     Hillplanet[n]=a*((1./3)*mplanet[n])^(1./3.)
     Rplanet[n]=sqrt(xplanet[n]^2+yplanet[n]^2)
     for j=0,nsec-1 do begin
        Alpha[j]=((j+0.5)/nsec)*2*!PI
        xHill[n,j]=Rplanet[n]*cos(ThetaPlanet[n])+Hillplanet[n]*cos(alpha[j]);Rplanet[n]
        yHill[n,j]=Rplanet[n]*sin(ThetaPlanet[n])+Hillplanet[n]*sin(alpha[j])
        HillplanetTheta[n,j]=atan(yHill[n,j],xHill[n,j])
        if (HillplanetTheta[n,j] < 0.) then HillplanetTheta[n,j]=$
           2.*!PI+HillplanetTheta[n,j]
;	HillplanetTheta[n,j]+=!PI
;        if HillplanetTheta[n,j] gt 2*!PI then HillplanetTheta[n,j]-=2*!PI
        HillplanetR[n,j]=sqrt(yHill[n,j]^2+xHill[n,j]^2)
     endfor
;  endfor
;  gapBorder[0]=5.
;  gapBorder[1]=20.
     if n eq 0 then gapBorder[0]=Rplanet[n]-(2*Hillplanet[n])
     if n eq nPlanets-1 then gapBorder[1]=Rplanet[n]+(2*Hillplanet[n])
     for i=0,nrad-1 do begin
        if (Rplanet[n] gt rad[i]) and (Rplanet[n] lt rad[i+1]) then begin
;print,ThetaPlanet[n],Rplanet[n]
           ;; for j=0,nsec-1 do begin
           ;;    Alpha[j]=((j+0.5)/nsec)*2*!PI
           ;;    if (ThetaPlanet[n] gt Alpha[j]) then print,i,j
           ;; endfor
        endif
     endfor
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Coloca la densidad en los arreglos
  dens=dblarr(nsec,nrad)

  dens_=strcompress(dataFolder+'/gasdens'+string(output)+'.dat'$
                    ,/remove_all)
  openr,1,dens_
  readu,1,dens
  close,1
  dens=dens*1.86  ;to recover to units of M_solar/AU^2
;Ciclos para calculo de densidad media
print,dens(0,*)
  densMed=dblarr(nrad)
  tot=dblarr(nrad)
  for r=0,nrad-1 do begin
     tot[r]=0.0
     for l=0,nsec-1 do begin
        tot[r]=dens[l,r]+tot[r]
     endfor
  endfor
  for r=0,nrad-1 do begin
     densMed[r]=tot[r]/nsec
  endfor

;Ciclos para calculo de densidad maxima
  wakel=dblarr(nrad)  
  waker=dblarr(nrad)

  Maxima=0.
  for r=0,nrad-1 do begin
     Maxima = max(dens(*,r))
     for l=0,nsec-1 do begin
        if (dens[l,r] eq Maxima) and (output gt 0) and (r ne 255) then begin
           densMaxima[r]=dens[l,r]
           wakel[r]=l*(2.*!PI)/(nsec-1)
           waker[r]=r*(rmax-rmin)/(nrad-1)+rmin
        endif
     endfor
  endfor

;Ciclos para calculo de arreglos para calculos de corrientes
  for r=0,nrad-1 do begin
     for l=0,nsec-1 do begin
        if (dens[l,r] gt densMed[r]) and (output gt 0) and (r ne 255) then begin
           CurrentsArray[l,r]=dens[l,r]
           CurrentsWakel[l]=l*(2.*!PI)/(nsec-1)
           CurrentsWaker[l,r]=r*(rmax-rmin)/(nrad-1)+rmin
        endif
     endfor
  endfor


;Coloca la velocidad radial en los arreglos
  vrad=dblarr(nsec,nrad)
  vrad_=strcompress(dataFolder+'/gasvy'+string(output)+'.dat'$
                    ,/remove_all)
  openr,1,vrad_
  readu,1,vrad
  close,1

  for i=0,nrad-1 do begin
     tot[i]=0.
     for j=0,nsec-1 do begin
        acret[j,i]=(2.*!PI/nsec)*rad[i]*dens[j,i]*vrad[j,i]
        if (acret[j,i] < 0.) then begin
           tot[i]=acret[j,i]+tot[i]
        endif
     endfor
     acretion[i]=tot[i]
  endfor

;Coloca la velocidad azimutal en los arreglos
  vtheta=dblarr(nsec,nrad)
  vtheta_=strcompress(dataFolder+'/gasvx'+string(output)+$
                      '.dat',/remove_all)
  openr,1,vtheta_
  readu,1,vtheta
  close,1
rmed = (rad(1:nrad)-rad(0:nrad-1))/2+rad(0:nrad-1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Transporta celdas, de coordenadas cartesianas 2D a 
;coordenadas cilíndricas 2D
  x=dblarr(nsec,nrad)
  y=dblarr(nsec,nrad)
  for i=0,nrad-1 do begin
     for j=0,nsec-1 do begin
        alpha[j]=((j+0.5)/nsec)*2*!PI
        if (j eq 0) then  alpha[j]=0.001
        if (j eq nsec-1) then  alpha[j]=2*!PI
        x[j,i]=rmed[i]*cos(alpha[j]);rad[i]
        y[j,i]=rmed[i]*sin(alpha[j])
     endfor
  endfor

  densMedTotal=total(dens)/(float(nsec)*float(nrad))
  densmax=max(dens)
  densmin=min(dens)
  acretmax=max(acret)
  acretmin=min(acret)

;-----------------
;0) ECCENTRICITY
;-----------------
;	Additional arrays declarations
  Ac = dblarr(nsec)
  Bc = dblarr(nsec)
  alpha = dblarr(nsec)
  vrmod = dblarr(nrad)
  e = dblarr(nrad)
  vr = total(vrad,1)/nsec
  vt = total(vtheta,1)/nsec
  for j=0,nrad-1 do begin
     for i=0,nsec-1 do begin
        alpha[i]=((i+0.5)/nsec)*2*!PI
        Ac[i]=vrad[i,j]*cos(alpha[i])*(alpha[2]-alpha[1])/!PI
        Bc[i]=vrad[i,j]*sin(alpha[i])*(alpha[2]-alpha[1])/!PI
     endfor
     vrmod[j]=total(sqrt((Ac*Ac)+(Bc*Bc)))
     e[j]=vrmod[j]/vt[j]
  endfor



  rmed = (rad(1:nrad)-rad(0:nrad-1))/2+rad(0:nrad-1)
  omegaframe = planet(9,output)
  vphi = vtheta
  density = dens
;--------------
;1) VORTICITY: curl(v) and VORTENSITY: curl(v)/Sigma
;--------------
;	Additional arrays declarations
  drrvphi = dblarr(nsec,nrad)      ;\partial{r v_phi} / \partial r
  dphivr = dblarr(nsec,nrad)       ;\partial{v_r}\includegraphics[width=0.3\textwidth]{/home/ramiro/Dropbox/visualization_tools/images/Currents/f10_p2/output300/TauTransition.png} / \partial phi
  rotv = dblarr(nsec,nrad)         ;rotv = (1/r) * (drrvphi - dphivr)
;	We will need to correct the azimuthal velocity from the frame's
;	angular velocity
  dphi = 2.0*!dpi/(nsec+0.0)
  vphi_t = vphi + omegaframe[0]*rmed##replicate(1.0,nsec)
;	We first calculate drrvphi
  for j=0, nsec-1 do begin
     for i=1, nrad-1 do begin
        drrvphi(j,i) = ( rmed(i)*vphi_t(j,i) - rmed(i-1)*vphi_t(j,i-1) ) / (rmed(i) - rmed(i-1) )
     endfor
     drrvphi(j,0) = drrvphi(j,1)
  endfor
;	Then we calculate dphivr
  for j=0, nsec-1 do begin
     if (j EQ 0) then jm1 = nsec-1
     if (j GT 0) then jm1 = j-1
     for i=0, nrad-1 do begin
        dphivr(j,i) = (vrad(j,i)-vrad(jm1,i))/dphi
     endfor
  endfor
;	We deduce the gas vorticity
  rotv(0:nsec-1,0:nrad-1) = (drrvphi(0:nsec-1,0:nrad-1)-dphivr(0:nsec-1,0:nrad-1)) / (rad(0:nrad-1)##replicate(1.0,nsec))
;	Centering in cells :
  vorticity=interpolate(rotv,indgen(nsec)+0.5,indgen(nrad)+0.5,/grid)
;	If one wishes to calculate the gas vortensity at same output number
  vortensity = vorticity/density
; --------------
; 2) Bernouilli
; --------------
  azimuths=(dindgen(nsec))/nsec*2*!dpi
  az=replicate(1,nrad)##azimuths
  rmedt=rmed##replicate(1,nsec)
  smooth = 0.6
  h = 0.05
  flare = 0.25
  r_array=rmed##replicate(1,nsec)
  temperature=h*h*r_array^(1+flare-1.5)
  pressure=temperature*density
  entropy=pressure/(density^1.4)
  soft=h*smooth
  dist = sqrt((xplanet[0]-rmedt*cos(az))^2. + (yplanet[0]-rmedt*sin(az))^2. + soft*soft)
  mp=mplanet[0]
  potp=-1.*mp/dist
  pot=potp-1./rmedt
  enthalpy=1.4/0.4*temperature
  indpot=mp*cos(az)*rmedt
  vradbig = [vrad[*,*],vrad[0,*]]
  vthetabig = [vtheta[*,*],vtheta[0,*]]
  vrCentered=interpolate(vradbig,dindgen(nsec),dindgen(nrad)+0.5,/grid)
  vpCentered=interpolate(vthetabig,dindgen(nsec)+0.5,dindgen(nrad),/grid)
  bernouilli=(vrCentered^2+vpCentered^2)/2.+pot+indpot+enthalpy-(rmedt*omegaframe[0])^2/2

; --------------
; 3) Angular momentum
; --------------   
  rmedt = rmed##replicate(1,nsec)
  vphi_t = vphi + omegaframe[0]*rmed##replicate(1.0,nsec)
  angmom = dens*rmedt*vphi_t
                                ;Centering in cells :
  angmom=interpolate(angmom,indgen(nsec)+0.5,indgen(nrad),/grid)

; --------------
; 4) Adiabatic Invariant G
; --------------
  dr=rad(1:nrad-1)-rad(0:nrad-2)
  dr=dr##replicate(1,nsec)
  drs = dblarr(nsec,nrad)
  drb = dblarr(nsec,nrad)
  adinv = dblarr(nsec,nrad)
  drs=(entropy(*,1:nrad-1)-entropy(*,0:nrad-2))/dr
  drb=(bernouilli(*,1:nrad-1)-bernouilli(*,0:nrad-2))/dr
  adinv=(1-temperature*drs/drb)/vortensity



  Theta=dblarr(nsec)

  R=1
  ResidualVtheta1=dblarr(nsec)
  VelRadial1=dblarr(nsec)
  for i=0,nrad-1 do begin
     if rmed[i] lt R then begin
        for j=0,nsec-1 do begin
           ResidualVtheta1[j] = (VpCentered[j,i]-sqrt(1/rmed[i]))/sqrt(1/rmed[i])
           VelRadial1[j] = VrCentered[j,i]
           Theta[j] = ((j+0.5)/nsec)*2*!PI
        endfor
     endif
  endfor

  R=2
  ResidualVtheta3=dblarr(nsec)
  VelRadial3=dblarr(nsec)
  for i=0,nrad-1 do begin
     if rmed[i] lt R then begin
        for j=0,nsec-1 do begin
           ResidualVtheta3[j] = (VpCentered[j,i]-sqrt(1/rmed[i]))/sqrt(1/rmed[i])
           VelRadial3[j] = VrCentered[j,i]
           Theta[j] = ((j+0.5)/nsec)*2*!PI
        endfor
     endif
  endfor

  R=3
  ResidualVtheta5=dblarr(nsec)
  VelRadial5=dblarr(nsec)
  for i=0,nrad-1 do begin
     if rmed[i] lt R then begin
        for j=0,nsec-1 do begin
           ResidualVtheta5[j] = (VpCentered[j,i]-sqrt(1/rmed[i]))/sqrt(1/rmed[i])
           VelRadial5[j] = VrCentered[j,i]
           Theta[j] = ((j+0.5)/nsec)*2*!PI
        endfor
     endif
  endfor

;Define color table\includegraphics[width=0.3\textwidth]{/home/ramiro/Dropbox/visualization_tools/images/Currents/f10_p2/output300/TauTransition.png}
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
;  window,0,xsize=600,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  plot,Theta,ResidualVtheta1,background=255,$                            
       color=0,xtitle='angulo',yrange=[-.08,.08],psym=0,xstyle=1,$
       ytitle='( Mom.Angular (gas) - Mom.Angular (Vk) ) / Mom.Angular (Vk)',charsize=1.5
  xyouts,500,530,'1 UA',charsize=1.5,color=0,/device
  oplot,Theta,ResidualVtheta3,color=100,psym=0
  xyouts,500,500,'2 UA',charsize=1.5,color=100,/device
  oplot,Theta,ResidualVtheta5,color=200,psym=0
  xyouts,500,470,'3 UA',charsize=1.5,color=200,/device
  write_png,strcompress(dischargeImagesFolder+'/ResidualMomentoAngular.png',/remove_all),tvrd(0,0,600,600,0,true=1)

;Define color table
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
;  window,1,xsize=600,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  plot,Theta,VelRadial1,background=255,$                            
       color=0,xtitle='angulo',yrange=[-0.15,0.15],psym=0,xstyle=1,$
       ytitle='Velocidad Radial',charsize=1.5
  xyouts,500,530,'1 UA',charsize=1.5,color=0,/device
  oplot,Theta,VelRadial3,color=100,psym=0
  xyouts,500,500,'2 UA',charsize=1.5,color=100,/device
  oplot,Theta,VelRadial5,color=200,psym=0
  xyouts,500,470,'3 UA',charsize=1.5,color=200,/device
  write_png,strcompress(dischargeImagesFolder+'/ResidualRadialVel.png',/remove_all),tvrd(0,0,600,600,0,true=1)



;---------------
;Cuantities on x-axis
  xrange_prev = rmed
;xrange_prev = alog10(e)

  xtitle_prev = 'Semieje mayor (UA)'
;xtitle_prev = 'Excentricidad'

;---------------
;Cuantities on y-axis
;yrange_prev = rmed
  yrange_prev = e


;ytitle_prev = 'Semieje mayor (UA)'
  ytitle_prev = 'Excentricidad'

;Define color table
  device,get_decompose=old_decomposed,decomposed=0
  loadct,0
;; window,1,xsize=600,ysize=600
;; plot,xrange_prev,yrange_prev,background=255,$                            
;; color=0,xtitle=xtitle_prev,psym=0,xstyle=1,$
;; ytitle=ytitle_prev,/ylog,$            
;; charsize=2.0

;Define tabla de colores
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3

  accretion = (2.*!PI/nsec)*dens*rmedt*VrCentered
  densCGS = dens*2./2.25*1.e7
  tau=9.6*densCGS

;zrange_prev = accretion         ;Msolar/year
;zrange_prev = vortensity
;zrange_prev = vorticity
zrange_prev = densCGS
;zrange_prev = alog10(angmom)
;  zrange_prev = alog10(tau)   ;Adimentional\includegraphics[width=0.3\textwidth]{/home/ramiro/Dropbox/visualization_tools/images/Currents/f10_p2/output300/TauTransition.png}
;zrange_prev = temperature
;zrange_prev = pressure
;zrange_prev = alog10(densCGS)
;zrange_prev = alog10(accretion)
;zrange_prev = alog10(angmom)
;zrange_prev = VrCentered        ;UA/year
;zrange_prev = VpCentered        ;UA/year
;zrange_prev = alog10(VpCentered)        ;UA/year
;zrange_prev = VthetaResidual

;Str_zrange_prev  = 'zrange_prev'
;print,Str_zrange_prev


;  window,2,xsize=600,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  contour,alog10(zrange_prev),x,y,/fill,charsize=1.5,$
          xrange=[-60,60],yrange=[-60,60],nlevels=300,$
	/xstyle,/ystyle,xtitle='AU',ytitle='AU',min_value=-5.,max_value=0

;Accretion of gas (minus sign)
  device,get_decompose=old_decomposed,decomposed=0
  loadct,33
  contour,accretion,x,y,levels=[-1e-9,-6.6e-10,-3.3e-10,-1e-10,-6.6e-11,-3.3e-11,-1e-11],c_colors=[50,80,110,140,170,200,230],/overplot
  xyouts,360,550,'Acc = -1.0e-11',charsize=1.5,color=230,/device
  xyouts,360,530,'Acc = -3.3e-11',charsize=1.5,color=200,/device
  xyouts,360,510,'Acc = -6.6e-11',charsize=1.5,color=170,/device
  xyouts,360,490,'Acc = -1.0e-10',charsize=1.5,color=140,/device
  xyouts,360,470,'Acc = -3.3e-10',charsize=1.5,color=110,/device
  xyouts,360,450,'Acc = -6.6e-10',charsize=1.5,color=80,/device
  xyouts,360,430,'Acc = -1.0e-9',charsize=1.5,color=50,/device

;Accretion of gas (positive sign)
;; device,get_decompose=old_decomposed,decomposed=0
;; loadct,33
;; contour,accretion,x,y,levels=[1e-11,3.3e-11,6.6e-11,1e-10,3.3e-10,6.6e-10,1e-9],c_colors=[50,80,110,140,170,200,230],/overplot
;; xyouts,440,550,'Acc = 1e-9',charsize=1.5,color=230,/device
;; xyouts,440,530,'Acc = 6.6e-10',charsize=1.5,color=200,/device
;; xyouts,440,510,'Acc = 3.3e-10',charsize=1.5,color=170,/device
;; xyouts,440,490,'Acc = 1e-10',charsize=1.5,color=140,/device
;; xyouts,440,470,'Acc = 6.6e-11',charsize=1.5,color=110,/device
;; xyouts,440,450,'Acc = 3.3e-10',charsize=1.5,color=80,/device
;; xyouts,440,430,'Acc = 1e-11',charsize=1.5,color=50,/device


  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
  oplot,xcero,ycero,psym=2
;  oplot,-xplanet,-yplanet,psym=1
  for n=0,nPlanets-1 do begin
     oplot,xHill(n,*),yHill(n,*),color=255
  endfor
  colscale=dindgen(101)-50
  mm=[min(zrange_prev),max(zrange_prev)]
  contour,transpose(replicate(1,2))##colscale,colscale,dindgen(2),nlev=101,levels=colscale,/fill,pos=[[0.3,0.95],[0.8,0.98]],/noerase, xticklen=0.1, xstyle=1,yticks=1,yminor=1,ystyle=1,xcharsize=0.001,ycharsize=0.001,font=-1
  xyouts,0.295,0.96,strcompress(string(mm(0))),/normal,align=1
  xyouts,0.8,0.96,strcompress(string(mm(1))),/normal,align=0.
  write_png,strcompress(dischargeImagesFolder+'/Accretion-Alog10(tau).png',/remove_all),tvrd(0,0,600,600,0,true=1)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Define tabla de colores
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
;produce gráfica (Currents)
;  window,3,xsize=600,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  contour,alog10(dens*2./2.25*1.e7),x,y,/fill,charsize=1.5,$
          xrange=[-60,60],yrange=[-60,60],nlevels=300,$
          /xstyle,/ystyle,xtitle='UA',ytitle='UA',$
          min_value=-3.,max_value=3.
  oplot,xcero,ycero,psym=1.5
  ;; oplot,xplanet,yplanet,psym=1
  ;; for n=0,nPlanets-1 do begin
  ;;    oplot,xHill(n,*),yHill(n,*),color=255
  ;; endfor
  ;; xyouts,100,100,strcompress('t = '+string(abs(output*outputinterval))+' a!Z(00f1)os'),charsize=1.5,color=0,/device
  ;; device,get_decompose=old_decomposed,decomposed=0
  ;; loadct,39
  ;; contour,alog10(dens),x,y,levels=[-3.,-2.,-1,-0.602059,-0.301029,-0.124938,0.,1.,2.],/overplot,$
  ;;         xrange=[-30,30],yrange=[-30,30],$
  ;;         /xstyle,/ystyle,xtitle='AU',ytitle='UA',$
  ;;         min_value=-3.,max_value=2.4,c_colors=[50,75,100,125,150,175,200,225,250]
  ;; xyouts,440,550,'!N!4R!3!N=100 g cm!E-2!N',charsize=1.5,color=250,/device
  ;; xyouts,440,530,'!N!4R!3!N=10 g cm!E-2!N',charsize=1.5,color=225,/device
  ;; xyouts,440,510,'!N!4R!3!N=1 g cm!E-2!N',charsize=1.5,color=200,/device
  ;; xyouts,440,490,'!N!4R!3!N=0.75 g cm!E-2!N',charsize=1.5,color=175,/device
  ;; xyouts,440,470,'!N!4R!3!N=0.5 g cm!E-2!N',charsize=1.5,color=150,/device
  ;; xyouts,440,450,'!N!4R!3!N=0.25 g cm!E-2!N',charsize=1.5,color=125,/device
  ;; xyouts,440,430,'!N!4R!3!N=0.1 g cm!E-2!N',charsize=1.5,color=100,/device
  ;; xyouts,440,410,'!N!4R!3!N=0.01 g cm!E-2!N',charsize=1.5,color=75,/device
  ;; xyouts,440,390,'!N!4R!3!N=0.001 g cm!E-2!N',charsize=1.5,color=50,/device
  write_png,strcompress(dischargeImagesFolder+'/Currents.png',/remove_all),tvrd(0,0,600,600,0,true=1)


  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;; produce gráfica (Acretion)
;; window,1,xsize=600,ysize=600
;; !P.MULTI = [0,1,1,0,0]
;; contour,acret,x,y,nlevels=2000,xrange=[-30,30],$
;; yrange=[-10,10],/xstyle,/ystyle,xtitle='AU',ytitle='UA',$
;; c_colors=[0]
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=-9.99999999e-8,max_value=-1.e-8,c_colors=[30],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=-9.99999999e-9,max_value=-1.e-9,c_colors=[60],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=-9.99999999e-10,max_value=-1.e-10,c_colors=[90],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=-9.99999999e-11,max_value=-1.e-11,c_colors=[120],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=-9.99999999e-12,max_value=-1.e-12,c_colors=[135],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=1.e-8,max_value=9.99999999e-8,c_colors=[150],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=1.e-9,max_value=9.99999999e-9,c_colors=[180],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=1.e-10,max_value=9.99999999e-10,c_colors=[210],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=1.e-11,max_value=9.99999999e-11,c_colors=[240],/overplot
;; contour,acret,x,y,/cell_fill,nlevels=5000,$
;; min_value=1.e-12,max_value=9.99999999e-12,c_colors=[250],/overplot
;; xyouts,10,500,'-e-8',charsize=1.,color=30,/device
;; xyouts,10,480,'-e-9',charsize=1.,color=60,/device
;; xyouts,10,460,'-e-10',charsize=1.,color=90,/device
;; xyouts,10,440,'-e-11',charsize=1.,color=120,/device
;; xyouts,10,420,'-e-12',charsize=1.,color=135,/device
;; xyouts,10,400,'e-8',charsize=1.,color=150,/device
;; xyouts,10,380,'e-9',charsize=1.,color=180,/device
;; xyouts,10,360,'e-10',charsize=1.,color=210,/device
;; xyouts,10,340,'e-11',charsize=1.,color=240,/device
;; xyouts,10,320,'e-12',charsize=1.,color=250,/device
;; oplot,xcero,ycero,psym=2
;; oplot,xplanet,yplanet,psym=1
;; for n=0,planets-1 do begin
;;   if (n eq 0) then begin
;;     oplot,xHill(n,*),yHill(n,*),color=255
;;   endif
;;   if (n eq 1) then begin
;;     oplot,xHill(n,*),yHill(n,*),color=255
;;   endif
;; endfor
;; write_png,strcompress('out'+string(k)+'/out'+string(k)+$
;; 'output'+string(m)+'Acretion.png',/remove_all),tvrd(0,0,600,600,0,true=1)

;Radius of circunference
;; r=dblarr(10)
;; r[0]=10. 
;; r[1]=15.
;; r[2]=30.

  rx=dblarr(10,nsec)
  ry=dblarr(10,nsec)
  for i=0,1 do begin
     for j=0,nsec-1 do begin
        alpha[j]=((j+0.5)/nsec)*2*!PI
        if (j eq 0) then  alpha[j]=0.001
        if (j eq nsec-1) then  alpha[j]=2*!PI
        rx[i,j]=gapBorder[i]*cos(alpha[j])
        ry[i,j]=gapBorder[i]*sin(alpha[j])
     endfor
  endfor



  loadct,0
;Calculo de profundidad optica con Xi=9.6 cm^2 g^⁻1
  tau=9.6*dens*2./2.25*1.e7
;print,min(alog10(tau)),max(alog10(tau)),max(dens*2./2.25*1.e7),min(dens*2./2.25*1.e7)
;Gráficas de profundidad óptica
;  window,4,xsize=600,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  contour,alog10(tau),x,y,charsize=1.5,levels=[-0.30103,0.30103],$
          c_color=[100,150],$
          xrange=[-50,50],yrange=[-50,50],$
          /xstyle,/fill,/ystyle,xtitle='AU',ytitle='AU' ;,$
  nPlanets=strtrim(nPlanets,2)
;  xyouts,150,530,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,520,nPlanets+' planets',charsize=1.5,color=0,/device
;  xyouts,150,100,strcompress(string(abs(output*1000.))+' years'),charsize=1.5,color=0,/device
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
  xyouts,390,530,'       !4s!X < 0.5',charsize=1.5,color=70,/device ;ñ=!Z(00f1)os
  xyouts,390,500,'0.5 < !4s!X < 2.0',charsize=1.5,color=160,/device
  xyouts,390,470,'2.0 < !4s!X',charsize=1.5,color=240,/device
  device,get_decompose=old_decomposed,decomposed=0
  loadct,0
  oplot,xcero,ycero,psym=2
   oplot,-xplanet,-yplanet,psym=1
  ;; for n=0,nPlanets-1 do begin
  ;;    if (n eq 0) then begin
  ;;       oplot,xHill(n,*),yHill(n,*),color=255
  ;;    endif
  ;;    if (n eq 1) then begin
  ;;       oplot,xHill(n,*),yHill(n,*),color=255
  ;;    endif
  ;;    if (n eq 2) then begin  
  ;;       oplot,xHill(n,*),yHill(n,*),color=255
  ;;    endif
  ;;    if (n eq 3) then begin
  ;;       oplot,xHill(n,*),yHill(n,*),color=255
  ;;    endif
  ;; endfor
 ;;  for i=0,2 do begin
   ;;   oplot,rx(i,*),ry(i,*),color=255
 ;;  endfor
  write_png,strcompress(dischargeImagesFolder+'/TauTransition.png',/remove_all),tvrd(0,0,600,600,0,true=1)
;filename=strcompress(dischargeImagesFolder+'/TauTransition.eps',/remove_all)






;set_plot,'ps'
;DEVICE, /ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, $
;   FILENAME=filename, XSIZE=600, $
;   YSIZE=600
;device,/close

;;;;;;;;;;;;;;;;;;
; Analisis of Optical depth calculation 
;;;;;;;;;;;;;;;;;;
  TauRanges=6
  SigmaTau=dblarr(nrad,TauRanges)
  AccTau=dblarr(nrad,TauRanges)
  Counter=dblarr(nrad,Tauranges)
  rad=dblarr(nrad+1)  
  DensTau=dblarr(nsec,nrad,TauRanges)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Coloca los espacios radiales en forma de arreglos
  rad_=strcompress(dataFolder+'/used_rad.dat',/remove_all)
  openr,1,rad_
  readf,1,rad
  close,1
  rmin=rad(0)
  rmax=rad(nrad)

  loadct,3
;Analisis of Optical depth calculation with Xi=9.6 cm^2 g^⁻1
  tau=9.6*dens;*2./2.25*1.e7
;  print,min(alog10(tau)),max(alog10(tau)),max(dens*2./2.25*1.e7),min(dens*2./2.25*1.e7)
;Gráficas de profundidad óptica
;  window,5,xsize=600,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  contour,alog10(dens*2./2.25*1.e7),x,y,/fill,charsize=1.5,nlevels=300,$
          xrange=[-50,50],yrange=[-50,50],$
          /xstyle,/ystyle,xtitle='UA',ytitle='UA',$
          min_value=-3.,max_value=2.4
;  xyouts,100,100,strcompress('t = '+string(abs(output*1000.))+' a!Z(00f1)os'),charsize=1.5,color=0,/device
  loadct,39
;  contour,alog10(tau),x,y,levels=[-2.,-1.,-0.30103,0.,1.,2.,3],c_colors=[30,70,100,130,170,200,230],/overplot
  contour,alog10(tau),x,y,levels=[-0.3010299957,0.3010299957],c_colors=[50,100],/overplot
;; min_value=0.,max_value=100000.,c_colors=[250,210,180,130]
;  xyouts,470,530,'!N!4s!3!N=0.01',charsize=1.5,color=30,/device
;  xyouts,470,510,'!N!4s!3!N=0.1',charsize=1.5,color=70,/device
  xyouts,470,490,'!N!4s!3!N=0.5',charsize=1.5,color=50,/device
  xyouts,470,470,'!N!4s!3!N=2.0',charsize=1.5,color=100,/device
;  xyouts,470,450,'!N!4s!3!N=10',charsize=1.5,color=170,/device
;  xyouts,470,430,'!N!4s!3!N=100',charsize=1.5,color=200,/device
;  xyouts,470,410,'!N!4s!3!N=1000',charsize=1.5,color=230,/device
  loadct,0
  oplot,xcero,ycero,psym=2
  oplot,-xplanet,-yplanet,psym=1
  for n=0,nPlanets-1 do begin
     if (n eq 0) then begin
        oplot,xHill(n,*),yHill(n,*),color=255
     endif
     if (n eq 1) then begin
        oplot,xHill(n,*),yHill(n,*),color=255
     endif
     if (n eq 2) then begin
        oplot,xHill(n,*),yHill(n,*),color=255
     endif
     if (n eq 3) then begin
        oplot,xHill(n,*),yHill(n,*),color=255
     endif
;     print,xplanet,yplanet
  endfor
  write_png,strcompress(dischargeImagesFolder+'/Tau.png',/remove_all),tvrd(0,0,600,600,0,true=1)



;; window,3,xsize=600,ysize=600
;; !P.MULTI = [0,1,1,0,0]
;; contour,tau,x,y,xrange=[-10,10],yrange=[-10,10],$
;; /xstyle,/ystyle,xtitle='UA',ytitle='UA',levels=[1.,10.,100.,1000.],$
;; min_value=0.,max_value=100000.,c_colors=[130,96,80,0]
;; xyouts,470,530,'!N!4s!3!N=1',charsize=1.5,color=130,/device
;; xyouts,470,510,'!N!4s!3!N=10',charsize=1.5,color=96,/device
;; xyouts,470,490,'!N!4s!3!N=100',charsize=1.5,color=80,/device
;; xyouts,470,470,'!N!4s!3!N=1000',charsize=1.5,color=0,/device
;; oplot,xcero,ycero,psym=2
;; oplot,xplanet,yplanet,psym=1
;; for n=0,planets-1 do begin
;;   if (n eq 0) then begin
;;     oplot,xHill(n,*),yHill(n,*),color=255
;;   endif
;;   if (n eq 1) then begin
;;     oplot,xHill(n,*),yHill(n,*),color=255
;;   endif
;; endfor
;; write_png,strcompress('out'+string(k)+'/out'+string(k)+$
;; 'output'+string(m)+'TauGrande_to_10_micras.png',/remove_all),tvrd(0,0,600,600,0,true=1)

;nrad=256

;Construye sitema coordenado (r,phi) 
  xx=dblarr(nsec,nrad)
  yy=dblarr(nsec,nrad)
  for i=0,nrad-1 do begin
     for j=0,nsec-1 do begin
        dens[j,i]=dens[j,i]
        xx[j,i]=j*2*!PI/(long(nsec)-1) 
        yy[j,i]=i
        radii[j,i]=rad[i]
     endfor
  endfor
;Imprime mapa de densidad
  loadct,3
;  window,6,xsize=1000,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[1000,600],Set_Pixel_Depth=24
  contour,alog10(dens*2./2.25*1.e7),xx,radii,nlevels=300,/fill,/xstyle,/ystyle,$
          xrange=[0,2*!PI],/irregular,/ylog,ytitle='AU',xtitle='!N!4u!3!N',$
          charsize=2,min_value=min(alog10(dens*2./2.25*1.e7)),max_value=max(alog10(dens*2./2.25*1.e7)),color=0,BACKGROUND=255,thick=3
;  xyouts,100,100,strcompress('t = '+string(abs(output*outputinterval))+' a!Z(00f1)os'),charsize=1.5,color=255,/device
;  device,get_decompose=old_decomposed,decomposed=0
;  loadct,9
;  contour,accretion,xx,radii,levels=[-1e-11,-1e-12],c_colors=[0,0],c_linestyle=[1,2],/overplot
;levels=[-1e-9,-6.6e-10,-3.3e-10,-1e-10,-6.6e-11,-3.3e-11,-1e-11],c_colors=[50,80,110,140,170,200,230],/overplot
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
  oplot,HillplanetTheta,HillplanetR,color=255,psym=3
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
  contour,alog10(9.6*dens*2./2.25*1.e7),xx,radii,levels=[-0.3010299957,0.,0.3010299957],c_color=[40,70,100],thick=3,/overplot
  device,get_decompose=old_decomposed,decomposed=0
  loadct,0
  contour,alog10(9.6*dens*2./2.25*1.e7),xx,radii,levels=[-2.0,-1.0],c_color=[80,120],thick=3,/overplot
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
  oplot,HillplanetTheta,HillplanetR,color=255,psym=3
  colscale=dindgen(101)-50
  mm=[min(dens*2./2.25*1.e7),max(dens*2./2.25*1.e7)]
  contour,transpose(replicate(1,2))##colscale,colscale,dindgen(2),nlev=101,levels=colscale,/fill,pos=[[0.3,0.95],[0.8,0.98]],/noerase, xticklen=0.1, xstyle=1,yticks=1,yminor=1,ystyle=1,xcharsize=0.001,ycharsize=0.001,font=-1,color=0
  xyouts,0.295,0.96,string(mm(0))+'g/cm^2',/normal,align=1,color=0,charsize=1.5
  xyouts,0.8,0.96,strcompress(string(mm(1))+'g/cm^2'),/normal,align=0.,color=0,charsize=1.5
  write_png,strcompress(dischargeImagesFolder+'/Density.png',/remove_all),tvrd(0,0,1000,600,0,true=1)



  Set_plot,'z'
  Device,Set_resolution=[1000,600],Set_Pixel_Depth=24
  contour,alog10(dens*2./2.25*1.e7),xx,radii,nlevels=300,/fill,/xstyle,/ystyle,$
          xrange=[0,2*!PI],/irregular,/ylog,ytitle='AU',xtitle='!N!4u!3!N',$
          charsize=2,min_value=min(alog10(dens*2./2.25*1.e7)),max_value=max(alog10(dens*2./2.25*1.e7)),color=0,BACKGROUND=255
;  xyouts,100,100,strcompress('t = '+string(abs(output*outputinterval))+' a!Z(00f1)os'),charsize=1.5,color=255,/device
;  device,get_decompose=old_decomposed,decomposed=0
;  loadct,9
  contour,accretion,xx,radii,levels=[-1e-13],c_colors=[0],c_linestyle=[2],/overplot
;levels=[-1e-9,-6.6e-10,-3.3e-10,-1e-10,-6.6e-11,-3.3e-11,-1e-11],c_colors=[50,80,110,140,170,200,230],/overplot
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
  oplot,HillplanetTheta,HillplanetR,color=255,psym=3
  colscale=dindgen(101)-50
  mm=[min(dens*2./2.25*1.e7),max(dens*2./2.25*1.e7)]
  contour,transpose(replicate(1,2))##colscale,colscale,dindgen(2),nlev=101,levels=colscale,/fill,pos=[[0.3,0.95],[0.8,0.98]],/noerase, xticklen=0.1, xstyle=1,yticks=1,yminor=1,ystyle=1,xcharsize=0.001,ycharsize=0.001,font=-1,color=0
  xyouts,0.295,0.96,string(mm(0))+'g/cm^2',/normal,align=1,color=0,charsize=1.5
  xyouts,0.8,0.96,strcompress(string(mm(1))+'g/cm^2'),/normal,align=0.,color=0,charsize=1.5
  write_png,strcompress(dischargeImagesFolder+'/Density-AccretionIsocontour.png',/remove_all),tvrd(0,0,1000,600,0,true=1)
;  endfor


  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
;Calculo de profundidad optica con Xi=9.6 cm^2 g^⁻1
  tau=9.6*dens*2./2.25*1.e7
  Set_plot,'z'
  Device,Set_resolution=[600,600],Set_Pixel_Depth=24
  !P.MULTI = [0,1,1,0,0]
  contour,alog10(tau),x,y,charsize=1.5,levels=[alog10(0.1),alog10(0.3),alog10(0.5),alog10(0.75),alog10(1.0),alog10(1.25),alog10(1.5),alog10(1.75),alog10(2.0)],$
          c_color=[50,60,120,130,140,150,160,170,240],$
;          c_color=[50,100],$
          xrange=[-50,50],yrange=[-50,50],$
          /xstyle,/fill,/ystyle,xtitle='AU',ytitle='AU',color=0,BACKGROUND=255
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
  write_png,strcompress(dischargeImagesFolder+'/TauTransition<50.png',/remove_all),tvrd(0,0,600,600,0,true=1)



  dataFolderStructures='~/Dropbox/visualization_tools/images/V1247_Ori-alpha_visc-seven_planets/StructuresQuantification/'

  gasdens_temp=dblarr(nsec,nrad)
  gasdens_l_m=dblarr(nsec,nrad,TauRanges,15)
  mfin=dblarr(Tauranges)
  mfin[0]=7
  mfin[1]=4
  mfin[2]=3
  mfin[3]=2
  mfin[4]=2
  mfin[5]=1  
  for l=0,TauRanges-1 do begin
    for m=1,mfin[l] do begin
      gasdens_temp_=strcompress(dataFolderStructures+'/gasdens_l='+string(l)+'_m='+string(m),/remove_all)
      openr,1,gasdens_temp_
      readf,1,gasdens_temp
      close,1
      gasdens_l_m[*,*,l,m]=gasdens_temp
    endfor
  endfor
;print,gasdens_l_m(*,*,2,8)
;stop

;Construye sitema coordenado (r,phi) 
  xx=dblarr(nsec,nrad)
  yy=dblarr(nsec,nrad)
  for i=0,nrad-1 do begin
     for j=0,nsec-1 do begin
        dens[j,i]=dens[j,i]
        xx[j,i]=j*2*!PI/(long(nsec)-1) 
        yy[j,i]=i
        radii[j,i]=rad[i]
     endfor
  endfor
;Imprime mapa de densidad
  loadct,0
;  window,6,xsize=1000,ysize=600
  Set_plot,'z'
  Device,Set_resolution=[1000,600],Set_Pixel_Depth=24
;  contour,alog10(dens*2./2.25*1.e7),xx,radii,/cell_fill,/xstyle,/ystyle,$;levels=[-0.30103,0.30103],$
  gasdens_l_m(*,*,0,0)=-6.
;  contour,gasdens_l_m(*,*,0,0),xx,radii,/fill,/xstyle,/ystyle,$
  contour,alog10(9.6*dens*2./2.25*1.e7),xx,radii,/xstyle,/ystyle,/fill,levels=[-6.],$;/fill,$
          xrange=[0,2*!PI],ytitle='UA',xtitle='!N!4u!3!N',c_color=[100],$
          charsize=2,/ylog,yrange=[1.,100.],thick=2,color=0,BACKGROUND=255
;  xyouts,100,100,strcompress('t = '+string(abs(output*1000.))+' a!Z(00f1)os'),charsize=1.5,color=0,/device
;  oplot,HillplanetTheta,HillplanetR,color=255,psym=3
;  contour,alog10(tau),xx,radii,levels=[-0.30103,0.30103],c_colors=[100,200],/overplot

  for l=0,TauRanges-1 do begin
    if l eq 0 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,0
      for m=1,mfin[l] do begin
       contour,alog10(9.6*gasdens_l_m(*,*,l,m)*2./2.25*1.e7),xx,radii,/cell_fill,levels=[-6.],c_color=[80],/irregular,/overplot
      endfor
    endif
    if l eq 1 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,0
      for m=1,mfin[l] do begin
        contour,alog10(9.6*gasdens_l_m(*,*,l,m)*2./2.25*1.e7),xx,radii,/cell_fill,levels=[-6.],c_color=[100],/irregular,/overplot
      endfor
    endif
    if l eq 2 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,0
      for m=1,mfin[l] do begin
        contour,alog10(9.6*gasdens_l_m(*,*,l,m)*2./2.25*1.e7),xx,radii,/cell_fill,levels=[-6.],c_color=[120],/irregular,/overplot
      endfor
    endif
    if l eq 3 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      for m=1,mfin[l] do begin
        contour,alog10(9.6*gasdens_l_m(*,*,l,m)*2./2.25*1.e7),xx,radii,/cell_fill,levels=[-10.],c_color=[40],/irregular,/overplot;39,70;9,160;39,250
      endfor
    endif

    if l eq 4 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      for m=1,mfin[l] do begin
        contour,alog10(9.6*gasdens_l_m(*,*,l,m)*2./2.25*1.e7),xx,radii,/cell_fill,levels=[-10.],c_color=[60],/irregular,/overplot;39,70;9,160;39,250
      endfor
    endif
    if l eq 5 then begin 
      device,get_decompose=old_decomposed,decomposed=0
      loadct,39
      for m=1,mfin[l] do begin
        contour,alog10(9.6*gasdens_l_m(*,*,l,m)*2./2.25*1.e7),xx,radii,/cell_fill,levels=[-10.],c_color=[100],/irregular,/overplot;39,70;9,160;39,250
      endfor
    endif
  endfor
;  contour,alog10(gasdens_l_m(*,*,2,2)*2./2.25*1.e7),xx,radii,levels=[0.30103],c_colors=[100],/overplot
;  contour,alog10(gasdens_l_m(*,*,2,1)*2./2.25*1.e7),xx,radii,/cell_fill,c_color=[250],/overplot  
  device,get_decompose=old_decomposed,decomposed=0
  loadct,3
  oplot,HillplanetTheta,HillplanetR,color=255,psym=3
  gapBorderIn=dblarr(nsec)
  gapBorderIn(*)=gapBorder[0]
  gapBorderOut=dblarr(nsec)
  gapBorderOut(*)=gapBorder[1]
print,gapborder[1]
;  oplot,xx(*,0),gapBorderIn,color=255,psym=0,thick=2
;  oplot,xx(*,0),gapBorderOut,color=255,psym=0,thick=2
  write_png,strcompress(dischargeImagesFolder+'/Structures.png',/remove_all),tvrd(0,0,1000,600,0,true=1)





;Imprime mapa de corrientes
;; loadct,3
;; window,4,xsize=1000,ysize=600
;; contour,vrad,xx,radii,nlevels=300,/fill,/xstyle,/ystyle,$
;; xrange=[0,2*!PI],/irregular,/ylog,ytitle='UA',xtitle='!N!4u!3!N',$
;; charsize=1.5,min_value=-1.,max_value=1.
;; loadct,18
;; contour,vrad,xx,radii,nlevels=300,/fill,/xstyle,/ystyle,$
;; xrange=[0,2*!PI],/irregular,/ylog,$
;; charsize=1.5,min_value=-1.,max_value=1.,/overplot
;; loadct,3
;; xyouts,100,100,strcompress('t = '+string(abs(m*1000.))+' a!Z(00f1)os'),charsize=1.5,color=0,/device
;; oplot,HillplanetTheta,HillplanetR,color=255,psym=3
;; write_png,strcompress('out'+string(k)+'/imagenes/out'+string(k)+$
;; 'output'+string(m)+'Vrad.png',/remove_all),tvrd(0,0,1000,600,0,true=1)

;acret=radii*dens*vrad*2.*!PI/364
;MinAcret=min(acret)
;MaxAcret=max(acret)
;print,vrad
;print,acret,min(acret),max(acret)
;; for i=0,nrad-1 do begin
;; print,i,total(acret(*,i))
;; endfor
;; loadct,0
;; window,4,xsize=1000,ysize=600
;; contour,acret,xx,radii,nlevels=300,/fill,/xstyle,/ystyle,$
;; xrange=[0,2*!PI],/irregular,/ylog,ytitle='UA',xtitle='!N!4u!3!N',$
;; charsize=1.5,min_value=MinAcret,max_value=MaxAcret
;; loadct,3
;; contour,acret,xx,radii,nlevels=300,/cell_fill,/xstyle,/ystyle,$
;; xrange=[0,2*!PI],/irregular,/ylog,$
;; charsize=1.5,min_value=MinAcret,max_value=0.,/overplot
;; loadct,39
;; contour,acret,xx,radii,levels=[-1e-8,-1e-9,-1e-10,-1e-11,-1e-12],c_color=[50,100,150,200,250],/overplot
;; xyouts,100,230,'dM/dt=-1e-8',charsize=1.5,color=50,/device
;; xyouts,100,210,'dM/dt=-1e-9',charsize=1.5,color=100,/device
;; xyouts,100,190,'dM/dt=-1e-10',charsize=1.5,color=150,/device
;; xyouts,100,170,'dM/dt=-1e-11',charsize=1.5,color=200,/device
;; xyouts,100,150,'dM/dt=-1e-12',charsize=1.5,color=250,/device
;; loadct,0
;; xyouts,100,100,strcompress('t = '+string(abs(m*1000.))+' a!Z(00f1)os'),charsize=1.5,color=0,/device
;; oplot,HillplanetTheta,HillplanetR,color=255,psym=3
;; write_png,strcompress('out'+string(k)+'/imagenes/out'+string(k)+$
;; 'output'+string(m)+'Acret.png',/remove_all),tvrd(0,0,1000,600,0,true=1)

;  endfor

  RETURN
end

