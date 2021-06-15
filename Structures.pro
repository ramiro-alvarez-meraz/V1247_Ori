Pro Structures

;  f=1                           ;Accretion parameter
;  nPlanets=2                   ;Number of planets
;  output_=300                    ;The output

;for output=300,output_,100 do begin
;  dischargeImagesFolder='~/Dropbox/visualization_tools/images/StructuresQuantification/f'+string(f)+'_p'+string(nPlanets)+'/output'+string(output)
;  dataFolder='~/ram_vogto/fargo/out'+string(out)
  nPlanets=7                     ;Number of planets
  output=30                     ;The output

  if nPlanets eq 5 then npStr='five'
  if nPlanets eq 6 then npStr='six'
  if nPlanets eq 7 then npStr='seven'
  diskcase='alpha_visc-'+npStr+'_planets'		;The case
  dataFolder='~/fargo3d_outputs/V1247_Ori-'+diskcase
  file_mkdir,strcompress('~/Dropbox/visualization_tools/images/V1247_Ori-'+diskcase,/remove_all)
  dischargeImagesFolder='~/Dropbox/visualization_tools/images/V1247_Ori-'+diskcase+'/StructuresQuantification'
  file_mkdir,strcompress(dischargeImagesFolder,/remove_all)
  dims=dblarr(8)
  dims_=strcompress(dataFolder+'/dims.dat',/remove_all)
  openr,1,dims_
  readf,1,dims
  close,1
  nsec = uint(dims(7))
  nrad = uint(dims(6))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; planet_=strcompress(dataFolder+'/out'+string(out)+$
  ;;                     'PlanetSemieje.dat',/remove_all)
  ;; openr,1,planet_

  ;; Semieje=dblarr(nPlanets)
  ;; for l=0,1100 do begin
  ;;    readf,1,Semieje
  ;;    if (l eq output) then break
  ;; endfor
  ;; close,1
  ;; SemiejeMax=max(Semieje)
  ;; SemiejeMin=min(Semieje)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  xPlanet=dblarr(nPlanets)
  yPlanet=dblarr(nPlanets)
  vxPlanet=dblarr(nPlanets)
  vyPlanet=dblarr(nPlanets)
  mPlanet=dblarr(nPlanets)
  ThetaPlanet=dblarr(nPlanets)
  Rplanet=dblarr(nPlanets)
  HillPlanet=dblarr(nPlanets)
  Semieje=dblarr(nPlanets)
  for i=0,nPlanets-1 do begin
     planetAsString = STRTRIM(i, 2)
;     outAsString = STRTRIM(out, 2)
     planet_=strcompress(dataFolder+'/planet'+planetAsString+'.dat',/remove_all)
     openr,1,planet_
     s=dblarr(9)
     planet=dblarr(9,output+1)
     oldS=0
     WHILE ~ EOF(1) DO BEGIN
        readf,1,s 
        if ((s[0] eq oldS) and (s[0] ne 0)) then continue 
        nn=abs(s[0])
        planet[*,nn]=s
        oldS=s[0]
        if (s[0] eq output) then goto,jump
     ENDWHILE
     jump:
     xPlanet[i]=planet[1,output]
     yPlanet[i]=planet[2,output]
     vxPlanet[i]=planet[3,output]
     vyPlanet[i]=planet[4,output] 
     mPlanet[i]=planet[5,output]
     close,1
  endfor
  h=xPlanet*vyPlanet-yPlanet*vxPlanet
  d=sqrt(xPlanet*xPlanet+yPlanet*yPlanet)
  Ax=xPlanet*vyPlanet*vyPlanet-yPlanet*vxPlanet*vyPlanet-xPlanet/d
  Ay=yPlanet*vxPlanet*vxPlanet-xPlanet*vxPlanet*vyPlanet-yPlanet/d
  e=sqrt(Ax*Ax+Ay*Ay)
  aa=h*h/(1-e*e)
  SemiejeMax=max(aa)
  SemiejeMin=min(aa)
  for i=0,nPlanets-1 do begin
     Rplanet[i]=sqrt(xPlanet[i]^2+yPlanet[i]^2)
     ThetaPlanet[i]=atan(yPlanet[i],xPlanet[i])+(1.0/nsec*2*!PI)
     if (ThetaPlanet[i] < 0.) then ThetaPlanet[i]=2.*!PI+ThetaPlanet[i]+(1.0/nsec*2*!PI)
     xplanet[i]=Rplanet[i]*cos(ThetaPlanet[i])
     yplanet[i]=Rplanet[i]*sin(ThetaPlanet[i])
  endfor
  h=xPlanet*vyPlanet-yPlanet*vxPlanet
  d=sqrt(xPlanet*xPlanet+yPlanet*yPlanet)
  Ax=xPlanet*vyPlanet*vyPlanet-yPlanet*vxPlanet*vyPlanet-xPlanet/d
  Ay=yPlanet*vxPlanet*vxPlanet-xPlanet*vxPlanet*vyPlanet-yPlanet/d
  e=sqrt(Ax*Ax+Ay*Ay)
  aa=h*h/(1-e*e)
  SemiejeMax=max(aa)
  SemiejeMin=min(aa)
  Hillplanet=aa(1-e)*((1./3)*mplanet)^(1./3.)
  distPlanetCell=dblarr(nsec,nrad)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Optical depth ranges:
  TauRanges=6
  tau_flag=dblarr(TauRanges-1)
  tau_flag[0]=0.01
  tau_flag[1]=0.1
  tau_flag[2]=0.5
  tau_flag[3]=1.0
  tau_flag[4]=2.0

;Maximun number of structures per tau range:
  Struct_max=50 

;Used arrays:
  rad=dblarr(nrad+1)                  ;Radial array
  dens=dblarr(nsec,nrad)              ;Density array 
  vrad=dblarr(nsec,nrad)              ;Radial velocity array
  vtheta=dblarr(nsec,nrad)            ;Theta velocity array
  SigmaTau=dblarr(nrad,TauRanges)     ;Mean azimutally density per tau range
  Counter=dblarr(nrad,Tauranges)      ;Count the cells radially per tau range
  DensTau=dblarr(nsec,nrad,TauRanges) ;Density array per tau range
  i_=dblarr(TauRanges)                ;i index for initial structure
  j_=dblarr(Tauranges)                ;j index for initial structure
  Block_prev=dblarr(TauRanges)        ;Array for count structures
  Block=dblarr(nsec,nrad,TauRanges,Struct_max) ;Structure array per tau range
  A=dblarr(nsec,nrad)                 ;Area Cell array
  MassTau=dblarr(TauRanges)
  ATau=dblarr(TauRanges)
  MeanFillingFactorAtGap=dblarr(TauRanges)

;Put the radial space on array form:
  rad_=strcompress(dataFolder+$
                   '/used_rad.dat',/remove_all)
  openr,1,rad_
  readf,1,rad
  close,1
  rmin=rad(0)
  rmax=rad(nrad)

;Put the density space on array form:
  dens_=strcompress(dataFolder+$
                    '/gasdens'+string(output)+'.dat'$
                    ,/remove_all)
  openr,1,dens_
  readu,1,dens
  close,1
  dens=dens*1.86 ;Because M_*=1.86M_solar

;Put the radial velocity on array form:
  vrad_=strcompress(dataFolder+$
                    '/gasvy'+string(output)+'.dat'$
                    ,/remove_all)
  openr,1,vrad_
  readu,1,vrad
  close,1

;Put the azimutal velocity on array form:
  vtheta_=strcompress(dataFolder+$
                      '/gasvx'+string(output)+$
                      '.dat',/remove_all)
  openr,1,vtheta_
  readu,1,vtheta
  close,1

;Put the area space on arrray form
  for i=0,nrad do begin
     if (i gt 0) then begin
        for j=0,nsec-1 do begin
           A[j,i-1]=2*!PI*(rad[i]^2-rad[i-1]^2)/nsec
        endfor
     endif
  endfor

  

  for i=0,nrad-1 do begin
     if ((rad[i] ge (Rplanet[0]-(2*HillPlanet[0]))) and (rad[i] le (Rplanet[nPlanets-1]-(2*HillPlanet[nPlanets-1])))) then begin
        for j=0,nsec-1 do begin
              tau=9.6*dens[j,i]*2./2.25*1.e7 ;Optical depth Troilita,etc.              
;Optically thin:
;              if (tau lt 0.5) then begin
;                 l=0
              if (tau lt tau_flag[0]) then begin
                 l=0
                 MassTau[l] += dens[j,i]*A[j,i]
                 ATau[l] += A[j,i]
              endif
;Optically intermediate:
;              if (tau ge 0.5) and (tau le 2.0) then begin
;                 l=1
              for ii=1,TauRanges-2 do begin
                if (tau ge tau_flag[ii-1]) and (tau le tau_flag[ii]) then begin
	           l=ii
                   MassTau[l] += dens[j,i]*A[j,i]
                   ATau[l] += A[j,i]
                endif
              endfor
;Optically thick:
;              if (tau gt 2.0) then begin
;                 l=2
              if (tau gt tau_flag[TauRanges-2]) then begin
                 l=TauRanges-1
                 MassTau[l] += dens[j,i]*A[j,i]
                 ATau[l] += A[j,i]
              endif
           endfor
     endif
  endfor
  for i=0,TauRanges-1 do begin
  MeanFillingFactorAtGap[i]=Atau[i]/total(ATau(*))
  endfor
  
;;;;;;;;;;;;;;;;;;
;;;
;;;     Analisis of Optical depth: 
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  rmed=dblarr(nrad)
  rmedt=dblarr(nsec,nrad)
  AccTauBlock=dblarr(nsec,nrad,TauRanges,Struct_max)
  AccRadPerBlock=dblarr(nrad,TauRanges,Struct_max)
  AccTotalPerBlock=dblarr(TauRanges,Struct_max)
  rmed = (rad(1:nrad)-rad(0:nrad-1))/2+rad(0:nrad-1)
  rmedt=rmed##replicate(1,nsec)
  vradbig = [vrad[*,*],vrad[0,*]]
  vrCentered=interpolate(vradbig,dindgen(nsec),$
                         dindgen(nrad)+0.5,/grid)
@StructuresDefineOpticalDepth.pro
  DensTau_real=DensTau
@StructuresBlocks.pro 

;We need to recover the original dens array
  total_dens=dblarr(nsec,nrad)
  total_dens_=strcompress(dataFolder+$
                          '/gasdens'+string(output)+'.dat'$
                          ,/remove_all)
  openr,1,total_dens_
  readu,1,total_dens
  close,1
  total_dens=1.86*total_dens
  for l=0,TauRanges-1 do begin
     for m=1,Block_prev[l] do begin
        openw,1,strcompress(dischargeImagesFolder+'/l='+string(l)+$
                            '_m='+string(m),/remove_all)
        printf,1,Block(*,*,l,m)
        close,1
        new_discharge=dblarr(nsec,nrad) 
        for i=0,nrad-1 do begin
           for j=0,nsec-1 do begin
              if (Block[j,i,l,m] gt 0) then begin
                 new_discharge[j,i]=total_dens[j,i]
              endif
           endfor
        endfor
        openw,2,strcompress(dischargeImagesFolder+'/gasdens_l='+string(l)+$
                            '_m='+string(m),/remove_all)
        printf,2,new_discharge(*,*)
	close,2
     endfor
  endfor

  MassRadPerBlock=dblarr(nrad,TauRanges,Struct_max) ;Structure Mass array per tau range
  AreaRadPerBlock=dblarr(nrad,TauRanges,Struct_max) ;Structure Area array per tau range
  BlockRad_ini=dblarr(TauRanges,Struct_max)
  BlockRad_fin=dblarr(TauRanges,Struct_max)
  AreaTotalPerBlock=dblarr(TauRanges,Struct_max)
  Structures=dblarr(TauRanges,Struct_max)
  AccTauRadPerBlock=dblarr(TauRanges,Struct_max,nPlanets+1)
  MassTauRadPerBlock=dblarr(TauRanges,Struct_max,nPlanets+1)
  r=dblarr(nPlanets+1)
  r[0]=10
  r[1]=15
  r[2]=30

;Escribe las estructuras en archivos y manipula datos en estructuras:
  for l=0,TauRanges-1 do begin
     for m=1,Block_prev[l] do begin
        for i=0,nrad-1 do begin
           for j=0,nsec-1 do begin
              if (Block[j,i,l,m] eq m) then begin
                 BlockRad_ini[l,m]=rad[i]
                 goto,jump7          
              endif
           endfor
        endfor
        jump7:
        for i=nrad-1,0,-1 do begin
           for j=nsec-1,0,-1 do begin
              if (Block[j,i,l,m] eq m) then begin
                 BlockRad_fin[l,m]=rad[i]
                 goto,jump8          
              endif 
           endfor
        endfor
        ;; for i=0,nrad-1 do begin
        ;;    if (rad[i] gt BlockRad_ini[l,m]) then begin
        ;;       contador=0
        ;;       for j=0,nsec-1 do begin
        ;;          if (Block[j,i,l,m] eq m) then begin
        ;;             contador += 1
        ;;          endif
        ;;       endfor
        ;;       if (contador eq 0) then begin
        ;;          BlockRad_fin[l,m]=rad[i]
        ;;          goto,jump8
        ;;       endif
        ;;    endif
        ;; endfor
        jump8:
        if (BlockRad_ini[l,m] eq rad[nrad-1]) $
        then BlockRad_fin[l,m]=BlockRad_ini[l,m]+1
        for i=0,nrad-1 do begin
           MassRadPerBlock_prov=0
           AreaRadPerBlock_prov=0
           AccRadPerBlock_prov=0
           for j=0,nsec-1 do begin
              if (Block[j,i,l,m] eq m) then begin
                 MassRadPerBlock_prov += dens[j,i]*A[j,i]
                 AreaRadPerBlock_prov += A[j,i]
                 AccRadPerBlock_prov += (2.*!PI/nsec)*dens[j,i]*$
                                        rmedt[j,i]*VrCentered[j,i]
;      AccTau[j,i,l,m] += (2.*!PI/nsec)*dens[j,i]*rmedt[j,i]*VrCentered[j,i]
              endif
           endfor
           MassRadPerBlock[i,l,m] = MassRadPerBlock_prov
           AreaRadPerBlock[i,l,m] = AreaRadPerBlock_prov
           AccRadPerBlock[i,l,m] = AccRadPerBlock_prov
        endfor
        AreaTotalPerBlock[l,m]=total(AreaRadPerBlock(*,l,m))
        AccTotalPerBlock[l,m]=total(AccRadPerBlock(*,l,m))
        for i=0,nrad-1 do begin
           for n=0,nPlanets do begin
              if (r[n] gt rad[i]) $
                 and (r[n] lt rad[i+1]) then begin
                 for j=0,nsec-1 do begin
                    if (Block[j,i,l,m] eq m) then begin                       
                       AccTauRadPerBlock[l,m,n]+=(2.*!PI/nsec)*dens[j,i]*$
                                                 rmedt[j,i]*VrCentered[j,i]
                       MassTauRadPerBlock[l,m,n]+=dens[j,i]*A[j,i]
                    endif
                 endfor
                 if AccTauRadPerBlock[l,m,n] gt 0. then AccTauRadPerBlock[l,m,n]=0.
              endif
           endfor
        endfor
     endfor
  endfor

;Graphics:
  device,get_decompose=old_decomposed,decomposed=0
  loadct,14
;loadct,39
@StructuresMassRadPerBlock.pro
@StructuresAreaRadPerBlock.pro
@StructuresAreaTotalPerBlock.pro  
@StructuresArea.pro
@StructuresSigmaTau.pro
@StructuresAccTau.pro
;@StructuresFillingFactor.pro
@StructuresFillingFactorGap.pro
@StructuresBorders.pro
@StructuresMass.pro
@StructuresMassPerRadii.pro
@StructuresAccPerRadii.pro
;f=STRTRIM(f, 2)
;nPlanets=STRTRIM(nPlanets,2)
;output=STRTRIM(output,2)
;  spawn,'rm ~/Dropbox/visualization_tools/images/StructuresQuantification/f'+f+'_p'+nPlanets+'/output'+output+'/l=*_m=*' 

;endfor
 
 RETURN
end
