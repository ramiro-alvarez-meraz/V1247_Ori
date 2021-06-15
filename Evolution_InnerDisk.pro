Pro Evolution_InnerDisk

  nPlanets=7                     ;Number of planets
  outputs=30                     ;The outputs

  if nPlanets eq 3 then npStr='three'
  if nPlanets eq 5 then npStr='five'
  if nPlanets eq 6 then npStr='six'
  if nPlanets eq 7 then npStr='seven'
  diskcase='alpha_visc-'+npStr+'_planets-no_initial_gap'		;The case
  diskcase='alpha_visc-'+npStr+'_planets'		;The case
  dataFolder='~/fargo3d_outputs/V1247_Ori-'+diskcase
  dischargeImagesFolder='~/Dropbox/visualization_tools/images/V1247_Ori-'+diskcase
  file_mkdir,strcompress(dischargeImagesFolder,/remove_all)

;Celdas del disco a considerar
  dims=dblarr(8)
  dims_=strcompress(dataFolder+'/dims.dat',/remove_all)
  openr,1,dims_
  readf,1,dims
  close,1
  nsec = uint(dims(7))
  nrad = uint(dims(6))
  Ninterm=20000
  Ninterm=long(Ninterm)

;Arreglos de cantidades físicas
  dens=dblarr(nsec,nrad)
  vrad=dblarr(nsec,nrad)
  vtheta=dblarr(nsec,nrad)
  densMed=dblarr(nrad,outputs+1)
  DensMaxima=dblarr(nrad,outputs+1)
  vradMed=dblarr(nrad,outputs+1)
  vradInner=dblarr(nsec,outputs+1)
  vthetaInner=dblarr(nsec,outputs+1)
  densInner=dblarr(nsec,outputs+1)
  CellMass=dblarr(nsec,outputs+1)
  InnerMass=dblarr(outputs+1)
  MassCondition1=dblarr(outputs+1)
  MassCondition2=dblarr(outputs+1)
  MassCondition3=dblarr(outputs+1)
  MassCondition4=dblarr(outputs+1)
  MassCondition5=dblarr(outputs+1)
  MassCondition6=dblarr(outputs+1)
  MassDensMayor=dblarr(outputs+1)
  MassDensMenor=dblarr(outputs+1)
  AcretMass=dblarr(nsec,nrad,outputs+1)
  TotalAcretMass=dblarr(nrad,outputs+1)
  StellarAcretMass=dblarr(outputs+1)
  tot=dblarr(nrad)
  rad=dblarr(nrad+1)
  areacell=dblarr(nrad)
  planet=dblarr(9,outputs+1)
  MdotStellar=dblarr(outputs+1)
  acretion=dblarr(nrad,outputs+1)
  acret=dblarr(nsec,nrad)
  vrmax=1.e-5
  px=dblarr(nPlanets,2)
  py=dblarr(nPlanets,2)
;  PlanetSemieje=dblarr(planets,outputs+1)
;  RplanetsMed=dblarr(outputs+1)
;  orbit=dblarr(6,Ninterm*outputs+1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Coloca los espacios radiales en forma de arreglos
  rad_=strcompress(dataFolder+'/used_rad.dat',/remove_all)
  openr,1,rad_
  readf,1,rad
  close,1
  rad=rad;/1.49597871e13
  rmin=rad(0)
  rmax=rad(nrad)
print,rad
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;comienza loop para barrer datos de outputs de cada corrida
  for m=0,outputs,5 do begin
    if (m eq 0) then begin
      for nb=0,nPlanets-1 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Para leer datos de masa acretada por el borde interno, del planeta 0
        planet_=strcompress(dataFolder+'/planet'+string(nb)+$
        '.dat',/remove_all)
        openr,1,planet_
        readf,1,planet 
        close,1 
        if (nb eq 0) then begin
          MdotStellar[0]=0.0
          for o=1,outputs do begin
            ;Evolucion de acreción 
           MdotStellar[o]=(planet[6,o]-planet[6,o-1])/(Ninterm/20)
          endfor
        endif
        px[nb,0]=planet[1,0]
        px[nb,1]=planet[1,0]
      endfor
    endif 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Coloca la densidad en los arreglos
    dens_=strcompress(dataFolder+'/gasdens'+string(m)+$
    '.dat'$
    ,/remove_all)
    openr,1,dens_
    readu,1,dens
    close,1
;print,dens
    dens=dens*2./2.25*1.e7*1.86
;print,dens

    for l=0,nsec-1 do begin
      densInner[l,m]=dens(l,2)
    endfor

;Ciclos para calculo de densidad media
    for r=0,nrad-1 do begin
      tot[r]=0.0
      for l=0,nsec-1 do begin
        tot[r]=dens[l,r]+tot[r]
      endfor
    endfor
    for r=0,nrad-1 do begin
      densMed[r,m]=tot[r]/nsec
    endfor
;print,densMed(*,m)
;Ciclos para calculo de densidad maxima
    Maxima=0.
    for r=0,nrad-1 do begin
      Maxima = max(dens(*,r))
      for l=0,nsec-1 do begin
        if (dens[l,r] eq Maxima) and (m gt 0) then begin;and (r ne 63) then begin
          densMaxima[r,m]=dens[l,r]
        endif
      endfor
    endfor

;Coloca la velocidad radial en los arreglos
    vrad_=strcompress(dataFolder+'/gasvy'+string(m)+$
    '.dat',/remove_all)
    openr,1,vrad_
    readu,1,vrad
    close,1

;Calcula tasas de acrecion radiales y azimutales
    for i=0,nrad-1 do begin
      tot[i]=0.
      for j=0,nsec-1 do begin
        acret[j,i]=(2.*!PI/nsec)*rad[i]*dens[j,i]*vrad[j,i]
        if (acret[j,i] < 0.) then begin
          tot[i]=acret[j,i]+tot[i]
        endif
      endfor
      acretion[i,m]=tot[i]
    endfor

;Cálculo de velocidades radiales medias
    for l=0,nsec-1 do begin
      vradInner[l,m]=vrad(l,2);Radial velocity at radial cell number 2
    endfor
    for r=0,nrad-1 do begin
      VradMed[r,m]=0.0
      for l=0,nsec-1 do begin
        VradMed[r,m]=vrad[l,r]+VradMed[r,m]
      endfor
      VradMed[r,m]=VradMed[r,m]/nsec
    endfor

;Coloca la velocidad azimutal en los arreglos
    vtheta_=strcompress(dataFolder+'/gasvx'+string(m)+$
    '.dat',/remove_all)
    openr,1,vtheta_
    readu,1,vtheta
    close,1
    for l=0,nsec-1 do begin
      vthetaInner[l,m]=vtheta(l,2)
    endfor
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Define un arreglo del número de celdas azimutales
    Ns=findgen(nsec)
    Ns=Ns+1.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  for l=0,nsec-1 do begin
;    vrvis[l] = vrmax
;  endfor
    densInnerMed=total(densInner(*,m))/nsec
;determina la densidad media de las celdas internas

    cellmass1=0.0
    cellmass2=0.0
    cellmass3=0.0
    cellmass4=0.0
    cellmass5=0.0
    cellmass6=0.0
    cellmass7=0.0
    cellmass8=0.0


    for l=0,nsec-1 do begin
      CellMass[l,m]=densInner(l,m)*AreaCell[0]
;determina la masa de la celda interna
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;valores acotados por densidad media
; y velocidad radial máxima

      if (vrad(l,2) lt vrmax) and (densInner[l,m] lt densMed[2,m]) $
      then begin
       cellmass1=CellMass[l,m]+cellmass1
      endif

      if (vrad(l,2) gt vrmax) and (vrad(l,2) lt 0.0) and $
        (densInner[l,m] lt densMed[2,m]) then begin
        cellmass2=CellMass[l,m]+cellmass2
      endif

      if (vrad(l,2) gt 0.0) and (densInner[l,m] lt densMed[2,m]) then begin
        cellmass3=CellMass[l,m]+cellmass3
      endif
 
      if (vrad(l,2) lt vrmax) and (densInner[l,m] gt densMed[2,m]) $
      then begin
        cellmass4=CellMass[l,m]+cellmass4
      endif

      if (vrad(l,2) gt vrmax) and (vrad(l,2) lt 0.0) and $
        (densInner[l,m] gt densMed[2,m]) then begin
        cellmass5=CellMass[l,m]+cellmass5
      endif

      if (vrad(l,2) gt 0.0) and (densInner[l,m] gt densMed[2,m]) then begin
        cellmass6=CellMass[l,m]+cellmass6
      endif

      if (densInner[l,m] gt densMed[2,m]) then begin
        cellmass7=CellMass[l,m]+cellmass7
      endif

      if (densInner[l,m] lt densMed[2,m]) then begin
        cellmass8=CellMass[l,m]+cellmass8
      endif
    endfor
    InnerMass[m]=total(CellMass(*,m))
    MassCondition1[m]=cellmass1
    MassCondition2[m]=cellmass2  
    MassCondition3[m]=cellmass3  
    MassCondition4[m]=cellmass4  
    MassCondition5[m]=cellmass5
    MassCondition6[m]=cellmass6
    MassDensMayor[m]=cellmass7
    MassDensMenor[m]=cellmass8

    for r=0,nrad-1 do begin
      TotalAcretMass[r,m]=0.0
      for l=0,nsec-1 do begin
        if (vrad(l,r) lt 0.0) then begin
          AcretMass[l,r,m]=(2.*!PI/nsec)*rad[r]*dens[l,r]*vrad[l,r]
          TotalAcretMass[r,m]=AcretMass[l,r,m]+TotalAcretMass[r,m]
        endif
      endfor
    endfor
    StellarAcretMass[m]=TotalAcretMass[2,m]
print,StellarAcretMass[m]
;    RplanetsMed[m]=(PlanetSemieje[0,m]+PlanetSemieje[1,m])/planets

    for i=1,nrad-1 do begin
;      if RplanetsMed[m] gt rad[i-1] and RplanetsMed[m] lt rad[i] then begin
if i eq 20 then begin
          totalneg=0.
          totalpos=0.
;        for i=i-5,i+5 do begin
          maxdensi=max(dens(*,i))
          totaldensi=total(dens(*,i))
          ;MeddensMedMax=(densMed[i,m]+maxdensi)/2
    ;; for r=0,nrad-1 do begin
    ;;   tot[r]=0.0
    ;;   for l=0,nsec-1 do begin
    ;;     tot[r]=dens[l,r]+tot[r]
    ;;   endfor
    ;; endfor
    ;; for r=0,nrad-1 do begin
    ;;   densMed[r,m]=tot[r]/nsec
    ;; endfor



          for j=0,nsec-1 do begin
            if (vrad[j,i] lt 0.0) then begin
;              totalneg=2*!PI/nsec*vrad[j,i]*rad[i]*dens[j,i]+totalneg
              totalneg=dens[j,i]+totalneg
            endif else begin
;              totalpos=2*!PI/nsec*vrad[j,i]*rad[i]*dens[j,i]+totalpos  
              totalpos=dens[j,i]+totalpos
            endelse

                ;else begin
            ;endelse

            ;; if (dens[j,i] ge densMed[i,m]) then begin
            ;;   print,m,j,i,dens[j,i] 
            ;; endif else begin
            ;; endelse
          endfor
;          print,totalneg,totalpos,totalneg/totalpos
;print,m,i,rad[i],densMed[i]
;        endfor
;          print,totalneg,totalpos,totalneg/totalpos
      endif
    endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  endfor         
;termina loop para barrer los distintos datos (densidad,
;velocidad, etc) de cada carpeta

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Detecta datos máximos y mínimos
vradInnerMax=max(vradInner)
vradInnerMin=min(vradInner)
vthetaInnerMax=max(vthetaInner)
vthetaInnerMin=min(vthetaInner) 
densMedMax=max(densMed)
densMedMin=min(densMed)
densInnerMax=max(densInner)
densInnerMin=min(densInner)
CellMassMax=max(CellMass)
CellMassMin=min(CellMass)
cellmass1Max=max(cellmass1)
cellmass1Min=min(cellmass1)
cellmass2Max=max(cellmass2)
cellmass2Min=min(cellmass2)
cellmass3Max=max(cellmass3)
cellmass3Min=min(cellmass3)
cellmass4Max=max(cellmass4)
cellmass4Min=min(cellmass4)
cellmass5Max=max(cellmass5)
cellmass5Min=min(cellmass5)
StellarAcretMassMax=max(abs(StellarAcretMass))
StellarAcretMassMin=min(abs(StellarAcretMass))
TotalAcretMassMax=max(abs(TotalAcretMass))
TotalAcretMassMin=min(abs(TotalAcretMass))

MassMax=max(InnerMass)
time=findgen(m)*Ninterm/20

;Define tabla de colores
device,get_decompose=old_decomposed,decomposed=0
loadct,39
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;produce una gráfica (Ns,vrad)
;; window,0,xsize=600,ysize=600
;; plot,ns,vradInner(*,0),background=255,xrange=[1,nsec],$
;; yrange=[vradInnerMin,vradInnerMax],xstyle=1,ystyle=1,$
;; color=0,xtitle='Celda Azimutal',ytitle='Velocidad radial (ua/a!Z(00f1)o)',$
;; psym=0,charsize=2.0
;; for m=1,outputs do begin
;; oplot,ns,vradInner(*,m),color=m;255*m/outputs
;; endfor
;; write_png,strcompress('out'+string(k)+'/imagenes/out'+string(k)+$
;; 'Ns-vrad.png',/remove_all),tvrd(0,0,600,600,1,true=1) 

;produce una gráfica (Rmed,densMed)
;window,1,xsize=600,ysize=600
;dimensionless to cgs 2./2.25*1.e7
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
plot,rad,densMed(*,0),background=255,$
xrange=[rmin,rmax],xstyle=1,ystyle=1,$
yrange=[0.0002,10.],/xlog,/ylog,$
color=0,xtitle='Radius ( AU )',ytitle='<  !N!4R!3!N (g cm!E-2!N)  >',$
psym=0,charsize=1.5
;        f=STRTRIM(f, 2)
  nPlanetsStr=strtrim(nPlanets,2)
;xyouts,260,540,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,520,nPlanetsStr+' planets with 2:1 mutual resonance',charsize=1.5,color=0,/device
;  xyouts,150,490,'(2:1 mutual resonance)',charsize=1.5,color=0,/device
;  xyouts,150,465,'alpha visc disk case',charsize=1.5,color=0,/device
;for m=1,outputs do begin
;oplot,rad,2./2.25*1.e7*densMed(*,m),color=m*0.84;255*m/outputs
;endfor
;-----------------------------------------------------------------
;-----------------------------------------------------------------
; oplot,rad,densMed(*,2),color=40
; xyouts,430,540,'2.e0 years',charsize=1.5,color=40,/device
; oplot,rad,densMed(*,4),color=80
; xyouts,430,515,'4.e0 years',charsize=1.5,color=80,/device
; oplot,rad,densMed(*,6),color=120
; xyouts,430,490,'6.e0 years',charsize=1.5,color=120,/device
;  oplot,rad,densMed(*,8),color=160
;  xyouts,430,465,'8.e0 years',charsize=1.5,color=160,/device
;  oplot,rad,densMed(*,1),color=200
;  xyouts,430,440,'1.e1 years',charsize=1.5,color=200,/device
;;  oplot,rad,densMed(*,5),color=240
;;  xyouts,430,415,'6.e1 years',charsize=1.5,color=240,/device
;----------------------------------------------------------------
;----------------------------------------------------------------
;oplot,rad,densMed(*,6),color=20
;xyouts,480,540,'1.0e5 years',charsize=1.0,color=20,/device
;oplot,rad,densMed(*,7),color=40
;xyouts,480,525,'1.2e5 years',charsize=1.0,color=40,/device
;oplot,rad,2./2.25*1.e7*densMed(*,140),color=60
;xyouts,480,510,'1.4e5 years',charsize=1.0,color=60,/device
;oplot,rad,densMed(*,5),color=80
;xyouts,480,495,'1.6e5 years',charsize=1.0,color=80,/device
;oplot,rad,2./2.25*1.e7*densMed(*,180),color=100
;xyouts,480,480,'1.8e5 years',charsize=1.0,color=100,/device
;oplot,rad,densMed(*,10),color=120
;xyouts,480,465,'1.0e4 years',charsize=1.0,color=120,/device
;oplot,rad,2./2.25*1.e7*densMed(*,220),color=140
;xyouts,480,450,'2.2e5 years',charsize=1.0,color=140,/device
;oplot,rad,densMed(*,20),color=160
;xyouts,480,435,'2.e4 years',charsize=1.0,color=160,/device
;oplot,rad,2./2.25*1.e7*densMed(*,260),color=180
;xyouts,480,420,'2.6e5 years',charsize=1.0,color=180,/device
;oplot,rad,densMed(*,30),color=180
;xyouts,480,405,'3.e4 years',charsize=1.0,color=180,/device
;oplot,rad,densMed(*,8),color=220
;xyouts,480,390,'8.0e4 years',charsize=1.0,color=220,/device
;for i=10,outputs,10 do begin
;  oplot,rad,densMed(*,i),color=i/4.
;endfor
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
  oplot,rad,densMed(*,10),color=80,thick=3,linestyle=2
  device,get_decompose=old_decomposed,decomposed=0
  loadct,9
  oplot,rad,densMed(*,20),color=160,thick=3,linestyle=3
  device,get_decompose=old_decomposed,decomposed=0
  loadct,39
  oplot,rad,densMed(*,30),color=240,thick=3,linestyle=4

;----------------------------------------------------------------
;oplot,rad,densMed(*,),color=240
;xyouts,480,375,'1.e5 years',charsize=1.0,color=240,/device
;oplot,rad,2./2.25*1.e7*densMed(*,220),color=250
;xyouts,480,375,'2.2e5 years',charsize=1.0,color=250,/device
  nsecStr=strtrim(nsec,2)
  nradStr=strtrim(nrad,2)
  rmaxStr=strtrim(round(rmax),2)
  rminStr=strtrim(round(rmin),2)
  py(*,0)=0.0001
  py(*,1)=10.
print,px,py
  for i=0,nplanets-1 do begin
    oplot,px(i,*),py(i,*),color=0,linestyle=1
  endfor
;  oplot,px(),py(),psym=1,color=0
write_png,strcompress(dischargeImagesFolder+$
'/Rmed-densMed'+rminStr+'<r<'+rmaxStr+'AU_'+nsecStr+'-'+nradStr+'.png',/remove_all),tvrd(0,0,600,600,1,true=1)


;Define tabla de colores
device,get_decompose=old_decomposed,decomposed=0
loadct,39
     Set_plot,'z'
     Device,Set_resolution=[600,600],Set_Pixel_Depth=24
plot,rad,densMed(*,5)*9.6,background=255,$
xrange=[rmin,rmax],xstyle=1,ystyle=1,$
yrange=[0.02,900.],/xlog,/ylog,$
color=0,xtitle='Radius ( AU )',ytitle='< tau >',$
psym=0,charsize=1.5
;        f=STRTRIM(f, 2)
  nPlanetsStr=strtrim(nPlanets,2)
;xyouts,260,540,'f='+f,charsize=1.5,color=0,/device
  xyouts,150,515,nPlanetsStr+' planets',charsize=1.5,color=0,/device
  xyouts,150,490,'(2:1 mutual resonance)',charsize=1.5,color=0,/device
  xyouts,150,465,'invisc disk case',charsize=1.5,color=0,/device
  oplot,rad,densMed(*,10)*9.6,color=200
  xyouts,430,440,'1.e1 years',charsize=1.5,color=200,/device
;write_png,strcompress(dischargeImagesFolder+$
;'/Rmed-tauMed.png',/remove_all),tvrd(0,0,600,600,1,true=1)
;stop

;produce una gráfica (ns,densInner)
;loadct,0
;window,2,xsize=600,ysize=600
;plot,ns,2./2.25*1.e7*densInner(*,0),background=255,$
;xrange=[1,nsec],$
;yrange=[2./2.25*1.e7*densInnerMin,10.*2./2.25*1.e7*densInnerMax],/ylog,$
;color=0,xtitle='ns',ytitle='densidad de celda interna',$
;psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,ns,2./2.25*1.e7*densInner(*,m),color=round(255*m/outputs)
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'ns-densInner.gif',/remove_all),tvrd(0,0,600,600,1)

;produce una gráfica (Rmed,Mdot/(nu*3PI))
;loadct,0
;window,3,xsize=600,ysize=600
;plot,Rmed,-2./2.25*1.e7*(1.e-8)/(3.*!PI*(2./3.)*Rmed*vrad(0,*)),background=255,$
;xrange=[rmin,RmedMax],$
;yrange=[2./2.25*1.e7*densMedMin,10.*2./2.25*1.e7*densMedMax],/xlog,/ylog,$
;color=0,xtitle='Rmed',ytitle='Mdot/(nu*3PI)',$
;psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,Rmed,2./2.25*1.e7*(2./3.*1.e-8*Rmed*vrad(0,*)),color=round(255*m/outputs)
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'Rmed-Mdotbynu3PI.gif',/remove_all),tvrd(0,0,600,600,1)

;produce una gráfica (Ns,vtheta)
;window,4,xsize=600,ysize=600
;plot,ns,vthetaInner(*,0),background=255,xrange=[1,nsec],$
;yrange=[vthetaInnerMin,vthetaInnerMax],$
;color=0,xtitle='Celda Azimutal',ytitle='Velocidad Azimutal',$
;psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,ns,vthetaInner(*,m),color=255*m/outputs
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'Ns-vtheta.gif',/remove_all),tvrd(0,0,600,600,1) 

;produce una gráfica (Ns,Mass)
;window,5,xsize=600,ysize=600
;plot,ns,CellMass(*,0),background=255,xrange=[1,nsec],$
;yrange=[CellMassMin,CellMassMax],$
;color=0,xtitle='Celda Azimutal',ytitle='Masa de celda interior',$
;psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,ns,CellMass(*,m),color=255*m/outputs
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'Ns-Mass.gif',/remove_all),tvrd(0,0,600,600,1)

;produce una gráfica (time,InnerMass)
;window,6,xsize=600,ysize=600
;plot,time,InnerMass,background=255,xrange=[0,outputs],/ylog,$
;yrange=[1.e-7,1.e-5],$
;color=0,xtitle='Tiempo',ytitle='Masas de anillo interior',$
;psym=0,charsize=1.0
;xyouts,100,170,'Masa interior total',charsize=1.0,$
;color=0,/device
;loadct,10
;oplot,time,MassCondition1,color=40
;xyouts,100,150,'Masa interior vr<vrmax y rho<rhomed',$
;charsize=1.0,color=40,/device
;oplot,time,MassCondition2,color=80
;xyouts,100,130,'Masa interior vr>vrmax, vr<0.0 y rho<rhomed',$
;charsize=1.0,color=80,/device
;oplot,time,MassCondition3,color=120
;xyouts,100,110,'Masa interior vr>0.0 y rho<rhomed',$
;charsize=1.0,color=120,/device
;oplot,time,MassCondition4,color=160
;xyouts,100,90,'Masa interior vr<vrmax y rho>rhomed',$
;charsize=1.0,color=160,/device
;oplot,time,MassCondition5,color=200
;xyouts,100,70,'Masa interior vr>vrmax, vr<0.0 y rho>rhomed'$
;,charsize=1.0,color=200,/device
;oplot,time,MassCondition6,color=240
;xyouts,100,50,'Masa interior vr>0.0 y rho>rhomed',$
;charsize=1.0,color=240,/device
;oplot,time,MassDensMayor,color=200
;xyouts,100,70,'Masa interior con rho>rhomed',$
;charsize=1.0,color=200,/device
;oplot,time,MassDensMenor,color=240
;xyouts,100,50,'Masa interior con rho<rhomed',$
;charsize=1.0,color=240,/device
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'time-InnerMass.gif',/remove_all),tvrd(0,0,600,600,1)

;produce gráfica (Time,AcretionMass)
;window,7,xsize=600,ysize=600
;plot,time,abs(StellarAcretMass),background=255,$
;xrange=[0,outputs],yrange=[1.e-10,1.e-5],xstyle=1,ystyle=1,$
;color=0,xtitle='tiempo (a!Z(00f1)os)',/ylog,$
;ytitle='tasa de acrecion ()',psym=0,charsize=1.0
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'time-StellarAcretMass.gif',/remove_all),tvrd(0,0,600,600,1)

;produce gráfica (Time,AcretionMass)
;window,8,xsize=600,ysize=600
;plot,ns,abs(AcretMass(*,2,0)),background=255,$
;xrange=[1,nsec],yrange=[1.e-12,1.e-7],$
;color=0,xtitle='ns',/ylog,$
;ytitle='tasa de acrecion',psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,ns,abs(AcretMass(*,2,m)),color=255*m/outputs
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'ns-AzimutalStellarAcretMass.gif',/remove_all),tvrd(0,0,600,600,1)

;produce gráfica (Time,AcretionMass)
;window,9,xsize=600,ysize=600
;plot,Rmed,abs(TotalAcretMass(*,0)),background=255,$
;xrange=[0,RmedMax],yrange=[1.e-11,1.e-5],$;*!PI*TotalAcretMassMax],$
;color=0,xtitle='Radio (UA)',/ylog,$
;ytitle='tasa de acrecion',psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,Rmed,abs(TotalAcretmass(*,m)),color=255*m/outputs
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'Rmed-TotalAcretMass.gif',/remove_all),tvrd(0,0,600,600,1)

;produce gráfica (Time,AcretionMass)
;window,10,xsize=600,ysize=600
;plot,rad,vradMed(*,0),background=255,$
;xrange=[rmin,rmax],yrange=[-.1,.1],/xlog,$
;yrange=[vradInnerMin,vradInnerMax],/xlog,$
;color=0,xtitle='Radio (UA)',$
;ytitle='vradMed',psym=0,charsize=1.0
;for m=1,outputs do begin
;oplot,rad,vradMed(*,m),color=255*m/outputs
;endfor
;write_gif,strcompress('out'+string(k)+'/out'+string(k)+$
;'Rmed-vradMed.gif',/remove_all),tvrd(0,0,600,600,1)

;produce gráfica (Time,MdotStellar)
;window,11,xsize=600,ysize=600
;plot,time,MdotStellar,background=255,xstyle=1,$
;xrange=[1,outputs*Ninterm/20],yrange=[MdotStellar[outputs],1.e-7],$
;color=0,xtitle='tiempo (a!Z(00f1)os)',/ylog,$
;ytitle='dM/dt  (M!L!9n!3!N/a!Z(00f1)o)',psym=0,charsize=1.0
;oplot,time,MdotStellar,color=50
;write_png,strcompress('out'+string(k)+'/imagenes/out'+string(k)+$
;'time-MdotStellar.png',/remove_all),tvrd(0,0,600,600,1,true=1)
;endfor
    RETURN
end

