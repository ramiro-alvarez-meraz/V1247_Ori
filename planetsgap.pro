Pro planetsgap

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Constants
np=7    ;Planets number
msource=1.86   ;Source mass
tsystem=7.4e6  ;Time of system
rInDisk=0.18   ;Protoplanetary inner radius
ringap=2.0     ;Inner gap radius
routgap=60.    ;Outer gap radius 44
rOutDisk=100.  ;Outer disk simulated (we consider this source have 500 AU)

;Array quantities
rp=dblarr(np)
mp=dblarr(np)
a=dblarr(np)
alpha=dblarr(np-1)
mult_alpha_inv=dblarr(np-1)
;The folder to discharge the generated images
;  dischargeImageFolder='~/Dropbox/visualization_tools/images/V1247_Ori'
  dischargeImageFolder='/media/ramiro/MiniStation/Dropbox/visualization_tools/images/V1247_Ori'

factor=2.      ;Factor of 0.01 M_star considerated in the protoplanetary disk
Sigma0=52.*factor    ;The surface density at AU multiplied by a factor
cgsToAU_grbycm2=(1.49597871e13*1.49597871e13)/1.9891e33 ;Conversion CGS surface density to Astronomical units
MnoDepleted=2.*!pi*Sigma0*(routgap-rinDisk)*cgsToAU_grbycm2 ;Total mass considering neither gap or inner disk 
Mgap=!pi*((routgap^2)-(ringap^2))*0.1*cgsToAU_grbycm2      ;Mass contained at gap, considering constant surface density
MinDisk=2.*!pi*Sigma0*(ringap-rinDisk)*cgsToAU_grbycm2/factor*0.01   ;Mass of inner disk (depleted by factor*500 in order to have an optically thickdisk at 2 AU)
Mdepleted=MnoDepleted-Mgap-MinDisk     ;Mass depleted (we do not include)
;print,MnoDepleted,Mgap,MinDisk,Mdepleted
Mtot_planets=0.9*Mdepleted		;Total mas of planets
Macc_star=0.1*Mdepleted			;Stellar accreted mass
MaccRate=Macc_star/tsystem		;Mean stellar accretion rate 
print,MaccRate
  for i=0,np-2 do begin      ;Loop to determine the inverse multiplicative alpha factors
    alpha[i]=0.7+(i*(0.2/(np-2.)))
    if i eq 0 then mult_alpha_inv[i]=1./alpha[i]
    if i eq 1 then mult_alpha_inv[i]=1./alpha[i]/alpha[i-1]
    if i eq 2 then mult_alpha_inv[i]=1./alpha[i]/alpha[i-1]/alpha[i-2]
    if i eq 3 then mult_alpha_inv[i]=1./alpha[i]/alpha[i-1]/alpha[i-2]/alpha[i-3]
    if i eq 4 then mult_alpha_inv[i]=1./alpha[i]/alpha[i-1]/alpha[i-2]/alpha[i-3]/alpha[i-4]
    if i eq 5 then mult_alpha_inv[i]=1./alpha[i]/alpha[i-1]/alpha[i-2]/alpha[i-3]/alpha[i-4]/alpha[i-5]
  endfor
sum_mult_alpha_inv=total(mult_alpha_inv)
mp[0]=Mtot_planets/(1.+sum_mult_alpha_inv)   ;Mass of inner planet
  x=mp[0]      ;Change of variable
  rgap=routgap-ringap   ;Gap size
  for i=1,np-1 do begin
    mp[i]=x*mult_alpha_inv[i-1]          ;Mass of planets (not include the inner planet)
  endfor
;  a[0]=4.*(x/(3.*msource))^(1./3.)   ;Distance from gap inner border to inner planet
;  rp[0]=ringap/(1.-a[0])            ;Inner planet position
  a[np-1]=8.*(mp[np-1]/(3.*msource))^(1./3.)   ;
  rp[np-1]=routgap/(1.+a[np-1])                 ;Position of extern planet
  for i=np-2,0,-1 do begin
;     rp[i]=(2.)^(2./3.)*rp[i-1]     ;Planets position (not include the inner planet)
     rp[i]=rp[i+1]/((2.)^(2./3.))    
  endfor
;  outer=rp[np-1]+(8.*(mp[np-1]/(3.*msource))^(1./3.)*rp[np-1])  ;Probe of gap open by planets
  print,rp,total(mp),mp,mp/1.87		;case separation in resonance

;Hereafter we make the calculation if is considered a separation between planets by mutual Hill Radius

  ;; a[0]=4.*(x/(3.*msource))^(1./3.)
  ;; a[1]=4.*(x/(3.*alpha[0]*msource))^(1./3.)
  ;; a[2]=4.*(x/(3.*alpha[0]*alpha[1]*msource))^(1./3.)
  ;; a[3]=4.*(x/(3.*alpha[0]*alpha[1]*alpha[2]*msource))^(1./3.)
  ;; rp[0]=ringap/(1.-a[0])
  ;; rp[1]=(ringap+(a[0]*rp[0]))/(1.-a[1])
  ;; rp[2]=((1.+a[1])/(1.-a[2]))*rp[1]
  ;; rp[3]=((1.+a[2])/(1.-a[3]))*rp[2]
  ;; probe_gap=2*((rp[0]*a[0])+(rp[1]*a[1])+(rp[2]*a[2])+(rp[3]*a[3]))
;  print,rp,mp,probe_gap,rgap

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  a[0]=2.*(x/(3.*msource))^(1./3.)
;  a[1]=2.*(x/(3.*alpha[0]*msource))^(1./3.)
;  a[2]=2.*(x/(3.*alpha[0]*alpha[1]*msource))^(1./3.)
;  a[3]=2.*(x/(3.*alpha[0]*alpha[1]*alpha[2]*msource))^(1./3.)
;  rp[0]=ringap/(1.-a[0])
;  rp[1]=(2.)^(2./3.)*rp[0]
;  rp[2]=(2.)^(2./3.)*rp[1]
;  rp[3]=(2.)^(2./3.)*rp[2]
;  distin1=rp[0]-ringap
;  dist12=rp[1]-rp[0]
;  dist23=rp[2]-rp[1]
;  dist34=rp[3]-rp[2]
;  dist4out=routgap-rp[3]
;  print,rp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  b=routgap/(ringap*0.25)
;  probe=(1.+a[3])/(1.-a[0])
;print,probe,b
;  for i=1,np-1 do begin
;     mp[i]=mp[i-1]/alpha[i-1]
;  endfor

;  rp[0]=ringap/(1.-(5.*((mp[0]/(3.*msource))^(1./3.))))
;  rp[np-1]=routgap/(1.+(5.*((mp[np-1]/(3.*msource))^(1./3.))))

;print,mp
;Routine for calculate the mass and radius of hte planet, considering
;that all planets have the same mass
  ;; for N=1,np do begin
  ;;    new_sum=0.
  ;;    for i=1,N do begin
  ;;       rp[i,N]=ringap+(i*rgap/(N+1))
  ;;       new_sum+=rp[i,N]
  ;;    endfor
  ;;    mp[N]=(((rgap/(10.*new_sum)))^3.)*3.*msource
  ;; endfor
  ;; rp(0,*)=mp(*)
  ;; print,rp

  ;; b=rgap/(10.*(3.*msource)^(-1./3.))
  ;; f12=1.2                         ;f=(5./3.)^(1./3.)
  ;; m12=0.02
  ;; m22=f12^(3.)*m12
  ;; r12=ringap/(1.-(5.*((m12/(3.*msource))^(1./3.))))
  ;; r22=routgap/(1.+(5.*((m22/(3.*msource))^(1./3.))))
  ;; probegap=10.*(((m12/(3.*msource))^(1./3.)*r12)+((m22/(3.*msource))^(1./3.)*r22))


;  r22=((b/m12)^(1./3.)-r12)/f

 ; print,r12,r22,m12,m22,rgap,probegap
  ;; for N=1,np do begin
  ;;    for i=1,N do begin
        
  ;;    endfor
  ;; endfor
  
  RETURN
END
