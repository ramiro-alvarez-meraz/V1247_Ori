  for i=0,nrad-1 do begin
     for j=0,nsec-1 do begin
        tau=9.6*dens[j,i]*2./2.25*1.e7 ;Optical depth Troilita,etc.
        
;Optically thin:
;        if (tau lt tauInner) then begin
;           l=0
;           Counter[i,l] += 1
;           SigmaTau[i,l] += dens[j,i]
;           DensTau[j,i,l]=dens[j,i]
;        endif
;        
;Optically transitional:
;        if (tau ge tauInner) and (tau le tauOuter) then begin
;           l=1
;           Counter[i,l] += 1
;           SigmaTau[i,l]+=dens[j,i]
;           DensTau[j,i,l]=dens[j,i]
;        endif
;        
;Optically thick:
;        if (tau gt tauOuter) then begin
;           l=2
;           Counter[i,l] += 1
;           SigmaTau[i,l]+=dens[j,i]
;           DensTau[j,i,l]=dens[j,i]
;        endif

;New optical depth ranges
        if (tau lt tau_flag[0]) then begin
           l=0
           Counter[i,l] += 1
           SigmaTau[i,l] += dens[j,i]
           DensTau[j,i,l]=dens[j,i]
        endif
        for ii=1,TauRanges-2 do begin
          if (tau ge tau_flag[ii-1]) and (tau le tau_flag[ii]) then begin
	     l=ii
             Counter[i,l] += 1
             SigmaTau[i,l]+=dens[j,i]
             DensTau[j,i,l]=dens[j,i]
	  endif
        endfor
        if (tau gt tau_flag[TauRanges-2]) then begin
           l=TauRanges-1
           Counter[i,l] += 1
           SigmaTau[i,l]+=dens[j,i]
           DensTau[j,i,l]=dens[j,i]
        endif
     endfor
  endfor

