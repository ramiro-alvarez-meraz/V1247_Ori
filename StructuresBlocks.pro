
;Begin the count of blocks:
  for l=0,TauRanges-1 do begin
        openw,1,strcompress(dischargeImagesFolder+'/DensTau_l='+string(l),/remove_all)
        printf,1,DensTau(*,*,l)
        close,1

     tot=total(DensTau[*,*,l])
     if tot eq 0 then begin
        m=0
        Block_prev[l]=m
        goto,jump0
     endif else begin
        m=1
     endelse
     jump6:

;Encuentra la celda de inicio de la estructura
     for i=0,nrad-1 do begin
        for j=0,nsec-1 do begin
           if DensTau[j,i,l] gt 0 then begin
              j_[l]=j
              i_[l]=i
              Block[j,i,l,m]=m
              Block_prev[l]=Block[j,i,l,m]
              goto,jump1
           endif
        endfor
     endfor
     jump1:
     print,'l    =',l,'           m      =',uint(m)
;The structure is found in the first ring
        for i=0,nrad-1 do begin
           if (i eq i_[l]) then begin
              for j=0,nsec-1 do begin
                 if (j gt j_[l]) and (j le nsec-1) then begin
                    if (DensTau[j,i,l] gt 0) and (j lt nsec-1) then begin
                       Block[j,i,l,m]=Block_prev[l]
                       Block_prev[l]=Block[j,i,l,m]
                    endif
                    if (DensTau[j,i,l] gt 0) and (j eq nsec-1) then begin
                       Block[j,i,l,m]=Block[j-1,i,l,m]
                       Block_prev[l]=Block[j,i,l,m]
                    endif
                    if (DensTau[j,i,l] eq 0) then goto,jump2
                 endif
              endfor     
              jump2:
              for j=nsec-1,0,-1 do begin
                 if j eq nsec-1 then begin
                    if (Denstau[j,i,l] gt 0) $
                       and (Block[0,i,l,m] eq Block_prev[l]) then begin
                       Block[j,i,l,m]=Block_prev[l]
                       Block_prev[l]=Block[j,i,l,m]
                    endif
                 endif
                 if (j lt nsec-1) and (j gt 0) then begin
                    if (Denstau[j,i,l] gt 0) $
                       and (Block[j+1,i,l,m] eq Block_prev[l]) then begin
                       Block[j,i,l,m]=Block[j+1,i,l,m]
                       Block_prev[l]=Block[j,i,l,m]
                       if (Denstau[j,i,l] eq 0) then goto,jump3
                    endif
                 endif
              endfor
           endif
        endfor
        jump3:

;We try to find the structure squeleton 
        m=Block_prev[l]
        for i=0,nrad-1 do begin
           if (i gt i_[l]) then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (DensTau[j,i,l] gt 0) and $
                       ((Block[nsec-1,i-1,l,m] gt 0) $
                        or (Block[j,i-1,l,m] gt 0) $
                        or (Block[j+1,i-1,l,m] gt 0)) then begin
                       Block[j,i,l,m]=m
                       Block_prev[l]=Block[j,i,l,m]
                    endif                    
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (DensTau[j,i,l] gt 0) and $
                       ((Block[j-1,i-1,l,m] gt 0) $
                        or (Block[j,i-1,l,m] gt 0) $
                        or (Block[j+1,i-1,l,m] gt 0)) then begin                
                       Block[j,i,l,m]=m
                       Block_prev[l]=Block[j,i,l,m]
                    endif
                 endif
                 if (j eq nsec-1) then begin
                    if (DensTau[j,i,l] gt 0) and $
                       ((Block[j-1,i-1,l,m] gt 0) $
                        or (Block[j,i-1,l,m] gt 0) $
                        or (Block[0,i-1,l,m] gt 0)) then begin
                       Block[j,i,l,m]=m
                       Block_prev[l]=Block[j,i,l,m]
                    endif                    
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                    if (DensTau[j,i,l] gt 0) and $
                       ((Block[j-1,i-1,l,m] gt 0) $
                        or (Block[j,i-1,l,m] gt 0) $
                        or (Block[0,i-1,l,m] gt 0)) then begin
                       Block[j,i,l,m]=m
                       Block_prev[l]=Block[j,i,l,m]
                    endif                    
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (DensTau[j,i,l] gt 0) and $
                       ((Block[j-1,i-1,l,m] gt 0) $
                        or (Block[j,i-1,l,m] gt 0) $
                        or (Block[j+1,i-1,l,m] gt 0)) then begin                
                       Block[j,i,l,m]=m
                       Block_prev[l]=Block[j,i,l,m]
                    endif
                 endif
                 if (j eq 0) then begin
                    if (DensTau[j,i,l] gt 0) and $
                       ((Block[nsec-1,i-1,l,m] gt 0) $
                        or (Block[j,i-1,l,m] gt 0) $
                        or (Block[j+1,i-1,l,m] gt 0)) then begin
                       Block[j,i,l,m]=m
                       Block_prev[l]=Block[j,i,l,m]
                    endif                    
                 endif
              endfor
           endif
        endfor
        for i=nrad-2,0,-1 do begin
           for j=0,nsec-1 do begin
              if (j eq 0) then begin
                 if (DensTau[j,i,l] gt 0) and $
                    ((Block[nsec-1,i+1,l,m] gt 0) $
                     or (Block[j,i+1,l,m] gt 0) $
                     or (Block[j+1,i+1,l,m] gt 0)) then begin
                    Block[j,i,l,m]=m
                    Block_prev[l]=Block[j,i,l,m]
                 endif                    
              endif
              if (j gt 0) and (j lt nsec-1) then begin
                 if (DensTau[j,i,l] gt 0) and $
                    ((Block[j-1,i+1,l,m] gt 0) $
                     or (Block[j,i+1,l,m] gt 0) $
                     or (Block[j+1,i+1,l,m] gt 0)) then begin                
                    Block[j,i,l,m]=m
                    Block_prev[l]=Block[j,i,l,m]
                 endif
              endif
              if (j eq nsec-1) then begin
                 if (DensTau[j,i,l] gt 0) and $
                    ((Block[j-1,i+1,l,m] gt 0) $
                     or (Block[j,i+1,l,m] gt 0) $
                     or (Block[0,i+1,l,m] gt 0)) then begin
                    Block[j,i,l,m]=m
                    Block_prev[l]=Block[j,i,l,m]
                 endif                    
              endif
           endfor
           for j=nsec-1,0,-1 do begin
              if (j eq nsec-1) then begin
                 if (DensTau[j,i,l] gt 0) and $
                    ((Block[j-1,i+1,l,m] gt 0) $
                     or (Block[j,i+1,l,m] gt 0) $
                     or (Block[0,i+1,l,m] gt 0)) then begin
                    Block[j,i,l,m]=m
                    Block_prev[l]=Block[j,i,l,m]
                 endif                    
              endif
              if (j gt 0) and (j lt nsec-1) then begin
                 if (DensTau[j,i,l] gt 0) and $
                    ((Block[j-1,i+1,l,m] gt 0) $
                     or (Block[j,i+1,l,m] gt 0) $
                     or (Block[j+1,i+1,l,m] gt 0)) then begin                
                    Block[j,i,l,m]=m
                    Block_prev[l]=Block[j,i,l,m]
                 endif
              endif
              if (j eq 0) then begin
                 if (DensTau[j,i,l] gt 0) and $
                    ((Block[nsec-1,i+1,l,m] gt 0) $
                     or (Block[j,i+1,l,m] gt 0) $
                     or (Block[j+1,i+1,l,m] gt 0)) then begin
                    Block[j,i,l,m]=m
                    Block_prev[l]=Block[j,i,l,m]
                 endif                    
              endif
           endfor
        endfor
;To determine the lack of material of each structure
        for i=0,nrad-1 do begin 
           for j=0,nsec-1 do begin
              if (DensTau[j,i,l] gt 0) and $
                 (Block[j,i,l,m] ne Block_prev[l]) $
              then Block[j,i,l,m]=999
           endfor
        endfor

        jump4:
;The path of 0 to nrad-1 
        for i=0,nrad-1 do begin
           if i eq 0 then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]                    
                 endif
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l] 
                 endif
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
           endif
           if i gt 0 and i lt nrad-1 then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
           endif
           if i eq nrad-1 then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
           endif
        endfor

;The path of nrad-1 to 0
        for i=nrad-1,0,-1 do begin
           if i eq nrad-1 then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                   if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
           endif
           if i gt 0 and i lt nrad-1 then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[0,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[j-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m) $
                            or (Block[nsec-1,i-1,l,m] eq m) $
                            or (Block[j,i-1,l,m] eq m) $
                            or (Block[j+1,i-1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
           endif
           if i eq 0 then begin
              for j=0,nsec-1 do begin
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
              for j=nsec-1,0,-1 do begin
                 if (j eq nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[0,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[0,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[j-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[j-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
                 if (j eq 0) then begin
                    if (Block[j,i,l,m] eq 999) $
                       and ((Block[nsec-1,i,l,m] eq m) $
                            or (Block[j+1,i,l,m] eq m) $
                            or (Block[nsec-1,i+1,l,m] eq m) $
                            or (Block[j,i+1,l,m] eq m) $
                            or (Block[j+1,i+1,l,m] eq m)) $
                    then Block[j,i,l,m]=Block_prev[l]
                 endif
              endfor
           endif
        endfor
        
;Here begin the search of the missing cells
        for i=0,nrad-1 do begin
           if i eq 0 then begin
              for j=0,nsec-1 do begin 
                 if j eq 0 then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[nsec-1,i,l,m] eq m) $
                             or (Block[j+1,i,l,m] eq m) $
                             or (Block[nsec-1,i+1,l,m] eq m) $
                             or (Block[j,i+1,l,m] eq m) $
                             or (Block[j+1,i+1,l,m] eq m)) $
                     then goto,jump4                      
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[j-1,i,l,m] eq m) $
                             or (Block[j+1,i,l,m] eq m) $
                             or (Block[j-1,i+1,l,m] eq m) $
                             or (Block[j,i+1,l,m] eq m) $
                             or (Block[j+1,i+1,l,m] eq m)) $
                     then goto,jump4
                 endif
                 if j eq nsec-1 then begin
                    if (Block[j,i,l,m] eq 999) $
                        and ((Block[j-1,i,l,m] eq m) $
                             or (Block[0,i,l,m] eq m) $
                             or (Block[j-1,i+1,l,m] eq m) $
                             or (Block[j,i+1,l,m] eq m) $
                             or (Block[0,i+1,l,m] eq m)) $
                     then goto,jump4                     
                 endif
              endfor
           endif
           if (i gt 0) and (i lt nrad-1) then begin
              for j=0,nsec-1 do begin 
                 if j eq 0 then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[nsec-1,i,l,m] eq m) $
                             or (Block[j+1,i,l,m] eq m) $
                             or (Block[nsec-1,i+1,l,m] eq m) $
                             or (Block[j,i+1,l,m] eq m) $
                             or (Block[j+1,i+1,l,m] eq m) $
                             or (Block[j+1,i-1,l,m] eq m) $
                             or (Block[j,i-1,l,m] eq m) $
                             or (Block[nsec-1,i-1,l,m] eq m)) $
                     then goto,jump4                    
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[j-1,i,l,m] eq m) $
                             or (Block[j+1,i,l,m] eq m) $
                             or (Block[j-1,i+1,l,m] eq m) $
                             or (Block[j,i+1,l,m] eq m) $
                             or (Block[j+1,i+1,l,m] eq m) $
                             or (Block[j+1,i-1,l,m] eq m) $
                             or (Block[j,i-1,l,m] eq m) $
                             or (Block[j-1,i-1,l,m] eq m)) $
                     then goto,jump4
                 endif
                 if j eq nsec-1 then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[j-1,i,l,m] eq m) $
                             or (Block[0,i,l,m] eq m) $
                             or (Block[j-1,i+1,l,m] eq m) $
                             or (Block[j,i+1,l,m] eq m) $
                             or (Block[0,i+1,l,m] eq m) $
                             or (Block[0,i-1,l,m] eq m) $
                             or (Block[j,i-1,l,m] eq m) $
                             or (Block[j-1,i-1,l,m] eq m)) $
                     then goto,jump4                    
                 endif
              endfor
           endif
           if j eq nrad-1 then begin
              for j=0,nsec-1 do begin 
                 if (j eq 0) then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[nsec-1,i,l,m] eq m) $
                             or (Block[j+1,i,l,m] eq m) $
                             or (Block[j+1,i-1,l,m] eq m) $
                             or (Block[j,i-1,l,m] eq m) $
                             or (Block[nsec-1,i-1,l,m] eq m)) $
                     then goto,jump4                    
                 endif
                 if (j gt 0) and (j lt nsec-1) then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[j-1,i,l,m] eq m) $
                             or (Block[j+1,i,l,m] eq m) $
                             or (Block[j+1,i-1,l,m] eq m) $
                             or (Block[j,i-1,l,m] eq m) $
                             or (Block[j-1,i-1,l,m] eq m)) $
                     then goto,jump4
                 endif
                 if (j eq nsec-1) then begin
                     if (Block[j,i,l,m] eq 999) $
                        and ((Block[j-1,i,l,m] eq m) $
                             or (Block[0,i,l,m] eq m) $
                             or (Block[0,i-1,l,m] eq m) $
                             or (Block[j,i-1,l,m] eq m) $
                             or (Block[j-1,i-1,l,m] eq m)) $
                     then goto,jump4                    
                 endif
              endfor
           endif
        endfor
     jump0:
     Block_prov1=Block
     m_prev=Block_prev[l]

;Aqui termina la busqueda de los elementos de la estructura
     for i=0,nrad-1 do begin
        for j=0,nsec-1 do begin
           if (Block[j,i,l,m] eq 999) $
           then Block[j,i,l,m]=0
        endfor
     endfor
     Block_prov2=Block

;Para determinar existencia de otra estructura
     for i=0,nrad-1 do begin
        for j=0,nsec-1 do begin
           if (Block_prov1[j,i,l,m] eq 999) then begin
              Block_prev[l]+=1
              goto,jump5
           endif
        endfor
     endfor
     jump5:

;Hace un salto hacia la estructura, en caso de haberse encontrado
     m=Block_prev[l]
     if (m gt m_prev) then begin
        for i=0,nrad-1 do begin
           for j=0,nsec-1 do begin
              if Block_prov2[j,i,l,m_prev] gt 0 $
              then DensTau[j,i,l]=0
           endfor
        endfor
        goto,jump6
     endif
  endfor
