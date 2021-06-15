Pro spectrum

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Constants^  
  parsec=double(3.08572d18)
  cc=double(2.99792458d10)      ; Light speed  [cm/s]
  SEDcgsToMKS=1.e-3

;Source parameters
  rsource=double(2.3)           ;Source Radi in Solar Radii 
  distance=double(385.0)             ;Distance to source in parsec
;  temp=double(4350.0)            ;Source Temperature
;  temp=double(5780.0)		;Solar Temperature
  temp=double(7250.0)		;Solar Temperature
  msource=double(1.86)

;The data folder, and the folder for discharge the generated images
  dataFolder='~/Dropbox/radmc-3d/version_0.27/examples/run_ppdisk_gui/spectrum.out'
;  dataFolder='/media/ramiro/MiniStation/Dropbox/radmc-3d/version_0.27/examples/run_ppdisk_gui/spectrum.out'
  dischargeImageFolder='~/Dropbox/visualization_tools/images/V1247_Ori'
;  dischargeImageFolder='/media/ramiro/MiniStation/Dropbox/visualization_tools/images/V1247_Ori'
;stop  
  spectrum=dblarr(2,100)
  spectrum_=dataFolder
  openr,1,spectrum_
  header=strarr(3)
  readf,1,header
  readf,1,spectrum
  close,1

;Stellar spectrum
  msun=double(1.9889200d33)
  rsun=double(6.96d10)          ;Solar Radii
  mstr=msource*msun
  rstr=rsource*rsun
  distance_cgs=distance*parsec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  temp=double(5780.)           ;Solar Temperature
  nu=dblarr(100)
  bplanck=dblarr(100)
  factor=dblarr(100)
  sstr=dblarr(100)
  for i=0,99 do begin
    nu[i]=cc*(1.0d4)/spectrum[0,i]
    bplanck[i]=double(1.47455d-47)*nu[i]*nu[i]*nu[i]/(exp(double(4.7989d-11)*nu[i]/temp)-1.0)+double(1.d-290)
    factor[i]=rstr*rstr*!pi/(distance_cgs*distance_cgs);(4.*!pi*mstr)
    sstr[i]=factor[i]*bplanck[i]
  endfor

;
  device,get_decompose=old_decomposed,decomposed=0
  loadct,5 
  window,0,xsize=600,ysize=600
;       ThisDevice=!D.name
;        Set_Plot,'z'
;        Device,Set_Resolution=[600,600],Set_Pixel_Depth=24
  plot,spectrum(0,*),spectrum(1,*)*nu(*)/distance/distance*SEDcgsToMKS,background=255,$
;  plot,spectrum(0,*),spectrum(1,*)*nu(*)/distance/distance*SEDcgsToMKS,background=255,$
       color=0,xtitle='!4k!X [!4l!Xm]',xstyle=1,ystyle=1,$
;       ytitle='!4m!XF!I!4m!X!N [erg cm!E-2!N s!E-1!N]',$
       ytitle='!4k!XF!I!4k!X!N [W m!E-2!N]',$
       charsize=1.5,/ylog,/xlog,xrange=[2.e-1,2.5e2],yrange=[2.e-14,5.e-12]
  write_png,strcompress(dischargeImageFolder+'/V1247_Ori_spectrum.png',/remove_all),$
            tvrd(0,0,600,600,0,true=1)


RETURN
END
