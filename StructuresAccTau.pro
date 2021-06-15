  ;; window,8,xsize=600,ysize=600
  ;; plot,rad,abs(AccTau(*,0)),background=255,$                            
  ;;      color=0,xtitle='Semieje mayor',yrange=[1.e-12,1.e-5],psym=0,xstyle=1,$
  ;;      ytitle='AccTau',charsize=1.5,/xlog,/ylog
  ;; xyouts,420,530,'       tau < 0.5',charsize=1.5,color=0,/device
  ;; oplot,rad,abs(AccTau(*,1)),color=100,psym=0
  ;; xyouts,420,500,'0.5 < tau < 3.0',charsize=1.5,color=100,/device
  ;; oplot,rad,abs(AccTau(*,2)),color=200,psym=0
  ;; xyouts,420,470,'3.0 < tau',charsize=1.5,color=200,/device
  ;; write_png,strcompress('out'+string(k)+'/imagenes/out'+string(k)+$
  ;;                       'output'+string(m)+'AccTau-SemiejeMayor.png',/remove_all),tvrd(0,0,600,600,0,true=1)
