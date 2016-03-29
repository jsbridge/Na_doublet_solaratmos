;+
;
; Joanna Bridge, 2/2012
; 
; This program calculates the opacity of an atmosphere.  Put lambda in Angstroms
;
;-
pro opacity, lambda, temp, Pe

Hbf = neutralHbf(lambda, temp, Pe)
Hff = neutralHff(lambda, temp, Pe)
Hminbf = Hminusbf(lambda, temp, Pe)
Hminff = Hminusff(lambda, temp, Pe)

opacity = Hbf + Hff + Hminbf*1d-18 + Hminff

set_plot, 'ps'

;device, file='a.ps'
;plot, lambda, opacity/(1d-26*Pe), title=textoidl('Opacity vs. Wavelength'), xtitle=textoidl('\lambda (Angstroms)'), ytitle=textoidl('\kappa_\nu/P_e cm^2/H atom per dyne/cm^2, unit=10^{-26}'), linestyle=0, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, xrange=[2000, 20000], xstyle=1, xcharsize=.87
;oplot, lambda, (Hff+Hbf)/(1d-26*Pe), linestyle=2, thick=4
;oplot, lambda, Hminff/(1d-26*Pe), linestyle=3, thick=4
;oplot, lambda, Hminbf*1d-18/(1d-26*Pe), linestyle=5, thick=4
;legend, [textoidl('\kappa_{total}'), textoidl('\kappa_H'), textoidl('\kappa_{H^-_{ff}}'), textoidl('\kappa_{H^-_{bf}}')], linestyle=[0, 2, 3, 5], thick=4, charthick=5, charsize=1.2, bthick=5, pos=[1.4d4, 7.6]

;device, file='b.ps'
;plot, lambda, opacity/(1d-26*Pe), title=textoidl('Opacity vs. Wavelength'), xtitle=textoidl('\lambda (Angstroms)'), ytitle=textoidl('\kappa_\nu/P_e cm^2/H atom per dyne/cm^2, unit=10^{-26}'), linestyle=0, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, xrange=[2000, 20000], xstyle=1, xcharsize=.87
;oplot, lambda, (Hff+Hbf)/(1d-26*Pe), linestyle=2, thick=4
;oplot, lambda, Hminff/(1d-26*Pe), linestyle=3, thick=4
;oplot, lambda, Hminbf*1d-18/(1d-26*Pe), linestyle=5, thick=4
;legend, [textoidl('\kappa_{total}'), textoidl('\kappa_H'), textoidl('\kappa_{H^-_{ff}}'), textoidl('\kappa_{H^-_{bf}}')], linestyle=[0, 2, 3, 5], thick=4, charthick=5, charsize=1.2, bthick=5, pos=[1.4d4, 3.8]

;device, file='c.ps'
;plot, lambda, opacity/(1d-26*Pe), title=textoidl('Opacity vs. Wavelength'), xtitle=textoidl('\lambda (Angstroms)'), ytitle=textoidl('\kappa_\nu/P_e cm^2/H atom per dyne/cm^2, unit=10^{-26}'), linestyle=0, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, xrange=[2000, 20000], xstyle=1, xcharsize=.87, yrange=[0,5.2], ystyle=1
;oplot, lambda, (Hff+Hbf)/(1d-26*Pe), linestyle=2, thick=4
;oplot, lambda, Hminff/(1d-26*Pe), linestyle=3, thick=4
;oplot, lambda, Hminbf*1d-18/(1d-26*Pe), linestyle=5, thick=4
;legend, [textoidl('\kappa_{total}'), textoidl('\kappa_H'), textoidl('\kappa_{H^-_{ff}}'), textoidl('\kappa_{H^-_{bf}}')], linestyle=[0, 2, 3, 5], thick=4, charthick=5, charsize=1.2, bthick=5, pos=[1.4d4, 4.72]

device, file='d.ps'
plot, lambda, opacity/(1d-26*Pe), title=textoidl('Opacity vs. Wavelength'), xtitle=textoidl('\lambda (Angstroms)'), ytitle=textoidl('\kappa_\nu/P_e cm^2/H atom per dyne/cm^2, unit=10^{-26}'), linestyle=0, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, xrange=[2000, 20000], xstyle=1, xcharsize=.87
oplot, lambda, (Hff+Hbf)/(1d-26*Pe), linestyle=2, thick=4
oplot, lambda, Hminff/(1d-26*Pe), linestyle=3, thick=4
oplot, lambda, Hminbf*1d-18/(1d-26*Pe), linestyle=5, thick=4
legend, [textoidl('\kappa_{total}'), textoidl('\kappa_H'), textoidl('\kappa_{H^-_{ff}}'), textoidl('\kappa_{H^-_{bf}}')], linestyle=[0, 2, 3, 5], thick=4, charthick=5, charsize=1.2, bthick=5, pos=[3.1d3, 480]

device, /close
set_plot, 'x'
end
