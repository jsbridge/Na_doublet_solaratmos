;+
;
; Plot ALL the sodium lines
;
;-
pro sodium_plot, temp, Pgas, Pe, vel

lambda = interpol([5835d, 5955d], 1000d)
kappa1 = dblarr(n_elements(lambda))
kappa2 = dblarr(n_elements(lambda))
opc = dblarr(n_elements(lambda))

for i = 0, n_elements(lambda)-1 do begin
   kappa1[i] = sodium_D2(lambda[i], temp, Pgas, Pe, vel)
   kappa2[i] = sodium_D1(lambda[i], temp, Pgas, Pe, vel)
   opc[i] = total_opacity(lambda[i], temp, Pe, Pgas)
endfor

set_plot, 'ps'
device, file = 'opacityNa.ps'
plot, lambda, kappa1+kappa2, /ylog, linestyle=4, xrange=[5835d, 5955d], xstyle=1, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('Na I D Doublet Opacity (D_1 = 5896 '+String(197B)+', D_2 = 5890 '+String(197B)+')'), xtitle=textoidl('\lambda (')+String(197B)+')', ytitle=textoidl('\kappa (cm^2/g)')
oplot, lambda, opc+kappa1+kappa2, thick=4
legend, [textoidl('\kappa(continuum)+\kappa(Na I)'), textoidl('\kappa(Na I)')], linestyle=[0, 4], thick=4, charthick=5, charsize=.94, bthick=5
device, /close


set_plot, 'x'


end
