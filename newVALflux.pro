;+
;
; Calculates the emergent flux from the VALIIIC atmosphere
;
;-
pro newVALflux

;restore, 'important.sav'
;defsysv, '!matchabund', matchabund
;defsysv, '!VALdata', VALdata
;defsysv, '!partit', partit
;defsysv, '!solarabund', solarabund
;defsysv, '!ioniztn', ioniztn

Pgas = !VALdata.Ptotal * !VALdata.Pgas_Ptotal
lambda = interpol([5885d, 5900d], 1d3)
Pe = !VALdata.n_e * 1.3806503d-16 * !VALdata.T
nu = 2.99792458d18/lambda ; lambda in angstroms

opacityD1 = dblarr(n_elements(lambda), n_elements(!VALdata.T))
opacityD2 = dblarr(n_elements(lambda), n_elements(!VALdata.T))
cont_op = dblarr(n_elements(lambda), n_elements(!VALdata.T))

for i = 0, n_elements(lambda)-1 do begin
   for j = 0, n_elements(!VALdata.T)-1 do begin
      opacityD1[i, j] = sodium_D1(lambda[i], !VALdata.T[j], Pgas[j], Pe[j], !VALdata.v[j])
      opacityD2[i, j] = sodium_D2(lambda[i], !VALdata.T[j], Pgas[j], Pe[j], !VALdata.v[j])
      cont_op[i, j] = total_opacity(lambda[i], !VALdata.T[j], Pe[j], Pgas[j])
   endfor
endfor

opacity = opacityD1 + opacityD2 + cont_op

tau = dblarr(n_elements(lambda), n_elements(!VALdata.T))

for i = 0, n_elements(lambda)-1 do begin
   for j = 1, n_elements(!VALdata.T)-1 do begin
      tau[i, j] = abs(tsum_d(!VALdata.H[0:j] * 1d5, opacity[i, 0:j] * !VALdata.rho[0:j]))
   endfor
endfor

phi = newVALphi(tau, lambda)

flux = !pi * VALphi(tau, lambda)
newflux = !pi * phi.val


set_plot, 'ps'  
device, file='niftydoublet.ps'
plot, lambda[findel(lambda,5889.0017):findel(lambda,5891.0017)], newflux[findel(lambda,5895):findel(lambda,5897)]*1.3, xrange=[5889,5891], xstyle=1, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, ytitle=textoidl('Flux (erg cm^{-2} s^{-1})'), xtitle=textoidl('\Delta\lambda (')+String(197B)+')', title='Na I Lines, VALIIIC Atmosphere', xtickname=['-1.0', '-0.5', '0.0', '0.5', '1.0']
oplot,lambda[findel(lambda,5889):findel(lambda,5891)], newflux[findel(lambda,5889):findel(lambda,5891)]*1.3 , thick=4
oplot, lambda[findel(lambda,5889):findel(lambda,5891)], flux[findel(lambda,5889):findel(lambda,5891)], linestyle=4, thick=4
oplot, lambda[findel(lambda,5889):findel(lambda,5891)], flux[findel(lambda,5895):findel(lambda,5897)], linestyle=3, thick=4
legend, ['LTE D1', 'LTE D2', 'NLTE'], linestyle=[3,4,0], thick=4, charthick=5, charsize=1, bthick=5, pos=[5889.1, 2e-5]

;device, file='SvB.ps'
;plot, !VALdata.h, phi.S, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, ytitle=textoidl('B, S'), xtitle=textoidl('h (km)'), title='B and S vs. Height', linestyle=4, /ylog
;oplot, !VALdata.h, phi.B, thick=4, /ylog 
;legend, ['S', 'B'], linestyle=[4, 0], thick=4, charthick=5, charsize=1, bthick=5

device, /close


stop
end
