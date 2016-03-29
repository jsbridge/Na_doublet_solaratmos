;+
;
; Calculates the emergent flux from the VALIIIC atmosphere
;
;-
pro VALflux

;restore, 'important.sav'
;defsysv, '!matchabund', matchabund
;defsysv, '!VALdata', VALdata
;defsysv, '!partit', partit
;defsysv, '!solarabund', solarabund
;defsysv, '!ioniztn', ioniztn

Pgas = !VALdata.Ptotal * !VALdata.Pgas_Ptotal
lambda = interpol([5885d, 5900d], 300d)
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

flux = !pi * VALphi(tau, lambda)

set_plot, 'ps'
device, file = 'sodiumdoublet.ps'
plot, lambda, flux, /ylog, xrange=[5885d, 5900d], xstyle=1, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('Na I D Doublet Flux (D_1 = 5896 '+String(197B)+', D_2 = 5890 '+String(197B)+')'), xtitle=textoidl('\lambda (')+String(197B)+')', ytitle=textoidl('Flux (erg cm^{-2} s^{-1})')

device, file='tauvstau.ps'
plot, !VALdata.tau_500, tau[100,*], /ylog, /xlog, xrange=[1d-8, 100], xstyle=1,$
 xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, $
title=textoidl('\tau_\nu vs. \tau_{500} for \lambda = 5890 '+String(197B)), $
xtitle=textoidl('\tau_{500}'), ytitle=textoidl('\tau_\nu')

device, file='blownupD2.ps'
plot, lambda, flux, xrange=[5889,5891], xstyle=1, yrange=[0, 9d-5], ystyle=1, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('Na I D2 Flux, D_2 = 5890 '+String(197B)), xtitle=textoidl('\lambda (')+String(197B)+')', ytitle=textoidl('Flux (erg cm^{-2} s^{-1})')

device, file='tauvsnu.ps'
plot, nu, tau[*,48], /ylog, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('\tau_\nu = 0 Surface vs. \nu'), xtitle=textoidl('\nu (Hz)'), ytitle=textoidl('\tau_\nu')

device, /close
set_plot, 'x'


stop
end
