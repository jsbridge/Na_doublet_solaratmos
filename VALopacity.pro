;+
;
; Calculates total opacity for the VALIIIC atmosphere
; 
; 
;
;-
pro VALopacity

restore, 'important.sav'
defsysv, '!matchabund', matchabund
defsysv, '!VALdata', VALdata
defsysv, '!partit', partit
defsysv, '!solarabund', solarabund
defsysv, '!ioniztn', ioniztn

lambda = 5000d
Pgas = double(!VALdata.Pgas_Ptotal * !VALdata.Ptotal)
Pe = Pe_VALIIIC()
;lambda = interpol([0d, 20000d], 1000d)
opacity = dblarr(n_elements(!VALdata.T))
elec = 0.6648d-24*total(!matchabund.abund)*(Pe)/((Pgas - Pe)*total(!matchabund.abund*!matchabund.weight*1.6606d-24))

for i = 0, n_elements(!VALdata.T)-1 do begin
   Hbf = neutralHbf(lambda, !VALdata.T[i], Pe[i])
   Hff = neutralHff(lambda, !VALdata.T[i], Pe[i])
   Hminbf = Hminusbf(lambda,!VALdata.T[i] , Pe[i])
   Hminff = Hminusff(lambda, !VALdata.T[i], Pe[i])
   opacity[i] = Hbf + Hff + Hminbf + Hminff + elec[i]
endfor

dtau = dblarr(51)
dPgas = dblarr(51)

for i = 0, n_elements(!VALdata.tau_500)-2 do begin
   dtau[i] = !VALdata.tau_500[i+1] - !VALdata.tau_500[i]
   dPgas[i] = Pgas[i+1] - Pgas[i]
endfor

HSEopacity = 10^4.4377d * dtau / dPgas

tau = !VALdata.tau_500[1:*]

set_plot, 'ps'
device, file='tau_v_opacity.ps'
plot, !VALdata.tau_500, opacity, /xlog, /ylog, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('\kappa_{500} vs. \tau_{500}'), xtitle=textoidl('\tau_{500}'), ytitle=textoidl('\kappa_{500} (cm^2/g)'), psym=-6,  xrange=[3d-8, 10], xstyle=1, symsize=0.75
oplot, tau, HSEopacity, psym=-2, thick=4
legend, [textoidl('HSE \kappa_{500}'), textoidl('Calculated \kappa_{500}')], psym=[2,6], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close
set_plot, 'x'

stop

end

