;+
;
; Joanna Bridge 2/2012
;
; Radiative equilibrium check
;
;-
pro rad_eq

tau = interpol([0,10d], 1000d)
wave = interpol([1d-6, 20d], 100d)  ;this is in inverse microns
nu = wave*2.99792458d10/1d-4   ;get into cm^-1 and mult. by cm/s
fnu = dblarr(n_elements(wave), n_elements(tau))
F_re = dblarr(n_elements(tau))

for i = 0, n_elements(wave)-1 do begin
   fnu[i,*] = !pi*phi(tau, nu[i])     ; curly F
endfor

F_eff = 5.6704e-5 * 8700d^4

for i = 0, n_elements(tau)-1 do begin
   F_re[i] = int_tabulated(nu, fnu[*,i], /double)
endfor

set_plot, 'ps'
device, file = 'flux_vs_tau.ps'
plot, tau, F_re, title=textoidl('Flux vs. \tau'), xthick=5, ythick=5, charthick=5,  charsize=1.2, ycharsize=.8, thick=4, xtitle=textoidl('\tau'), ytitle=textoidl('Flux (erg s^{-1} cm^{-2})'), yrange=[2.5d11,3.75d11]
hline, F_eff, thick=4

device, file='error.ps'
plot, tau, (F_eff-F_re)/F_eff, title='Relative Error in Flux', xthick=5, ythick=5, charthick=5, charsize=1.2, thick=4, ytitle='Relative Error', xtitle=textoidl('\tau')

device, /close

end
