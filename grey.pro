;+
;
; Joanna Bridge, 2/2012
; 
; Does grey atmosphere stuff
;
;-
pro grey

;tau=0
tau = interpol([0,8d], 100d)
wave = interpol([1d-6, 12d], 100d)  ;this is in inverse microns
nu = wave*2.99792458d10/1d-4   ;get into cm^-1 and mult. by cm/s
jnu = dblarr(n_elements(wave), n_elements(tau))
fnu = dblarr(n_elements(wave), n_elements(tau))
snu = dblarr(n_elements(wave), n_elements(tau))
temp = 8700d*(0.75*tau+0.5)^(1d/4d)
;print, wave[2], wave[80]

for i = 0, n_elements(wave)-1 do begin
;   jnu[i,*] = lambda(tau, nu[i])
   fnu[i,*] = !pi*phi(tau, nu[i])
;   snu[i,*] = bnu_func(nu[i], temp)
endfor

;second derivative wrt S at tau=1
x = findel(1, tau)
data = snu[*,x]
dx = deriv(tau, data)
dx2 = deriv(tau, dx)

set_plot, 'ps'
device, file='deriv.ps'
plot, wave, dx2, xthick=5, ythick=5, charthick=5,  charsize=1.2, ycharsize=.9 , thick=4, title=textoidl('d^2S_\nu/d\tau^2|_{\tau=1} vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('d^2S_\nu/d\tau^2|_{\tau=1}')

;E-B for J and F
J = (1d/2d)* bnu_func(nu, 8700d*(0.75*0.5+0.5)^(1d/4d))
F = !pi* bnu_func(nu, 8700d*(0.75*(2d/3d)+0.5)^(1d/4d))

J_diff = (jnu-J)
F_diff = (fnu-F)


;device, file = 'J_diff.ps'
;plot, wave, J_diff, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, title=textoidl('Difference in J_\nu(0) vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('Difference')

;device, file='F_diff.ps'
;plot, wave, F_diff, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, title=textoidl('Difference in F_\nu(0) vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('Difference')



;device, file='jnu_v_wave.ps', /color
loadct, 39
;plot, wave, jnu, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, title=textoidl('J_\nu(0) vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('J_\nu(0)')
;oplot, wave, J, thick=4, color='120'
;legend, ['E-B Approximation', 'Integrated'], linestyle=[0,0], colors=[120, 0], thick=4, charthick=5, charsize=1.2, bthick=5,  pos=[5,5e-5]

;device, file='fnu_v_wave.ps'
;plot, wave, fnu, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, title=textoidl('F_\nu(0) vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('F_\nu(0)')
;oplot, wave, F, thick=4, color='120'
;legend, ['E-B Approximation', 'Integrated'], thick=4, charthick=5, charsize=1.2, bthick=5,linestyle=[0,0], colors=[120, 0], pos=[5,3.2e-4]

device, file='fnu_v_tau(1).ps', /color
plot, tau, jnu[80,*], title=textoidl('F_\nu, J_\nu, S_\nu vs. \tau_\nu at High Frequency'), xtitle=textoidl('\tau_\nu'), xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, linestyle=5
oplot, tau, fnu[80,*], thick=4, linestyle=5, color='200'
oplot, tau, snu[80,*], thick=4, linestyle=5, color='120'
legend, [textoidl('J_\nu'), textoidl('F_\nu'), textoidl('S_\nu')], linestyle=[5, 5 , 5], colors=[0,200, 120],thick=4, charthick=5, charsize=1.2, bthick=5 

;device, file='fnu_v_tau(2).ps'
;plot, tau, jnu[50,*], title=textoidl('F_\nu, J_\nu, S_\nu vs. \tau_\nu at Mid Wavenumber'), xtitle=textoidl('\tau_\nu'),  xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4
;oplot, tau, fnu[50,*],  thick=4, linestyle=3
;oplot, tau, snu[50,*], thick=4, linestyle=5
;legend, [textoidl(' F_\nu'), textoidl('J_\nu'), textoidl('S_\nu')], linestyle=[3, 0, 5], thick=4, charthick=5, charsize=1.2, bthick=5 

device, file='fnu_v_tau(3).ps', /color
plot, tau, jnu[2,*], title=textoidl('F_\nu, J_\nu, S_\nu vs. \tau_\nu at Low Frequency'), xtitle=textoidl('\tau_\nu'),  xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, linestyle=5
oplot, tau, fnu[2,*],  thick=4, linestyle=5, color='200'
oplot, tau, snu[2,*], thick=4, linestyle=5, color='120'
legend, [textoidl(' J_\nu'), textoidl('F_\nu'), textoidl('S_\nu')],linestyle=[5, 5, 5],colors=[0,200, 120], thick=4, charthick=5, charsize=1.2, bthick=5  , pos=[5, 1e-5]
device, /close

loadct, 0
set_plot, 'x'

end
