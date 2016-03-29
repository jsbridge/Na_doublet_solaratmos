;+
;
; Joanna Bridge 2/2012
;
; This function uses snu.pro and lambda.pro and plots J vs. S. Tau and
; a arrays are given, but program can be modified to accept arbitrary
; tau and a's.
;
;-
pro plotJvS

a = [1, 0, 0]
tau = interpol([0d, 3d], 100d)

J = lambda(tau, a)
S = snu(tau, a)

set_plot, 'ps'
device, file = 'JvS.ps'
plot,  tau, J, xrange=[0, 3], yrange=[0, 4], linestyle=3, xthick=5, ythick=5, charthick=5,  charsize=1.2, thick=4, title=textoidl('J_\nu vs. S_\nu for Linear Source Function'), xtitle=textoidl('\tau_\nu')
oplot, tau, S, thick=4
legend, [textoidl('S_\nu'), textoidl( 'J_\nu')], linestyle=[0, 3], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close


end
