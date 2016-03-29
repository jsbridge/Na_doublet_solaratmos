;+
; Joanna Bridge, 1/2012
;
; This program compares the numerical integral of the emergent
; intensity with a quadratic source function to the Eddington-Barbier
; approximation of the same source function. The emergent intensity is
; evaluated at the surface and mu=1.
;
;-

pro emerg_int

tau = dindgen(50)
a_2 = dindgen(100)
val = dblarr(n_elements(a_2))

for i = 0, n_elements(a_2)-1 do begin                
   f = (1.+tau+a_2[i]*tau^2)*exp(-tau)            ; dependent array (a_0=a_1=1) for each a_2 value
   val[i] = int_tabulated(tau, f, /double, /sort) ; integrate over tau for each f value
endfor
                          
source = 1. + 1. + a_2  ; create array of I from E-B approximation where I=S(tau=mu=1)

set_plot, 'ps'
device, file='EB_v_integI.ps'
plot, a_2, val, title=textoidl('Integrated I_\nu and Eddington-Barbier I_\nu v. a_2'), ytitle=textoidl('I_\nu'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5, psym=4, charsize=1.2, xcharsize=1.25, ycharsize=1.3, thick=2
oplot, a_2, source, psym=1, thick=2
legend, [textoidl(' Integrated I_\nu'), textoidl(' E-B I_\nu')], psym=[4, 1], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close

diff = (val-source)/((val+source)/2)   ; Calculate percent difference between E-B and integrated I

set_plot, 'ps'
device, file='percent_diff.ps'
plot, a_2, diff, title=textoidl('Percent Difference v. a_2'), ytitle=textoidl('Percent Difference'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5, thick=5, charsize=1.2, xcharsize=1.25, ycharsize=1.3, yrange=[0, 1], ystyle=1
device, /close


end
