;+
; Joanna Bridge 1/2012
;
; This program compares the percent difference between the
; int_tabulated function and the analytic solution for various sample sizes
;-
pro emerg_int_check

a_2 = dindgen(100)
tau1 = dindgen(500)/10
tau2 = dindgen(50)
val1 = dblarr(n_elements(a_2))
val2 = dblarr(n_elements(a_2))

for i = 0, n_elements(a_2)-1 do begin                
   f = (1.+tau1+a_2[i]*tau1^2)*exp(-tau1)            ; dependent array (a_0=a_1=1) for each a_2 value
   val1[i] = int_tabulated(tau1, f, /double, /sort) ; integrate over tau for each f value
endfor

for i = 0, n_elements(a_2)-1 do begin                
   f = (1.+tau2+a_2[i]*tau2^2)*exp(-tau2)            ; dependent array (a_0=a_1=1) for each a_2 value
   val2[i] = int_tabulated(tau2, f, /double, /sort) ; integrate over tau for each f value
endfor

val3 = 2.*(a_2 + 1.)

set_plot, 'ps'
device, file='percent_diff2.ps'
plot, a_2, (val2-val3)/((val2+val3)/2), linestyle=4, /ylog, yrange=[1d-9, 1d-1], ystyle=1, title=textoidl('Percent Difference v. a_2'), ytitle=textoidl('Percent Difference'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5, thick=5, charsize=1.2, xcharsize=1.25, ycharsize=1.3
oplot, a_2, (val1-val3)/((val1+val3)/2), thick=5
legend, ['50 elements', '500 elements'], linestyle=[4, 0], thick=4, charthick=5, charsize=1.2, bthick=5, pos=[52, .5d-7]
device, /close

tau3=dindgen(1000)/10
tau4 = dindgen(700)/10
f = (1.+tau3+5*tau3^2)*exp(-tau3)
f1 = (1.+tau4+5*tau4^2)*exp(-tau4)
b = int_tabulated(tau3, f, /double, /sort)
b2 = int_tabulated(tau4, f1, /double, /sort)
print, b, b2, val2[5], val3[5]
end
