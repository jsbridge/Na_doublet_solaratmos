;+
;
; Joanna Bridge, 2/2012
;
; This program evaluates integral of the exponential integral for n =
;1, 2, 3. It does it using both int_tabulated as well as tsum and
;plots the two relative errors as a function of sampling of the x array
;
;-
pro int_expint

b = interpol([10d, 500000d], 100d)
int1=dblarr(n_elements(b))
int2=dblarr(n_elements(b))
int3=dblarr(n_elements(b))
num=dblarr(n_elements(b))
int1t=dblarr(n_elements(b))
int2t=dblarr(n_elements(b))
int3t=dblarr(n_elements(b))

for i = 0, n_elements(b)-1 do begin
   x = interpol([0d, 100d], b[i])
   num[i] = n_elements(x)
   n1 = [(1-(expint(2, x[1], /double))), (expint(1, x[1:*], /double))]
   n2 = expint(2, x, /double)
   n3 = expint(3, x, /double)
   int1[i] = int_tabulated(x, n1, /double, /sort)
   int2[i] = int_tabulated(x, n2, /double, /sort)
   int3[i] = int_tabulated(x, n3, /double, /sort)
   int1t[i] = tsum_d(x, n1)
   int2t[i] = tsum_d(x, n2)
   int3t[i] = tsum_d(x, n3)
endfor

err1 = abs(1d0-int1)
err2 = abs(((1d/2d)-int2)/(1d/2d))
err3 = abs(((1d/3d)-int3)/(1d/3d))
err1t = abs(1d0-int1t)
err2t = abs(((1d/2d)-int2t)/(1d/2d))
err3t = abs(((1d/3d)-int3t)/(1d/3d))

set_plot, 'ps'
device, file='error1.ps'
plot, num, err1, /ylog, title=textoidl('Relative Error of E_1(x) v. X Array Size'), ytitle=textoidl('Relative Error of E_1(x)'), xtitle='X Array Size', xthick=5, ythick=5, charthick=5, psym=-4, charsize=1.2, xcharsize=1.15, ycharsize=1.25, thick=4
oplot, num, err1t, psym=-1, thick=4
legend, ['Trapezoid Rule', 'Newton-Coates'], psym=[1, 4], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close

set_plot, 'ps'
device, file='error2.ps'
plot, num, err2, /ylog, title=textoidl('Relative Error of E_2(x) v. X Array Size'), ytitle=textoidl('Relative Error of E_2(x)'), xtitle='X Array Size', xthick=5, ythick=5, charthick=5, psym=-4, charsize=1.2, xcharsize=1.15, ycharsize=1.25, thick=4
oplot, num, err2t, psym=-1, thick=4
legend, ['Trapezoid Rule', 'Newton-Coates'], psym=[1, 4], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close

set_plot, 'ps'
device, file='error3.ps'
plot, num, err3, /ylog, title=textoidl('Relative Error of E_3(x) v. X Array Size'), ytitle=textoidl('Relative Error of E_3(x)'), xtitle='X Array Size', xthick=5, ythick=5, charthick=5, psym=-4, charsize=1.2, xcharsize=1.15, ycharsize=1.25, thick=4
oplot, num, err3t, psym=-1, thick=4
legend, ['Trapezoid Rule', 'Newton-Coates'], psym=[1, 4], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close



end
