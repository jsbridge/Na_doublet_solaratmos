;+
;
; Joanna Bridge, 2/2012
;
; This program evaluates integral of the exponential integral for n =
;1, 2, 3. It does it using both int_tabulated as well as tsum and
;plots the two relative errors as a function of sampling of the x array
;
;-
pro int_expint_len


d = interpol([15, 60], 86)
d1 = interpol([3, 60], 86)
int1=dblarr(n_elements(d))
int2=dblarr(n_elements(d))
int3=dblarr(n_elements(d))
num=dblarr(n_elements(d))
num1 = dblarr(n_elements(d1))

for i = 0, n_elements(d)-1 do begin
   x = interpol([0d, d[i]], ((d[i]+1)*10000d))
   x1 = interpol([0d, d1[i]], ((d1[i]+1)*10000d))
   num[i] = max(x)
   num1[i] = max(x)
   n1 = [(1-(expint(2, x1[1], /double))), (expint(1, x1[1:*], /double))]
   n2 = expint(2, x, /double)
   n3 = expint(3, x, /double)
   int1[i] = int_tabulated(x1, n1, /double, /sort)
   int2[i] = int_tabulated(x, n2, /double, /sort)
   int3[i] = int_tabulated(x, n3, /double, /sort)
endfor

err1 = abs(1d0-int1)
err2 = abs(((1d/2d)-int2)/(1d/2d))
err3 = abs(((1d/3d)-int3)/(1d/3d))

set_plot, 'ps'
device, file='error1_ch.ps'
plot, num1, err1, title=textoidl('Relative Error of E_1(x) v. X Array Length'), ytitle=textoidl('Relative Error of E_1(x)'), xtitle='X Array Length', xthick=5, ythick=5, charthick=5, psym=-4, charsize=1.2, xcharsize=1.15, ycharsize=.85, thick=4, xrange=[15,60], xstyle=1, yrange=[0.0003, 0.00038], ystyle=1
device, /close

set_plot, 'ps'
device, file='error2_ch.ps'
plot, num, err2, /ylog, title=textoidl('Relative Error of E_2(x) v. X Array Length'), ytitle=textoidl('Relative Error of E_2(x)'), xtitle='X Array Length', xthick=5, ythick=5, charthick=5, psym=-4, charsize=1.2, xcharsize=1.15, ycharsize=1.25, thick=4, xrange=[15,60], xstyle=1
device, /close

set_plot, 'ps'
device, file='error3_ch.ps'
plot, num, err3, /ylog, title=textoidl('Relative Error of E_3(x) v. X Array Length'), ytitle=textoidl('Relative Error of E_3(x)'), xtitle='X Array Length', xthick=5, ythick=5, charthick=5, psym=-4, charsize=1.2, xcharsize=1.15, ycharsize=1.25, thick=4, xrange=[15,60], xstyle=1
device, /close



end
