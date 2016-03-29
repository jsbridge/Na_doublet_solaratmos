;+
;Joanna Bridge, 1/2012
;
;This program checks the validity of bnu_integ by creating different
;length arrays and plotting the result of the integral vs. the number 
;of values
; 
;-

pro bnu_check

result = dblarr(4)
value = lonarr(4)

for i=1,6 do begin
   wave = dindgen(10L^(i+1))/10+0.1  ; intervals of 0.1
   value[i-1] = n_elements(wave)
   integ = bnu_int(wave)
   result[i-1] = integ
endfor

result = result/(10.^10)

set_plot, 'ps'
device, file='result_v_num.ps'
plot, value, result, /xlog, thick=6, xthick=6, ythick=6, charthick=3, title='Result v. Size of Array', ytitle=textoidl('Result of Integral/10^{10} (erg cm^{-2})'), xtitle='Size of Wavenumber Array', yrange=[5.7090, 5.70905], ycharsize=.8, ystyle=1
device, /close


end
