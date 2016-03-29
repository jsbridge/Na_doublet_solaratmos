;+
;
; This programs reproduces Rutten Figure 8.8
;
;-
pro density

VALdata = VALdata('VALIIIC.txt')
a = min(where(VALdata.h lt 900))
Pe = VALdata.n_e[a:*] * 1.3806503d-16 * VALdata.T[a:*]
temp = VALdata.T[a:*]
b = size(where(VALdata.h lt 900))
saha = dblarr(b[1])

for i = 0, b[1]-1 do begin
   saha[i] = saha('H', temp[i], 1, Pe[i])
endfor

n_p = (saha * VALdata.n_H[a:*])/(Pe)


set_plot, 'ps'
device, file='Fig8.8_a.ps'
plot, VALdata.h[a:*], alog10(n_p), xrange=[800, 0], xstyle=1, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title='Log n vs. h', xtitle='h (km)', ytitle='Log n', linestyle=4, yrange=[8,14], ystyle=1
oplot, VALdata.h[a:*], alog10(VALdata.n_e[a:*]), thick=4
legend, [textoidl('n_e'), textoidl('n_p')], linestyle=[0,4], thick=4, charthick=5, charsize=1.2, bthick=5

device, /close
set_plot, 'x'

end
