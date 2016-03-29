;+
;
; This attempts to reproduce Figure 8.8b
;
;-
pro solarabund

readcol, 'solarabund.txt', name, abund, abund1, logabund, format='A,D,I,D', /silent
VALdata = VALdata('VALIIIC.txt')
a = min(where(VALdata.h lt 900))
Pe = VALdata.n_e[a:*] * 1.3806503d-16 * VALdata.T[a:*]
temp = VALdata.T[a:*]
b = size(where(VALdata.h lt 900))
N_Mg = dblarr(b[1])
N_Fe = dblarr(b[1])
N_Si = dblarr(b[1])
N_H = dblarr(b[1])

c = where(name eq 'Mg')
for i = 0, b[1]-1 do begin
   N_Mg[i] = (10^logabund[c]* (saha('Mg', temp[i], 1, Pe[i])/Pe[i])/(1d + (saha('Mg', temp[i], 1, Pe[i])/Pe[i])))
endfor

d = where(name eq 'Fe')
for i = 0, b[1]-1 do begin
   N_Fe[i] = (10^logabund[d]* (saha('Fe', temp[i], 1, Pe[i])/Pe[i])/(1d + (saha('Fe', temp[i], 1, Pe[i])/Pe[i])))
endfor

f = where(name eq 'Si')
for i = 0, b[1]-1 do begin
   N_Si[i] = (10^logabund[f]* (saha('Si', temp[i], 1, Pe[i])/Pe[i])/(1d + (saha('Si', temp[i], 1, Pe[i])/Pe[i])))
endfor


g = where(name eq 'H')
for i = 0, b[1]-1 do begin
   N_H[i] = (10^logabund[g]* (saha('H', temp[i], 1, Pe[i])/Pe[i])/(1d + (saha('H', temp[i], 1, Pe[i])/Pe[i])))
endfor

e_Mg = VALdata.n_h[a:*] * N_Mg/VALdata.n_e[a:*]
e_Fe = VALdata.n_h[a:*] * N_Fe/VALdata.n_e[a:*]
e_Si = VALdata.n_h[a:*] * N_Si/VALdata.n_e[a:*]
e_H = VALdata.n_h[a:*] * N_H/VALdata.n_e[a:*]

set_plot, 'ps'
device, file='Fig8.8_b.ps'
plot, VALdata.h[a:*], e_Mg, xrange=[800, 0], /ylog, xstyle=1, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('Contribution to n_e vs. h'), xtitle='h (km)', ytitle=textoidl('Contribution to n_e'), yrange=[0.1, 1], ystyle=1
oplot, VALdata.h[a:*], e_Si, linestyle=1, thick=4
oplot, VALdata.h[a:*], e_Fe, linestyle=3, thick=4
oplot, VALdata.h[a:*], e_H, linestyle=5, thick=4
legend, ['H', 'Mg', 'Fe','Si'], linestyle=[5, 0, 3, 1], thick=4, charthick=5, charsize=1, bthick=5, pos=[785,.0275]
device, /close
set_plot, 'x'

end
