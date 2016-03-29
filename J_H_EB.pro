;+
; Joanna Bridge, 1/2012
;
; This program compares the numerical integral of the J and H  with a 
; quadratic source function to the Eddington-Barbier
; approximation of the same source function.
;
;-
pro J_H_EB

tauJ = dindgen(500)/10+1d-8
tauH = dindgen(500)/10
a_2 = dindgen(100)
valJ = dblarr(n_elements(a_2))
valH = dblarr(n_elements(a_2))
f = dblarr(n_elements(tauJ))
f1 = dblarr(n_elements(tauH))

for i = 0, n_elements(a_2)-1 do begin
   for j = 0, n_elements(tauJ)-1 do begin
      E = expint(1, tauJ[j], /double)
      f[j] = 0.5*(1d0+tauJ[j]+a_2[i]*tauJ[j]^2)*E ;dependent array (a_0=a_1=1) for each a_2
   endfor
   valJ[i] = int_tabulated(tauJ, f, /double, /sort)
endfor

for i = 0, n_elements(a_2)-1 do begin
   for j = 0, n_elements(tauH)-1 do begin
      E = expint(2, tauH[j], /double)
      f1[j] = 0.5*(1d0+tauH[j]+a_2[i]*tauH[j]^2)*E
   endfor
   valH[i] = int_tabulated(tauH, f1, /double, /sort)
endfor         
               
J = 0.5 + 0.25 + a_2[1:*]*0.125  ; E-B approx. for J
H = (1d + (2d/3d) + a_2*(4d/9d))/(4d); E-B approx. for H


set_plot, 'ps'
device, file='EBJ_v_integral.ps'
plot, a_2, valJ, title=textoidl('Integrated J_\nu and Eddington-Barbier J_\nu v. a_2'), ytitle=textoidl('J_\nu'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5, charsize=1.2, xcharsize=1.25, ycharsize=1.3, thick=5
oplot, a_2, J, linestyle=3, thick=5
legend, [textoidl(' Integrated J_\nu'), textoidl(' E-B J_\nu')], linestyle=[0, 3], thick=5, charthick=5, charsize=1.2, bthick=5
device, /close

set_plot, 'ps'
device, file='EBH_v_integral.ps'
plot, a_2, valH, title=textoidl('Integrated H_\nu and Eddington-Barbier H_\nu v. a_2'), ytitle=textoidl('H_\nu'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5,  charsize=1.2, xcharsize=1.25, ycharsize=1.3, thick=5
oplot, a_2, H, linestyle=3, thick=5
legend, [textoidl(' Integrated H_\nu'), textoidl(' E-B H_\nu')], linestyle=[0, 3], thick=5, charthick=5, charsize=1.2, bthick=5
device, /close

diffJ = 100*(valJ-J)/valJ   ; Calculate percent differenceerror between E-B and integrated J and H
diffH = 100*(valH-H)/valH 

set_plot, 'ps'
device, file='percent_diffJ.ps'
plot, a_2, diffJ, title=textoidl('Percent Difference J v. a_2'), ytitle=textoidl('Percent Difference'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5, thick=5, charsize=1.2, xcharsize=1.25, ycharsize=1.3
device, /close

set_plot, 'ps'
device, file='percent_diffH.ps'
plot, a_2, diffH, title=textoidl('Percent Difference H v. a_2'), ytitle=textoidl('Percent Difference'), xtitle=textoidl('a_2'), xthick=5, ythick=5, charthick=5, thick=5, charsize=1.2, xcharsize=1.25, ycharsize=1.3
device, /close

end
