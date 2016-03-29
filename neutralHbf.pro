;+
;
; Joanna Bridge
;
; This calculates the opacity for neutral bound free Hydrogen. Lambda
; should be in angstroms (and start at 0) (and also have enough points
; (like 1/lambda or close))
;
;-
function neutralHbf, lambda, temp, Pe

a0 = 1.0449d-26               ; when lambda is in Angstroms!
R = 1.0968d-3                 ; Rydberg in Angstroms^-1
I = 13.5984d                  ; eV 
theta = 5040d/temp
loge = 0.43429d
;lam = interpol([0d, 20000d], 20000)
lam = dindgen(20000)
chi_lam = 1.2398d4/lam

for j = 1, 5 do begin
   if j eq 1 then begin
      a = max(where(lam le 912))
      val = 0
      for n = 1, 3 do begin
         gbf = 1d - (0.3456d/((lam[0:a]*R)^(1d/3d)))*(lam[0:a]*R/n^2 - 0.5d)
         chi = I*(1d - 1d/n^2)  ; eV
         val += gbf*(10^(-theta*chi))/n^3
      endfor
      
      chi3 = I*(1d - 1d/(j+3d)^2) 
      Hbf1 = a0*((lam[0:a])^3)*(val + (loge/(2d*theta*I))*(10^(-theta*chi3) - 10^(-theta*I)))
   endif

   if j eq 2 then begin
      a = max(where(lam le 912))
      b = max(where(lam lt 3746))
      val = 0
      for n = 2, 4 do begin
         gbf = 1d - (0.3456d/((lam[a+1:b]*R)^(1d/3d)))*(lam[a+1:b]*R/n^2 - 0.5d)
         chi = I*(1d - 1d/n^2)  ; eV
         val += gbf*(10^(-theta*chi))/n^3
      endfor
      
      chi3 = I*(1d - 1d/(j+3d)^2) 
      Hbf2 = a0*((lam[a+1:b])^3)*(val + (loge/(2d*theta*I))*(10^(-theta*chi3) - 10^(-theta*I)))
   endif

   if j eq 3 then begin
       a = max(where(lam le 3746))
       b = max(where(lam lt 8206))
       val = 0
       for n = 3, 5 do begin
         gbf = 1d - (0.3456d/((lam[a+1:b]*R)^(1d/3d)))*(lam[a+1:b]*R/n^2 - 0.5d)
         chi = I*(1d - 1d/n^2)  ; eV
         val += gbf*(10^(-theta*chi))/n^3
      endfor
      
      chi3 = I*(1d - 1d/(j+3d)^2)
      Hbf3 = a0*((lam[a+1:b])^3)*(val + (loge/(2d*theta*I))*(10^(-theta*chi3) - 10^(-theta*I)))
   endif

   if j eq 4 then begin
      a = max(where(lam le 8206))
      b = max(where(lam lt 14588))
      val = 0
      for n = 4, 6 do begin
         gbf = 1d - (0.3456d/((lam[a+1:b]*R)^(1d/3d)))*(lam[a+1:b]*R/n^2 - 0.5d)
         chi = I*(1d - 1d/n^2)  ; eV
         val += gbf*(10^(-theta*chi))/n^3
      endfor
      
      chi3 = I*(1d - 1d/(j+3d)^2) 
      Hbf4 = a0*((lam[a+1:b])^3)*(val + (loge/(2d*theta*I))*(10^(-theta*chi3) - 10^(-theta*I)))
   endif

   if j eq 5 then begin
      a = max(where(lam lt 14588))
      val = 0
      for n = 5, 7 do begin
         gbf = 1d - (0.3456d/((lam[a+1:*]*R)^(1d/3d)))*(lam[a+1:*]*R/n^2 - 0.5d)
         chi = I*(1d - 1d/n^2)  ; eV
         val += gbf*(10^(-theta*chi))/n^3
      endfor
      
      chi3 = I*(1d - 1d/(j+3d)^2)
      Hbf5 = a0*((lam[a+1:*])^3)*(val + (loge/(2d*theta*I))*(10^(-theta*chi3) - 10^(-theta*I)))
   endif

endfor

c =findel(lambda, lam)
Hbf_tot = [Hbf1, Hbf2, Hbf3, Hbf4, Hbf5]*(1d - 10^(-1d*chi_lam*theta))/(1d + (saha('H', temp, 1, Pe))[0])

Hbf = Hbf_tot[c]/total(!matchabund.abund*!matchabund.weight*1.6606d-24)

return, Hbf

end
