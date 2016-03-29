;+ 
;
; Joanna Bridge
;
; Calculates neutral Hydrogen free free opacity
;
;-
function neutralHff, lambda,  temp, Pe

a0 = 1.0449d-26               ; when lambda is in Angstroms!
R = 1.0968d-3                 ; Rydberg in Angstroms^-1
I = 13.5984d                  ; eV 
loge = 0.43429d
chi_lam = 1.2398d4/lambda
theta = 5040d/temp

gff = 1d + (0.3456d/((lambda*R)^(1d/3d)))*(loge/(theta*chi_lam) + 0.5d   )

Hff = a0*lambda^3*gff*(loge/(2*theta*I))*10^(-I*theta)

return, Hff*(1d - 10^(-1d*chi_lam*theta))/(1d + (saha('H', temp, 1, Pe)))/total(!matchabund.abund*!matchabund.weight*1.6606d-24)

end
