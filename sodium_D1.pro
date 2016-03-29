;+
;
; Models the opacity of the NaI D1 lines (lambda = 5896)
; Question: when doing two lines, do you sum the opacities?
;
;-
function sodium_D1, lambda, temp, Pgas, Pe, vel

;restore, 'important.sav'
;defsysv, '!matchabund', matchabund
;defsysv, '!VALdata', VALdata
;defsysv, '!partit', partit
;defsysv, '!solarabund', solarabund
;defsysv, '!ioniztn', ioniztn

chi_lam = 1.2398d4/lambda ; lambda in angstroms
abund_sum = total(!matchabund.abund*!matchabund.weight)
theta = 5040d/temp
m_Na = 22.99d * 1.66053886d-24 ; grams
lam_core = 5896d   ; angstroms
e = 4.803242d-10   ; cgs
c = 2.99792458d10  ; cm/s
k = 1.3806503d-16  ; ergs/K
g_l = 2d  
g_u = 2d
A_ul = 1.23d8/g_u

x = where(!matchabund.name eq 'Na')
Q = (!matchabund.abund)[x]

gam4 = 10^(19d + ((-15.33d * 2d)/3d) + alog10(Pe) - ((5d * alog10(temp))/6d))

y = where(!ioniztn.name eq 'Na')
chi = 0 ;2.1028523  ; eV
C6 = 0.3d-30 * ((1d/(!ioniztn.ion1[y] - chi - chi_lam)^2d) - (1d/(!ioniztn.ion1[y] - chi)^2d))
gam6 = 10^(20d + 0.4 * alog10(C6) + alog10(Pgas) - 0.7d * alog10(temp))

gam_up = A_ul ; doesn't include the 4 pi b/c NIST seems to include 4pi already
;gam_low = 
gam_nat = gam_up; + gam_low

gamma = gam4 + gam6 + gam_nat

lam_D = (lam_core/c) * ((2d * k *temp/m_Na) + (vel * 1d5)^2d)^(0.5d)
a = 2.65d-20 * (lam_core)^2d * gamma/lam_D
;a = (gamma * (lam_core*1d-8)^2d)/(4d * !pi * lam_D * c)
u = (lambda - lam_core)/lam_D
H = voigt(a, u)

;boltz = g_l * exp(-E/(k * temp))/partition('Na', temp)
boltz = g_l /partition('Na', temp)
saha = saha('Na', temp, 1, Pe)
N_NE = (1d/(1d + saha)) * boltz

f = 1.884d-15 * (lam_core^2) * (g_u/g_l) * A_ul /(4d * !pi)

;kappa = ((!pi^0.5d * e^2d)/(m_Na * c)) * (H/lam_D) * (A/abund_sum) *
;((lam_core*1d-8)^2d) * f * (N_NE) * (1d - 10^(-chi_lam*theta))

kappa = 4.995d-21 * (H/lam_D) * (Q/(abund_sum*1.6606d-24)) * ((lam_core)^2d) * f * (N_NE) * (1d - 10^(-chi_lam*theta))

return, kappa

end
