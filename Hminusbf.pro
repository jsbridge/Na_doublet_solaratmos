;+
;
; Joanna Bridge
;
; Calculates opacity for H minus bound free
;
;-
function Hminusbf, lambda, temp, Pe

theta = 5040d/temp
chi_lam = 1.2398d4/lambda

alphabfHminus = 1.99654d0 - 1.18267d-5*lambda + 2.64243d-6*lambda^2 - 4.40524d-10*lambda^3 + 3.23992d-14*lambda^4 - 1.39568d-18*lambda^5 + 2.78701d-23*lambda^6

Hminusbf = (4.158d-10)*Pe*(theta^(5d/2d))*(10^(.754d*theta))*alphabfHminus

; This section is for producing Gray's plots, don't need generally??
;x = min(where(Hminusbf lt 0))
;Hminusbf = [Hminusbf[0:x-1], dblarr(x)]

return, Hminusbf*(1d - 10^(-1d*chi_lam*theta))/(1d + saha('H', temp, 1, Pe))*10^(-18d)/total(!matchabund.abund*!matchabund.weight*1.6606d-24)

end
