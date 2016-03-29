;+
;
; Joanna Bridge
;
; Calculates opacity for H minus free free
;
;-
function Hminusff, lambda, temp, Pe

theta = 5040d/temp

f0 = -2.2763d0-1.6850d0*alog10(lambda)+0.76661d0*(alog10(lambda))^2-0.053346d0*(alog10(lambda))^3
f1 = 15.2827d0-9.2846d0*alog10(lambda)+1.99281d0*(alog10(lambda))^2-0.142631d0*(alog10(lambda))^3
f2 = -197.789d0+190.266d0*alog10(lambda)-67.9775d0*(alog10(lambda))^2+10.6913d0*(alog10(lambda))^3-0.625151d0*(alog10(lambda))^4
alphaffHminus = (10^(-26d))*10^(f0+f1*alog10(theta)+f2*(alog10(theta))^2)

Hminusff = alphaffHminus*Pe

return, Hminusff/(1d + saha('H', temp, 1, Pe))/total(!matchabund.abund*!matchabund.weight*1.6606d-24)

end
