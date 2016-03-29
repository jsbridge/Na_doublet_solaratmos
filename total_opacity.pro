;+
;
; Calculates total opacity (plus electron scattering) given
; wavelength, gas pressure, and Pe. Returns in cm^2/g
; 
; 
;
;-
function total_opacity, lambda, temp, Pe, Pgas

elec = 0.6648d-24*total(!matchabund.abund)*(Pe)/((Pgas - Pe)*total(!matchabund.abund*!matchabund.weight*1.6606d-24))

Hbf = neutralHbf(lambda, temp, Pe)
Hff = neutralHff(lambda, temp, Pe)
Hminbf = Hminusbf(lambda, temp, Pe)
Hminff = Hminusff(lambda, temp, Pe)

opacity = Hbf + Hff + Hminbf + Hminff + elec

print, 'Hbf', Hbf
print, 'Hff', Hff
print, 'Hminbf', Hminbf
print, 'Hminff', Hminff
print, 'e', elec

return, opacity

end
