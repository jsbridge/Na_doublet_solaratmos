;+
;Joanna Bridge, 1/2012
;
;This program takes in an array of wavenumbers and returns the
;integral of the Planck function for those wavenumbers at T=7500 K in
;ergs/cm^2.  Wavenumbers should be given in microns.
; 
;-

function bnu_int, wave

nu = wave*2.99792458d10/1d-4
bb = bnu_func(nu, 7500.)
integ = int_tabulated(nu, bb, /double, /sort)

return, integ

end
