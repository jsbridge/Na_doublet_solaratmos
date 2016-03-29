;+
;Joanna Bridge, 1/2012
;
;This program calculates the Planck distribution given a temperature
;
;-
function bnu_func, nu, temp

bnu = (2d*6.62606876d-27*nu^3/2.99792458d10^2)*(1d/(exp((6.62606876d-27*nu)/(1.3806503d-16*temp))-1d))

return, bnu

end
