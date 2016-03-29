;+
;
; Joanna Bridge, 2/2012
;
; This function is a lambda operator.  It takes a source function and
; spits out the mean intensity, J_nu
; Tau out to 60, sampling, 500,000
;-
function lambda, tau, a
; use tau from 0 to 3

t = interpol([0.0001d, 50d], 500000d)
jnu = dblarr(n_elements(tau))

for i = 0, n_elements(tau)-1 do begin
   int = (1d/2d)*snu(tau[i], a)*expint(1, abs(t - tau[i]), /double)
   jnu[i] = int_tabulated( t, int, /double)
endfor 

return, jnu

end
