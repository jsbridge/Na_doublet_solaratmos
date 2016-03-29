;+
;
; Joanna Bridge, 2/2012
;
; This function is a lambda operator.  It takes a source function and
; and tau array, and assumes that tau starts at 0.
; 
;-
function lambda, tau, nu

t = interpol([0d, 11d], 10000d)+1d-6
val = dblarr(n_elements(tau))
temp = 8700d*(0.75*t+0.5)^(1d/4d)

for i = 0, n_elements(tau)-1 do begin
   int = (1d/2d)* bnu_func(nu, temp)* expint(1, abs(t - tau[i]), /double)
   val[i] = int_tabulated( t, int, /double)
endfor 

return, val

end
;[(func[0]*expint(2, t[1], /double)),
