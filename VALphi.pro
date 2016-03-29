;+
;
; 
; 
; Phi operator for VALIIIC
;
;-
function VALphi, tau, lambda

nu = 2.99792458d18/lambda ; lambda in angstroms
val = dblarr(n_elements(lambda))
temp = interpol(!VALdata.T, 1d6)


for i = 0, n_elements(nu)-1 do begin
   t1 = interpol(tau[i,*], 1d6)
   t2 = interpol([max(tau[i,*]), max(tau[i,*])+100d], 1d6)
   t = [t1, t2]
   temp_g = !VALdata.T[51]*(0.75*t2+0.5)^(1d/4d)
   temp_arr = [temp, temp_g]
   bb = bnu_func(nu[i], temp_arr)
   int = 2d * bb * expint(2, abs(t), /double)
   val[i] = tsum_d( t, int) 

endfor

return, val

end
