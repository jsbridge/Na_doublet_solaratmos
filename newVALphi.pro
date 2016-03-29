;+
;
; 
; 
; Phi operator for VALIIIC
;
;-
function newVALphi, tau, lambda

readcol, 'source.txt', height, intensity, format='D,D'

nu = 2.99792458d18/lambda ; lambda in angstroms
val = dblarr(n_elements(lambda))
temp = interpol(!VALdata.T, 1d6)


for i = 0, n_elements(nu)-1 do begin
   t1 = interpol(tau[i,*], 1d4)
   t2 = interpol([max(tau[i,*]), max(tau[i,*])+100d], 1d4)
   t = [t1, t2]
   temp_g = !VALdata.T[51]*(0.75*t2+0.5)^(1d/4d)
   temp_arr = [temp, temp_g]
   bb1 = bnu_func(nu[i], temp_arr)
   bb = [interpol(intensity, 1d4)*1d-10, bb1]
   int = 2d * bb * expint(2, abs(t), /double)
   val[i] = tsum_d( t, int) 

endfor

phi = create_struct('val', val, 'B', bb1, 'S', bb)

return, phi

end
