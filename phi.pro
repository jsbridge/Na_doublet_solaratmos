;+
;
; Joanna Bridge, 2/2012
; 
; Phi operator
;
;-
function phi, tau, nu

;t = interpol([1d-6, (max(tau)+10)], 1d4)
val = dblarr(n_elements(tau))
int = dblarr(n_elements(tau))
;temp = 8700d*(0.75*t+0.5)^(1d/4d)

for i = 0, n_elements(tau)-1 do begin
   t1 = interpol([1d-6, tau[i]], 1d4)
   t2 = interpol([tau[i], (tau[i]+10)], 1d4)
   temp1 = 8700d*(0.75*t1+0.5)^(1d/4d)
   temp2 = 8700d*(0.75*t2+0.5)^(1d/4d)
;   b = findel(tau[i], t)
;   if b gt 1 then begin
   if tau[i] gt 0 then begin
      pos = 2d * bnu_func(nu, temp2) * expint(2, ( t2-tau[i] ), /double)
      neg = 2d * bnu_func(nu, temp1) * expint(2, ( tau[i]-t1 ), /double)
      val[i] = int_tabulated( t2, pos, /double) - int_tabulated( t1, neg, /double)
   endif else begin 
      int =  2d * bnu_func(nu, [temp1,temp2[1:*]]) * expint(2, abs([t1,t2[1:*]]-tau[i]), /double)
      val[i] = int_tabulated([t1,t2[1:*]] , int, /double)
   endelse
endfor 

return, val

end
