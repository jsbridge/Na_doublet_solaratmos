;+
;
; This program uses e_pressure to calculate the Pe at the range of
; VALIIIC temperatures and Pgas.
;
;-
function Pe_VALIIIC

Pgas = double(!VALdata.Pgas_Ptotal * !VALdata.Ptotal)
Pe = dblarr(n_elements(!VALdata.T))
abund = dblarr(n_elements(!VALdata.T))

for i = 0, n_elements(!VALdata.T)-1 do begin
   Pe[i] = e_pressure(Pgas[i], double(!VALdata.T[i]))
endfor

return, Pe

end
