;+
;
; Joanna Bridge, 2/2012
;
; This function calculates the log of the Saha equation for any given T
; and species.  If you're talking about H- to HI, input 'H-' for symb.
; Note: this returns PHI/Pe
;
;-
function saha, symb, temp, state, Pe

;readcol, 'ioniz.txt', num, name, weight, ion1, ion2, ion3, format='I,A,D,D,D,D', /silent
theta = 5040d/temp

if symb eq 'H-' then begin
   logN1_N2Pe = -0.1762d + (2.5d*alog10(temp))-(theta*0.755d)+alog10(partition('H', temp)/partition('H-', temp))
endif else begin
   a = where(symb eq !ioniztn.name)

   if state eq 1 then begin     ; first ionization
      U = partition(symb, temp)
      U1 = partition(symb+'+', temp)
      logN1_N2Pe = -0.1762d + (2.5d*alog10(temp))-(theta*!ioniztn.ion1[a])+alog10(U1/U)
   endif

   if state eq 2 then begin     ; second ionization
      U1 = partition(symb+'+', temp)
      U2 = partition(symb+'++', temp)
      logN1_N2Pe = -0.1762d + (2.5d*alog10(temp))-(theta*!ioniztn.ion2[a])+alog10(U2/U1)
   endif
endelse

saha = 10^(logN1_N2Pe)/Pe

return, saha

end
