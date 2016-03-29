;+
;
; Joanna Bridge, 2/2012
;
; This program calculates the partition function for any species at an
; arbitrary temparture.
;
;-
function partition, symb, temp

if symb eq 'HII' then return, 1d else $
if symb eq 'H-' then return, 1d else $
if symb eq 'H+' then return, 1d else $

;readcol, 'names.txt', name, format='A', /silent
;readcol, 'partition.txt', th0, th1, th2, th3, th4, th5, th6, th7, th8, th9, log_g0, format='D,D,D,D,D,D,D,D,D,D,D', /silent


th0 = double(10^!partit.th0)
th1 = double(10^!partit.th1)
th2 = double(10^!partit.th2)
th3 = double(10^!partit.th3)
th4 = double(10^!partit.th4)
th5 = double(10^!partit.th5)
th6 = double(10^!partit.th6)
th7 = double(10^!partit.th7)
th8 = double(10^!partit.th8)
th9 = double(10^!partit.th9)

x = where(symb eq !partit.name)   ; index of symbol - x-th row
theta = 5040d/temp

if (theta ge 0.2) and (theta le 2.0) then begin

   title = double([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
   a = max(where(theta gt title)) ; interpolate btwn a and a+1
 ;  if a lt 9 then begin
   num = interpol([title[a], title[a+1]], 3000d)
   b = findel(num, theta)
  ; endif

   case a of
      0: if th0[x] eq -1 then print, 'Outisde the range of this function; choose a lower temperature.' else return, (interpol([th0[x], th1[x]], 3000d))[b]
      1: return, (interpol([th1[x], th2[x]], 3000d))[b]
      2: return, (interpol([th2[x], th3[x]], 3000d))[b]
      3: return, (interpol([th3[x], th4[x]], 3000d))[b]
      4: return, (interpol([th4[x], th5[x]], 3000d))[b]
      5: return, (interpol([th5[x], th6[x]], 3000d))[b]
      6: return, (interpol([th6[x], th7[x]], 3000d))[b]
      7: return, (interpol([th7[x], th8[x]], 3000d))[b]
      8: return, (interpol([th8[x], th9[x]], 3000d))[b]
    ;  9: print, 'Outisde the range of this function; choose a lower temperature.'
   endcase

endif

if theta gt 2.0 then return, 10^(!partit.log_g0[x])

if theta lt 0.2 then print, 'Outisde the range of this function; choose a lower temperature.'

end
