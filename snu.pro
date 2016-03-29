;+
;
; Joanna Bridge 2/2012
;
; This function will generate a source function given a set of taus
; and a's. If tau has more than one element, than the array
; returned will have each column be for a new tau value.
;
;-
function snu, tau, a

S = 0

for i = 0, n_elements(a)-1 do begin
   S += a[i] * tau^i
endfor
   
return, S

end
