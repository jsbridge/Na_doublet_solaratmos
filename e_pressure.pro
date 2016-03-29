;+
;
; This function calculates Pe for some temperature and gas pressure
;
;-
function e_pressure, Pgas, temp

;This first chunk of code is because there is not agreement between
;ioniz.txt, parition.txt, and solarabund.txt in the species present,
;so they must be whittled down so all are accounted for.

narr = 'H'

for i = 0, n_elements(!solarabund.name)-1 do begin
   if where((!solarabund.name[i]+'+') eq !partit.name) ne -1 then begin
      narr = [narr,!solarabund.name[i]]
   endif
endfor

solar1 = 'b'
solar2 = 0

for i = 0, n_elements(narr)-1 do begin
   if min(where(!solarabund.name[i] eq narr)) ne -1 then begin
      solar1 = [solar1, !solarabund.name[i]]
      solar2 = [solar2, double(!solarabund.logabund[i])]
   endif
endfor

solar1 = solar1[1:*]
solar2 = solar2[1:*]

new = 0
for k = 0, n_elements(!ioniztn.name)-1 do begin
   if min(where(!ioniztn.name[k] eq solar1)) ne -1 then begin
      new = [new, double(!ioniztn.weight[k])]
   endif
endfor
new = new[1:*]

;light_sum = total(10^solar2[0:29])
;light_sum1 = total((10^solar2)*new[0:29])
;print, light_sum, light_sum1
;metal_sum = total(10^solar2[2:29])
;metal_sum1 = total((10^solar2[2:29])*new[2:29])
;print, metal_sum, metal_sum1

;c = create_struct('name', solar1,'abund', 10^solar2, 'weight', new)
;defsysv, '!matchabund', c

numer = 0
denom = 0
j = 1
Pe = dblarr(10000)

;For hot stars, Pe = 0.5*Pgas is a good guess; for cooler stars, more
;like Pe = sqrt(saha*Pgas)
if temp gt 25200d then return, 0.5d * Pgas else $

if temp lt 6500d then begin
   Pe[0] = 0.5d * Pgas
   Pe[1] = 0.5d * Pgas + 0.02d 
endif

if (temp ge 6500d) and (temp le 25200d) then begin
   Pe[0] = sqrt(Pgas * saha('H', temp, 1, 0.5d * Pgas))
   Pe[1] = Pe[0] + 0.02d
endif

while abs(Pe[j] - Pe[j-1]) gt 0.01d do begin

   for i = 0, n_elements(solar1)-1 do begin
      a = saha(solar1[i], temp, 1, Pe[j])
      num = (10^(solar2[i]))*(a/(1d + a))
      numer += num
      den = (10^(solar2[i]))*(1+(a/(1d + a)))
      denom += den
   endfor
  
   j += 1
   Pe[j] = Pgas*numer/denom
  ; print, Pe[j]
  ; print, j-1 
endwhile

return, Pe[j]

end
