;+
;
; This function reads the VALIIIC.txt file and out puts the data for
; use. It returns an array.
;
;-
function VALdata, file

dat = read_ascii(file, comment_symbol='#')
data = create_struct('h', dat.field01[0,*], 'm', dat.field01[1,*], 'tau_500', dat.field01[2,*], 'T', dat.field01[3,*], 'V', dat.field01[4,*], 'n_H', dat.field01[5,*], 'n_e', dat.field01[6,*], 'Ptotal', dat.field01[7,*], 'Pgas_Ptotal', dat.field01[8,*], 'rho', dat.field01[9,*])

return, data

end
