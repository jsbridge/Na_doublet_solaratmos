;+
;
; Does stuff with Pe
;
;-
pro Pe_height

Pe = Pe_VALIIIC()
h = !VALdata.h
Pe1 = !VALdata.n_e * 1.3806503d-16 * !VALdata.T

set_plot, 'ps'
device, file='Pe_height.ps'
plot, h, Pe, xthick=5, ythick=5, charthick=5,  charsize=1, thick=4, title=textoidl('P_e vs. h'), xtitle='h (km)', ytitle=textoidl('P_e'), linestyle=4, xrange=[-75, 2540], xstyle=1, /ylog, ycharsize=0.85
oplot, h, Pe1, thick=4
legend, ['Ideal Gas Law', 'Gray Eqn. (9.8)'], linestyle=[0,4], thick=4, charthick=5, charsize=1.2, bthick=5
device, /close
set_plot, 'x'

end
