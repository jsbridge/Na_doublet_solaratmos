;+
;Joanna Bridge, 1/2012
;
;This program calculates and plots the Planck distribution
;vs. wavenumber in microns^-1, log(B_nu) vs. wavenumber, and log(B_nu)
;vs. log(wavenumber).  Each plot is saved to a postscript file.
;
;It calls the function bnu_func.
;
;Three different temperatures are used: 10000 K, 7000 K, and 3000 K
;
;-
pro bnu

wave = dindgen(1200)/100       ;this is in inverse microns
nu = wave*2.99792458d10/1d-4   ;get into cm^-1 and mult. by cm/s
bb1 = bnu_func(nu, 3000.)
bb2 = bnu_func(nu, 7000.)
bb3 = bnu_func(nu, 10000.)

set_plot, 'ps'
device, file='B_v_wave.ps'
plot, wave, bb3, title=textoidl('B_\nu vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('B_\nu (erg cm^{-2})'), ycharsize=.8, thick=6, xthick=6, ythick=6, charthick=3
oplot, wave, bb1, linestyle=3, thick=6
oplot, wave, bb2,linestyle=2, thick=6
legend, ['10,000K', '7,000 K', '3,000 K'], linestyle=[0, 2,3], pos=[8.5, 1.8e-4], thick=6, charthick=3
device, /close

set_plot, 'ps'
device, file='logB_v_wave.ps'
plot, wave, bb3,/ylog, title=textoidl('log(B_\nu) vs. Wavenumber'), xtitle=textoidl('Wavenumber (\mum^{-1})'), ytitle=textoidl('log(B_\nu) (erg cm^{-2})'), thick=6, xthick=6, ythick=6, charthick=3
oplot, wave, bb1, linestyle=3, thick=6
oplot, wave, bb2, linestyle=2, thick=6
legend, ['10,000K', '7,000 K', '3,000 K'], linestyle=[0, 2,3], pos=[8.5, 10^(-3.3)], thick=6, charthick=3
device, /close

set_plot, 'ps'
device, file='logB_v_logwave.ps'
plot,  wave, bb3,/xlog, /ylog,  title=textoidl('log(B_\nu) vs. log(wavenumber)'), xtitle=textoidl('log(wavenumber) (\mum^{-1})'), ytitle=textoidl('log(B_\nu) (erg cm^{-2})'), xrange=[0.1, 15.85], xstyle=1,  thick=6, xthick=6, ythick=6, charthick=3
oplot,  wave, bb1, linestyle=3, thick=6
oplot, wave,  bb2,  linestyle=2, thick=6
legend, ['10,000K', '7,000 K', '3,000 K'], linestyle=[0, 2, 3], pos=[.12,10^(-3.3)], thick=6, charthick=3
device, /close

end
