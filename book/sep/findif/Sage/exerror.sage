a,k=var('a,k')
f(a,k)=exp(-a*k^2)
tay(a,k)=1 - a*k^2 + (a^2*k^4)/2
ptay=plot(tay(2/3,k)-f(2/3,k),(k,0,pi),legend_label='a=2/3',thickness=2)+plot(tay(4/3,k)-f(4/3,k),(k,0,pi),legend_label='a=4/3',thickness=2,linestyle='--',color='black')
ptay.set_legend_options(title='Explicit Error')
ptay.save(filename='junk_sage.pdf',frame=True,gridlines=True,axes_labels=['wavenumber',''])