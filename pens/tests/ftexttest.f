	real xp, yp, xpath, ypath, xup, yup
	integer i
	character*10 filename
	data NP/18/

	call vp_init()
c
c  choose standard style for vplot
c
	call vp_style(0)
	call vp_scale(.7,.7)
	call vp_fat(2)
c
	call vp_message('This is a message!')
c
	do i=1,NP-1
		xp = 1.8 + .9 * (i-1)
		yp = .4 + .5 * (i-1)
		xpath = .2 + (7. - i) / 9.
		ypath = .2
		xup = -.15
		yup = .35
		call vp_color(mod(i,6)+2)
		call vp_tfont((i-1),2,3)
		call vp_ugtext
     1         (xp,yp,xpath,ypath,xup,yup,'> SEP Geophysics <')
	enddo
c
c  end of session
c
	stop
	end
