        integer offset, xpix, ypix, blast
	real    xll, yll, xur, yur, ppi
	character *1 array(25,25)
	integer  ii, jj

	offset = 0
	xpix = 25
	ypix = 25
	xll = 0.
	yll = 0.
	xur = 10.
	yur = 10.
	ppi = 0
	blast = 0

	call vp_init()

	do 100 ii = 1 , 25
	do 100 jj = 1 , 25
	array(jj, ii) = char(7)
100	continue

	do 200 ii = 1 , 25
	array(ii, ii) = char(0)
200	continue

	do 300 ii = 1 , 25
	array(ii, 12) = char(4)
300	continue

	do 400 ii = 1 , 25
	array(12, ii) = char(2)
400	continue

	call vp_raster (array, .false., offset, xpix, ypix, 
     c       xll, yll, xur, yur, 1)

        stop
	end
