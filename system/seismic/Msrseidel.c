/* Amplitude balancing of a 2-D amplitude map. */
/*
  Copyright (C) 1994 The Board of Trustees of Stanford University

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
  Use Gauss-Seidel algorithm to approach the problem of amplitude
  balancing of a 2-D amplitude map. The amplitude map should have
  the geological trend removed and should be normalized by the rms
  value. Here we assume an amplitude model of the form:

  A(total) = A(receiver) * A(source) * A(earth) * A(noise)

  where A() is the amplitude.

  Once A(earth) is removed, we can solve the following optimization
  problem:

  min(s,h) || A(total) - A(s) * A(h) * A(y) * A(r) ||^2

  s:  source
  h:  offset
  y:  midpoint
  r:  receiver

  I use the Gauss-Seidel method to get the value of the source-,
  receiver-, and midpoint-consistent balancing coefficients.
*/

# include <rsf.h>
# include <assert.h>


int main(int argc, char* argv[])
{
    float *tabamp                            ;
    float *tabsrc, *taboff, *tabmid, *tabrcv ;
    float *numo, *nums, *numy, *numr         ;
    float *deno, *dens, *deny, *denr         ;
    float *synthso, *synthrv                 ;
    float  o, s, y, r, d                     ;
    float  oh, os, oy, or, dh, ds, dy, dr    ;
    float  sloc, hloc, yloc, rloc            ;
    int    nh, ns, ny, nr, nhs               ;
    int    N = 1, niter                      ;
    int    i, j, k, l                        ;
    sf_file in, off, src, mid, rcv, so, rv;

    sf_init(argc, argv);
    in = sf_input("in");
    off = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint  (in, "n1", &nh)) sf_error("No n1= in input");
    if (!sf_histint  (in, "n2", &ns)) sf_error("No n2= in input");
    if (!sf_histfloat(in, "o1", &oh)) sf_error("No o1= in input");
    if (!sf_histfloat(in, "o2", &os)) sf_error("No o2= in input");
    if (!sf_histfloat(in, "d1", &dh)) sf_error("No d1= in input");
    if (!sf_histfloat(in, "d2", &ds)) sf_error("No d2= in input");

    ny = 2 * (ns - 1) + nh ;
    oy = oh/2              ;
    dy = dh/2              ;

    nr = ns - 1 + nh ;
    or = os + oh     ;
    dr = dh          ;

    if (!sf_getint( "niter", &niter )) niter = 1 ;
    /* number of iterations */

    tabamp = sf_floatalloc( nh * ns);

    taboff = sf_floatalloc( nh);
    tabsrc = sf_floatalloc( ns);
    tabmid = sf_floatalloc( ny);
    tabrcv = sf_floatalloc( nr);

    sf_putint(off,"n1", nh); 
    sf_putint(off,"n2", N );
    sf_putfloat(off,"o1", oh);
    sf_putfloat(off,"d1", dh);
    sf_putstring(off,"label1", "Offset"   );
    sf_putstring(off,"label2", "Amplitude");

    src = sf_output("src"); /* source-consistent coefficients */
    sf_putint(src,"n1", ns);
    sf_putint(src,"n2", N );
    sf_putfloat(src,"o1", os);
    sf_putfloat(src,"d1", ds);
    sf_putstring(src,"label1", "Shot"     );
    sf_putstring(src,"label2", "Amplitude");

    mid = sf_output("mid"); /* midpoint-consistent coefficients */
    sf_putint(mid,"n1", ny);
    sf_putint(mid,"n2", N );
    sf_putfloat(mid,"o1", oy);
    sf_putfloat(mid,"d1", dy);
    sf_putstring(mid,"label1", "Midpoint" );
    sf_putstring(mid,"label2", "Amplitude");

    rcv = sf_output("rcv"); /* receiver-consistent coefficients */
    sf_putint(rcv,"n1", nr);
    sf_putint(rcv,"n2", N );
    sf_putfloat(rcv,"o1", or);
    sf_putfloat(rcv,"d1", dr);
    sf_putstring(rcv,"label1", "Receiver" );
    sf_putstring(rcv,"label2", "Amplitude");

    so = sf_output("so");
    rv = sf_output("rv");

    nhs = nh * ns;

    sf_floatread(tabamp, nhs,in);

    numo = sf_floatalloc( nh);
    deno = sf_floatalloc( nh);

    nums = sf_floatalloc( ns);
    dens = sf_floatalloc( ns);

    numy = sf_floatalloc( ny);
    deny = sf_floatalloc( ny);

    numr = sf_floatalloc( nr);
    denr = sf_floatalloc( nr);

    synthso = sf_floatalloc(nhs);
    synthrv = sf_floatalloc(nhs);

    for ( i = 0 ; i < nh ; i ++ )
    {
	numo[i] = 0.0 ; deno[i] = 0.0 ; taboff[i] = 1.0 ;
    }

    for ( j = 0 ; j < ns ; j ++ )
    {
	nums[j] = 0.0 ; dens[j] = 0.0 ; tabsrc[j] = 1.0 ;
    }

    for ( k = 0 ; k < ny ; k ++ )
    {
	numy[k] = 0.0 ; deny[k] = 0.0 ; tabmid[k] = 1.0 ;
    }

    for ( l = 0 ; l < nr ; l ++ )
    {
	numr[l] = 0.0 ; denr[l] = 0.0 ; tabrcv[l] = 1.0 ;
    }

    while ( niter != 0 ) 
    {
	for ( i = 0 ; i < nh ; i ++ )
        {
	    for ( j = 0 ; j < ns ; j ++ )
	    {
		hloc = oh + i * dh   ;
		sloc = os + j * ds   ;
		yloc = sloc + hloc/2 ;
		rloc = sloc + hloc   ;

		k = (int) floorf( (yloc - oy)/dy + 0.5 );
		l = (int) floorf( (rloc - or)/dr + 0.5 );

		assert( k >= 0 && k < ny );
		assert( l >= 0 && l < nr );

		d = tabamp[nh * j + i] ;
		o = taboff[i]          ;
		s = tabsrc[j]          ;
		y = tabmid[k]          ;
		r = tabrcv[l]          ;

		numo[i] += d * s * y * r ;
		nums[j] += d * o * y * r ;
		numy[k] += d * o * s * r ;
		numr[l] += d * o * s * y ;

		deno[i] += s * y * r * s * y * r ;
		dens[j] += o * y * r * o * y * r ;
		deny[k] += o * s * r * o * s * r ;
		denr[l] += o * s * y * o * s * y ;
	    }
        }

	for ( i = 0 ; i < nh ; i ++ ) taboff[i] = numo[i]/deno[i] ;
	for ( j = 0 ; j < ns ; j ++ ) tabsrc[j] = nums[j]/dens[j] ;
	for ( k = 0 ; k < ny ; k ++ ) tabmid[k] = numy[k]/deny[k] ;
	for ( l = 0 ; l < nr ; l ++ ) tabrcv[l] = numr[l]/denr[l] ;

	niter -- ;
    }

    sf_floatwrite(taboff, nh,off);
    sf_floatwrite(tabsrc, ns,src);
    sf_floatwrite(tabmid, ny,mid);
    sf_floatwrite(tabrcv, nr,rcv);

    for ( i = 0 ; i < nh ; i ++ )
    {
	for ( j = 0 ; j < ns ; j ++ )
	{
	    hloc = oh + i * dh;
	    sloc = os + j * ds;
	    rloc = sloc + hloc;

	    l = (int) floorf( (rloc - or)/dr + 0.5 );

	    assert( l >= 0 && l < nr );

	    synthso[nh * j + i] = taboff[i] * tabsrc[j];
	    synthrv[nh * j + i] = tabrcv[l]            ;
	}
    }

    sf_floatwrite(synthso, nhs, so) ;
    sf_floatwrite(synthrv, nhs, rv) ;

    exit(0);
}
