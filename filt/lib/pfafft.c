/* Copyright (c) Colorado School of Mines, 2000.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
PFAFFT - Functions to perform Prime Factor (PFA) FFT's, in place

npfa		return valid n for complex-to-complex PFA
npfar		return valid n for real-to-complex/complex-to-real PFA
npfao		return optimal n for complex-to-complex PFA
npfaro		return optimal n for real-to-complex/complex-to-real PFA
pfacc		1D PFA complex to complex
pfacr		1D PFA complex to real
pfarc		1D PFA real to complex
pfamcc		multiple PFA complex to real
pfa2cc		2D PFA complex to complex
pfa2cr		2D PFA complex to real
pfa2rc		2D PFA real to complex

*****************************************************************************
Function Prototypes:
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, complex z[]);
void pfacr (int isign, int n, complex cz[], float rz[]);
void pfarc (int isign, int n, float rz[], complex cz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);
void pfa2cc (int isign, int idim, int n1, int n2, complex z[]);
void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], float rz[]);
void pfa2rc (int isign, int idim, int n1, int n2, float rz[], complex cz[]);

*****************************************************************************
npfa:
Input:
nmin		lower bound on returned value (see notes below)

Returned:	valid n for prime factor fft

******************************************************************************
npfao
Input:
nmin		lower bound on returned value (see notes below)
nmax		desired (but not guaranteed) upper bound on returned value

Returned:	valid n for prime factor fft

******************************************************************************
npfar
Input:
nmin		lower bound on returned value

Returned:	valid n for real-to-complex/complex-to-real prime factor fft

*****************************************************************************
npfaro:
Input:
nmin		lower bound on returned value
nmax		desired (but not guaranteed) upper bound on returned value

Returned:	valid n for real-to-complex/complex-to-real prime factor fft

******************************************************************************
pfacc:
Input:
isign		sign of isign is the sign of exponent in fourier kernel
n		length of transform (see notes below)
z		array[n] of complex numbers to be transformed in place

Output:
z		array[n] of complex numbers transformed

******************************************************************************
pfacr:
Input:
isign       sign of isign is the sign of exponent in fourier kernel
n           length of transform (see notes below)
cz          array[n/2+1] of complex values (may be equivalenced to rz)

Output:
rz          array[n] of real values (may be equivalenced to cz)

******************************************************************************
pfarc:
Input:
isign       sign of isign is the sign of exponent in fourier kernel
n           length of transform; must be even (see notes below)
rz          array[n] of real values (may be equivalenced to cz)

Output:
cz          array[n/2+1] of complex values (may be equivalenced to rz)

******************************************************************************
pfamcc:
Input:
isign       	sign of isign is the sign of exponent in fourier kernel
n           	number of complex elements per transform (see notes below)
nt          	number of transforms
k           	stride in complex elements within transforms
kt          	stride in complex elements between transforms
z           	array of complex elements to be transformed in place

Output:
z		array of complex elements transformed

******************************************************************************
pfa2cc:
Input:
isign       	sign of isign is the sign of exponent in fourier kernel
idim        	dimension to transform, either 1 or 2 (see notes)
n1          	1st (fast) dimension of array to be transformed (see notes)
n2          	2nd (slow) dimension of array to be transformed (see notes)
z           	array[n2][n1] of complex elements to be transformed in place

Output:
z		array[n2][n1] of complex elements transformed

******************************************************************************
pfa2cr:
Input:
isign       sign of isign is the sign of exponent in fourier kernel
idim        dimension to transform, which must be either 1 or 2 (see notes)
n1          1st (fast) dimension of array to be transformed (see notes)
n2          2nd (slow) dimension of array to be transformed (see notes)
cz          array of complex values (may be equivalenced to rz)

Output:
rz          array of real values (may be equivalenced to cz)

******************************************************************************
pfa2rc:
Input:
isign       sign of isign is the sign of exponent in fourier kernel
idim        dimension to transform, which must be either 1 or 2 (see notes)
n1          1st (fast) dimension of array to be transformed (see notes)
n2          2nd (slow) dimension of array to be transformed (see notes)
rz          array of real values (may be equivalenced to cz)

Output:
cz          array of complex values (may be equivalenced to rz)

******************************************************************************
Notes:
Table of valid n and cost for prime factor fft.  For each n, cost
was estimated to be the inverse of the number of ffts done in 1 sec
on an IBM RISC System/6000 Model 320H, by Dave Hale, 08/04/91.
(Redone by Jack Cohen for 15 sec to rebuild NTAB table on advice of
David and Gregory Chudnovsky, 05/03/94).
Cost estimates are least accurate for very small n.  An alternative method
for estimating cost would be to count multiplies and adds, but this
method fails to account for the overlapping of multiplies and adds
that is possible on some computers, such as the IBM RS/6000 family.

npfa:
The returned n will be composed of mutually prime factors from
the set {2,3,4,5,7,8,9,11,13,16}.  Because n cannot exceed
720720 = 5*7*9*11*13*16, 720720 is returned if nmin exceeds 720720.

npfao:
The returned n will be composed of mutually prime factors from
the set {2,3,4,5,7,8,9,11,13,16}.  Because n cannot exceed
720720 = 5*7*9*11*13*16, 720720 is returned if nmin exceeds 720720.
If nmin does not exceed 720720, then the returned n will not be 
less than nmin.  The optimal n is chosen to minimize the estimated
cost of performing the fft, while satisfying the constraint, if
possible, that n not exceed nmax.

npfar and npfaro:
Current implemenations of real-to-complex and complex-to-real prime 
factor ffts require that the transform length n be even and that n/2 
be a valid length for a complex-to-complex prime factor fft.  The 
value returned by npfar satisfies these conditions.  Also, see notes 
for npfa.

pfacc:
n must be factorable into mutually prime factors taken 
from the set {2,3,4,5,7,8,9,11,13,16}.  in other words,
	n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u
where
	0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1
is required for pfa to yield meaningful results.  this
restriction implies that n is restricted to the range
	1 <= n <= 720720 (= 5*7*9*11*13*16)

pfacr:
Because pfacr uses pfacc to do most of the work, n must be even 
and n/2 must be a valid length for pfacc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.

pfarc:
Because pfarc uses pfacc to do most of the work, n must be even 
and n/2 must be a valid length for pfacc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.

pfamcc:
To perform a two-dimensional transform of an n1 by n2 complex array 
(assuming that both n1 and n2 are valid "n"), stored with n1 fast 
and n2 slow:
    pfamcc(isign,n1,n2,1,n1,z); (to transform 1st dimension)
    pfamcc(isign,n2,n1,n1,1,z); (to transform 2nd dimension)

pfa2cc:
Only one (either the 1st or 2nd) dimension of the 2-D array is transformed.

If idim equals 1, then n2 transforms of n1 complex elements are performed; 
else, if idim equals 2, then n1 transforms of n2 complex elements are 
performed.

Although z appears in the argument list as a one-dimensional array,
z may be viewed as an n1 by n2 two-dimensional array:  z[n2][n1].

Valid n is computed via the "np" subroutines.

To perform a two-dimensional transform of an n1 by n2 complex array 
(assuming that both n1 and n2 are valid "n"), stored with n1 fast 
and n2 slow:  pfa2cc(isign,1,n1,n2,z);  pfa2cc(isign,2,n1,n2,z);

pfa2cr:
If idim equals 1, then n2 transforms of n1/2+1 complex elements to n1 real 
elements are performed; else, if idim equals 2, then n1 transforms of n2/2+1 
complex elements to n2 real elements are performed.

Although rz appears in the argument list as a one-dimensional array,
rz may be viewed as an n1 by n2 two-dimensional array:  rz[n2][n1].  
Likewise, depending on idim, cz may be viewed as either an n1/2+1 by 
n2 or an n1 by n2/2+1 two-dimensional array of complex elements.

Let n denote the transform length, either n1 or n2, depending on idim.
Because pfa2rc uses pfa2cc to do most of the work, n must be even 
and n/2 must be a valid length for pfa2cc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.

pfa2rc:
If idim equals 1, then n2 transforms of n1 real elements to n1/2+1 complex 
elements are performed; else, if idim equals 2, then n1 transforms of n2 
real elements to n2/2+1 complex elements are performed.

Although rz appears in the argument list as a one-dimensional array,
rz may be viewed as an n1 by n2 two-dimensional array:  rz[n2][n1].  
Likewise, depending on idim, cz may be viewed as either an n1/2+1 by 
n2 or an n1 by n2/2+1 two-dimensional array of complex elements.

Let n denote the transform length, either n1 or n2, depending on idim.
Because pfa2rc uses pfa2cc to do most of the work, n must be even 
and n/2 must be a valid length for pfa2cc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.

******************************************************************************
References:  
Temperton, C., 1985, Implementation of a self-sorting
in-place prime factor fft algorithm:  Journal of
Computational Physics, v. 58, p. 283-299.

Temperton, C., 1988, A new set of minimum-add rotated
rotated dft modules: Journal of Computational Physics,
v. 75, p. 190-198.

Press et al, 1988, Numerical Recipes in C, p. 417.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/27/89
*****************************************************************************/
/**************** end self doc ********************************/
#include <math.h>
#include <stdlib.h>

#include "pfafft.h"

#define NTAB 240
static struct {
int n;  float c;
} nctab[NTAB] = {
{       1, 0.000052 },
{       2, 0.000061 },
{       3, 0.000030 },
{       4, 0.000053 },
{       5, 0.000066 },
{       6, 0.000067 },
{       7, 0.000071 },
{       8, 0.000062 },
{       9, 0.000079 },
{      10, 0.000080 },
{      11, 0.000052 },
{      12, 0.000069 },
{      13, 0.000103 },
{      14, 0.000123 },
{      15, 0.000050 },
{      16, 0.000086 },
{      18, 0.000108 },
{      20, 0.000101 },
{      21, 0.000098 },
{      22, 0.000135 },
{      24, 0.000090 },
{      26, 0.000165 },
{      28, 0.000084 },
{      30, 0.000132 },
{      33, 0.000158 },
{      35, 0.000138 },
{      36, 0.000147 },
{      39, 0.000207 },
{      40, 0.000156 },
{      42, 0.000158 },
{      44, 0.000176 },
{      45, 0.000171 },
{      48, 0.000185 },
{      52, 0.000227 },
{      55, 0.000242 },
{      56, 0.000194 },
{      60, 0.000215 },
{      63, 0.000233 },
{      65, 0.000288 },
{      66, 0.000271 },
{      70, 0.000248 },
{      72, 0.000247 },
{      77, 0.000285 },
{      78, 0.000395 },
{      80, 0.000285 },
{      84, 0.000209 },
{      88, 0.000332 },
{      90, 0.000321 },
{      91, 0.000372 },
{      99, 0.000400 },
{     104, 0.000391 },
{     105, 0.000358 },
{     110, 0.000440 },
{     112, 0.000367 },
{     117, 0.000494 },
{     120, 0.000413 },
{     126, 0.000424 },
{     130, 0.000549 },
{     132, 0.000480 },
{     140, 0.000450 },
{     143, 0.000637 },
{     144, 0.000497 },
{     154, 0.000590 },
{     156, 0.000626 },
{     165, 0.000654 },
{     168, 0.000536 },
{     176, 0.000656 },
{     180, 0.000611 },
{     182, 0.000730 },
{     195, 0.000839 },
{     198, 0.000786 },
{     208, 0.000835 },
{     210, 0.000751 },
{     220, 0.000826 },
{     231, 0.000926 },
{     234, 0.000991 },
{     240, 0.000852 },
{     252, 0.000820 },
{     260, 0.001053 },
{     264, 0.000987 },
{     273, 0.001152 },
{     280, 0.000952 },
{     286, 0.001299 },
{     308, 0.001155 },
{     312, 0.001270 },
{     315, 0.001156 },
{     330, 0.001397 },
{     336, 0.001173 },
{     360, 0.001259 },
{     364, 0.001471 },
{     385, 0.001569 },
{     390, 0.001767 },
{     396, 0.001552 },
{     420, 0.001516 },
{     429, 0.002015 },
{     440, 0.001748 },
{     455, 0.001988 },
{     462, 0.001921 },
{     468, 0.001956 },
{     495, 0.002106 },
{     504, 0.001769 },
{     520, 0.002196 },
{     528, 0.002127 },
{     546, 0.002454 },
{     560, 0.002099 },
{     572, 0.002632 },
{     585, 0.002665 },
{     616, 0.002397 },
{     624, 0.002711 },
{     630, 0.002496 },
{     660, 0.002812 },
{     693, 0.002949 },
{     715, 0.003571 },
{     720, 0.002783 },
{     728, 0.003060 },
{     770, 0.003392 },
{     780, 0.003553 },
{     792, 0.003198 },
{     819, 0.003726 },
{     840, 0.003234 },
{     858, 0.004354 },
{     880, 0.003800 },
{     910, 0.004304 },
{     924, 0.003975 },
{     936, 0.004123 },
{     990, 0.004517 },
{    1001, 0.005066 },
{    1008, 0.003902 },
{    1040, 0.004785 },
{    1092, 0.005017 },
{    1144, 0.005599 },
{    1155, 0.005380 },
{    1170, 0.005730 },
{    1232, 0.005323 },
{    1260, 0.005112 },
{    1287, 0.006658 },
{    1320, 0.005974 },
{    1365, 0.006781 },
{    1386, 0.006413 },
{    1430, 0.007622 },
{    1456, 0.006679 },
{    1540, 0.007032 },
{    1560, 0.007538 },
{    1584, 0.007126 },
{    1638, 0.007979 },
{    1680, 0.007225 },
{    1716, 0.008961 },
{    1820, 0.008818 },
{    1848, 0.008427 },
{    1872, 0.009004 },
{    1980, 0.009398 },
{    2002, 0.010830 },
{    2145, 0.012010 },
{    2184, 0.010586 },
{    2288, 0.012058 },
{    2310, 0.011673 },
{    2340, 0.011700 },
{    2520, 0.011062 },
{    2574, 0.014313 },
{    2640, 0.013021 },
{    2730, 0.014606 },
{    2772, 0.013216 },
{    2860, 0.015789 },
{    3003, 0.016988 },
{    3080, 0.014911 },
{    3120, 0.016393 },
{    3276, 0.016741 },
{    3432, 0.018821 },
{    3465, 0.018138 },
{    3640, 0.018892 },
{    3696, 0.018634 },
{    3960, 0.020216 },
{    4004, 0.022455 },
{    4095, 0.022523 },
{    4290, 0.026087 },
{    4368, 0.023474 },
{    4620, 0.024590 },
{    4680, 0.025641 },
{    5005, 0.030303 },
{    5040, 0.025253 },
{    5148, 0.030364 },
{    5460, 0.031250 },
{    5544, 0.029412 },
{    5720, 0.034404 },
{    6006, 0.037500 },
{    6160, 0.034091 },
{    6435, 0.040214 },
{    6552, 0.037221 },
{    6864, 0.042735 },
{    6930, 0.040214 },
{    7280, 0.042980 },
{    7920, 0.045872 },
{    8008, 0.049505 },
{    8190, 0.049834 },
{    8580, 0.055762 },
{    9009, 0.057034 },
{    9240, 0.054945 },
{    9360, 0.056818 },
{   10010, 0.066667 },
{   10296, 0.065502 },
{   10920, 0.068182 },
{   11088, 0.065217 },
{   11440, 0.075000 },
{   12012, 0.078534 },
{   12870, 0.087719 },
{   13104, 0.081081 },
{   13860, 0.084270 },
{   15015, 0.102740 },
{   16016, 0.106383 },
{   16380, 0.105634 },
{   17160, 0.119048 },
{   18018, 0.123967 },
{   18480, 0.119048 },
{   20020, 0.137615 },
{   20592, 0.140187 },
{   21840, 0.154639 },
{   24024, 0.168539 },
{   25740, 0.180723 },
{   27720, 0.180723 },
{   30030, 0.220588 },
{   32760, 0.241935 },
{   34320, 0.254237 },
{   36036, 0.254237 },
{   40040, 0.288462 },
{   45045, 0.357143 },
{   48048, 0.357143 },
{   51480, 0.384615 },
{   55440, 0.384615 },
{   60060, 0.454545 },
{   65520, 0.517241 },
{   72072, 0.576923 },
{   80080, 0.625000 },
{   90090, 0.833333 },
{  102960, 0.789474 },
{  120120, 1.153846 },
{  144144, 1.153846 },
{  180180, 1.875000 },
{  240240, 2.500000 },
{  360360, 3.750000 },
{  720720, 7.500000 },
};

int sf_npfa (int nmin)
/*****************************************************************************
Return smallest valid n not less than nmin for prime factor fft.
******************************************************************************
Input:
nmin		lower bound on returned value (see notes below)

Returned:	valid n for prime factor fft
******************************************************************************
Notes:
The returned n will be composed of mutually prime factors from
the set {2,3,4,5,7,8,9,11,13,16}.  Because n cannot exceed
720720 = 5*7*9*11*13*16, 720720 is returned if nmin exceeds 720720.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/28/89
Modified:  Dave Hale, Colorado School of Mines, 08/05/91
	For efficiency, use pre-computed table of valid n and costs.
*****************************************************************************/
{
	int i;
	for (i=0; i<NTAB-1 && nctab[i].n<nmin; ++i);
	return nctab[i].n;
}

int sf_npfao (int nmin, int nmax)
/*****************************************************************************
Return optimal n between nmin and nmax for prime factor fft.
******************************************************************************
Input:
nmin		lower bound on returned value (see notes below)
nmax		desired (but not guaranteed) upper bound on returned value

Returned:	valid n for prime factor fft
******************************************************************************
Notes:
The returned n will be composed of mutually prime factors from
the set {2,3,4,5,7,8,9,11,13,16}.  Because n cannot exceed
720720 = 5*7*9*11*13*16, 720720 is returned if nmin exceeds 720720.
If nmin does not exceed 720720, then the returned n will not be 
less than nmin.  The optimal n is chosen to minimize the estimated
cost of performing the fft, while satisfying the constraint, if
possible, that n not exceed nmax.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/13/89
Modified:  Dave Hale, Colorado School of Mines, 08/05/91
	For efficiency, use pre-computed table of valid n and costs.
*****************************************************************************/
{
	int i,j;
	for (i=0; i<NTAB-1 && nctab[i].n<nmin; ++i);
	for (j=i+1; j<NTAB-1 && nctab[j].n<=nmax; ++j)
		if (nctab[j].c<nctab[i].c) i = j;
	return nctab[i].n;
}

int sf_npfar (int nmin)
/*****************************************************************************
Return smallest valid n not less than nmin for real-to-complex or 
complex-to-real prime factor ffts.
******************************************************************************
Input:
nmin		lower bound on returned value

Returned:	valid n for real-to-complex/complex-to-real prime factor fft
******************************************************************************
Notes:
Current implemenations of real-to-complex and complex-to-real prime 
factor ffts require that the transform length n be even and that n/2 
be a valid length for a complex-to-complex prime factor fft.  The 
value returned by npfar satisfies these conditions.  Also, see notes 
for npfa.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/16/89
*****************************************************************************/
{
    return 2*sf_npfa((nmin+1)/2);
}

int sf_npfaro (int nmin, int nmax)
/*****************************************************************************
Return optimal n between nmin and nmax for real-to-complex or 
complex-to-real prime factor ffts
******************************************************************************
Input:
nmin		lower bound on returned value
nmax		desired (but not guaranteed) upper bound on returned value

Returned:	valid n for real-to-complex/complex-to-real prime factor fft
******************************************************************************
Notes:
Current implemenations of real-to-complex and complex-to-real prime 
factor ffts require that the transform length n be even and that n/2 
be a valid length for a complex-to-complex prime factor fft.  The 
value returned by npfaro satisfies these conditions.  Also, see notes 
for npfao.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/16/89
*****************************************************************************/
{
    return 2*sf_npfao((nmin+1)/2,(nmax+1)/2);
}

#define P120 0.120536680
#define P142 0.142314838
#define P173 0.173648178
#define P222 0.222520934
#define P239 0.239315664
#define P281 0.281732557
#define P342 0.342020143
#define P354 0.354604887
#define P382 0.382683432
#define P415 0.415415013
#define P433 0.433883739
#define P464 0.464723172
#define P540 0.540640817
#define P559 0.559016994
#define P568 0.568064747
#define P587 0.587785252
#define P623 0.623489802
#define P642 0.642787610
#define P654 0.654860734
#define P663 0.663122658
#define P707 0.707106781
#define P748 0.748510748
#define P755 0.755749574
#define P766 0.766044443
#define P781 0.781831482
#define P822 0.822983866
#define P841 0.841253533
#define P866 0.866025404
#define P885 0.885456026
#define P900 0.900968868
#define P909 0.909631995
#define P923 0.923879533
#define P935 0.935016243
#define P939 0.939692621
#define P951 0.951056516
#define P959 0.959492974
#define P970 0.970941817
#define P974 0.974927912
#define P984 0.984807753
#define P989 0.989821442
#define P992 0.992708874
#define NFAX 10

void sf_pfacc (int isign, int n, float complex cz[])
/*****************************************************************************
Prime factor fft:  complex to complex transform, in place
******************************************************************************
Input:
isign		sign of isign is the sign of exponent in fourier kernel
n		length of transform (see notes below)
z		array[n] of complex numbers to be transformed in place

Output:
z		array[n] of complex numbers transformed
******************************************************************************
Notes:
n must be factorable into mutually prime factors taken 
from the set {2,3,4,5,7,8,9,11,13,16}.  in other words,
	n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u
where
	0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1
is required for pfa to yield meaningful results.  this
restriction implies that n is restricted to the range
	1 <= n <= 720720 (= 5*7*9*11*13*16)
******************************************************************************
References:  
Temperton, C., 1985, Implementation of a self-sorting
in-place prime factor fft algorithm:  Journal of
Computational Physics, v. 58, p. 283-299.

Temperton, C., 1988, A new set of minimum-add rotated
rotated dft modules: Journal of Computational Physics,
v. 75, p. 190-198.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/27/89
*****************************************************************************/
{
	static int kfax[] = { 16,13,11,9,8,7,5,4,3,2 };
	register float *z=(float*)cz;
	register int j00,j01,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,jt;
	int nleft,jfax,ifac,jfac,jinc,jmax,ndiv,m,mm=0,mu=0,l;
	float t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,
		t6r,t6i,t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,
		t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,
		t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t20r,t20i,
		t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,
		t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t30r,t30i,
		t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,
		t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t40r,t40i,
		t41r,t41i,t42r,t42i,
		y1r,y1i,y2r,y2i,y3r,y3i,y4r,y4i,y5r,y5i,
		y6r,y6i,y7r,y7i,y8r,y8i,y9r,y9i,y10r,y10i,
		y11r,y11i,y12r,y12i,y13r,y13i,y14r,y14i,y15r,y15i,
		c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

	/* keep track of n left after dividing by factors */
	nleft = n;

	/* begin loop over possible factors (from biggest to smallest) */
	for (jfax=0; jfax<NFAX; jfax++) {

		/* skip if not a mutually prime factor of n */
        ifac = kfax[jfax];
        ndiv = nleft/ifac;
        if (ndiv*ifac!=nleft) continue;
 
		/* update n left and determine n divided by factor */
        nleft = ndiv;
        m = n/ifac;
 
		/* determine rotation factor mu and stride mm */
        for (jfac=1; jfac<=ifac; jfac++) {
			mu = jfac;
			mm = jfac*m;
			if (mm%ifac==1) break;
		}
 
		/* adjust rotation factor for sign of transform */
        if (isign<0) mu = ifac-mu;
 
		/* compute stride, limit, and pointers */
        jinc = 2*mm;
		jmax = 2*n;
        j00 = 0;
        j01 = j00+jinc;

		/* if factor is 2 */
        if (ifac==2) {
			for (l=0; l<m; l++) {
				t1r = z[j00]-z[j01];
				t1i = z[j00+1]-z[j01+1];
				z[j00] = z[j00]+z[j01];
				z[j00+1] = z[j00+1]+z[j01+1];
				z[j01] = t1r;
				z[j01+1] = t1i;
				jt = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
        j2 = j01+jinc;
        if (j2>=jmax) j2 = j2-jmax;

		/* if factor is 3 */
        if (ifac==3) {
			if (mu==1)
				c1 = P866;
			else
				c1 = -P866;
			for (l=0; l<m; l++) {
				t1r = z[j01]+z[j2];
				t1i = z[j01+1]+z[j2+1];
				y1r = z[j00]-0.5*t1r;
				y1i = z[j00+1]-0.5*t1i;
				y2r = c1*(z[j01]-z[j2]);
				y2i = c1*(z[j01+1]-z[j2+1]);
				z[j00] = z[j00]+t1r;
				z[j00+1] = z[j00+1]+t1i;
				z[j01] = y1r-y2i;
				z[j01+1] = y1i+y2r;
				z[j2] = y1r+y2i;
				z[j2+1] = y1i-y2r;
				jt = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j3 = j2+jinc;
		if (j3>=jmax) j3 = j3-jmax;

		/* if factor is 4 */
		if (ifac==4) {
			if (mu==1)
				c1 = 1.0;
			else
				c1 = -1.0;
			for (l=0; l<m; l++) {
				t1r = z[j00]+z[j2];
				t1i = z[j00+1]+z[j2+1];
				t2r = z[j01]+z[j3];
				t2i = z[j01+1]+z[j3+1];
				y1r = z[j00]-z[j2];
				y1i = z[j00+1]-z[j2+1];
				y3r = c1*(z[j01]-z[j3]);
				y3i = c1*(z[j01+1]-z[j3+1]);
				z[j00] = t1r+t2r;
				z[j00+1] = t1i+t2i;
				z[j01] = y1r-y3i;
				z[j01+1] = y1i+y3r;
				z[j2] = t1r-t2r;
				z[j2+1] = t1i-t2i;
				z[j3] = y1r+y3i;
				z[j3+1] = y1i-y3r;
				jt = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j4 = j3+jinc;
		if (j4>=jmax) j4 = j4-jmax;

		/* if factor is 5 */
		if (ifac==5) {
			if (mu==1) {
				c1 = P559;
				c2 = P951;
				c3 = P587;
			} else if (mu==2) {
				c1 = -P559;
				c2 = P587;
				c3 = -P951;
			} else if (mu==3) {
				c1 = -P559;
				c2 = -P587;
				c3 = P951;
			} else { 
				c1 = P559;
				c2 = -P951;
				c3 = -P587;
			}
			for (l=0; l<m; l++) {
				t1r = z[j01]+z[j4];
				t1i = z[j01+1]+z[j4+1];
				t2r = z[j2]+z[j3];
				t2i = z[j2+1]+z[j3+1];
				t3r = z[j01]-z[j4];
				t3i = z[j01+1]-z[j4+1];
				t4r = z[j2]-z[j3];
				t4i = z[j2+1]-z[j3+1];
				t5r = t1r+t2r;
				t5i = t1i+t2i;
				t6r = c1*(t1r-t2r);
				t6i = c1*(t1i-t2i);
				t7r = z[j00]-0.25*t5r;
				t7i = z[j00+1]-0.25*t5i;
				y1r = t7r+t6r;
				y1i = t7i+t6i;
				y2r = t7r-t6r;
				y2i = t7i-t6i;
				y3r = c3*t3r-c2*t4r;
				y3i = c3*t3i-c2*t4i;
				y4r = c2*t3r+c3*t4r;
				y4i = c2*t3i+c3*t4i;
				z[j00] = z[j00]+t5r;
				z[j00+1] = z[j00+1]+t5i;
				z[j01] = y1r-y4i;
				z[j01+1] = y1i+y4r;
				z[j2] = y2r-y3i;
				z[j2+1] = y2i+y3r;
				z[j3] = y2r+y3i;
				z[j3+1] = y2i-y3r;
				z[j4] = y1r+y4i;
				z[j4+1] = y1i-y4r;
				jt = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j5 = j4+jinc;
		if (j5>=jmax) j5 = j5-jmax;
		j6 = j5+jinc;
		if (j6>=jmax) j6 = j6-jmax;

		/* if factor is 7 */
		if (ifac==7) {
			if (mu==1) {
				c1 = P623;
				c2 = -P222;
				c3 = -P900;
				c4 = P781;
				c5 = P974;
				c6 = P433;
			} else if (mu==2) {
				c1 = -P222;
				c2 = -P900;
				c3 = P623;
				c4 = P974;
				c5 = -P433;
				c6 = -P781;
			} else if (mu==3) {
				c1 = -P900;
				c2 = P623;
				c3 = -P222;
				c4 = P433;
				c5 = -P781;
				c6 = P974;
			} else if (mu==4) {
				c1 = -P900;
				c2 = P623;
				c3 = -P222;
				c4 = -P433;
				c5 = P781;
				c6 = -P974;
			} else if (mu==5) {
				c1 = -P222;
				c2 = -P900;
				c3 = P623;
				c4 = -P974;
				c5 = P433;
				c6 = P781;
			} else {
				c1 = P623;
				c2 = -P222;
				c3 = -P900;
				c4 = -P781;
				c5 = -P974;
				c6 = -P433;
			}
			for (l=0; l<m; l++) {
				t1r = z[j01]+z[j6];
				t1i = z[j01+1]+z[j6+1];
				t2r = z[j2]+z[j5];
				t2i = z[j2+1]+z[j5+1];
				t3r = z[j3]+z[j4];
				t3i = z[j3+1]+z[j4+1];
				t4r = z[j01]-z[j6];
				t4i = z[j01+1]-z[j6+1];
				t5r = z[j2]-z[j5];
				t5i = z[j2+1]-z[j5+1];
				t6r = z[j3]-z[j4];
				t6i = z[j3+1]-z[j4+1];
				t7r = z[j00]-0.5*t3r;
				t7i = z[j00+1]-0.5*t3i;
				t8r = t1r-t3r;
				t8i = t1i-t3i;
				t9r = t2r-t3r;
				t9i = t2i-t3i;
				y1r = t7r+c1*t8r+c2*t9r;
				y1i = t7i+c1*t8i+c2*t9i;
				y2r = t7r+c2*t8r+c3*t9r;
				y2i = t7i+c2*t8i+c3*t9i;
				y3r = t7r+c3*t8r+c1*t9r;
				y3i = t7i+c3*t8i+c1*t9i;
				y4r = c6*t4r-c4*t5r+c5*t6r;
				y4i = c6*t4i-c4*t5i+c5*t6i;
				y5r = c5*t4r-c6*t5r-c4*t6r;
				y5i = c5*t4i-c6*t5i-c4*t6i;
				y6r = c4*t4r+c5*t5r+c6*t6r;
				y6i = c4*t4i+c5*t5i+c6*t6i;
				z[j00] = z[j00]+t1r+t2r+t3r;
				z[j00+1] = z[j00+1]+t1i+t2i+t3i;
				z[j01] = y1r-y6i;
				z[j01+1] = y1i+y6r;
				z[j2] = y2r-y5i;
				z[j2+1] = y2i+y5r;
				z[j3] = y3r-y4i;
				z[j3+1] = y3i+y4r;
				z[j4] = y3r+y4i;
				z[j4+1] = y3i-y4r;
				z[j5] = y2r+y5i;
				z[j5+1] = y2i-y5r;
				z[j6] = y1r+y6i;
				z[j6+1] = y1i-y6r;
				jt = j6+2;
				j6 = j5+2;
				j5 = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j7 = j6+jinc;
		if (j7>=jmax) j7 = j7-jmax;

		/* if factor is 8 */
		if (ifac==8) {
			if (mu==1) {
				c1 = 1.0;
				c2 = P707;
			} else if (mu==3) {
				c1 = -1.0;
				c2 = -P707;
			} else if (mu==5) {
				c1 = 1.0;
				c2 = -P707;
			} else {
				c1 = -1.0;
				c2 = P707;
			}
			c3 = c1*c2;
			for (l=0; l<m; l++) {
				t1r = z[j00]+z[j4];
				t1i = z[j00+1]+z[j4+1];
				t2r = z[j00]-z[j4];
				t2i = z[j00+1]-z[j4+1];
				t3r = z[j01]+z[j5];
				t3i = z[j01+1]+z[j5+1];
				t4r = z[j01]-z[j5];
				t4i = z[j01+1]-z[j5+1];
				t5r = z[j2]+z[j6];
				t5i = z[j2+1]+z[j6+1];
				t6r = c1*(z[j2]-z[j6]);
				t6i = c1*(z[j2+1]-z[j6+1]);
				t7r = z[j3]+z[j7];
				t7i = z[j3+1]+z[j7+1];
				t8r = z[j3]-z[j7];
				t8i = z[j3+1]-z[j7+1];
				t9r = t1r+t5r;
				t9i = t1i+t5i;
				t10r = t3r+t7r;
				t10i = t3i+t7i;
				t11r = c2*(t4r-t8r);
				t11i = c2*(t4i-t8i);
				t12r = c3*(t4r+t8r);
				t12i = c3*(t4i+t8i);
				y1r = t2r+t11r;
				y1i = t2i+t11i;
				y2r = t1r-t5r;
				y2i = t1i-t5i;
				y3r = t2r-t11r;
				y3i = t2i-t11i;
				y5r = t12r-t6r;
				y5i = t12i-t6i;
				y6r = c1*(t3r-t7r);
				y6i = c1*(t3i-t7i);
				y7r = t12r+t6r;
				y7i = t12i+t6i;
				z[j00] = t9r+t10r;
				z[j00+1] = t9i+t10i;
				z[j01] = y1r-y7i;
				z[j01+1] = y1i+y7r;
				z[j2] = y2r-y6i;
				z[j2+1] = y2i+y6r;
				z[j3] = y3r-y5i;
				z[j3+1] = y3i+y5r;
				z[j4] = t9r-t10r;
				z[j4+1] = t9i-t10i;
				z[j5] = y3r+y5i;
				z[j5+1] = y3i-y5r;
				z[j6] = y2r+y6i;
				z[j6+1] = y2i-y6r;
				z[j7] = y1r+y7i;
				z[j7+1] = y1i-y7r;
				jt = j7+2;
				j7 = j6+2;
				j6 = j5+2;
				j5 = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j8 = j7+jinc;
		if (j8>=jmax) j8 = j8-jmax;

		/* if factor is 9 */
		if (ifac==9) {
			if (mu==1) {
				c1 = P866;
				c2 = P766;
				c3 = P642;
				c4 = P173;
				c5 = P984;
			} else if (mu==2) {
				c1 = -P866;
				c2 = P173;
				c3 = P984;
				c4 = -P939;
				c5 = P342;
			} else if (mu==4) {
				c1 = P866;
				c2 = -P939;
				c3 = P342;
				c4 = P766;
				c5 = -P642;
			} else if (mu==5) {
				c1 = -P866;
				c2 = -P939;
				c3 = -P342;
				c4 = P766;
				c5 = P642;
			} else if (mu==7) {
				c1 = P866;
				c2 = P173;
				c3 = -P984;
				c4 = -P939;
				c5 = -P342;
			} else {
				c1 = -P866;
				c2 = P766;
				c3 = -P642;
				c4 = P173;
				c5 = -P984;
			}
			c6 = c1*c2;
			c7 = c1*c3;
			c8 = c1*c4;
			c9 = c1*c5;
			for (l=0; l<m; l++) {
				t1r = z[j3]+z[j6];
				t1i = z[j3+1]+z[j6+1];
				t2r = z[j00]-0.5*t1r;
				t2i = z[j00+1]-0.5*t1i;
				t3r = c1*(z[j3]-z[j6]);
				t3i = c1*(z[j3+1]-z[j6+1]);
				t4r = z[j00]+t1r;
				t4i = z[j00+1]+t1i;
				t5r = z[j4]+z[j7];
				t5i = z[j4+1]+z[j7+1];
				t6r = z[j01]-0.5*t5r;
				t6i = z[j01+1]-0.5*t5i;
				t7r = z[j4]-z[j7];
				t7i = z[j4+1]-z[j7+1];
				t8r = z[j01]+t5r;
				t8i = z[j01+1]+t5i;
				t9r = z[j2]+z[j5];
				t9i = z[j2+1]+z[j5+1];
				t10r = z[j8]-0.5*t9r;
				t10i = z[j8+1]-0.5*t9i;
				t11r = z[j2]-z[j5];
				t11i = z[j2+1]-z[j5+1];
				t12r = z[j8]+t9r;
				t12i = z[j8+1]+t9i;
				t13r = t8r+t12r;
				t13i = t8i+t12i;
				t14r = t6r+t10r;
				t14i = t6i+t10i;
				t15r = t6r-t10r;
				t15i = t6i-t10i;
				t16r = t7r+t11r;
				t16i = t7i+t11i;
				t17r = t7r-t11r;
				t17i = t7i-t11i;
				t18r = c2*t14r-c7*t17r;
				t18i = c2*t14i-c7*t17i;
				t19r = c4*t14r+c9*t17r;
				t19i = c4*t14i+c9*t17i;
				t20r = c3*t15r+c6*t16r;
				t20i = c3*t15i+c6*t16i;
				t21r = c5*t15r-c8*t16r;
				t21i = c5*t15i-c8*t16i;
				t22r = t18r+t19r;
				t22i = t18i+t19i;
				t23r = t20r-t21r;
				t23i = t20i-t21i;
				y1r = t2r+t18r;
				y1i = t2i+t18i;
				y2r = t2r+t19r;
				y2i = t2i+t19i;
				y3r = t4r-0.5*t13r;
				y3i = t4i-0.5*t13i;
				y4r = t2r-t22r;
				y4i = t2i-t22i;
				y5r = t3r-t23r;
				y5i = t3i-t23i;
				y6r = c1*(t8r-t12r);
				y6i = c1*(t8i-t12i);
				y7r = t21r-t3r;
				y7i = t21i-t3i;
				y8r = t3r+t20r;
				y8i = t3i+t20i;
				z[j00] = t4r+t13r;
				z[j00+1] = t4i+t13i;
				z[j01] = y1r-y8i;
				z[j01+1] = y1i+y8r;
				z[j2] = y2r-y7i;
				z[j2+1] = y2i+y7r;
				z[j3] = y3r-y6i;
				z[j3+1] = y3i+y6r;
				z[j4] = y4r-y5i;
				z[j4+1] = y4i+y5r;
				z[j5] = y4r+y5i;
				z[j5+1] = y4i-y5r;
				z[j6] = y3r+y6i;
				z[j6+1] = y3i-y6r;
				z[j7] = y2r+y7i;
				z[j7+1] = y2i-y7r;
				z[j8] = y1r+y8i;
				z[j8+1] = y1i-y8r;
				jt = j8+2;
				j8 = j7+2;
				j7 = j6+2;
				j6 = j5+2;
				j5 = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j9 = j8+jinc;
		if (j9>=jmax) j9 = j9-jmax;
		j10 = j9+jinc;
		if (j10>=jmax) j10 = j10-jmax;

		/* if factor is 11 */
		if (ifac==11) {
			if (mu==1) {
				c1 = P841;
				c2 = P415;
				c3 = -P142;
				c4 = -P654;
				c5 = -P959;
				c6 = P540;
				c7 = P909;
				c8 = P989;
				c9 = P755;
				c10 = P281;
			} else if (mu==2) {
				c1 = P415;
				c2 = -P654;
				c3 = -P959;
				c4 = -P142;
				c5 = P841;
				c6 = P909;
				c7 = P755;
				c8 = -P281;
				c9 = -P989;
				c10 = -P540;
			} else if (mu==3) {
				c1 = -P142;
				c2 = -P959;
				c3 = P415;
				c4 = P841;
				c5 = -P654;
				c6 = P989;
				c7 = -P281;
				c8 = -P909;
				c9 = P540;
				c10 = P755;
			} else if (mu==4) {
				c1 = -P654;
				c2 = -P142;
				c3 = P841;
				c4 = -P959;
				c5 = P415;
				c6 = P755;
				c7 = -P989;
				c8 = P540;
				c9 = P281;
				c10 = -P909;
			} else if (mu==5) {
				c1 = -P959;
				c2 = P841;
				c3 = -P654;
				c4 = P415;
				c5 = -P142;
				c6 = P281;
				c7 = -P540;
				c8 = P755;
				c9 = -P909;
				c10 = P989;
			} else if (mu==6) {
				c1 = -P959;
				c2 = P841;
				c3 = -P654;
				c4 = P415;
				c5 = -P142;
				c6 = -P281;
				c7 = P540;
				c8 = -P755;
				c9 = P909;
				c10 = -P989;
			} else if (mu==7) {
				c1 = -P654;
				c2 = -P142;
				c3 = P841;
				c4 = -P959;
				c5 = P415;
				c6 = -P755;
				c7 = P989;
				c8 = -P540;
				c9 = -P281;
				c10 = P909;
			} else if (mu==8) {
				c1 = -P142;
				c2 = -P959;
				c3 = P415;
				c4 = P841;
				c5 = -P654;
				c6 = -P989;
				c7 = P281;
				c8 = P909;
				c9 = -P540;
				c10 = -P755;
			} else if (mu==9) {
				c1 = P415;
				c2 = -P654;
				c3 = -P959;
				c4 = -P142;
				c5 = P841;
				c6 = -P909;
				c7 = -P755;
				c8 = P281;
				c9 = P989;
				c10 = P540;
			} else {
				c1 = P841;
				c2 = P415;
				c3 = -P142;
				c4 = -P654;
				c5 = -P959;
				c6 = -P540;
				c7 = -P909;
				c8 = -P989;
				c9 = -P755;
				c10 = -P281;
			}
			for (l=0; l<m; l++) {
				t1r = z[j01]+z[j10];
				t1i = z[j01+1]+z[j10+1];
				t2r = z[j2]+z[j9];
				t2i = z[j2+1]+z[j9+1];
				t3r = z[j3]+z[j8];
				t3i = z[j3+1]+z[j8+1];
				t4r = z[j4]+z[j7];
				t4i = z[j4+1]+z[j7+1];
				t5r = z[j5]+z[j6];
				t5i = z[j5+1]+z[j6+1];
				t6r = z[j01]-z[j10];
				t6i = z[j01+1]-z[j10+1];
				t7r = z[j2]-z[j9];
				t7i = z[j2+1]-z[j9+1];
				t8r = z[j3]-z[j8];
				t8i = z[j3+1]-z[j8+1];
				t9r = z[j4]-z[j7];
				t9i = z[j4+1]-z[j7+1];
				t10r = z[j5]-z[j6];
				t10i = z[j5+1]-z[j6+1];
				t11r = z[j00]-0.5*t5r;
				t11i = z[j00+1]-0.5*t5i;
				t12r = t1r-t5r;
				t12i = t1i-t5i;
				t13r = t2r-t5r;
				t13i = t2i-t5i;
				t14r = t3r-t5r;
				t14i = t3i-t5i;
				t15r = t4r-t5r;
				t15i = t4i-t5i;
				y1r = t11r+c1*t12r+c2*t13r+c3*t14r+c4*t15r;
				y1i = t11i+c1*t12i+c2*t13i+c3*t14i+c4*t15i;
				y2r = t11r+c2*t12r+c4*t13r+c5*t14r+c3*t15r;
				y2i = t11i+c2*t12i+c4*t13i+c5*t14i+c3*t15i;
				y3r = t11r+c3*t12r+c5*t13r+c2*t14r+c1*t15r;
				y3i = t11i+c3*t12i+c5*t13i+c2*t14i+c1*t15i;
				y4r = t11r+c4*t12r+c3*t13r+c1*t14r+c5*t15r;
				y4i = t11i+c4*t12i+c3*t13i+c1*t14i+c5*t15i;
				y5r = t11r+c5*t12r+c1*t13r+c4*t14r+c2*t15r;
				y5i = t11i+c5*t12i+c1*t13i+c4*t14i+c2*t15i;
				y6r = c10*t6r-c6*t7r+c9*t8r-c7*t9r+c8*t10r;
				y6i = c10*t6i-c6*t7i+c9*t8i-c7*t9i+c8*t10i;
				y7r = c9*t6r-c8*t7r+c6*t8r+c10*t9r-c7*t10r;
				y7i = c9*t6i-c8*t7i+c6*t8i+c10*t9i-c7*t10i;
				y8r = c8*t6r-c10*t7r-c7*t8r+c6*t9r+c9*t10r;
				y8i = c8*t6i-c10*t7i-c7*t8i+c6*t9i+c9*t10i;
				y9r = c7*t6r+c9*t7r-c10*t8r-c8*t9r-c6*t10r;
				y9i = c7*t6i+c9*t7i-c10*t8i-c8*t9i-c6*t10i;
				y10r = c6*t6r+c7*t7r+c8*t8r+c9*t9r+c10*t10r;
				y10i = c6*t6i+c7*t7i+c8*t8i+c9*t9i+c10*t10i;
				z[j00] = z[j00]+t1r+t2r+t3r+t4r+t5r;
				z[j00+1] = z[j00+1]+t1i+t2i+t3i+t4i+t5i;
				z[j01] = y1r-y10i;
				z[j01+1] = y1i+y10r;
				z[j2] = y2r-y9i;
				z[j2+1] = y2i+y9r;
				z[j3] = y3r-y8i;
				z[j3+1] = y3i+y8r;
				z[j4] = y4r-y7i;
				z[j4+1] = y4i+y7r;
				z[j5] = y5r-y6i;
				z[j5+1] = y5i+y6r;
				z[j6] = y5r+y6i;
				z[j6+1] = y5i-y6r;
				z[j7] = y4r+y7i;
				z[j7+1] = y4i-y7r;
				z[j8] = y3r+y8i;
				z[j8+1] = y3i-y8r;
				z[j9] = y2r+y9i;
				z[j9+1] = y2i-y9r;
				z[j10] = y1r+y10i;
				z[j10+1] = y1i-y10r;
				jt = j10+2;
				j10 = j9+2;
				j9 = j8+2;
				j8 = j7+2;
				j7 = j6+2;
				j6 = j5+2;
				j5 = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j11 = j10+jinc;
		if (j11>=jmax) j11 = j11-jmax;
		j12 = j11+jinc;
		if (j12>=jmax) j12 = j12-jmax;

		/* if factor is 13 */
		if (ifac==13) {
			if (mu==1) {
				c1 = P885;
				c2 = P568;
				c3 = P120;
				c4 = -P354;
				c5 = -P748;
				c6 = -P970;
				c7 = P464;
				c8 = P822;
				c9 = P992;
				c10 = P935;
				c11 = P663;
				c12 = P239;
			} else if (mu==2) {
				c1 = P568;
				c2 = -P354;
				c3 = -P970;
				c4 = -P748;
				c5 = P120;
				c6 = P885;
				c7 = P822;
				c8 = P935;
				c9 = P239;
				c10 = -P663;
				c11 = -P992;
				c12 = -P464;
			} else if (mu==3) {
				c1 = P120;
				c2 = -P970;
				c3 = -P354;
				c4 = P885;
				c5 = P568;
				c6 = -P748;
				c7 = P992;
				c8 = P239;
				c9 = -P935;
				c10 = -P464;
				c11 = P822;
				c12 = P663;
			} else if (mu==4) {
				c1 = -P354;
				c2 = -P748;
				c3 = P885;
				c4 = P120;
				c5 = -P970;
				c6 = P568;
				c7 = P935;
				c8 = -P663;
				c9 = -P464;
				c10 = P992;
				c11 = -P239;
				c12 = -P822;
			} else if (mu==5) {
				c1 = -P748;
				c2 = P120;
				c3 = P568;
				c4 = -P970;
				c5 = P885;
				c6 = -P354;
				c7 = P663;
				c8 = -P992;
				c9 = P822;
				c10 = -P239;
				c11 = -P464;
				c12 = P935;
			} else if (mu==6) {
				c1 = -P970;
				c2 = P885;
				c3 = -P748;
				c4 = P568;
				c5 = -P354;
				c6 = P120;
				c7 = P239;
				c8 = -P464;
				c9 = P663;
				c10 = -P822;
				c11 = P935;
				c12 = -P992;
			} else if (mu==7) {
				c1 = -P970;
				c2 = P885;
				c3 = -P748;
				c4 = P568;
				c5 = -P354;
				c6 = P120;
				c7 = -P239;
				c8 = P464;
				c9 = -P663;
				c10 = P822;
				c11 = -P935;
				c12 = P992;
			} else if (mu==8) {
				c1 = -P748;
				c2 = P120;
				c3 = P568;
				c4 = -P970;
				c5 = P885;
				c6 = -P354;
				c7 = -P663;
				c8 = P992;
				c9 = -P822;
				c10 = P239;
				c11 = P464;
				c12 = -P935;
			} else if (mu==9) {
				c1 = -P354;
				c2 = -P748;
				c3 = P885;
				c4 = P120;
				c5 = -P970;
				c6 = P568;
				c7 = -P935;
				c8 = P663;
				c9 = P464;
				c10 = -P992;
				c11 = P239;
				c12 = P822;
			} else if (mu==10) {
				c1 = P120;
				c2 = -P970;
				c3 = -P354;
				c4 = P885;
				c5 = P568;
				c6 = -P748;
				c7 = -P992;
				c8 = -P239;
				c9 = P935;
				c10 = P464;
				c11 = -P822;
				c12 = -P663;
			} else if (mu==11) {
				c1 = P568;
				c2 = -P354;
				c3 = -P970;
				c4 = -P748;
				c5 = P120;
				c6 = P885;
				c7 = -P822;
				c8 = -P935;
				c9 = -P239;
				c10 = P663;
				c11 = P992;
				c12 = P464;
			} else {
				c1 = P885;
				c2 = P568;
				c3 = P120;
				c4 = -P354;
				c5 = -P748;
				c6 = -P970;
				c7 = -P464;
				c8 = -P822;
				c9 = -P992;
				c10 = -P935;
				c11 = -P663;
				c12 = -P239;
			}
			for (l=0; l<m; l++) {
				t1r = z[j01]+z[j12];
				t1i = z[j01+1]+z[j12+1];
				t2r = z[j2]+z[j11];
				t2i = z[j2+1]+z[j11+1];
				t3r = z[j3]+z[j10];
				t3i = z[j3+1]+z[j10+1];
				t4r = z[j4]+z[j9];
				t4i = z[j4+1]+z[j9+1];
				t5r = z[j5]+z[j8];
				t5i = z[j5+1]+z[j8+1];
				t6r = z[j6]+z[j7];
				t6i = z[j6+1]+z[j7+1];
				t7r = z[j01]-z[j12];
				t7i = z[j01+1]-z[j12+1];
				t8r = z[j2]-z[j11];
				t8i = z[j2+1]-z[j11+1];
				t9r = z[j3]-z[j10];
				t9i = z[j3+1]-z[j10+1];
				t10r = z[j4]-z[j9];
				t10i = z[j4+1]-z[j9+1];
				t11r = z[j5]-z[j8];
				t11i = z[j5+1]-z[j8+1];
				t12r = z[j6]-z[j7];
				t12i = z[j6+1]-z[j7+1];
				t13r = z[j00]-0.5*t6r;
				t13i = z[j00+1]-0.5*t6i;
				t14r = t1r-t6r;
				t14i = t1i-t6i;
				t15r = t2r-t6r;
				t15i = t2i-t6i;
				t16r = t3r-t6r;
				t16i = t3i-t6i;
				t17r = t4r-t6r;
				t17i = t4i-t6i;
				t18r = t5r-t6r;
				t18i = t5i-t6i;
				y1r = t13r+c1*t14r+c2*t15r+c3*t16r+c4*t17r+c5*t18r;
				y1i = t13i+c1*t14i+c2*t15i+c3*t16i+c4*t17i+c5*t18i;
				y2r = t13r+c2*t14r+c4*t15r+c6*t16r+c5*t17r+c3*t18r;
				y2i = t13i+c2*t14i+c4*t15i+c6*t16i+c5*t17i+c3*t18i;
				y3r = t13r+c3*t14r+c6*t15r+c4*t16r+c1*t17r+c2*t18r;
				y3i = t13i+c3*t14i+c6*t15i+c4*t16i+c1*t17i+c2*t18i;
				y4r = t13r+c4*t14r+c5*t15r+c1*t16r+c3*t17r+c6*t18r;
				y4i = t13i+c4*t14i+c5*t15i+c1*t16i+c3*t17i+c6*t18i;
				y5r = t13r+c5*t14r+c3*t15r+c2*t16r+c6*t17r+c1*t18r;
				y5i = t13i+c5*t14i+c3*t15i+c2*t16i+c6*t17i+c1*t18i;
				y6r = t13r+c6*t14r+c1*t15r+c5*t16r+c2*t17r+c4*t18r;
				y6i = t13i+c6*t14i+c1*t15i+c5*t16i+c2*t17i+c4*t18i;
				y7r = c12*t7r-c7*t8r+c11*t9r-c8*t10r+c10*t11r-c9*t12r;
				y7i = c12*t7i-c7*t8i+c11*t9i-c8*t10i+c10*t11i-c9*t12i;
				y8r = c11*t7r-c9*t8r+c8*t9r-c12*t10r-c7*t11r+c10*t12r;
				y8i = c11*t7i-c9*t8i+c8*t9i-c12*t10i-c7*t11i+c10*t12i;
				y9r = c10*t7r-c11*t8r-c7*t9r+c9*t10r-c12*t11r-c8*t12r;
				y9i = c10*t7i-c11*t8i-c7*t9i+c9*t10i-c12*t11i-c8*t12i;
				y10r = c9*t7r+c12*t8r-c10*t9r-c7*t10r+c8*t11r+c11*t12r;
				y10i = c9*t7i+c12*t8i-c10*t9i-c7*t10i+c8*t11i+c11*t12i;
				y11r = c8*t7r+c10*t8r+c12*t9r-c11*t10r-c9*t11r-c7*t12r;
				y11i = c8*t7i+c10*t8i+c12*t9i-c11*t10i-c9*t11i-c7*t12i;
				y12r = c7*t7r+c8*t8r+c9*t9r+c10*t10r+c11*t11r+c12*t12r;
				y12i = c7*t7i+c8*t8i+c9*t9i+c10*t10i+c11*t11i+c12*t12i;
				z[j00] = z[j00]+t1r+t2r+t3r+t4r+t5r+t6r;
				z[j00+1] = z[j00+1]+t1i+t2i+t3i+t4i+t5i+t6i;
				z[j01] = y1r-y12i;
				z[j01+1] = y1i+y12r;
				z[j2] = y2r-y11i;
				z[j2+1] = y2i+y11r;
				z[j3] = y3r-y10i;
				z[j3+1] = y3i+y10r;
				z[j4] = y4r-y9i;
				z[j4+1] = y4i+y9r;
				z[j5] = y5r-y8i;
				z[j5+1] = y5i+y8r;
				z[j6] = y6r-y7i;
				z[j6+1] = y6i+y7r;
				z[j7] = y6r+y7i;
				z[j7+1] = y6i-y7r;
				z[j8] = y5r+y8i;
				z[j8+1] = y5i-y8r;
				z[j9] = y4r+y9i;
				z[j9+1] = y4i-y9r;
				z[j10] = y3r+y10i;
				z[j10+1] = y3i-y10r;
				z[j11] = y2r+y11i;
				z[j11+1] = y2i-y11r;
				z[j12] = y1r+y12i;
				z[j12+1] = y1i-y12r;
				jt = j12+2;
				j12 = j11+2;
				j11 = j10+2;
				j10 = j9+2;
				j9 = j8+2;
				j8 = j7+2;
				j7 = j6+2;
				j6 = j5+2;
				j5 = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
		j13 = j12+jinc;
		if (j13>=jmax) j13 = j13-jmax;
		j14 = j13+jinc;
		if (j14>=jmax) j14 = j14-jmax;
		j15 = j14+jinc;
		if (j15>=jmax) j15 = j15-jmax;

		/* if factor is 16 */
		if (ifac==16) {
			if (mu==1) {
				c1 = 1.0;
				c2 = P923;
				c3 = P382;
				c4 = P707;
			} else if (mu==3) {
				c1 = -1.0;
				c2 = P382;
				c3 = P923;
				c4 = -P707;
			} else if (mu==5) {
				c1 = 1.0;
				c2 = -P382;
				c3 = P923;
				c4 = -P707;
			} else if (mu==7) {
				c1 = -1.0;
				c2 = -P923;
				c3 = P382;
				c4 = P707;
			} else if (mu==9) {
				c1 = 1.0;
				c2 = -P923;
				c3 = -P382;
				c4 = P707;
			} else if (mu==11) {
				c1 = -1.0;
				c2 = -P382;
				c3 = -P923;
				c4 = -P707;
			} else if (mu==13) {
				c1 = 1.0;
				c2 = P382;
				c3 = -P923;
				c4 = -P707;
			} else {
				c1 = -1.0;
				c2 = P923;
				c3 = -P382;
				c4 = P707;
			}
			c5 = c1*c4;
			c6 = c1*c3;
			c7 = c1*c2;
			for (l=0; l<m; l++) {
				t1r = z[j00]+z[j8];
				t1i = z[j00+1]+z[j8+1];
				t2r = z[j4]+z[j12];
				t2i = z[j4+1]+z[j12+1];
				t3r = z[j00]-z[j8];
				t3i = z[j00+1]-z[j8+1];
				t4r = c1*(z[j4]-z[j12]);
				t4i = c1*(z[j4+1]-z[j12+1]);
				t5r = t1r+t2r;
				t5i = t1i+t2i;
				t6r = t1r-t2r;
				t6i = t1i-t2i;
				t7r = z[j01]+z[j9];
				t7i = z[j01+1]+z[j9+1];
				t8r = z[j5]+z[j13];
				t8i = z[j5+1]+z[j13+1];
				t9r = z[j01]-z[j9];
				t9i = z[j01+1]-z[j9+1];
				t10r = z[j5]-z[j13];
				t10i = z[j5+1]-z[j13+1];
				t11r = t7r+t8r;
				t11i = t7i+t8i;
				t12r = t7r-t8r;
				t12i = t7i-t8i;
				t13r = z[j2]+z[j10];
				t13i = z[j2+1]+z[j10+1];
				t14r = z[j6]+z[j14];
				t14i = z[j6+1]+z[j14+1];
				t15r = z[j2]-z[j10];
				t15i = z[j2+1]-z[j10+1];
				t16r = z[j6]-z[j14];
				t16i = z[j6+1]-z[j14+1];
				t17r = t13r+t14r;
				t17i = t13i+t14i;
				t18r = c4*(t15r-t16r);
				t18i = c4*(t15i-t16i);
				t19r = c5*(t15r+t16r);
				t19i = c5*(t15i+t16i);
				t20r = c1*(t13r-t14r);
				t20i = c1*(t13i-t14i);
				t21r = z[j3]+z[j11];
				t21i = z[j3+1]+z[j11+1];
				t22r = z[j7]+z[j15];
				t22i = z[j7+1]+z[j15+1];
				t23r = z[j3]-z[j11];
				t23i = z[j3+1]-z[j11+1];
				t24r = z[j7]-z[j15];
				t24i = z[j7+1]-z[j15+1];
				t25r = t21r+t22r;
				t25i = t21i+t22i;
				t26r = t21r-t22r;
				t26i = t21i-t22i;
				t27r = t9r+t24r;
				t27i = t9i+t24i;
				t28r = t10r+t23r;
				t28i = t10i+t23i;
				t29r = t9r-t24r;
				t29i = t9i-t24i;
				t30r = t10r-t23r;
				t30i = t10i-t23i;
				t31r = t5r+t17r;
				t31i = t5i+t17i;
				t32r = t11r+t25r;
				t32i = t11i+t25i;
				t33r = t3r+t18r;
				t33i = t3i+t18i;
				t34r = c2*t29r-c6*t30r;
				t34i = c2*t29i-c6*t30i;
				t35r = t3r-t18r;
				t35i = t3i-t18i;
				t36r = c7*t27r-c3*t28r;
				t36i = c7*t27i-c3*t28i;
				t37r = t4r+t19r;
				t37i = t4i+t19i;
				t38r = c3*t27r+c7*t28r;
				t38i = c3*t27i+c7*t28i;
				t39r = t4r-t19r;
				t39i = t4i-t19i;
				t40r = c6*t29r+c2*t30r;
				t40i = c6*t29i+c2*t30i;
				t41r = c4*(t12r-t26r);
				t41i = c4*(t12i-t26i);
				t42r = c5*(t12r+t26r);
				t42i = c5*(t12i+t26i);
				y1r = t33r+t34r;
				y1i = t33i+t34i;
				y2r = t6r+t41r;
				y2i = t6i+t41i;
				y3r = t35r+t40r;
				y3i = t35i+t40i;
				y4r = t5r-t17r;
				y4i = t5i-t17i;
				y5r = t35r-t40r;
				y5i = t35i-t40i;
				y6r = t6r-t41r;
				y6i = t6i-t41i;
				y7r = t33r-t34r;
				y7i = t33i-t34i;
				y9r = t38r-t37r;
				y9i = t38i-t37i;
				y10r = t42r-t20r;
				y10i = t42i-t20i;
				y11r = t36r+t39r;
				y11i = t36i+t39i;
				y12r = c1*(t11r-t25r);
				y12i = c1*(t11i-t25i);
				y13r = t36r-t39r;
				y13i = t36i-t39i;
				y14r = t42r+t20r;
				y14i = t42i+t20i;
				y15r = t38r+t37r;
				y15i = t38i+t37i;
				z[j00] = t31r+t32r;
				z[j00+1] = t31i+t32i;
				z[j01] = y1r-y15i;
				z[j01+1] = y1i+y15r;
				z[j2] = y2r-y14i;
				z[j2+1] = y2i+y14r;
				z[j3] = y3r-y13i;
				z[j3+1] = y3i+y13r;
				z[j4] = y4r-y12i;
				z[j4+1] = y4i+y12r;
				z[j5] = y5r-y11i;
				z[j5+1] = y5i+y11r;
				z[j6] = y6r-y10i;
				z[j6+1] = y6i+y10r;
				z[j7] = y7r-y9i;
				z[j7+1] = y7i+y9r;
				z[j8] = t31r-t32r;
				z[j8+1] = t31i-t32i;
				z[j9] = y7r+y9i;
				z[j9+1] = y7i-y9r;
				z[j10] = y6r+y10i;
				z[j10+1] = y6i-y10r;
				z[j11] = y5r+y11i;
				z[j11+1] = y5i-y11r;
				z[j12] = y4r+y12i;
				z[j12+1] = y4i-y12r;
				z[j13] = y3r+y13i;
				z[j13+1] = y3i-y13r;
				z[j14] = y2r+y14i;
				z[j14+1] = y2i-y14r;
				z[j15] = y1r+y15i;
				z[j15+1] = y1i-y15r;
				jt = j15+2;
				j15 = j14+2;
				j14 = j13+2;
				j13 = j12+2;
				j12 = j11+2;
				j11 = j10+2;
				j10 = j9+2;
				j9 = j8+2;
				j8 = j7+2;
				j7 = j6+2;
				j6 = j5+2;
				j5 = j4+2;
				j4 = j3+2;
				j3 = j2+2;
				j2 = j01+2;
				j01 = j00+2;
				j00 = jt;
			}
			continue;
		}
	}
}

void sf_pfacr (int isign, int n, float complex cz[], float rz[])
/*****************************************************************************
Prime factor fft:  complex to real transform
******************************************************************************
Input:
isign       sign of isign is the sign of exponent in fourier kernel
n           length of transform (see notes below)
cz          array[n/2+1] of complex values (may be equivalenced to rz)

Output:
rz          array[n] of real values (may be equivalenced to cz)
******************************************************************************
Notes:
Because pfacr uses pfacc to do most of the work, n must be even 
and n/2 must be a valid length for pfacc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.
******************************************************************************
References:  
Press et al, 1988, Numerical Recipes in C, p. 417.

Also, see notes and references for function pfacc.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/13/89
*****************************************************************************/
{
    int i,ir,ii,jr,ji,no2;
    float *z,tempr,tempi,sumr,sumi,difr,difi;
    double wr,wi,wpr,wpi,wtemp,theta;

    /* copy input to output and fix dc and nyquist */
    z = (float*)cz;
    for (i=2; i<n; i++)
        rz[i] = z[i];
    rz[1] = z[0]-z[n];
    rz[0] = z[0]+z[n];
    z = rz;

    /* initialize cosine-sine recurrence */
    theta = 2.0*SF_PI/(double)n;
    if (isign>0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle */
    no2 = n/2;
    for (ir=2,ii=3,jr=n-2,ji=n-1; ir<=no2; ir+=2,ii+=2,jr-=2,ji-=2) {
        sumr = z[ir]+z[jr];
        sumi = z[ii]+z[ji];
        difr = z[ir]-z[jr];
        difi = z[ii]-z[ji];
        tempr = wi*difr-wr*sumi;
        tempi = wi*sumi+wr*difr;
        z[ir] = sumr+tempr;
        z[ii] = difi+tempi;
        z[jr] = sumr-tempr;
        z[ji] = tempi-difi;
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }

    /* do complex to complex transform */
    sf_pfacc(isign,n/2,(float complex*)z);
}

void sf_pfarc (int isign, int n, float rz[], float complex cz[])
/*****************************************************************************
Prime factor fft:  real to complex transform
******************************************************************************
Input:
isign       sign of isign is the sign of exponent in fourier kernel
n           length of transform; must be even (see notes below)
rz          array[n] of real values (may be equivalenced to cz)

Output:
cz          array[n/2+1] of complex values (may be equivalenced to rz)
******************************************************************************
Notes:
Because pfarc uses pfacc to do most of the work, n must be even 
and n/2 must be a valid length for pfacc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.
******************************************************************************
References:  
Press et al, 1988, Numerical Recipes in C, p. 417.

Also, see notes and references for function pfacc.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/13/89
*****************************************************************************/
{
    int i,ir,ii,jr,ji,no2;
    float *z,tempr,tempi,sumr,sumi,difr,difi;
    double wr,wi,wpr,wpi,wtemp,theta;

    /* copy input to output while scaling */
    z = (float*)cz;
    for (i=0; i<n; i++)
        z[i] = 0.5*rz[i];

    /* do complex to complex transform */
    sf_pfacc(isign,n/2,cz);

    /* fix dc and nyquist */
    z[n] = 2.0*(z[0]-z[1]);
    z[0] = 2.0*(z[0]+z[1]);
    z[n+1] = 0.0;
    z[1] = 0.0;

    /* initialize cosine-sine recurrence */
    theta = 2.0*SF_PI/(double)n;
    if (isign<0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle */
    no2 = n/2;
    for (ir=2,ii=3,jr=n-2,ji=n-1; ir<=no2; ir+=2,ii+=2,jr-=2,ji-=2) {
        sumr = z[ir]+z[jr];
        sumi = z[ii]+z[ji];
        difr = z[ir]-z[jr];
        difi = z[ii]-z[ji];
        tempr = wi*difr+wr*sumi;
        tempi = wi*sumi-wr*difr;
        z[ir] = sumr+tempr;
        z[ii] = difi+tempi;
        z[jr] = sumr-tempr;
        z[ji] = tempi-difi;
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }
}

void sf_pfamcc (int isign, int n, int nt, int k, int kt, float complex cz[])
/*****************************************************************************
Prime factor fft:  multiple complex to complex transforms, in place
******************************************************************************
Input:
isign       	sign of isign is the sign of exponent in fourier kernel
n           	number of complex elements per transform (see notes below)
nt          	number of transforms
k           	stride in complex elements within transforms
kt          	stride in complex elements between transforms
z           	array of complex elements to be transformed in place

Output:
z		array of complex elements transformed
******************************************************************************
Notes:
n must be factorable into mutually prime factors taken 
from the set {2,3,4,5,7,8,9,11,13,16}.  in other words,
    n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u
where
    0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1
is required for pfamcc to yield meaningful results.  this
restriction implies that n is restricted to the range
    1 <= n <= 720720 (= 5*7*9*11*13*16)

To perform a two-dimensional transform of an n1 by n2 complex array 
(assuming that both n1 and n2 are valid "n"), stored with n1 fast 
and n2 slow:
    pfamcc(isign,n1,n2,1,n1,z); (to transform 1st dimension)
    pfamcc(isign,n2,n1,n1,1,z); (to transform 2nd dimension)
******************************************************************************
References:  
Temperton, C., 1985, Implementation of a self-sorting
in-place prime factor fft algorithm:  Journal of
Computational Physics, v. 58, p. 283-299.

Temperton, C., 1988, A new set of minimum-add rotated
rotated dft modules: Journal of Computational Physics,
v. 75, p. 190-198.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/15/89
*****************************************************************************/
{
    static int kfax[] = { 16,13,11,9,8,7,5,4,3,2 };
    register float *z=(float*)cz;
    register int j00,j01,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15;
    int nleft,jfax,ifac,jfac,iinc,imax,ndiv,m,mm=0,mu=0,l,istep,jstep,
        jt,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,it;
    float t1r,t1i,t2r,t2i,t3r,t3i,t4r,t4i,t5r,t5i,
        t6r,t6i,t7r,t7i,t8r,t8i,t9r,t9i,t10r,t10i,
        t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,
        t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t20r,t20i,
        t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,
        t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t30r,t30i,
        t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,
        t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t40r,t40i,
        t41r,t41i,t42r,t42i,
        y1r,y1i,y2r,y2i,y3r,y3i,y4r,y4i,y5r,y5i,
        y6r,y6i,y7r,y7i,y8r,y8i,y9r,y9i,y10r,y10i,
        y11r,y11i,y12r,y12i,y13r,y13i,y14r,y14i,y15r,y15i,
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    /* determine step within and between transforms */
    istep = 2*k;
    jstep = 2*kt;

    /* keep track of n left after dividing by factors */
    nleft = n;

    /* begin loop over possible factors (from biggest to smallest) */
    for (jfax=0; jfax<NFAX; jfax++) {

        /* skip if not a mutually prime factor of n */
        ifac = kfax[jfax];
        ndiv = nleft/ifac;
        if (ndiv*ifac!=nleft) continue;
 
        /* update n left and determine n divided by factor */
        nleft = ndiv;
        m = n/ifac;
 
        /* determine rotation factor mu and stride mm */
        for (jfac=1; jfac<=ifac; jfac++) {
            mu = jfac;
            mm = jfac*m;
            if (mm%ifac==1) break;
        }
 
        /* adjust rotation factor for sign of transform */
        if (isign<0) mu = ifac-mu;
 
        /* compute stride, limit, and pointers */
        iinc = istep*mm;
        imax = istep*n;
        i0 = 0;
        i1 = i0+iinc;

        /* if factor is 2 */
        if (ifac==2) {
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j00]-z[j01];
                    t1i = z[j00+1]-z[j01+1];
                    z[j00] = z[j00]+z[j01];
                    z[j00+1] = z[j00+1]+z[j01+1];
                    z[j01] = t1r;
                    z[j01+1] = t1i;
                    j00 += jstep;
                    j01 += jstep;
                }
                it = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i2 = i1+iinc;
        if (i2>=imax) i2 = i2-imax;

        /* if factor is 3 */
        if (ifac==3) {
            if (mu==1)
                c1 = P866;
            else
                c1 = -P866;
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j01]+z[j2];
                    t1i = z[j01+1]+z[j2+1];
                    y1r = z[j00]-0.5*t1r;
                    y1i = z[j00+1]-0.5*t1i;
                    y2r = c1*(z[j01]-z[j2]);
                    y2i = c1*(z[j01+1]-z[j2+1]);
                    z[j00] = z[j00]+t1r;
                    z[j00+1] = z[j00+1]+t1i;
                    z[j01] = y1r-y2i;
                    z[j01+1] = y1i+y2r;
                    z[j2] = y1r+y2i;
                    z[j2+1] = y1i-y2r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                }
                it = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i3 = i2+iinc;
        if (i3>=imax) i3 = i3-imax;

        /* if factor is 4 */
        if (ifac==4) {
            if (mu==1)
                c1 = 1.0;
            else
                c1 = -1.0;
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j00]+z[j2];
                    t1i = z[j00+1]+z[j2+1];
                    t2r = z[j01]+z[j3];
                    t2i = z[j01+1]+z[j3+1];
                    y1r = z[j00]-z[j2];
                    y1i = z[j00+1]-z[j2+1];
                    y3r = c1*(z[j01]-z[j3]);
                    y3i = c1*(z[j01+1]-z[j3+1]);
                    z[j00] = t1r+t2r;
                    z[j00+1] = t1i+t2i;
                    z[j01] = y1r-y3i;
                    z[j01+1] = y1i+y3r;
                    z[j2] = t1r-t2r;
                    z[j2+1] = t1i-t2i;
                    z[j3] = y1r+y3i;
                    z[j3+1] = y1i-y3r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                }
                it = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i4 = i3+iinc;
        if (i4>=imax) i4 = i4-imax;

        /* if factor is 5 */
        if (ifac==5) {
            if (mu==1) {
                c1 = P559;
                c2 = P951;
                c3 = P587;
            } else if (mu==2) {
                c1 = -P559;
                c2 = P587;
                c3 = -P951;
            } else if (mu==3) {
                c1 = -P559;
                c2 = -P587;
                c3 = P951;
            } else { 
                c1 = P559;
                c2 = -P951;
                c3 = -P587;
            }
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j01]+z[j4];
                    t1i = z[j01+1]+z[j4+1];
                    t2r = z[j2]+z[j3];
                    t2i = z[j2+1]+z[j3+1];
                    t3r = z[j01]-z[j4];
                    t3i = z[j01+1]-z[j4+1];
                    t4r = z[j2]-z[j3];
                    t4i = z[j2+1]-z[j3+1];
                    t5r = t1r+t2r;
                    t5i = t1i+t2i;
                    t6r = c1*(t1r-t2r);
                    t6i = c1*(t1i-t2i);
                    t7r = z[j00]-0.25*t5r;
                    t7i = z[j00+1]-0.25*t5i;
                    y1r = t7r+t6r;
                    y1i = t7i+t6i;
                    y2r = t7r-t6r;
                    y2i = t7i-t6i;
                    y3r = c3*t3r-c2*t4r;
                    y3i = c3*t3i-c2*t4i;
                    y4r = c2*t3r+c3*t4r;
                    y4i = c2*t3i+c3*t4i;
                    z[j00] = z[j00]+t5r;
                    z[j00+1] = z[j00+1]+t5i;
                    z[j01] = y1r-y4i;
                    z[j01+1] = y1i+y4r;
                    z[j2] = y2r-y3i;
                    z[j2+1] = y2i+y3r;
                    z[j3] = y2r+y3i;
                    z[j3+1] = y2i-y3r;
                    z[j4] = y1r+y4i;
                    z[j4+1] = y1i-y4r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                }
                it = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i5 = i4+iinc;
        if (i5>=imax) i5 = i5-imax;
        i6 = i5+iinc;
        if (i6>=imax) i6 = i6-imax;

        /* if factor is 7 */
        if (ifac==7) {
            if (mu==1) {
                c1 = P623;
                c2 = -P222;
                c3 = -P900;
                c4 = P781;
                c5 = P974;
                c6 = P433;
            } else if (mu==2) {
                c1 = -P222;
                c2 = -P900;
                c3 = P623;
                c4 = P974;
                c5 = -P433;
                c6 = -P781;
            } else if (mu==3) {
                c1 = -P900;
                c2 = P623;
                c3 = -P222;
                c4 = P433;
                c5 = -P781;
                c6 = P974;
            } else if (mu==4) {
                c1 = -P900;
                c2 = P623;
                c3 = -P222;
                c4 = -P433;
                c5 = P781;
                c6 = -P974;
            } else if (mu==5) {
                c1 = -P222;
                c2 = -P900;
                c3 = P623;
                c4 = -P974;
                c5 = P433;
                c6 = P781;
            } else {
                c1 = P623;
                c2 = -P222;
                c3 = -P900;
                c4 = -P781;
                c5 = -P974;
                c6 = -P433;
            }
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j01]+z[j6];
                    t1i = z[j01+1]+z[j6+1];
                    t2r = z[j2]+z[j5];
                    t2i = z[j2+1]+z[j5+1];
                    t3r = z[j3]+z[j4];
                    t3i = z[j3+1]+z[j4+1];
                    t4r = z[j01]-z[j6];
                    t4i = z[j01+1]-z[j6+1];
                    t5r = z[j2]-z[j5];
                    t5i = z[j2+1]-z[j5+1];
                    t6r = z[j3]-z[j4];
                    t6i = z[j3+1]-z[j4+1];
                    t7r = z[j00]-0.5*t3r;
                    t7i = z[j00+1]-0.5*t3i;
                    t8r = t1r-t3r;
                    t8i = t1i-t3i;
                    t9r = t2r-t3r;
                    t9i = t2i-t3i;
                    y1r = t7r+c1*t8r+c2*t9r;
                    y1i = t7i+c1*t8i+c2*t9i;
                    y2r = t7r+c2*t8r+c3*t9r;
                    y2i = t7i+c2*t8i+c3*t9i;
                    y3r = t7r+c3*t8r+c1*t9r;
                    y3i = t7i+c3*t8i+c1*t9i;
                    y4r = c6*t4r-c4*t5r+c5*t6r;
                    y4i = c6*t4i-c4*t5i+c5*t6i;
                    y5r = c5*t4r-c6*t5r-c4*t6r;
                    y5i = c5*t4i-c6*t5i-c4*t6i;
                    y6r = c4*t4r+c5*t5r+c6*t6r;
                    y6i = c4*t4i+c5*t5i+c6*t6i;
                    z[j00] = z[j00]+t1r+t2r+t3r;
                    z[j00+1] = z[j00+1]+t1i+t2i+t3i;
                    z[j01] = y1r-y6i;
                    z[j01+1] = y1i+y6r;
                    z[j2] = y2r-y5i;
                    z[j2+1] = y2i+y5r;
                    z[j3] = y3r-y4i;
                    z[j3+1] = y3i+y4r;
                    z[j4] = y3r+y4i;
                    z[j4+1] = y3i-y4r;
                    z[j5] = y2r+y5i;
                    z[j5+1] = y2i-y5r;
                    z[j6] = y1r+y6i;
                    z[j6+1] = y1i-y6r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                }
                it = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i7 = i6+iinc;
        if (i7>=imax) i7 = i7-imax;

        /* if factor is 8 */
        if (ifac==8) {
            if (mu==1) {
                c1 = 1.0;
                c2 = P707;
            } else if (mu==3) {
                c1 = -1.0;
                c2 = -P707;
            } else if (mu==5) {
                c1 = 1.0;
                c2 = -P707;
            } else {
                c1 = -1.0;
                c2 = P707;
            }
            c3 = c1*c2;
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j00]+z[j4];
                    t1i = z[j00+1]+z[j4+1];
                    t2r = z[j00]-z[j4];
                    t2i = z[j00+1]-z[j4+1];
                    t3r = z[j01]+z[j5];
                    t3i = z[j01+1]+z[j5+1];
                    t4r = z[j01]-z[j5];
                    t4i = z[j01+1]-z[j5+1];
                    t5r = z[j2]+z[j6];
                    t5i = z[j2+1]+z[j6+1];
                    t6r = c1*(z[j2]-z[j6]);
                    t6i = c1*(z[j2+1]-z[j6+1]);
                    t7r = z[j3]+z[j7];
                    t7i = z[j3+1]+z[j7+1];
                    t8r = z[j3]-z[j7];
                    t8i = z[j3+1]-z[j7+1];
                    t9r = t1r+t5r;
                    t9i = t1i+t5i;
                    t10r = t3r+t7r;
                    t10i = t3i+t7i;
                    t11r = c2*(t4r-t8r);
                    t11i = c2*(t4i-t8i);
                    t12r = c3*(t4r+t8r);
                    t12i = c3*(t4i+t8i);
                    y1r = t2r+t11r;
                    y1i = t2i+t11i;
                    y2r = t1r-t5r;
                    y2i = t1i-t5i;
                    y3r = t2r-t11r;
                    y3i = t2i-t11i;
                    y5r = t12r-t6r;
                    y5i = t12i-t6i;
                    y6r = c1*(t3r-t7r);
                    y6i = c1*(t3i-t7i);
                    y7r = t12r+t6r;
                    y7i = t12i+t6i;
                    z[j00] = t9r+t10r;
                    z[j00+1] = t9i+t10i;
                    z[j01] = y1r-y7i;
                    z[j01+1] = y1i+y7r;
                    z[j2] = y2r-y6i;
                    z[j2+1] = y2i+y6r;
                    z[j3] = y3r-y5i;
                    z[j3+1] = y3i+y5r;
                    z[j4] = t9r-t10r;
                    z[j4+1] = t9i-t10i;
                    z[j5] = y3r+y5i;
                    z[j5+1] = y3i-y5r;
                    z[j6] = y2r+y6i;
                    z[j6+1] = y2i-y6r;
                    z[j7] = y1r+y7i;
                    z[j7+1] = y1i-y7r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                }
                it = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i8 = i7+iinc;
        if (i8>=imax) i8 = i8-imax;

        /* if factor is 9 */
        if (ifac==9) {
            if (mu==1) {
                c1 = P866;
                c2 = P766;
                c3 = P642;
                c4 = P173;
                c5 = P984;
            } else if (mu==2) {
                c1 = -P866;
                c2 = P173;
                c3 = P984;
                c4 = -P939;
                c5 = P342;
            } else if (mu==4) {
                c1 = P866;
                c2 = -P939;
                c3 = P342;
                c4 = P766;
                c5 = -P642;
            } else if (mu==5) {
                c1 = -P866;
                c2 = -P939;
                c3 = -P342;
                c4 = P766;
                c5 = P642;
            } else if (mu==7) {
                c1 = P866;
                c2 = P173;
                c3 = -P984;
                c4 = -P939;
                c5 = -P342;
            } else {
                c1 = -P866;
                c2 = P766;
                c3 = -P642;
                c4 = P173;
                c5 = -P984;
            }
            c6 = c1*c2;
            c7 = c1*c3;
            c8 = c1*c4;
            c9 = c1*c5;
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j3]+z[j6];
                    t1i = z[j3+1]+z[j6+1];
                    t2r = z[j00]-0.5*t1r;
                    t2i = z[j00+1]-0.5*t1i;
                    t3r = c1*(z[j3]-z[j6]);
                    t3i = c1*(z[j3+1]-z[j6+1]);
                    t4r = z[j00]+t1r;
                    t4i = z[j00+1]+t1i;
                    t5r = z[j4]+z[j7];
                    t5i = z[j4+1]+z[j7+1];
                    t6r = z[j01]-0.5*t5r;
                    t6i = z[j01+1]-0.5*t5i;
                    t7r = z[j4]-z[j7];
                    t7i = z[j4+1]-z[j7+1];
                    t8r = z[j01]+t5r;
                    t8i = z[j01+1]+t5i;
                    t9r = z[j2]+z[j5];
                    t9i = z[j2+1]+z[j5+1];
                    t10r = z[j8]-0.5*t9r;
                    t10i = z[j8+1]-0.5*t9i;
                    t11r = z[j2]-z[j5];
                    t11i = z[j2+1]-z[j5+1];
                    t12r = z[j8]+t9r;
                    t12i = z[j8+1]+t9i;
                    t13r = t8r+t12r;
                    t13i = t8i+t12i;
                    t14r = t6r+t10r;
                    t14i = t6i+t10i;
                    t15r = t6r-t10r;
                    t15i = t6i-t10i;
                    t16r = t7r+t11r;
                    t16i = t7i+t11i;
                    t17r = t7r-t11r;
                    t17i = t7i-t11i;
                    t18r = c2*t14r-c7*t17r;
                    t18i = c2*t14i-c7*t17i;
                    t19r = c4*t14r+c9*t17r;
                    t19i = c4*t14i+c9*t17i;
                    t20r = c3*t15r+c6*t16r;
                    t20i = c3*t15i+c6*t16i;
                    t21r = c5*t15r-c8*t16r;
                    t21i = c5*t15i-c8*t16i;
                    t22r = t18r+t19r;
                    t22i = t18i+t19i;
                    t23r = t20r-t21r;
                    t23i = t20i-t21i;
                    y1r = t2r+t18r;
                    y1i = t2i+t18i;
                    y2r = t2r+t19r;
                    y2i = t2i+t19i;
                    y3r = t4r-0.5*t13r;
                    y3i = t4i-0.5*t13i;
                    y4r = t2r-t22r;
                    y4i = t2i-t22i;
                    y5r = t3r-t23r;
                    y5i = t3i-t23i;
                    y6r = c1*(t8r-t12r);
                    y6i = c1*(t8i-t12i);
                    y7r = t21r-t3r;
                    y7i = t21i-t3i;
                    y8r = t3r+t20r;
                    y8i = t3i+t20i;
                    z[j00] = t4r+t13r;
                    z[j00+1] = t4i+t13i;
                    z[j01] = y1r-y8i;
                    z[j01+1] = y1i+y8r;
                    z[j2] = y2r-y7i;
                    z[j2+1] = y2i+y7r;
                    z[j3] = y3r-y6i;
                    z[j3+1] = y3i+y6r;
                    z[j4] = y4r-y5i;
                    z[j4+1] = y4i+y5r;
                    z[j5] = y4r+y5i;
                    z[j5+1] = y4i-y5r;
                    z[j6] = y3r+y6i;
                    z[j6+1] = y3i-y6r;
                    z[j7] = y2r+y7i;
                    z[j7+1] = y2i-y7r;
                    z[j8] = y1r+y8i;
                    z[j8+1] = y1i-y8r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                }
                it = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i9 = i8+iinc;
        if (i9>=imax) i9 = i9-imax;
        i10 = i9+iinc;
        if (i10>=imax) i10 = i10-imax;

        /* if factor is 11 */
        if (ifac==11) {
            if (mu==1) {
                c1 = P841;
                c2 = P415;
                c3 = -P142;
                c4 = -P654;
                c5 = -P959;
                c6 = P540;
                c7 = P909;
                c8 = P989;
                c9 = P755;
                c10 = P281;
            } else if (mu==2) {
                c1 = P415;
                c2 = -P654;
                c3 = -P959;
                c4 = -P142;
                c5 = P841;
                c6 = P909;
                c7 = P755;
                c8 = -P281;
                c9 = -P989;
                c10 = -P540;
            } else if (mu==3) {
                c1 = -P142;
                c2 = -P959;
                c3 = P415;
                c4 = P841;
                c5 = -P654;
                c6 = P989;
                c7 = -P281;
                c8 = -P909;
                c9 = P540;
                c10 = P755;
            } else if (mu==4) {
                c1 = -P654;
                c2 = -P142;
                c3 = P841;
                c4 = -P959;
                c5 = P415;
                c6 = P755;
                c7 = -P989;
                c8 = P540;
                c9 = P281;
                c10 = -P909;
            } else if (mu==5) {
                c1 = -P959;
                c2 = P841;
                c3 = -P654;
                c4 = P415;
                c5 = -P142;
                c6 = P281;
                c7 = -P540;
                c8 = P755;
                c9 = -P909;
                c10 = P989;
            } else if (mu==6) {
                c1 = -P959;
                c2 = P841;
                c3 = -P654;
                c4 = P415;
                c5 = -P142;
                c6 = -P281;
                c7 = P540;
                c8 = -P755;
                c9 = P909;
                c10 = -P989;
            } else if (mu==7) {
                c1 = -P654;
                c2 = -P142;
                c3 = P841;
                c4 = -P959;
                c5 = P415;
                c6 = -P755;
                c7 = P989;
                c8 = -P540;
                c9 = -P281;
                c10 = P909;
            } else if (mu==8) {
                c1 = -P142;
                c2 = -P959;
                c3 = P415;
                c4 = P841;
                c5 = -P654;
                c6 = -P989;
                c7 = P281;
                c8 = P909;
                c9 = -P540;
                c10 = -P755;
            } else if (mu==9) {
                c1 = P415;
                c2 = -P654;
                c3 = -P959;
                c4 = -P142;
                c5 = P841;
                c6 = -P909;
                c7 = -P755;
                c8 = P281;
                c9 = P989;
                c10 = P540;
            } else {
                c1 = P841;
                c2 = P415;
                c3 = -P142;
                c4 = -P654;
                c5 = -P959;
                c6 = -P540;
                c7 = -P909;
                c8 = -P989;
                c9 = -P755;
                c10 = -P281;
            }
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j01]+z[j10];
                    t1i = z[j01+1]+z[j10+1];
                    t2r = z[j2]+z[j9];
                    t2i = z[j2+1]+z[j9+1];
                    t3r = z[j3]+z[j8];
                    t3i = z[j3+1]+z[j8+1];
                    t4r = z[j4]+z[j7];
                    t4i = z[j4+1]+z[j7+1];
                    t5r = z[j5]+z[j6];
                    t5i = z[j5+1]+z[j6+1];
                    t6r = z[j01]-z[j10];
                    t6i = z[j01+1]-z[j10+1];
                    t7r = z[j2]-z[j9];
                    t7i = z[j2+1]-z[j9+1];
                    t8r = z[j3]-z[j8];
                    t8i = z[j3+1]-z[j8+1];
                    t9r = z[j4]-z[j7];
                    t9i = z[j4+1]-z[j7+1];
                    t10r = z[j5]-z[j6];
                    t10i = z[j5+1]-z[j6+1];
                    t11r = z[j00]-0.5*t5r;
                    t11i = z[j00+1]-0.5*t5i;
                    t12r = t1r-t5r;
                    t12i = t1i-t5i;
                    t13r = t2r-t5r;
                    t13i = t2i-t5i;
                    t14r = t3r-t5r;
                    t14i = t3i-t5i;
                    t15r = t4r-t5r;
                    t15i = t4i-t5i;
                    y1r = t11r+c1*t12r+c2*t13r+c3*t14r+c4*t15r;
                    y1i = t11i+c1*t12i+c2*t13i+c3*t14i+c4*t15i;
                    y2r = t11r+c2*t12r+c4*t13r+c5*t14r+c3*t15r;
                    y2i = t11i+c2*t12i+c4*t13i+c5*t14i+c3*t15i;
                    y3r = t11r+c3*t12r+c5*t13r+c2*t14r+c1*t15r;
                    y3i = t11i+c3*t12i+c5*t13i+c2*t14i+c1*t15i;
                    y4r = t11r+c4*t12r+c3*t13r+c1*t14r+c5*t15r;
                    y4i = t11i+c4*t12i+c3*t13i+c1*t14i+c5*t15i;
                    y5r = t11r+c5*t12r+c1*t13r+c4*t14r+c2*t15r;
                    y5i = t11i+c5*t12i+c1*t13i+c4*t14i+c2*t15i;
                    y6r = c10*t6r-c6*t7r+c9*t8r-c7*t9r+c8*t10r;
                    y6i = c10*t6i-c6*t7i+c9*t8i-c7*t9i+c8*t10i;
                    y7r = c9*t6r-c8*t7r+c6*t8r+c10*t9r-c7*t10r;
                    y7i = c9*t6i-c8*t7i+c6*t8i+c10*t9i-c7*t10i;
                    y8r = c8*t6r-c10*t7r-c7*t8r+c6*t9r+c9*t10r;
                    y8i = c8*t6i-c10*t7i-c7*t8i+c6*t9i+c9*t10i;
                    y9r = c7*t6r+c9*t7r-c10*t8r-c8*t9r-c6*t10r;
                    y9i = c7*t6i+c9*t7i-c10*t8i-c8*t9i-c6*t10i;
                    y10r = c6*t6r+c7*t7r+c8*t8r+c9*t9r+c10*t10r;
                    y10i = c6*t6i+c7*t7i+c8*t8i+c9*t9i+c10*t10i;
                    z[j00] = z[j00]+t1r+t2r+t3r+t4r+t5r;
                    z[j00+1] = z[j00+1]+t1i+t2i+t3i+t4i+t5i;
                    z[j01] = y1r-y10i;
                    z[j01+1] = y1i+y10r;
                    z[j2] = y2r-y9i;
                    z[j2+1] = y2i+y9r;
                    z[j3] = y3r-y8i;
                    z[j3+1] = y3i+y8r;
                    z[j4] = y4r-y7i;
                    z[j4+1] = y4i+y7r;
                    z[j5] = y5r-y6i;
                    z[j5+1] = y5i+y6r;
                    z[j6] = y5r+y6i;
                    z[j6+1] = y5i-y6r;
                    z[j7] = y4r+y7i;
                    z[j7+1] = y4i-y7r;
                    z[j8] = y3r+y8i;
                    z[j8+1] = y3i-y8r;
                    z[j9] = y2r+y9i;
                    z[j9+1] = y2i-y9r;
                    z[j10] = y1r+y10i;
                    z[j10+1] = y1i-y10r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                }
                it = i10+istep;
                i10 = i9+istep;
                i9 = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i11 = i10+iinc;
        if (i11>=imax) i11 = i11-imax;
        i12 = i11+iinc;
        if (i12>=imax) i12 = i12-imax;

        /* if factor is 13 */
        if (ifac==13) {
            if (mu==1) {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = -P354;
                c5 = -P748;
                c6 = -P970;
                c7 = P464;
                c8 = P822;
                c9 = P992;
                c10 = P935;
                c11 = P663;
                c12 = P239;
            } else if (mu==2) {
                c1 = P568;
                c2 = -P354;
                c3 = -P970;
                c4 = -P748;
                c5 = P120;
                c6 = P885;
                c7 = P822;
                c8 = P935;
                c9 = P239;
                c10 = -P663;
                c11 = -P992;
                c12 = -P464;
            } else if (mu==3) {
                c1 = P120;
                c2 = -P970;
                c3 = -P354;
                c4 = P885;
                c5 = P568;
                c6 = -P748;
                c7 = P992;
                c8 = P239;
                c9 = -P935;
                c10 = -P464;
                c11 = P822;
                c12 = P663;
            } else if (mu==4) {
                c1 = -P354;
                c2 = -P748;
                c3 = P885;
                c4 = P120;
                c5 = -P970;
                c6 = P568;
                c7 = P935;
                c8 = -P663;
                c9 = -P464;
                c10 = P992;
                c11 = -P239;
                c12 = -P822;
            } else if (mu==5) {
                c1 = -P748;
                c2 = P120;
                c3 = P568;
                c4 = -P970;
                c5 = P885;
                c6 = -P354;
                c7 = P663;
                c8 = -P992;
                c9 = P822;
                c10 = -P239;
                c11 = -P464;
                c12 = P935;
            } else if (mu==6) {
                c1 = -P970;
                c2 = P885;
                c3 = -P748;
                c4 = P568;
                c5 = -P354;
                c6 = P120;
                c7 = P239;
                c8 = -P464;
                c9 = P663;
                c10 = -P822;
                c11 = P935;
                c12 = -P992;
            } else if (mu==7) {
                c1 = -P970;
                c2 = P885;
                c3 = -P748;
                c4 = P568;
                c5 = -P354;
                c6 = P120;
                c7 = -P239;
                c8 = P464;
                c9 = -P663;
                c10 = P822;
                c11 = -P935;
                c12 = P992;
            } else if (mu==8) {
                c1 = -P748;
                c2 = P120;
                c3 = P568;
                c4 = -P970;
                c5 = P885;
                c6 = -P354;
                c7 = -P663;
                c8 = P992;
                c9 = -P822;
                c10 = P239;
                c11 = P464;
                c12 = -P935;
            } else if (mu==9) {
                c1 = -P354;
                c2 = -P748;
                c3 = P885;
                c4 = P120;
                c5 = -P970;
                c6 = P568;
                c7 = -P935;
                c8 = P663;
                c9 = P464;
                c10 = -P992;
                c11 = P239;
                c12 = P822;
            } else if (mu==10) {
                c1 = P120;
                c2 = -P970;
                c3 = -P354;
                c4 = P885;
                c5 = P568;
                c6 = -P748;
                c7 = -P992;
                c8 = -P239;
                c9 = P935;
                c10 = P464;
                c11 = -P822;
                c12 = -P663;
            } else if (mu==11) {
                c1 = P568;
                c2 = -P354;
                c3 = -P970;
                c4 = -P748;
                c5 = P120;
                c6 = P885;
                c7 = -P822;
                c8 = -P935;
                c9 = -P239;
                c10 = P663;
                c11 = P992;
                c12 = P464;
            } else {
                c1 = P885;
                c2 = P568;
                c3 = P120;
                c4 = -P354;
                c5 = -P748;
                c6 = -P970;
                c7 = -P464;
                c8 = -P822;
                c9 = -P992;
                c10 = -P935;
                c11 = -P663;
                c12 = -P239;
            }
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                j11 = i11;
                j12 = i12;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j01]+z[j12];
                    t1i = z[j01+1]+z[j12+1];
                    t2r = z[j2]+z[j11];
                    t2i = z[j2+1]+z[j11+1];
                    t3r = z[j3]+z[j10];
                    t3i = z[j3+1]+z[j10+1];
                    t4r = z[j4]+z[j9];
                    t4i = z[j4+1]+z[j9+1];
                    t5r = z[j5]+z[j8];
                    t5i = z[j5+1]+z[j8+1];
                    t6r = z[j6]+z[j7];
                    t6i = z[j6+1]+z[j7+1];
                    t7r = z[j01]-z[j12];
                    t7i = z[j01+1]-z[j12+1];
                    t8r = z[j2]-z[j11];
                    t8i = z[j2+1]-z[j11+1];
                    t9r = z[j3]-z[j10];
                    t9i = z[j3+1]-z[j10+1];
                    t10r = z[j4]-z[j9];
                    t10i = z[j4+1]-z[j9+1];
                    t11r = z[j5]-z[j8];
                    t11i = z[j5+1]-z[j8+1];
                    t12r = z[j6]-z[j7];
                    t12i = z[j6+1]-z[j7+1];
                    t13r = z[j00]-0.5*t6r;
                    t13i = z[j00+1]-0.5*t6i;
                    t14r = t1r-t6r;
                    t14i = t1i-t6i;
                    t15r = t2r-t6r;
                    t15i = t2i-t6i;
                    t16r = t3r-t6r;
                    t16i = t3i-t6i;
                    t17r = t4r-t6r;
                    t17i = t4i-t6i;
                    t18r = t5r-t6r;
                    t18i = t5i-t6i;
                    y1r = t13r+c1*t14r+c2*t15r+c3*t16r+c4*t17r+c5*t18r;
                    y1i = t13i+c1*t14i+c2*t15i+c3*t16i+c4*t17i+c5*t18i;
                    y2r = t13r+c2*t14r+c4*t15r+c6*t16r+c5*t17r+c3*t18r;
                    y2i = t13i+c2*t14i+c4*t15i+c6*t16i+c5*t17i+c3*t18i;
                    y3r = t13r+c3*t14r+c6*t15r+c4*t16r+c1*t17r+c2*t18r;
                    y3i = t13i+c3*t14i+c6*t15i+c4*t16i+c1*t17i+c2*t18i;
                    y4r = t13r+c4*t14r+c5*t15r+c1*t16r+c3*t17r+c6*t18r;
                    y4i = t13i+c4*t14i+c5*t15i+c1*t16i+c3*t17i+c6*t18i;
                    y5r = t13r+c5*t14r+c3*t15r+c2*t16r+c6*t17r+c1*t18r;
                    y5i = t13i+c5*t14i+c3*t15i+c2*t16i+c6*t17i+c1*t18i;
                    y6r = t13r+c6*t14r+c1*t15r+c5*t16r+c2*t17r+c4*t18r;
                    y6i = t13i+c6*t14i+c1*t15i+c5*t16i+c2*t17i+c4*t18i;
                    y7r = c12*t7r-c7*t8r+c11*t9r-c8*t10r+c10*t11r-c9*t12r;
                    y7i = c12*t7i-c7*t8i+c11*t9i-c8*t10i+c10*t11i-c9*t12i;
                    y8r = c11*t7r-c9*t8r+c8*t9r-c12*t10r-c7*t11r+c10*t12r;
                    y8i = c11*t7i-c9*t8i+c8*t9i-c12*t10i-c7*t11i+c10*t12i;
                    y9r = c10*t7r-c11*t8r-c7*t9r+c9*t10r-c12*t11r-c8*t12r;
                    y9i = c10*t7i-c11*t8i-c7*t9i+c9*t10i-c12*t11i-c8*t12i;
                    y10r = c9*t7r+c12*t8r-c10*t9r-c7*t10r+c8*t11r+c11*t12r;
                    y10i = c9*t7i+c12*t8i-c10*t9i-c7*t10i+c8*t11i+c11*t12i;
                    y11r = c8*t7r+c10*t8r+c12*t9r-c11*t10r-c9*t11r-c7*t12r;
                    y11i = c8*t7i+c10*t8i+c12*t9i-c11*t10i-c9*t11i-c7*t12i;
                    y12r = c7*t7r+c8*t8r+c9*t9r+c10*t10r+c11*t11r+c12*t12r;
                    y12i = c7*t7i+c8*t8i+c9*t9i+c10*t10i+c11*t11i+c12*t12i;
                    z[j00] = z[j00]+t1r+t2r+t3r+t4r+t5r+t6r;
                    z[j00+1] = z[j00+1]+t1i+t2i+t3i+t4i+t5i+t6i;
                    z[j01] = y1r-y12i;
                    z[j01+1] = y1i+y12r;
                    z[j2] = y2r-y11i;
                    z[j2+1] = y2i+y11r;
                    z[j3] = y3r-y10i;
                    z[j3+1] = y3i+y10r;
                    z[j4] = y4r-y9i;
                    z[j4+1] = y4i+y9r;
                    z[j5] = y5r-y8i;
                    z[j5+1] = y5i+y8r;
                    z[j6] = y6r-y7i;
                    z[j6+1] = y6i+y7r;
                    z[j7] = y6r+y7i;
                    z[j7+1] = y6i-y7r;
                    z[j8] = y5r+y8i;
                    z[j8+1] = y5i-y8r;
                    z[j9] = y4r+y9i;
                    z[j9+1] = y4i-y9r;
                    z[j10] = y3r+y10i;
                    z[j10+1] = y3i-y10r;
                    z[j11] = y2r+y11i;
                    z[j11+1] = y2i-y11r;
                    z[j12] = y1r+y12i;
                    z[j12+1] = y1i-y12r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                    j11 += jstep;
                    j12 += jstep;
                }
                it = i12+istep;
                i12 = i11+istep;
                i11 = i10+istep;
                i10 = i9+istep;
                i9 = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
        i13 = i12+iinc;
        if (i13>=imax) i13 = i13-imax;
        i14 = i13+iinc;
        if (i14>=imax) i14 = i14-imax;
        i15 = i14+iinc;
        if (i15>=imax) i15 = i15-imax;

        /* if factor is 16 */
        if (ifac==16) {
            if (mu==1) {
                c1 = 1.0;
                c2 = P923;
                c3 = P382;
                c4 = P707;
            } else if (mu==3) {
                c1 = -1.0;
                c2 = P382;
                c3 = P923;
                c4 = -P707;
            } else if (mu==5) {
                c1 = 1.0;
                c2 = -P382;
                c3 = P923;
                c4 = -P707;
            } else if (mu==7) {
                c1 = -1.0;
                c2 = -P923;
                c3 = P382;
                c4 = P707;
            } else if (mu==9) {
                c1 = 1.0;
                c2 = -P923;
                c3 = -P382;
                c4 = P707;
            } else if (mu==11) {
                c1 = -1.0;
                c2 = -P382;
                c3 = -P923;
                c4 = -P707;
            } else if (mu==13) {
                c1 = 1.0;
                c2 = P382;
                c3 = -P923;
                c4 = -P707;
            } else {
                c1 = -1.0;
                c2 = P923;
                c3 = -P382;
                c4 = P707;
            }
            c5 = c1*c4;
            c6 = c1*c3;
            c7 = c1*c2;
            for (l=0; l<m; l++) {
                j00 = i0;
                j01 = i1;
                j2 = i2;
                j3 = i3;
                j4 = i4;
                j5 = i5;
                j6 = i6;
                j7 = i7;
                j8 = i8;
                j9 = i9;
                j10 = i10;
                j11 = i11;
                j12 = i12;
                j13 = i13;
                j14 = i14;
                j15 = i15;
                for (jt=0; jt<nt; jt++) {
                    t1r = z[j00]+z[j8];
                    t1i = z[j00+1]+z[j8+1];
                    t2r = z[j4]+z[j12];
                    t2i = z[j4+1]+z[j12+1];
                    t3r = z[j00]-z[j8];
                    t3i = z[j00+1]-z[j8+1];
                    t4r = c1*(z[j4]-z[j12]);
                    t4i = c1*(z[j4+1]-z[j12+1]);
                    t5r = t1r+t2r;
                    t5i = t1i+t2i;
                    t6r = t1r-t2r;
                    t6i = t1i-t2i;
                    t7r = z[j01]+z[j9];
                    t7i = z[j01+1]+z[j9+1];
                    t8r = z[j5]+z[j13];
                    t8i = z[j5+1]+z[j13+1];
                    t9r = z[j01]-z[j9];
                    t9i = z[j01+1]-z[j9+1];
                    t10r = z[j5]-z[j13];
                    t10i = z[j5+1]-z[j13+1];
                    t11r = t7r+t8r;
                    t11i = t7i+t8i;
                    t12r = t7r-t8r;
                    t12i = t7i-t8i;
                    t13r = z[j2]+z[j10];
                    t13i = z[j2+1]+z[j10+1];
                    t14r = z[j6]+z[j14];
                    t14i = z[j6+1]+z[j14+1];
                    t15r = z[j2]-z[j10];
                    t15i = z[j2+1]-z[j10+1];
                    t16r = z[j6]-z[j14];
                    t16i = z[j6+1]-z[j14+1];
                    t17r = t13r+t14r;
                    t17i = t13i+t14i;
                    t18r = c4*(t15r-t16r);
                    t18i = c4*(t15i-t16i);
                    t19r = c5*(t15r+t16r);
                    t19i = c5*(t15i+t16i);
                    t20r = c1*(t13r-t14r);
                    t20i = c1*(t13i-t14i);
                    t21r = z[j3]+z[j11];
                    t21i = z[j3+1]+z[j11+1];
                    t22r = z[j7]+z[j15];
                    t22i = z[j7+1]+z[j15+1];
                    t23r = z[j3]-z[j11];
                    t23i = z[j3+1]-z[j11+1];
                    t24r = z[j7]-z[j15];
                    t24i = z[j7+1]-z[j15+1];
                    t25r = t21r+t22r;
                    t25i = t21i+t22i;
                    t26r = t21r-t22r;
                    t26i = t21i-t22i;
                    t27r = t9r+t24r;
                    t27i = t9i+t24i;
                    t28r = t10r+t23r;
                    t28i = t10i+t23i;
                    t29r = t9r-t24r;
                    t29i = t9i-t24i;
                    t30r = t10r-t23r;
                    t30i = t10i-t23i;
                    t31r = t5r+t17r;
                    t31i = t5i+t17i;
                    t32r = t11r+t25r;
                    t32i = t11i+t25i;
                    t33r = t3r+t18r;
                    t33i = t3i+t18i;
                    t34r = c2*t29r-c6*t30r;
                    t34i = c2*t29i-c6*t30i;
                    t35r = t3r-t18r;
                    t35i = t3i-t18i;
                    t36r = c7*t27r-c3*t28r;
                    t36i = c7*t27i-c3*t28i;
                    t37r = t4r+t19r;
                    t37i = t4i+t19i;
                    t38r = c3*t27r+c7*t28r;
                    t38i = c3*t27i+c7*t28i;
                    t39r = t4r-t19r;
                    t39i = t4i-t19i;
                    t40r = c6*t29r+c2*t30r;
                    t40i = c6*t29i+c2*t30i;
                    t41r = c4*(t12r-t26r);
                    t41i = c4*(t12i-t26i);
                    t42r = c5*(t12r+t26r);
                    t42i = c5*(t12i+t26i);
                    y1r = t33r+t34r;
                    y1i = t33i+t34i;
                    y2r = t6r+t41r;
                    y2i = t6i+t41i;
                    y3r = t35r+t40r;
                    y3i = t35i+t40i;
                    y4r = t5r-t17r;
                    y4i = t5i-t17i;
                    y5r = t35r-t40r;
                    y5i = t35i-t40i;
                    y6r = t6r-t41r;
                    y6i = t6i-t41i;
                    y7r = t33r-t34r;
                    y7i = t33i-t34i;
                    y9r = t38r-t37r;
                    y9i = t38i-t37i;
                    y10r = t42r-t20r;
                    y10i = t42i-t20i;
                    y11r = t36r+t39r;
                    y11i = t36i+t39i;
                    y12r = c1*(t11r-t25r);
                    y12i = c1*(t11i-t25i);
                    y13r = t36r-t39r;
                    y13i = t36i-t39i;
                    y14r = t42r+t20r;
                    y14i = t42i+t20i;
                    y15r = t38r+t37r;
                    y15i = t38i+t37i;
                    z[j00] = t31r+t32r;
                    z[j00+1] = t31i+t32i;
                    z[j01] = y1r-y15i;
                    z[j01+1] = y1i+y15r;
                    z[j2] = y2r-y14i;
                    z[j2+1] = y2i+y14r;
                    z[j3] = y3r-y13i;
                    z[j3+1] = y3i+y13r;
                    z[j4] = y4r-y12i;
                    z[j4+1] = y4i+y12r;
                    z[j5] = y5r-y11i;
                    z[j5+1] = y5i+y11r;
                    z[j6] = y6r-y10i;
                    z[j6+1] = y6i+y10r;
                    z[j7] = y7r-y9i;
                    z[j7+1] = y7i+y9r;
                    z[j8] = t31r-t32r;
                    z[j8+1] = t31i-t32i;
                    z[j9] = y7r+y9i;
                    z[j9+1] = y7i-y9r;
                    z[j10] = y6r+y10i;
                    z[j10+1] = y6i-y10r;
                    z[j11] = y5r+y11i;
                    z[j11+1] = y5i-y11r;
                    z[j12] = y4r+y12i;
                    z[j12+1] = y4i-y12r;
                    z[j13] = y3r+y13i;
                    z[j13+1] = y3i-y13r;
                    z[j14] = y2r+y14i;
                    z[j14+1] = y2i-y14r;
                    z[j15] = y1r+y15i;
                    z[j15+1] = y1i-y15r;
                    j00 += jstep;
                    j01 += jstep;
                    j2 += jstep;
                    j3 += jstep;
                    j4 += jstep;
                    j5 += jstep;
                    j6 += jstep;
                    j7 += jstep;
                    j8 += jstep;
                    j9 += jstep;
                    j10 += jstep;
                    j11 += jstep;
                    j12 += jstep;
                    j13 += jstep;
                    j14 += jstep;
                    j15 += jstep;
                }
                it = i15+istep;
                i15 = i14+istep;
                i14 = i13+istep;
                i13 = i12+istep;
                i12 = i11+istep;
                i11 = i10+istep;
                i10 = i9+istep;
                i9 = i8+istep;
                i8 = i7+istep;
                i7 = i6+istep;
                i6 = i5+istep;
                i5 = i4+istep;
                i4 = i3+istep;
                i3 = i2+istep;
                i2 = i1+istep;
                i1 = i0+istep;
                i0 = it;
            }
            continue;
        }
    }
}


void sf_pfa2cc (int isign, int idim, int n1, int n2, float complex cz[])
/*****************************************************************************
Prime factor fft:  2-D complex to complex transforms, in place
******************************************************************************
Input:
isign       	sign of isign is the sign of exponent in fourier kernel
idim        	dimension to transform, either 1 or 2 (see notes)
n1          	1st (fast) dimension of array to be transformed (see notes)
n2          	2nd (slow) dimension of array to be transformed (see notes)
z           	array[n2][n1] of complex elements to be transformed in place

Output:
z		array[n2][n1] of complex elements transformed
******************************************************************************
Notes:
Only one (either the 1st or 2nd) dimension of the 2-D array is transformed.

If idim equals 1, then n2 transforms of n1 complex elements are performed; 
else, if idim equals 2, then n1 transforms of n2 complex elements are 
performed.

Although z appears in the argument list as a one-dimensional array,
z may be viewed as an n1 by n2 two-dimensional array:  z[n2][n1].

Let n denote the transform length, either n1 or n2, depending on idim.
Then, n must be factorable into mutually prime factors taken 
from the set {2,3,4,5,7,8,9,11,13,16}.  in other words,
    n = 2**p * 3**q * 5**r * 7**s * 11**t * 13**u
where
    0 <= p <= 4,  0 <= q <= 2,  0 <= r,s,t,u <= 1
is required for pfa2cc to yield meaningful results.  this
restriction implies that n is restricted to the range
    1 <= n <= 720720 (= 5*7*9*11*13*16)

To perform a two-dimensional transform of an n1 by n2 complex array 
(assuming that both n1 and n2 are valid "n"), stored with n1 fast 
and n2 slow:  pfa2cc(isign,1,n1,n2,z);  pfa2cc(isign,2,n1,n2,z);
******************************************************************************
References:  
Temperton, C., 1985, Implementation of a self-sorting
in-place prime factor fft algorithm:  Journal of
Computational Physics, v. 58, p. 283-299.

Temperton, C., 1988, A new set of minimum-add rotated
rotated dft modules: Journal of Computational Physics,
v. 75, p. 190-198.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/15/89
*****************************************************************************/
{
    int n,nt,k,kt;

    /* determine transform length, number of transforms, and strides */
    if (idim==1) {
        n = n1;
        nt = n2;
        k = 1;
        kt = n1;
    } else {
        n = n2;
        nt = n1;
        k = n1;
        kt = 1;
    }

    /* do multiple complex to complex transforms */
    sf_pfamcc(isign,n,nt,k,kt,cz);
}

void sf_pfa2cr (int isign, int idim, int n1, int n2, 
		float complex cz[], float rz[])
/*****************************************************************************
Prime factor fft:  2-D complex to real transforms
******************************************************************************
Input:
isign       sign of isign is the sign of exponent in fourier kernel
idim        dimension to transform, which must be either 1 or 2 (see notes)
n1          1st (fast) dimension of array to be transformed (see notes)
n2          2nd (slow) dimension of array to be transformed (see notes)
cz          array of complex values (may be equivalenced to rz)

Output:
rz          array of real values (may be equivalenced to cz)
******************************************************************************
Notes:
If idim equals 1, then n2 transforms of n1/2+1 complex elements to n1 real 
elements are performed; else, if idim equals 2, then n1 transforms of n2/2+1 
complex elements to n2 real elements are performed.

Although rz appears in the argument list as a one-dimensional array,
rz may be viewed as an n1 by n2 two-dimensional array:  rz[n2][n1].  
Likewise, depending on idim, cz may be viewed as either an n1/2+1 by 
n2 or an n1 by n2/2+1 two-dimensional array of complex elements.

Let n denote the transform length, either n1 or n2, depending on idim.
Because pfa2rc uses pfa2cc to do most of the work, n must be even 
and n/2 must be a valid length for pfa2cc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.
******************************************************************************
References:  
Press et al, 1988, Numerical Recipes in C, p. 417.

Also, see notes and references for function pfa2cc.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/13/89
Modified:  Dave Hale, Colorado School of Mines, 11/26/89
	Fixed bug that gave wrong answers for idim==2 and cz not 
	equivalenced to rz.  The bug was in copying the input
	to the output in the first section below - the j index
	was not being incremented along with i1.  Fixed by adding
	j++ after the i1++.
*****************************************************************************/
{
    int i1,i2,j,k,it,jt,kt,n,nt,itmul,itinc;
    float *z,*temp,tempr,tempi,sumr,sumi,difr,difi;
    double wr,wi,wpr,wpi,wtemp,theta;

    /* if transforming dimension 1 */
    if (idim==1) {

        /* copy input to output and fix dc and nyquist */
        z = (float*)cz;
        for (i2=0,jt=0,kt=0; i2<n2; i2++,jt+=n1,kt+=(n1+2)) {
            rz[jt+1] = z[kt]-z[kt+n1];
            rz[jt] = z[kt]+z[kt+n1];
            for (i1=2,j=jt+2,k=kt+2; i1<n1; i1++,j++,k++)
                rz[j] = z[k];
        }
        z = rz;

        /* set transform length, number of transforms and strides */
        n = n1;
        nt = n2;
        itmul = 1;
        itinc = n1;

    /* else, if transforming dimension 2 */
    } else {

        /* copy input to output and fix dc and nyquist */
        z = (float*)cz;
        for (i2=1; i2<n2/2; i2++) {
            for (i1=0,j=i2*n1*2; i1<n1*2; i1++,j++)
                rz[j] = z[j];
        }
        for (i1=0,j=n1*n2; i1<n1*2; i1+=2,j+=2) {
            rz[i1+1] = z[i1]-z[j];
            rz[i1] = z[i1]+z[j];
        }
        z = rz;

        /* set transform length, number of transforms and strides */
        n = n2;
        nt = n1;
        itmul = n1;
        itinc = 2;
    }

    /* initialize cosine-sine recurrence */
    theta = 2.0*SF_PI/(double)n;
    if (isign>0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle transforms simultaneously */
    for (j=2,k=n-2; j<=n/2; j+=2,k-=2) {
        jt = j*itmul;
        kt = k*itmul;
        for (it=0; it<nt; it++) {
            sumr = z[jt]+z[kt];
            sumi = z[jt+1]+z[kt+1];
            difr = z[jt]-z[kt];
            difi = z[jt+1]-z[kt+1];
            tempr = wi*difr-wr*sumi;
            tempi = wi*sumi+wr*difr;
            z[jt] = sumr+tempr;
            z[jt+1] = difi+tempi;
            z[kt] = sumr-tempr;
            z[kt+1] = tempi-difi;
            jt += itinc;
            kt += itinc;
        }
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }

    /* if transforming dimension 1 */
    if (idim==1) {

        /* transform as complex elements */
        sf_pfa2cc(isign,1,n1/2,n2,(float complex*)z);

    /* else, if transforming dimension 2 */
    } else {

        /* transform as complex elements */
        sf_pfa2cc(isign,2,n1,n2/2,(float complex*)z);

        /* unmerge even and odd vectors */
        temp = (float*)malloc(n1*sizeof(float));
        for (i2=0; i2<n2; i2+=2) {
            for (i1=0,j=i2*n1+1; i1<n1; i1++,j+=2)
                temp[i1] = z[j];
            for (i1=0,j=i2*n1,k=i2*n1; i1<n1; i1++,j+=2,k++)
                z[k] = z[j];
            for (i1=0,j=(i2+1)*n1; i1<n1; i1++,j++)
                z[j] = temp[i1];
        }
        free(temp);
    }
}

void sf_pfa2rc (int isign, int idim, int n1, int n2, 
		float rz[], float complex cz[])
/*****************************************************************************
Prime factor fft:  2-D real to complex transforms
******************************************************************************
Input:
isign       sign of isign is the sign of exponent in fourier kernel
idim        dimension to transform, which must be either 1 or 2 (see notes)
n1          1st (fast) dimension of array to be transformed (see notes)
n2          2nd (slow) dimension of array to be transformed (see notes)
rz          array of real values (may be equivalenced to cz)

Output:
cz          array of complex values (may be equivalenced to rz)
******************************************************************************
Notes:
If idim equals 1, then n2 transforms of n1 real elements to n1/2+1 complex 
elements are performed; else, if idim equals 2, then n1 transforms of n2 
real elements to n2/2+1 complex elements are performed.

Although rz appears in the argument list as a one-dimensional array,
rz may be viewed as an n1 by n2 two-dimensional array:  rz[n2][n1].  
Likewise, depending on idim, cz may be viewed as either an n1/2+1 by 
n2 or an n1 by n2/2+1 two-dimensional array of complex elements.

Let n denote the transform length, either n1 or n2, depending on idim.
Because pfa2rc uses pfa2cc to do most of the work, n must be even 
and n/2 must be a valid length for pfa2cc.  The simplest way to
obtain a valid n is via n = npfar(nmin).  A more optimal n can be 
obtained with npfaro.
******************************************************************************
References:  
Press et al, 1988, Numerical Recipes in C, p. 417.

Also, see notes and references for function pfa2cc.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/13/89
*****************************************************************************/
{
    int i1,i2,j,k,it,jt,kt,n,nt,itmul,itinc;
    float *z,*temp,tempr,tempi,sumr,sumi,difr,difi;
    double wr,wi,wpr,wpi,wtemp,theta;

    /* copy input to output while scaling */
    z = (float*)cz;
    for (i2=0,j=0; i2<n2; i2++)
        for (i1=0; i1<n1; i1++,j++)
            z[j] = 0.5*rz[j];

    /* if transforming dimension 1 */
    if (idim==1) {

        /* transform as complex elements */
        sf_pfa2cc(isign,1,n1/2,n2,cz);

        /* shift rows to make room for nyquist */
        z = (float*)cz;
        for (i2=n2-1; i2>0; i2--) {
            jt = i2*n1+n1-1;
            kt = jt+i2*2;
            for (i1=n1-1,j=jt,k=kt; i1>=0; i1--,j--,k--)
                z[k] = z[j];
        }

        /* set transform length, number of transforms and strides */
        n = n1;
        nt = n2;
        itmul = 1;
        itinc = n1+2;

    /* else, if transforming dimension 2 */
    } else {

        /* merge even and odd vectors */
        temp = z+n1*n2;
        for (i2=0; i2<n2; i2+=2) {
            for (i1=0,j=i2*n1; i1<n1; i1++,j++)
                temp[i1] = z[j];
            for (i1=0,j=(i2+1)*n1,k=i2*n1+1; i1<n1; i1++,j++,k+=2)
                z[k] = z[j];
            for (i1=0,j=i2*n1; i1<n1; i1++,j+=2)
                z[j] = temp[i1];
        }

        /* transform as complex elements */
        sf_pfa2cc(isign,2,n1,n2/2,cz);

        /* set transform length, number of transforms and strides */
        n = n2;
        nt = n1;
        itmul = n1;
        itinc = 2;
    }

    /* fix dc and nyquist for each transform */
    for (it=0,j=0,k=n*itmul; it<nt; it++,j+=itinc,k+=itinc) {
        z[k] = 2.0*(z[j]-z[j+1]);
        z[j] = 2.0*(z[j]+z[j+1]);
        z[k+1] = 0.0;
        z[j+1] = 0.0;
    }

    /* initialize cosine-sine recurrence */
    theta = 2.0*SF_PI/(double)n;
    if (isign<0) theta = -theta;
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0+wpr;
    wi = wpi;

    /* twiddle transforms simultaneously */
    for (j=2,k=n-2; j<=n/2; j+=2,k-=2) {
        jt = j*itmul;
        kt = k*itmul;
        for (it=0; it<nt; it++) {
            sumr = z[jt]+z[kt];
            sumi = z[jt+1]+z[kt+1];
            difr = z[jt]-z[kt];
            difi = z[jt+1]-z[kt+1];
            tempr = wi*difr+wr*sumi;
            tempi = wi*sumi-wr*difr;
            z[jt] = sumr+tempr;
            z[jt+1] = difi+tempi;
            z[kt] = sumr-tempr;
            z[kt+1] = tempi-difi;
            jt += itinc;
            kt += itinc;
        }
        wtemp = wr;
        wr += wr*wpr-wi*wpi;
        wi += wi*wpr+wtemp*wpi;
    }
}

/* 	$Id: pfafft.c,v 1.2 2003/09/29 14:34:55 fomels Exp $	 */

