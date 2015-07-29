/* Make traces with reverberations for testing deconvolution. */
/*
  Copyright (c) 2007, Colorado School of Mines,
  All rights reserved.
  
  
  Redistribution and use in source and binary forms, with or 
  without modification, are permitted provided that the following 
  conditions are met:
  
  *  Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
  *  Redistributions in binary form must reproduce the above 
  copyright notice, this list of conditions and the following 
  disclaimer in the documentation and/or other materials provided 
  with the distribution.
  *  Neither the name of the Colorado School of Mines nor the names of
  its contributors may be used to endorse or promote products 
  derived from this software without specific prior written permission.
  
  Warranty Disclaimer:
  THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
  COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
  POSSIBILITY OF SUCH DAMAGE.
  
  
  Export Restriction Disclaimer:
  We believe that CWP/SU: Seismic Un*x is a low technology product that does
  not appear on the Department of Commerce CCL list of restricted exports.
  Accordingly, we believe that our product meets the qualifications of
  an ECCN (export control classification number) of EAR99 and we believe
  it fits the qualifications of NRR (no restrictions required), and
  is thus not subject to export restrictions of any variety.
  
  Approved Reference Format:
  In publications, please refer to SU as per the following example:
  Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
  Release No. __: an open source software  package for seismic 
  research and processing, 
  Center for Wave Phenomena, Colorado School of Mines.
  
  Articles about SU in peer-reviewed journals:
  Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
  Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
  Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
  Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
  Acknowledgements:
  SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
  School of Mines, partially based on Stanford Exploration Project (SEP) 
  software.
*/

/* Make traces with 2 way reverberations--suggested by Ken */
/* Modified by Sergey Fomel for inclusion with Madagascar. */

#include <rsf.h>

#define ISODD(n) ((n) & 01)

int main(int argc, char* argv[])
{
    float x;	    /* trace value       */
    int ns=512;	    /* number of samples */
    int ntr=8;	    /* number of traces  */
    int dtime=16;   /* one way time in samples through "ocean" layer */
    float rbot=0.8; /* reflection strength of "ocean" bottom         */
    int h=100;	    /* location in samples of two way reverb train 	 */
    float amp=0.2;  /* strength of reflector */
    int loc=170;    /* location of reflector on trace 1 in time samples */
    int dip=12;	    /* dip of reflector in time samples */
    int i, j, k, sgn;
    float *trace;
    sf_file traces;

    sf_init(argc,argv);
    traces = sf_output("out");
    sf_setformat(traces,"native_float");
    sf_putint(traces,"n1",ns);
    sf_putint(traces,"n2",ntr);

    sf_putfloat(traces,"o1",0);
    sf_putfloat(traces,"d1",0.004);
    sf_putfloat(traces,"o2",1);
    sf_putfloat(traces,"d2",1);
    
    sf_putstring(traces,"label1","Time");
    sf_putstring(traces,"unit1","s");
    sf_putstring(traces,"label2","Trace");

    trace = sf_floatalloc(ns);
    
    for (j = 0; j < ntr; j++) {
	for (i = 0; i < ns; i++) {
	    if (i >= h && ((i-h) % dtime == 0)) {
		k = (i-h)/dtime;
		sgn = (ISODD(k) ? -1 : 1);
		x = sgn * (k+1) * powf(rbot, k);
	    } else {
		x = 0.0f;
	    }
	    if (i == loc + j*dip) x += amp;
	    trace[i] = x;
	}
	sf_floatwrite(trace,ns,traces);
    } 

    exit(0);
}
