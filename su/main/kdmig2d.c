/* 2-D Prestack Kirchhoff depth migration (SU version). */
/*
  Copyright � 2010, Colorado School of Mines,
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
  Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1),
                     data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
  Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
  Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
  Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
  Acknowledgements:
  SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
  School of Mines, partially based on Stanford Exploration Project (SEP) 
  software.
*/
/*
 * Author:  Zhenyue Liu, 03/01/95,  Colorado School of Mines 
 * Modifcations:
 *    Gary Billings, Talisman Energy, Sept 2005:  added absoff, limoff.
 */
/* Modified by Vladimir Bashkardin for inclusion with Madagascar. */

#include <math.h>

#include <rsf.h>

/*
 Notes:
 1. Traveltime tables exist on coarse grids, with dimension ns*nxt*nzt.
    In the migration process, traveltimes are interpolated into shot/gephone
    positions and output grids.
 2. Input seismic traces can be any type of gathers (common shot, common offset,
    common CDP, and so on). 
 3. Migrated traces are output in CDP gathers if velocity analysis
    is required, with dimension nxo*noff*nzo. 
 4. If the offset value of an input trace is not in the offset array
    of output, the nearest one in the array is chosen.
 5. Memory requirement for this program is about
       [ns*nxt*nzt+noff*nxo*nzo+4*nr*nzt+5*nxt*nzt+npa*(2*ns*nxt*nzt
        +noff*nxo*nzo+4*nxt*nzt)]*4 bytes
    where nr = 1+min(nxt*dxt,0.5*offmax+aperx)/dxo.
 6. Amplitudes are computed using the reference velocity profile, v(z),
    specified by the parameters v0= and dvz=.
 8. if limoff=0, input traces from outside the range defined by off0, doff,
    noff, will get migrated into the extremal offset bins/planes.  E.g. if
    absoff=0 and limoff=0, all traces with gx<sx will get migrated into the
    off0 bin.
*/

void resit (int nx, float fx, float dx, int nz, int nr, float dr,
            float **tb, float **t, float x0);
void interpx (int nxt, float fxt, float dxt, int nx, float fx, float dx,
              int nzt, float **tt, float **t);
void sum2 (int nx, int nz, float a1, float a2, float **t1, float **t2, float **t);
void timeb (int nr, int nz, float dr, float dz, float fz, float a,
            float v0, float **t, float **p, float **sig, float **ang);
void mig2d (float *trace, int nt, float ft, float dt,
            float sx, float gx, float **mig, float aperx,
            int nx, float fx, float dx, float nz, float fz, float dz,
            int ls, int mtmax, float dxm, float fmax, float angmax,
            float **tb, float **pb, float **cs0b, float **angb, int nr,
            float **tsum, int nzt, float fzt, float dzt, int nxt, float fxt, float dxt,
            int npv, float **cssum, float **tvsum, float **mig1);

static kiss_fftr_cfg forw = NULL, invs = NULL;

/* define */
#define RSCALE_KDMIG 1000.0
#define NINT(x) ((int)roundf(x))

int main (int argc, char **argv) {
    int nt;   /* number of time samples in input data          */
    int nzt;  /* number of z-values in traveltime table        */
    int nxt;  /* number of x-values in traveltime table        */
    int nzo;  /* number of z-values in output data             */
    int nxo;  /* number of x-values in output data             */
    int ns;   /* number of sources                             */
    int noff; /* number of offsets in output data              */
    int nr, nsi, ngi;
    int is, io, iis, iig, ixo, izo; /* counters */
    int ls, ntr, jtr, ktr, mtr, npv, mtmax;
    int absoff, limoff;
    float ft, fzt, fxt, fzo, fxo, fs, off0, offset,
          dt, dzt, dxt, dzo, dxo, ds, doff, dxm,
          ext, ezt, ezo, exo, es, s, scal;
    float fsi, fgi, dsi, dgi;
    float v0, dvz, fmax, angmax, offmax, rmax, aperx, sx = 0, gx = 0;
    float ***mig = NULL, ***ttab = NULL, **tb = NULL, **pb = NULL, *data = NULL;
    float **cs0b = NULL, **angb = NULL, **tsum = NULL, **tt = NULL;
    float **tvsum = NULL, ***mig1 = NULL, ***cs = NULL, ***tv = NULL, **cssum = NULL;
    float rscale; /* scaling factor for roundoff */
    char *ttfile, *tvfile, *csfile, *outfile1, *unit;
    sf_file infp, outfp, ttfp, tvfp = NULL, out1fp = NULL,  csfp = NULL;
    bool verb;
    float as,res;

    sf_init (argc, argv);

    infp = sf_input ("in");
    outfp = sf_output ("out");

    ttfile = sf_getstring ("ttfile");
    /* input traveltime tables */
    if (NULL == ttfile) sf_error ("Need ttfile=");
    ttfp = sf_input (ttfile);
    /* input traveltime tables */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (!sf_histint (infp, "n1", &nt)) sf_error ("No n1= in input");
    /* number of time samples per trace of input data */
    if (!sf_histfloat (infp, "o1", &ft)) sf_error ("No o1= in input");
    /* first time sample of input data */
    if (!sf_histfloat (infp, "d1", &dt)) sf_error ("No d1= in input");
    /* time sampling interval of input data */

    /* Determine shot/receiver geometry from header parameters;
       add irregular geometry support later (see ktmig.c for example) */
    if (!sf_histint (infp, "n2", &ngi)) sf_error ("No n2= in input");
    if (!sf_histfloat (infp, "o2", &fgi)) sf_error ("No o2= in input");
    if (!sf_histfloat (infp, "d2", &dgi)) sf_error ("No d2= in input");
    if (!sf_histint (infp, "n3", &nsi)) sf_error ("No n3= in input");
    if (!sf_histfloat (infp, "o3", &fsi)) sf_error ("No o3= in input");
    if (!sf_histfloat (infp, "d3", &dsi)) sf_error ("No d3= in input");

    if (!sf_histint (ttfp, "n1", &nzt)) sf_error ("No n1= in ttfile=");
    /* number of depth samples in traveltime table */
    if (!sf_histfloat (ttfp, "o1", &fzt)) sf_error ("No o1= in ttfile=");
    /* first depth sample in traveltime table */
    if (!sf_histfloat (ttfp, "d1", &dzt)) sf_error ("No d1= in ttfile=");
    /* depth interval in traveltime table */
    if (!sf_histint (ttfp, "n2", &nxt)) sf_error ("No n2= in ttfile=");
    /* number of lateral samples in traveltime table */
    if (!sf_histfloat (ttfp, "o2", &fxt)) sf_error ("No o2= in ttfile=");
    /* first lateral sample in traveltime table */
    if (!sf_histfloat (ttfp, "d2", &dxt)) sf_error ("No d2= in ttfile=");
    /* lateral interval in traveltime table */
    if (!sf_histint (ttfp, "n3", &ns)) sf_error ("No n3= in ttfile=");
    /* number of sources */
    if (!sf_histfloat (ttfp, "o3", &fs)) sf_error ("No o3= in ttfile=");
    /* x-coordinate of first source */
    if (!sf_histfloat (ttfp, "d3", &ds)) sf_error ("No d3= in ttfile=");
    /* x-coordinate increment of sources */

    if (!sf_getfloat ("dxm", &dxm)) dxm = 0.5*ds;
    /* sampling interval of midpoints */

    if (!sf_getfloat ("rscale", &rscale)) rscale = RSCALE_KDMIG;
    /* scaling for roundoff error suppression */

    ext = NINT(rscale*(fxt+(nxt-1)*dxt))/rscale;
    ezt = NINT(rscale*(fzt+(nzt-1)*dzt))/rscale;
    es = NINT(rscale*(fs+(ns-1)*ds))/rscale;

    if (!sf_getint ("nxo", &nxo)) nxo = (nxt-1)*2+1;
    /* number of output traces */
    if (!sf_getfloat ("fxo", &fxo)) fxo = fxt;
    /* x-coordinate of first output trace */
    if (!sf_getfloat ("dxo", &dxo)) dxo = dxt*0.5;
    /* horizontal spacing of output trace */
    if (!sf_getint ("nzo", &nzo)) nzo = (nzt-1)*5+1;
    /* number of points in output trace */
    if (!sf_getfloat ("fzo", &fzo)) fzo = fzt;
    /* z-coordinate of first point in output trace */
    if (!sf_getfloat ("dzo", &dzo)) dzo = dzt*0.2;
    /* vertical spacing of output trace */

    exo = NINT(rscale*(fxo+(nxo-1)*dxo))/rscale;
    ezo = NINT(rscale*(fzo+(nzo-1)*dzo))/rscale;

    if (verb) {
        sf_warning (" fxt=%f fxo=%f", fxt, fxo);
        sf_warning (" ext=%f exo=%f", ext, exo);
        sf_warning (" fzt=%f fzo=%f", fzt, fzo);
        sf_warning (" ezt=%f ezo=%f", ezt, ezo);
        sf_warning (" es=%f \n", es);
    }
    if (fxt > fxo || ext < exo || fzt > fzo || ezt < ezo) {
        sf_warning ("This condition must NOT be satisfied: fxt>fxo || ext<exo || fzt>fzo || ezt<ezo");
        sf_warning ("fxt=%f fxo=%f ext=%f exo=%f fzt=%f fzo=%f ezt=%f ezo=%f",
                    fxt, fxo, ext, exo, fzt, fzo, ezt, ezo);
        sf_error ("Migration output range is out of traveltime table");
    }

    if (!sf_getfloat ("v0", &v0)) v0 = 1.5;
    /* reference velocity value at surface */
    if (!sf_getfloat ("dvz", &dvz)) dvz = 0;
    /* reference velocity vertical gradient */
    if (!sf_getfloat ("angmax", &angmax)) angmax = 60.;
    /* migration angle aperature from vertical */
    if (angmax < 0.00001) sf_error ("angmax= must be positive");
    mtmax = 2*dxm*sin(angmax*SF_PI/180.)/(v0*dt);
    if (mtmax < 1) mtmax = 1;

    if (!sf_getfloat ("aperx", &aperx)) aperx = 0.5*nxt*dxt;
    /* migration lateral aperature */
    if (!sf_getfloat ("offmax", &offmax)) offmax = 3.0;
    /* maximum absolute offset allowed in migration */
    if (!sf_getfloat ("fmax", &fmax)) fmax = 0.25/dt;
    /* frequency-highcut for input traces */
    if (!sf_getint ("noff", &noff)) noff = 1;
    /* number of offsets in output */
    if (!sf_getfloat ("off0", &off0)) off0 = 0.;
    /* first offest in output */
    if (!sf_getfloat ("doff", &doff)) doff = 0.1;
    /* offset increment in output */

    if (!sf_getint ("ls", &ls)) ls = 1;
    /* flag for line source */
    if (!sf_getint ("absoff", &absoff)) absoff = 0;
    /* 1 - use absolute value of offset, 0 - use offset =gx-sx */
    if (!sf_getint ("limoff", &limoff)) limoff = 0;
    /* 1 - limit traces used by offset, 0 - use all traces */
    if (!sf_getint ("ntr", &ntr)) ntr = sf_leftsize (infp, 1);
    /* maximum number of input traces to be migrated */
    if (!sf_getint ("mtr", &mtr)) mtr = 100;
    /* print verbal information at every mtr traces */
    if (!sf_getint ("npv", &npv)) npv = 0;
    /* 1 - compute quantities for velocity analysis */

    if (npv) {
        tvfile = sf_getstring ("tvfile");
        /* input file of traveltime variation tables */
        if (NULL == tvfile) sf_error ("Need tvfile=");
        tvfp = sf_input (tvfile);
        /* input file of traveltime variation tables */
        csfile = sf_getstring ("csfile");
        /* input file of cosine tables */
        if (NULL == csfile) sf_error ("Need csfile=");
        csfp = sf_input (csfile);
        /* input file of cosine tables */
        outfile1 = sf_getstring ("outfile1");
        /* file containning additional migration output */
        if (NULL == outfile1) sf_error ("Need outfile1=");
        out1fp = sf_output (outfile1);
        /* file containning additional migration output */
    }

    if (verb) {
        sf_warning (" Migration parameters");
        sf_warning (" ================n");
        sf_warning (" nzt=%d fzt=%g dzt=%g", nzt, fzt, dzt);
        sf_warning (" nxt=%d fxt=%g dxt=%g", nxt, fxt, dxt);
        sf_warning (" ns=%d fs=%g ds=%g", ns, fs, ds);
        sf_warning ("");
        sf_warning (" nzo=%d fzo=%g dzo=%g", nzo, fzo, dzo);
        sf_warning (" nxo=%d fxo=%g dxo=%g", nxo, fxo, dxo);
    }

    /* compute reference traveltime and slowness  */
    rmax = SF_MAX(es-fxt,ext-fs);
    rmax = SF_MIN(rmax,0.5*offmax+aperx);
    nr = 2+(int)(rmax/dxo);
    tb = sf_floatalloc2 (nzt, nr);
    pb = sf_floatalloc2 (nzt, nr);
    cs0b = sf_floatalloc2 (nzt, nr);
    angb = sf_floatalloc2 (nzt, nr);
    timeb (nr, nzt, dxo, dzt, fzt, dvz, v0, tb, pb, cs0b, angb);

    if (verb) {
        sf_warning (" nt=%d ft=%g dt=%g", nt, ft, dt);
        sf_warning (" dxm=%g fmax=%g", dxm, fmax);
        sf_warning (" noff=%d off0=%g doff=%g", noff, off0, doff);
        sf_warning (" v0=%g dvz=%g", v0, dvz);
        sf_warning (" aperx=%g offmax=%g angmax=%g", aperx, offmax, angmax);
        sf_warning (" ntr=%d mtr=%d ls=%d npv=%d", ntr, mtr, ls, npv);
        sf_warning (" absoff=%d limoff=%d", absoff, limoff);
    }

    /* allocate space */
    mig = sf_floatalloc3 (nzo, nxo, noff);
    ttab = sf_floatalloc3 (nzt, nxt, ns);
    tt = sf_floatalloc2 (nzt, nxt);
    tsum = sf_floatalloc2 (nzt, nxt);
    if (npv){
        tv = sf_floatalloc3 (nzt, nxt, ns);
        tvsum = sf_floatalloc2 (nzt, nxt);
        cs = sf_floatalloc3 (nzt, nxt, ns);
        cssum = sf_floatalloc2 (nzt, nxt);
    }
    if (!npv) 
        mig1 = sf_floatalloc3 (1, 1, noff);
    else
        mig1 = sf_floatalloc3 (nzo, nxo, noff);

     memset ((void*)mig[0][0], 0, noff*nxo*nzo*sizeof(float)); 
     if (npv)
         memset ((void *)mig1[0][0], 0, noff*nxo*nzo*sizeof(float));

     if (verb)
         sf_warning (" input traveltime tables");
                         
    /* compute traveltime residual */
    for (is = 0; is < ns; ++is) {
        sf_floatread (ttab[is][0], nxt*nzt, ttfp);
        s = fs + is*ds;
        resit (nxt, fxt, dxt, nzt, nr, dxo, tb, ttab[is], s);
        if (npv) {
            sf_floatread (tv[is][0], nxt*nzt, tvfp);
            sf_floatread (cs[is][0], nxt*nzt, csfp);
        }
    }
    sf_fileclose (ttfp);
    if (npv) {
        sf_fileclose (tvfp);
        sf_fileclose (csfp);
    }

    /* Set up output dimensions */
    sf_putint (outfp, "n1", nzo);
    sf_putint (outfp,"n2", nxo);
    sf_putfloat (outfp, "o1", fzo);
    sf_putfloat (outfp, "d1", dzo);
    sf_putfloat (outfp,"o2", fxo);
    sf_putfloat (outfp, "d2", dxo);
    sf_putstring (outfp, "label1", "Depth");
    sf_putstring (outfp, "label2", "Lateral");
    unit = sf_histstring (infp, "unit2");
    if (NULL != unit) sf_putstring (outfp, "unit1", unit);
    sf_putint (outfp, "n3", noff);
    sf_putfloat (outfp ,"o3", off0);
    sf_putfloat (outfp,"d3", doff);
    sf_putstring (outfp, "label3", "Offset");
    if (NULL != unit) sf_putstring (outfp, "unit3", unit);
    if (npv) {
        sf_putint (out1fp, "n1", nzo);
        sf_putint (out1fp,"n2", nxo);
        sf_putfloat (out1fp, "o1", fzo);
        sf_putfloat (out1fp, "d1", dzo);
        sf_putfloat (out1fp,"o2", fxo);
        sf_putfloat (out1fp, "d2", dxo);
        sf_putstring (out1fp, "label1", "Depth");
        sf_putstring (out1fp, "label2", "Lateral");
        if (NULL != unit) sf_putstring (out1fp, "unit1", unit);
        sf_putint (out1fp, "n3", noff);
        sf_putfloat (out1fp ,"o3", off0);
        sf_putfloat (out1fp,"d3", doff);
        sf_putstring (out1fp, "label3", "Offset");
        if (NULL != unit) sf_putstring (out1fp, "unit3", unit);
    }

    if (verb) {
        sf_warning (" start migration ... ");
        sf_warning ("");
    }

    jtr = 1;
    ktr = 0;
    iis = iig = 0;

    if (verb)
      sf_warning (" fs=%g es=%g offmax=%g", fs, es, offmax);

    data = sf_floatalloc (nt);

    do {
        sf_floatread (data, nt, infp); 
        /* determine offset index */

        sx = fsi + iis*dsi;
        gx = sx + fgi + iig*dgi;

        offset = gx - sx;
        if (absoff && offset < 0) offset =-offset;
        io = (int)((offset-off0)/doff+0.5);
        if (limoff && (io < 0 || io >= noff)) continue;

        if (io < 0) io = 0;
        if (io >= noff) io = noff-1;

        if (SF_MIN(sx,gx) >= fs && SF_MAX(sx,gx) <= es && 
            SF_MAX(gx-sx,sx-gx) <= offmax) {
            /* migrate this trace */

            as = (sx-fs)/ds;
            is = (int)as;
            if (is == ns-1)is = ns - 2;
            res = as-is;
            if (res <= 0.01) res = 0.0;
            if (res >= 0.99) res = 1.0;
            sum2 (nxt, nzt, 1-res, res, ttab[is], ttab[is+1], tsum);
            if (npv) {
                sum2 (nxt, nzt, 1-res, res, tv[is], tv[is+1], tvsum);
                sum2 (nxt, nzt, 1-res, res, cs[is], cs[is+1], cssum);
            }

            as = (gx-fs)/ds;
            is = (int)as;
            if (is == ns-1) is = ns-2;
            res = as-is;
            if (res <= 0.01) res = 0.0;
            if (res >= 0.99) res = 1.0;
            sum2 (nxt, nzt, 1-res, res, ttab[is], ttab[is+1], tt);
            sum2 (nxt, nzt, 1, 1, tt, tsum, tsum);
            if (npv) {
                sum2 (nxt, nzt, 1-res, res, tv[is], tv[is+1], tt);
                sum2 (nxt, nzt, 1, 1, tt, tvsum, tvsum);
                sum2 (nxt, nzt, 1-res, res, cs[is], cs[is+1], tt);
                sum2 (nxt, nzt, 1, 1, tt, cssum, cssum);
            }

            mig2d (data, nt, ft, dt, sx, gx, mig[io], aperx,
                   nxo, fxo, dxo, nzo, fzo, dzo,
                   ls, mtmax, dxm, fmax, angmax,
                   tb, pb, cs0b, angb, nr, tsum, nzt, fzt, dzt, nxt, fxt, dxt,
                   npv, cssum, tvsum, mig1[io]);

            ktr++;
            if ((jtr - 1) % mtr == 0 && verb)
                sf_warning (" migrated trace %d", jtr);
        }
        /* Next trace */
        jtr++;
        iig++;
        if (ngi == iig) { /* Next shot */
            iig = 0;
            iis++;
        }
    } while (jtr < ntr);

    if (verb)
        sf_warning (" migrated %d traces in total", ktr);

    scal = 4.0/sqrt(SF_PI)*dxm/v0;
    for (io = 0; io < noff; io++) {
        for (ixo = 0; ixo < nxo; ixo++) {
            for (izo = 0; izo < nzo; ++izo)
                mig[io][ixo][izo] *=scal;
            /* write out */
            sf_floatwrite (mig[io][ixo], nzo, outfp);

            if (npv){
                for (izo = 0; izo < nzo; ++izo)
                    mig1[io][ixo][izo] *=scal;
                /* write out */
                sf_floatwrite (mig1[io][ixo], nzo, out1fp);
            }
        }
    }

    if (verb) {
        sf_warning ("");
        sf_warning (" output done");
    }

    free (data);
    free (tsum[0]); free (tsum);
    free (tt[0]); free (tt);
    free (pb[0]); free (pb);
    free (tb[0]); free (tb);
    free (cs0b[0]); free (cs0b);
    free (angb[0]); free (angb);
    free (ttab[0][0]); free (ttab[0]); free (ttab);
    free (mig[0][0]); free (mig[0]); free (mig);
    free (mig1[0][0]); free (mig1[0]); free (mig1);
    if (npv) {
       free (tv[0][0]); free (tv[0]); free (tv);
       free (cs[0][0]); free (cs[0]); free (cs);
       free (tvsum[0]); free (tvsum);
       free (cssum[0]); free (cssum);
    }
    if (forw)
        free (forw);
    if (invs)
        free (invs);
    return 0;
}

/* residual traveltime calculation based on reference time  */
void resit (int nx, float fx, float dx, int nz, int nr, float dr,
            float **tb, float **t, float x0) {
        int ix,iz,jr;
        float xi,ar,sr,sr0;

        for(ix=0; ix<nx; ++ix){
                xi = fx+ix*dx-x0;
                ar = abs(xi)/dr;
                jr = (int)ar;
                sr = ar-jr;
                sr0 = 1.0-sr;
                if(jr>nr-2) jr = nr-2;
                for(iz=0; iz<nz; ++iz)
                        t[ix][iz] -= sr0*tb[jr][iz]+sr*tb[jr+1][iz];
        }
} 

/* lateral interpolation */

/* sum of two tables */
void sum2 (int nx, int nz, float a1, float a2, float **t1, float **t2, float **t) {
        int ix,iz;

        for(ix=0; ix<nx; ++ix) 
                for(iz=0; iz<nz; ++iz)
                        t[ix][iz] = a1*t1[ix][iz]+a2*t2[ix][iz];
}
 
/* compute  reference traveltime and slowness */
void timeb (int nr, int nz, float dr, float dz, float fz, float a,
            float v0, float **t, float **p, float **cs0, float **ang) {
        int  ir,iz;
        float r,z,v,rc,oa,temp,rou,zc;


        if( a==0.0) {
                for(ir=0,r=0;ir<nr;++ir,r+=dr)
                        for(iz=0,z=fz;iz<nz;++iz,z+=dz){
                                rou = sqrt(r*r+z*z);
                                if(rou<dz) rou = dz;
                                t[ir][iz] = rou/v0;
                                p[ir][iz] = r/(rou*v0);
                                cs0[ir][iz] = z/rou;
                                ang[ir][iz] = asin(r/rou);
                        }
        } else {
                oa = 1.0/a;         zc = v0*oa;
                for(ir=0,r=0;ir<nr;++ir,r+=dr)
                        for(iz=0,z=fz+zc;iz<nz;++iz,z+=dz){
                                rou = sqrt(r*r+z*z);
                                v = v0+a*(z-zc);
                                if(ir==0){ 
                                        t[ir][iz] = log(v/v0)*oa;
                                        p[ir][iz] = 0.0;
                                        ang[ir][iz] = 0.0;
                                        cs0[ir][iz] = 1.0;
                                } else {
                                        rc = (r*r+z*z-zc*zc)/(2.0*r);
                                        rou = sqrt(zc*zc+rc*rc);
                                        t[ir][iz] = log((v*(rou+rc))
                                                /(v0*(rou+rc-r)))*oa;
                                        p[ir][iz] = sqrt(rou*rou-rc*rc)
                                                /(rou*v0);
                                        temp = v*p[ir][iz];
                                        if(temp>1.0) temp = 1.0;
                                        ang[ir][iz] = asin(temp);
                                        cs0[ir][iz] = rc/rou;
                                }
                        }
        }
}

void filt (float *trace, int nt, float dt, float fmax, int ls, int m, float *trf);

void mig2d (float *trace, int nt, float ft, float dt,
            float sx, float gx, float **mig, float aperx,
            int nx, float fx, float dx, float nz, float fz, float dz,
            int ls, int mtmax, float dxm, float fmax, float angmax,
            float **tb, float **pb, float **cs0b, float **angb, int nr,
            float **tsum, int nzt, float fzt, float dzt, int nxt, float fxt, float dxt,
            int npv, float **cssum, float **tvsum, float **mig1)
/*****************************************************************************
Migrate one trace 
******************************************************************************
Input:
*trace                   one seismic trace 
nt                       number of time samples in seismic trace
ft                       first time sample of seismic trace
dt                       time sampleing interval in seismic trace
sx,gx                    lateral coordinates of source and geophone 
aperx                    lateral aperature in migration
nx,fx,dx,nz,fz,dz        dimension parameters of migration region
ls                       =1 for line source; =0 for point source
mtmax                    number of time samples in triangle filter
dxm                      midpoint sampling interval
fmax                     frequency-highcut for input trace
angmax                   migration angle aperature from vertical
tb,pb,cs0b,angb          reference traveltime, lateral slowness, cosine of
                         incident angle, and emergent angle
nr                       number of lateral samples in reference quantities
tsum                     sum of residual traveltimes from shot and receiver
nxt,fxt,dxt,nzt,fzt,dzt  dimension parameters of traveltime table
npv=0                    flag of computing quantities for velocity analysis
cssume                   sum of cosine of emergence angles from shot and recerver
tvsum                    sum of  traveltime variations from shot and recerver

Output:
mig                      migrated section
mig1                     additional migrated section for velocity analysis if npv>0
*****************************************************************************/
{
        int nxf,nxe,nxtf,nxte,ix,iz,iz0,izt0,nzp,jrs,jrg,jz,jt,mt,jx;
        float xm,x,dis,rxz,ar,srs,srg,srs0,srg0,sigp,z0,rdz,ampd,res0,
                angs,angg,cs0s,cs0g,ax,ax0,pmin,
                odt=1.0/dt,pd,az,sz,sz0,at,td,res,temp;
        float **tmt,**ampt,**ampti,**ampt1=NULL,*tm,*amp,*ampi,*amp1=NULL,
                *tzt,*trf,*zpt;

        tmt = sf_floatalloc2(nzt,nxt);
        ampt = sf_floatalloc2(nzt,nxt);
        ampti = sf_floatalloc2(nzt,nxt);
        tm = sf_floatalloc(nzt);
        tzt = sf_floatalloc(nzt);
        amp = sf_floatalloc(nzt);
        ampi = sf_floatalloc(nzt);
        zpt = sf_floatalloc(nxt);
        trf = sf_floatalloc(nt+2*mtmax);
        if(npv) {
                ampt1 = sf_floatalloc2(nzt,nxt);
                amp1 = sf_floatalloc(nzt);
        }

        z0 = (fz-fzt)/dzt;
        rdz = dz/dzt;
        pmin = 1.0/(2.0*dxm*fmax);
        
        filt(trace,nt,dt,fmax,ls,mtmax,trf);
/*        for (int ii = 0; ii < nt; ii++)
            trf[ii] = trace[ii];*/

        xm = 0.5*(sx+gx);
        rxz = (angmax==90)?0.0:1.0/tan(angmax*SF_PI/180.);
        nxtf = (xm-aperx-fxt)/dxt;
        if(nxtf<0) nxtf = 0;
        nxte = (xm+aperx-fxt)/dxt+1;
        if(nxte>=nxt) nxte = nxt-1;

        /* compute amplitudes and filter length */
        for(ix=nxtf; ix<=nxte; ++ix){
                x = fxt+ix*dxt;
                dis = (xm>=x)?xm-x:x-xm;
                izt0 = ((dis-dxt)*rxz-fzt)/dzt-1;
                if(izt0<0) izt0 = 0;
                if(izt0>=nzt) izt0 = nzt-1;

                ar = (sx>=x)?(sx-x)/dx:(x-sx)/dx;
                jrs = (int)ar;
                if(jrs>nr-2) jrs = nr-2;
                srs = ar-jrs;
                srs0 = 1.0-srs;
                ar = (gx>=x)?(gx-x)/dx:(x-gx)/dx;
                jrg = (int)ar;
                if(jrg>nr-2) jrg = nr-2;
                srg = ar-jrg;
                srg0 = 1.0-srg;
                sigp = ((sx-x)*(gx-x)>0)?1.0:-1.0;
                zpt[ix] = fzt+(nzt-1)*dzt;

                for(iz=izt0; iz<nzt; ++iz){
                        angs = srs0*angb[jrs][iz]+srs*angb[jrs+1][iz]; 
                        angg = srg0*angb[jrg][iz]+srg*angb[jrg+1][iz]; 
                        cs0s = srs0*cs0b[jrs][iz]+srs*cs0b[jrs+1][iz]; 
                        cs0g = srg0*cs0b[jrg][iz]+srg*cs0b[jrg+1][iz]; 
                        ampd = (cs0s+cs0g)*cos(0.5*(angs-angg));
                        if(ampd<0.0) ampd = -ampd;
                        ampt[ix][iz] = ampd;

                        pd = srs0*pb[jrs][iz]+srs*pb[jrs+1][iz]+sigp 
                             *(srg0*pb[jrg][iz]+srg*pb[jrg+1][iz]);
                        if(pd<0.0) pd = -pd;
                        temp = pd*dxm*odt;
                        if(temp<1) temp = 1.0;
                        if(temp>mtmax) temp = mtmax;
                        ampti[ix][iz] = ampd/(temp*temp);
                        tmt[ix][iz] = temp;
                        if(pd<pmin && zpt[ix]>fzt+(nzt-1.1)*dzt) 
                                zpt[ix] = fzt+iz*dzt;

                    if(npv){
                        if(cssum[ix][iz]<1.0) 
                             ampt1[ix][iz] = 0; 
                        else 
                             ampt1[ix][iz] = tvsum[ix][iz]/cssum[ix][iz];
                    }
                }
        }

        nxf = (xm-aperx-fx)/dx+0.5;
        if(nxf<0) nxf = 0;
        nxe = (xm+aperx-fx)/dx+0.5;
        if(nxe>=nx) nxe = nx-1;
        
        /* interpolate amplitudes and filter length along lateral        */
        for(ix=nxf; ix<=nxe; ++ix){
                x = fx+ix*dx;
                dis = (xm>=x)?xm-x:x-xm;
                izt0 = (dis*rxz-fzt)/dzt;
                if(izt0<0) izt0 = 0;
                if(izt0>=nzt) izt0 = nzt-1;
                iz0 = (dis*rxz-fz)/dz;
                if(iz0<0) iz0 = 0;
                if(iz0>=nz) iz0 = nz-1;

                ax = (x-fxt)/dxt;
                jx = (int)ax;
                ax = ax-jx;
                if(ax<=0.01) ax = 0.;
                if(ax>=0.99) ax = 1.0;
                ax0 = 1.0-ax;
                if(jx>nxte-1) jx = nxte-1;
                if(jx<nxtf) jx = nxtf;

                ar = (sx>=x)?(sx-x)/dx:(x-sx)/dx;
                jrs = (int)ar;
                if(jrs>nr-2) jrs = nr-2;
                srs = ar-jrs;
                srs0 = 1.0-srs;
                ar = (gx>=x)?(gx-x)/dx:(x-gx)/dx;
                jrg = (int)ar;
                if(jrg>nr-2) jrg = nr-2;
                srg = ar-jrg;
                srg0 = 1.0-srg;

                for(iz=izt0; iz<nzt; ++iz){
                    tzt[iz] = ax0*tsum[jx][iz]+ax*tsum[jx+1][iz]
                                +srs0*tb[jrs][iz]+srs*tb[jrs+1][iz]
                                +srg0*tb[jrg][iz]+srg*tb[jrg+1][iz];

                    amp[iz] = ax0*ampt[jx][iz]+ax*ampt[jx+1][iz];
                    ampi[iz] = ax0*ampti[jx][iz]+ax*ampti[jx+1][iz];
                    tm[iz] = ax0*tmt[jx][iz]+ax*tmt[jx+1][iz];

                    if(npv) 
                            amp1[iz] = ax0*ampt1[jx][iz]+ax*ampt1[jx+1][iz];

                }

                nzp = (ax0*zpt[jx]+ax*zpt[jx+1]-fz)/dz+1.5;
                if(nzp<iz0) nzp = iz0;
                if(nzp>nz) nzp = nz;

                /* interpolate along depth if operater aliasing         */
                for(iz=iz0; iz<nzp; ++iz) {
                        az = z0+iz*rdz;
                        jz = (int)az;
                        if(jz>=nzt-1) jz = nzt-2;
                        sz = az-jz;
                        sz0 = 1.0-sz;
                        td = sz0*tzt[jz]+sz*tzt[jz+1];
                        at = (td-ft)*odt+mtmax;
                        jt = (int)at;
                        if(jt > mtmax && jt < nt+mtmax-1){
                            ampd = sz0*ampi[jz]+sz*ampi[jz+1];
                            mt = (int)(0.5+sz0*tm[jz]+sz*tm[jz+1]);
                            res = at-jt;
                            res0 = 1.0-res;
                             temp = (res0*(-trf[jt-mt]+2.0*trf[jt]-trf[jt+mt]) 
                                +res*(-trf[jt-mt+1]+2.0*trf[jt+1]
                                -trf[jt+mt+1]))*ampd;
                            mig[ix][iz] += temp;

                            if(npv) 
                                mig1[ix][iz]  += temp
                                        *(sz0*amp1[jz]+sz*amp1[jz+1]);
                        }
                }

                /* interpolate along depth if not operater aliasing         */
                for(iz=nzp; iz<nz; ++iz) {
                        az = z0+iz*rdz;
                        jz = (int)az;
                        if(jz>=nzt-1) jz = nzt-2;
                        sz = az-jz;
                        sz0 = 1.0-sz;
                        td = sz0*tzt[jz]+sz*tzt[jz+1];
                        at = (td-ft)*odt;
                        jt = (int)at;
                        if(jt > 0 && jt < nt-1){
                            ampd = sz0*amp[jz]+sz*amp[jz+1];
                            res = at-jt;
                            res0 = 1.0-res;
                             temp = (res0*trace[jt]+res*trace[jt+1])*ampd; 
                            mig[ix][iz] += temp;
                            if(npv) 
                                mig1[ix][iz]  += temp
                                        *(sz0*amp1[jz]+sz*amp1[jz+1]);
                        }
                }

        }

        free (ampt[0]); free (ampt);
        free (ampti[0]); free (ampti);
        free (tmt[0]); free (tmt);
        free (amp);
        free (ampi);
        free (zpt);
        free (tm);
        free (tzt);
        free (trf);
        if (npv) {
           free (amp1);
           free (ampt1[0]); free (ampt1);
        }
}

void filt (float *trace, int nt, float dt, float fmax, int ls, int m, float *trf)
/* Low-pass filter, integration and phase shift for input data
   input: 
    trace(nt)        single seismic trace
    fmax             high cut frequency
    ls               ls=1, line source; ls=0, point source
  output:
    trace(nt)        filtered and phase-shifted seismic trace
    tracei(nt)       filtered, integrated and phase-shifted seismic trace
 */
{
        static int nfft=0, itaper, nw, nwf;
        static float *taper, *amp, *ampi, dw;
        int  it, iw, itemp;
        float temp, ftaper, const2, *rt;
/*        sf_complex *ct;*/
        kiss_fft_cpx *ct;

        fmax *= 2.0*SF_PI;
        ftaper = 0.1*fmax;
        const2 = sqrt(2.0);

        if(nfft==0) {
                  /* Set up FFT parameters */
                  nfft = 2*kiss_fft_next_fast_size((nt+m+1)/2);
/*                  nfft = npfaro(nt+m, 2 * (nt+m));
                  if (nfft >= SU_NFLTS || nfft >= 720720)
                            err("Padded nt=%d -- too big", nfft);*/

                  nw = nfft/2 + 1;
                dw = 2.0*SF_PI/(nfft*dt);

                itaper = 0.5+ftaper/dw;
                taper = sf_floatalloc(2*itaper+1);
                for(iw=-itaper; iw<=itaper; ++iw){
                        temp = (float)iw/(1.0+itaper); 
                        taper[iw+itaper] = (1-temp)*(1-temp)*(temp+2)/4;
                }

                nwf = 0.5+fmax/dw;
                if(nwf>nw-itaper-1) nwf = nw-itaper-1;
                amp = sf_floatalloc(nwf+itaper+1);
                ampi = sf_floatalloc(nwf+itaper+1);
                amp[0] = ampi[0] = 0.;
                for(iw=1; iw<=nwf+itaper; ++iw){
                        amp[iw] = sqrt(dw*iw)/nfft;
                        ampi[iw] = 0.5/(1-cos(iw*dw*dt));
                }
                forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
                invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
        }

          /* Allocate fft arrays */
          rt   = sf_floatalloc(nfft);
/*          ct   = sf_complexalloc (nw);*/
          ct = (kiss_fft_cpx*) sf_complexalloc(nw);

          memcpy(rt, trace, nt*sizeof(float));
          memset((void *) (rt + nt), (int) '\0', (nfft-nt)*sizeof(float)); 
/*          pfarc(1, nfft, rt, ct);*/
          kiss_fftr(forw,rt,ct);

        for(iw=nwf-itaper;iw<=nwf+itaper;++iw){
                itemp = iw-(nwf-itaper);
                ct[iw].r = taper[itemp]*ct[iw].r; 
                ct[iw].i = taper[itemp]*ct[iw].i; 
        }
        for(iw=nwf+itaper+1;iw<nw;++iw){
                ct[iw].r = 0.; 
                ct[iw].i = 0.; 
        }

                 if(!ls){
                for(iw=0; iw<=nwf+itaper; ++iw){
                        /* phase shifts PI/4         */
                        temp = (ct[iw].r-ct[iw].i)*amp[iw]*const2;
                        ct[iw].i = (ct[iw].r+ct[iw].i)*amp[iw]*const2;
                        ct[iw].r = temp;
                    }
        } else {
                for(iw=0; iw<=nwf+itaper; ++iw){
                        ct[iw].i = ct[iw].i*amp[iw];
                        ct[iw].r = ct[iw].r*amp[iw];
                }
        }                  
/*          pfacr(-1, nfft, ct, rt);*/
          kiss_fftri(invs,ct,rt);
                
          /* Load traces back in */
        for (it=0; it<nt; ++it) trace[it] = rt[it];

          /* Integrate traces   */
        for(iw=0; iw<=nwf+itaper; ++iw){
                ct[iw].i = ct[iw].i*ampi[iw];
                ct[iw].r = ct[iw].r*ampi[iw];
        }
/*          pfacr(-1, nfft, ct, rt);*/
          kiss_fftri(invs,ct,rt);
          for (it=0; it<m; ++it)  trf[it] = rt[nfft-m+it];
          for (it=0; it<nt+m; ++it)  trf[it+m] = rt[it];

        free (rt);
        free (ct);
}

