/* Velocity Scan by 2D Parametric Development of Reflections */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// files
sf_file dataFile;
sf_file outFile;
sf_file auxFile;          // output file - semblance

int   tNum_;                 
float tStart_;
float tStep_;

int   recNum_;                 
float recStart_;
float recStep_;

int   shotNum_;                 
float shotStart_;
float shotStep_;

// stacking params
int   pNum_;                 
float pStart_;
float pStep_;

// scan parameters
int   vNum_;
float vStep_;
float vStart_;

int   wh_;

int main (int argc, char* argv[]) 
{
    //
    char* corUnit;
    char* unit;

    int zoSize, dataSize, tNumRed;
    float *zo, *zoSq, *semb, *data;
    int *count;
    int iv, is, ir, ip, forDataInd, it, tInd, dataInd, indZO;
    float vel, shotPos, curOffset, halfOffset, fabsOffset, offsetSq;
    float curPos, l0, forA, t0, a, t, forLim, limitLeft, limitRight;
    float bef, aft, sample, sampleSq, sqSample, curSemb;
    int ic, iw, ts, ind, totalCount, offset, vwhalf;

// Initialize RSF 
    sf_init (argc,argv);

// INPUT FILES
    dataFile = sf_input ("in");
    /* common-offset sections */
    outFile  = sf_output("out");
    /*  */

    if ( NULL != sf_getstring("aux") ) {
	/* output file containing semblance measure of CIGs stacking */ 
	auxFile  = sf_output ("aux");
    } else {
	sf_error ("Need output: partial zero-offset sections");
    }

    // data params

    if ( !sf_histint   (dataFile, "n1", &tNum_)   ) sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile, "d1", &tStep_)  ) sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile, "o1", &tStart_) ) sf_error ("Need o1= in input");

    if ( !sf_histint   (dataFile, "n2", &recNum_)   ) sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile, "d2", &recStep_)  ) sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile, "o2", &recStart_) ) sf_error ("Need o2= in input");

    if ( !sf_histint   (dataFile, "n3", &shotNum_)   ) sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile, "d3", &shotStep_)  ) sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile, "o3", &shotStart_) ) sf_error ("Need o3= in input");

//    if ( !sf_histint   (dataFile, "n4", &dp.yNum)   ) sf_error ("Need n4= in input");
//    if ( !sf_histfloat (dataFile, "d4", &dp.yStep)  ) sf_error ("Need d4= in input");
//    if ( !sf_histfloat (dataFile, "o4", &dp.yStart) ) sf_error ("Need o4= in input");

	
    // data
    // time - in ms
//    corUnit = (char*) "ms"; unit = sf_histstring (dataFile, "unit1"); if (!unit) sf_error ("unit1 in data file is not defined");
//    if ( strcmp (corUnit, unit) ) { dp.zStep *= 1000; dp.zStart *= 1000; }
    // receiver - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit2"); if (!unit) sf_error ("unit2 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { recStep_ *= 1000; recStart_ *= 1000; }
    // source - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit3"); if (!unit) sf_error ("unit3 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { shotStep_ *= 1000; shotStart_ *= 1000; }

    // offset - in m
    //   corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit4"); if (!unit) sf_error ("unit4 in data file is not defined");
//    if ( strcmp (corUnit, unit) ) { dp.hStep *= 1000; dp.hStart *= 1000; }

    if ( !sf_getfloat ("po",    &pStart_) ) pStart_ = shotStart_;
    /* start position in stack section */
    if ( !sf_getint ("pn", &pNum_) ) pNum_ = recNum_;
    /* number of positions in stack section */
    if (!pNum_) {sf_warning ("vn value is changed to 1"); pNum_ = recNum_;}
    if ( !sf_getfloat ("pd",    &pStep_) ) pStep_ = recStep_;
    /* increment of positions in stack section */
    if (!pStep_) {sf_warning ("pd value is changed to 50"); pStep_ = recStep_;}


    if ( !sf_getint ("vn",    &vNum_) ) vNum_ = 1;
    /* number of scanned velocities  */
    if (!vNum_) {sf_warning ("vn value is changed to 1"); vNum_ = 1;}
    if ( !sf_getfloat ("vo",    &vStart_) ) vStart_ = 1500;
    /* start velocity */
    if (!vStart_) {sf_warning ("vn value is changed to 1500"); vStart_ = 1500.0;}
    if ( !sf_getfloat ("vd",    &vStep_) ) vStep_ = 50;
    /* increment of velocities */
    if (!vStep_) {sf_warning ("vd value is changed to 50"); vStep_ = 50.f;}

    if ( !sf_getint ("wh",   &wh_) )   wh_ = 11;
    /* height of a vertical window for semblance calculation */
    if (!wh_) {sf_warning ("vertical window size is changed to 1"); wh_ = 1;}

    sf_putint    (outFile, "n1", tNum_);
    sf_putint    (outFile, "n2", pNum_);
    sf_putint    (outFile, "n3", vNum_);
    sf_putfloat  (outFile, "d1", tStep_); 
    sf_putfloat  (outFile, "d2", pStep_);
    sf_putint    (outFile, "d3", vStep_);
    sf_putfloat  (outFile, "o1", tStart_); 
    sf_putfloat  (outFile, "o2", pStart_);
    sf_putfloat  (outFile, "o3", vStart_);
    sf_putstring (outFile, "label1", "time"); sf_putstring(outFile, "label2", "inline"); sf_putstring(outFile, "label3", "velocity");
    sf_putstring(outFile, "unit3", "m/s");

    sf_putint    (auxFile, "n1", tNum_);
    sf_putint    (auxFile, "n2", pNum_);
    sf_putint    (auxFile, "n3", vNum_);
    sf_putfloat  (auxFile, "d1", tStep_); 
    sf_putfloat  (auxFile, "d2", pStep_);
    sf_putint    (auxFile, "d3", vStep_);
    sf_putfloat  (auxFile, "o1", tStart_); 
    sf_putfloat  (auxFile, "o2", pStart_);
    sf_putfloat  (auxFile, "o3", vStart_);
    sf_putstring (auxFile, "label1", "time"); sf_putstring (auxFile, "label2", "inline"); sf_putstring (auxFile, "label3", "velocity");
    sf_putstring (auxFile, "unit3", "m/s");

    zoSize = pNum_ * tNum_;
    dataSize = shotNum_ * recNum_ * tNum_;
    tNumRed = tNum_ - 2;
	
    zo    = sf_floatalloc (zoSize);
    zoSq  = sf_floatalloc (zoSize);
    semb  = sf_floatalloc (zoSize);
    count = sf_intalloc (zoSize);

    data = sf_floatalloc (dataSize);
    sf_floatread (data, dataSize, dataFile);

    for (iv = 0; iv < vNum_; ++iv) {

	vel = vStart_ + iv * vStep_;
    
	memset ( zo,    0, zoSize * sizeof (float) );   
	memset ( zoSq,  0, zoSize * sizeof (float) );   
	memset ( semb,  0, zoSize * sizeof (float) );   
	memset ( count, 0, zoSize * sizeof (int)   );   

	// loop over shots
	for (is = 0; is < shotNum_; ++is) {				
	    sf_warning ("scanning: velocity %d of %d - %g; shot %d of %d;", iv + 1, vNum_, vel, is + 1, shotNum_);
	    shotPos = shotStart_ + shotStep_ * is;
	    // loop over receivers
	    for (ir = 0; ir < recNum_; ++ir) {						
		curOffset = recStart_ + recStep_ * ir;
		halfOffset = curOffset / 2.f;
		fabsOffset = fabs (curOffset);
		offsetSq = curOffset * halfOffset;
		forDataInd = (is * recNum_ + ir) * tNum_;
#ifdef _OPENMP 
#pragma omp parallel for
#endif					
		for (ip = 0; ip < pNum_; ++ip) {
		    curPos = pStart_ + ip * pStep_;
		    l0 = curPos - shotPos;			
		    if (fabsOffset <= fabs (l0) || curOffset * l0 <= 0) continue;
		    forA = 4 * l0 * (curOffset - l0);
		    for (it = 0; it < tNum_; ++it) {	
			t0 = tStart_ + it * tStep_;

			// get time
			a = t0 * t0 / forA;
			t = fabsOffset * sqrt (a + 1 / pow (vel, 2) );

			// calc curve limits
			forLim = offsetSq / (vel * t);
			limitLeft  = halfOffset - forLim;
			limitRight = halfOffset + forLim;					
					
			if (l0 < limitLeft || l0 > limitRight) continue;

			tInd = (t - tStart_) / tStep_;
			if (tInd < 0 || tInd > tNumRed) continue; 

			bef = (t - tInd * tStep_) / tStep_;
			aft = 1.f - bef;

			dataInd = forDataInd + tInd;
			sample = data [dataInd] * aft + data [dataInd + 1] * bef;

			indZO   = ip * tNum_ + it;
			zo    [indZO] += sample;
			zoSq  [indZO] += sample*sample;
			count [indZO] += 1;									
		    }
		}
	    }
	}
	// semblance calculation
	vwhalf = wh_ / 2;
	for (ip = 0; ip < pNum_; ++ip) {
	    ts = ip * tNum_;
#ifdef _OPENMP 
#pragma omp parallel for
#endif					
	    for (it = 0; it < tNum_; ++it) {	

		sampleSq = 0.f;	
		sqSample = 0.f;
			
		/* ccount = 0; */
		totalCount = 0;

		for (ic = 0, iw = it - vwhalf; ic < wh_; ++ic, ++iw) {
		    if (iw < 0 || iw > tNumRed) continue;
		    ind = ts + iw;
		    sampleSq   += pow (zo [ind], 2);
		    sqSample   += zoSq [ind];
		    if (totalCount < count [ind]) totalCount = count [ind];
		}
		curSemb = sqSample && totalCount ? sampleSq / ( totalCount * sqSample ) : 0.f;
		semb [ts + it] = curSemb;
	    }
	}	

	offset = zoSize * iv * sizeof (float);
	sf_seek (outFile, offset, SEEK_SET);
	sf_floatwrite (semb, zoSize, outFile);
	sf_seek (auxFile, offset, SEEK_SET);
	sf_floatwrite (zo, zoSize, auxFile);
    }

    free (data);
    free (zo);
    free (zoSq);
    free (semb);
    free (count);

    sf_fileclose (dataFile);
    sf_fileclose (outFile);
    sf_fileclose (auxFile);

    return 0;
}
