/* 3D time scattering-angle Kirchhoff migration  */
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

#include <rsf.hh>

#include "support.hh"
#include "tmigratorBase.hh"
#include "tmigrator2D.hh"
#include "tmigrator3D.hh"
#include "sembler.hh"

RunParamsTmigda  rp; // migration (running) params
VolumeParams     dp; // data params
VolumeParams     vp; // velocity params
VolumeParams     ip; // image params
GatherParams     gp; // gather params

// files
sf_file dataFile;
sf_file velFile;

sf_file imageFile;
sf_file sembFile;
sf_file dagFile;
sf_file cigFile;

//  Causal integration of a trace[n]
void applyCasualIntegration (float *trace, int n) {

    for (int i = 1; i < n; ++i)
        trace[i] += trace[i - 1];
        
    return;      
}

// Anticausal integrations of a trace[n]
void applyAnticasualIntegration (float *trace, int n) {

    for (int i = n - 2; i >= 0; --i)
        trace[i] += trace[i + 1];

    return;
}

void readOffsetSection (int iOffset, float* offsetSection) {

    const int tracesNum = dp.xNum * dp.yNum;
    const int sectionSize = dp.zNum * tracesNum;
    memset (offsetSection, 0, sectionSize * sizeof (float));

    const size_t startPos = iOffset * sectionSize * sizeof(float);

    sf_seek (dataFile, startPos, SEEK_SET);
    sf_floatread (offsetSection, sectionSize, dataFile);

    if (rp.isAA) {
	float* ptrData = offsetSection;
	for (int it = 0; it < tracesNum; ++it, ptrData += dp.zNum) {
	    applyCasualIntegration (ptrData, dp.zNum);
	    applyAnticasualIntegration (ptrData, dp.zNum);
	}
    }
	
    return;
}

void checkImageParams () {

    const float zVelMax = vp.zStart + vp.zStep * (vp.zNum - 1);
    const float xVelMax = vp.xStart + vp.xStep * (vp.xNum - 1);
    const float yVelMax = vp.yStart + vp.yStep * (vp.yNum - 1);

    // time checking

    if (vp.zStart - ip.zStart > 1E-6) {
	float diff = vp.zStart - ip.zStart;
	ip.zNum -= (int) floorf(diff / ip.zStep + 1);
	ip.zStart += (diff / ip.zStep + 1) * ip.zStep;
	sf_warning ("first migrated time was changed to %g", ip.zStart);			
    }
    if (ip.zNum <= 0) sf_error ("wrong image time sample number");	
    float zImageMax = ip.zStart + ip.zStep * (ip.zNum - 1);
    if (zImageMax - zVelMax > 1E-6) {
	float diff = zImageMax - zVelMax; 
	ip.zNum -= (int) floorf(diff / ip.zStep + 1);
	sf_warning ("time sample number was changed to %d", ip.zNum);	
    }
    if (ip.zNum <= 0) sf_error ("wrong image time sample number");	

    // inline checking

    if (vp.xStart - ip.xStart > 1E-6) {
	float diff = vp.xStart - ip.xStart;
	ip.xNum -= (int) floorf(diff / ip.xStep + 1);
	ip.xStart += (diff / ip.xStep + 1) * ip.xStep;
	sf_warning ("first inline coord was changed to %g", ip.xStart);			
    }
    if (ip.xNum <= 0) sf_error ("wrong image inline number");	
    float xImageMax = ip.xStart + ip.xStep * (ip.xNum - 1);
    if (xImageMax - xVelMax > 1E-6) {
	float diff = xImageMax - xVelMax; 
	ip.xNum -= (int) floorf(diff / ip.xStep + 1);
	sf_warning ("inline number was changed to %d", ip.xNum);	
    }
    if (ip.xNum <= 0) sf_error ("wrong image inline number");	

    // crossline checking

	if (!rp.is3D)
		return; // 2D migration is runnig

    if (vp.yStart - ip.yStart > 1E-6) {
	float diff = vp.yStart - ip.yStart;
	ip.yNum -= (int) floorf(diff / ip.yStep + 1);
	ip.yStart += (diff / ip.yStep + 1) * ip.yStep;
	sf_warning ("first crossline coord was changed to %g", ip.yStart);			
    }
    if (ip.yNum <= 0) sf_error ("wrong image crossline number");	
    float yImageMax = ip.yStart + ip.yStep * (ip.yNum - 1);
    if (yImageMax - yVelMax > 1E-6) {
	float diff = yImageMax - yVelMax; 
	ip.yNum -= (int) floorf(diff / ip.yStep + 1);
	sf_warning ("crossline number was changed to %d", ip.yNum);	
    }
    if (ip.yNum <= 0) sf_error ("wrong image crossline number");	
	
    return;
}

void prepareVelocityTrace (int taskX, int taskY, float* velTrace) {

    const int velZNum = vp.zNum;

    // -- clean previous velocity trace	
    memset (velTrace, 0, velZNum * sizeof (float));   

    // -- read velocity trace
    const float geoX = ip.xStart + taskX * ip.xStep;
    const float geoY = ip.yStart + taskY * ip.yStep;	
	
    const int velX = (int) roundf((geoX - vp.xStart) / vp.xStep);
    const int velY = rp.is3D ? (int) roundf((geoY - vp.yStart) / vp.yStep) : 0;

    const size_t startPos = (velX + velY * vp.xNum) * vp.zNum * sizeof(float);

    sf_seek (velFile, startPos, SEEK_SET);
    sf_floatread (velTrace, vp.zNum, velFile);

    if (!rp.isVelMS) {	// velocity is in km/s - we reduce it to m/s
	float* ptrVel = velTrace;
	const int size = vp.zNum;
	for (int iv = 0; iv < size; ++iv, ++ptrVel)
	    *ptrVel *= 1000;
    }

    return;
}

int main (int argc, char* argv[]) {

// Initialize RSF 
    sf_init (argc,argv);

// INPUT FILES
    dataFile = sf_input ("in");
    /* common-offset sections */

    if ( NULL != sf_getstring("vel") ) {
	/* velocity model file (velocity in m/s) */ 
	velFile  = sf_input ("vel");
	float firstvel;
	sf_floatread (&firstvel, 1, velFile);
	rp.isVelMS = true;		
	if (firstvel < 15) {
	    sf_warning ("it seems that velocity is in km/s - will be multiplied on 1000");	
	    rp.isVelMS = false;					
	}			

    } else { sf_error ("Need input: velocity model"); }

// check that the input is float 
    if ( SF_FLOAT != sf_gettype (dataFile) ) sf_error ("Need float input: common-offset sections");
    if ( SF_FLOAT != sf_gettype (velFile)  ) sf_error ("Need float input: velocity model");

// OUTPUT FILES

    imageFile  = sf_output("out");
    /* migrated gathers in the dip-angle domain */

    if ( NULL != sf_getstring("semb") ) {
	/* output file containing semblance measure of CIGs stacking */ 
	sembFile  = sf_output ("semb"); rp.isSemb = true;
    } else { rp.isSemb = false; }

    if ( NULL != sf_getstring("dag") ) {
	/* output file containing CIGs in the dip-angle domain */ 
	dagFile  = sf_output ("dag"); rp.isDag = true;
    } else { rp.isDag = false; }

    if ( NULL != sf_getstring("cig") ) {
	/* output file containing CIGs in the surface-offset domain */ 
	cigFile  = sf_output ("cig"); rp.isCig = true;
    } else { rp.isCig = false; }

    // data params

    if ( !sf_histint   (dataFile, "n1", &dp.zNum) )   sf_error ("Need n1= in input");
    if ( !sf_histfloat (dataFile, "d1", &dp.zStep) )  sf_error ("Need d1= in input");
    if ( !sf_histfloat (dataFile, "o1", &dp.zStart) ) sf_error ("Need o1= in input");

    if ( !sf_histint   (dataFile, "n2", &dp.xNum) )   sf_error ("Need n2= in input");
    if ( !sf_histfloat (dataFile, "d2", &dp.xStep) )  sf_error ("Need d2= in input");
    if ( !sf_histfloat (dataFile, "o2", &dp.xStart) ) sf_error ("Need o2= in input");

    if ( !sf_histint   (dataFile, "n3", &dp.yNum) )   sf_error ("Need n3= in input");
    if ( !sf_histfloat (dataFile, "d3", &dp.yStep) )  sf_error ("Need d3= in input");
    if ( !sf_histfloat (dataFile, "o3", &dp.yStart) ) sf_error ("Need o3= in input");

    if ( !sf_histint   (dataFile, "n4", &dp.hNum) )   sf_error ("Need n4= in input");
    if ( !sf_histfloat (dataFile, "d4", &dp.hStep) )  sf_error ("Need d4= in input");
    if ( !sf_histfloat (dataFile, "o4", &dp.hStart) ) sf_error ("Need o4= in input");

    //
    char* corUnit;
    char* unit;
	
    // data
    // time - in ms
    corUnit = (char*) "ms"; unit = sf_histstring (dataFile, "unit1"); if (!unit) sf_error ("unit1 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.zStep *= 1000; dp.zStart *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit2"); if (!unit) sf_error ("unit2 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.xStep *= 1000; dp.xStart *= 1000; }
    // crossline - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit3"); if (!unit) sf_error ("unit3 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.yStep *= 1000; dp.yStart *= 1000; }
    // offset - in m
    corUnit = (char*) "m"; unit = sf_histstring (dataFile, "unit4"); if (!unit) sf_error ("unit4 in data file is not defined");
    if ( strcmp (corUnit, unit) ) { dp.hStep *= 1000; dp.hStart *= 1000; }
	
    // VELOCITY MODEL PARAMS

    if ( !sf_histint   (velFile, "n1", &vp.zNum) )   sf_error ("Need n1= in velocity file");
    if ( !sf_histfloat (velFile, "d1", &vp.zStep) )  sf_error ("Need d1= in velocity file");
    if ( !sf_histfloat (velFile, "o1", &vp.zStart) ) sf_error ("Need o1= in velocity file");

    if ( !sf_histint   (velFile, "n2", &vp.xNum) )   sf_error ("Need n2= in velocity file");
    if ( !sf_histfloat (velFile, "d2", &vp.xStep) )  sf_error ("Need d2= in velocity file");
    if ( !sf_histfloat (velFile, "o2", &vp.xStart) ) sf_error ("Need o2= in velocity file");

    if ( !sf_histint   (velFile, "n3", &vp.yNum) )   sf_error ("Need n3= in velocity file");
    if ( !sf_histfloat (velFile, "d3", &vp.yStep) )  sf_error ("Need d3= in velocity file");
    if ( !sf_histfloat (velFile, "o3", &vp.yStart) ) sf_error ("Need o3= in velocity file");

    // time - in ms
    corUnit = (char*) "ms"; unit = sf_histstring (velFile, "unit1"); if (!unit) sf_error ("unit1 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { vp.zStep *= 1000; vp.zStart *= 1000; }
    // inline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit2"); if (!unit) sf_error ("unit2 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { vp.xStep *= 1000; vp.xStart *= 1000; }
    // crossline - in m
    corUnit = (char*) "m"; unit = sf_histstring (velFile, "unit3"); if (!unit) sf_error ("unit3 in velocity model is not defined");
    if ( strcmp (corUnit, unit) ) { vp.yStep *= 1000; vp.yStart *= 1000; }

// Migration parameters
    if (!sf_getbool ("is3d",       &rp.is3D))       rp.is3D = false;
    /* if y, apply 3D migration */
	int axis2label (0);
    if ( !sf_getint ("axis2label", &axis2label) )   axis2label = 0;
	/* 0 - shot; 1 - cmp; 2 - receiver */
    if (!sf_getbool ("isAA",       &rp.isAA))       rp.isAA = true;
    /* if y, apply anti-aliasing */
    if (!sf_getbool ("isDipAz",    &rp.isDipAz))    rp.isDipAz = true;
    /* if y, apply dip/azimuth mode; if n, apply inline/crossline angle mode */
    if (!sf_getint  ("hmign",   &rp.hMigNum)) rp.hMigNum = dp.hNum;	
    /* number of migrated offsets */
    if (!sf_getint  ("sembWindow",   &rp.sembWindow)) rp.sembWindow = 11;	
    /* vertical window for semblance calculation (in samples) */
    if (!sf_getfloat  ("edgeTaper",   &rp.edgeTaper)) rp.edgeTaper = 5.f;	
    /* edge taper for dip-angle gathers (in degree) */

    // IMAGE PARAMS
    if (!sf_getint ("itn", &ip.zNum))        ip.zNum = dp.zNum;	
    /* number of imaged times */
    if (!sf_getint ("ixn", &ip.xNum))        ip.xNum = dp.xNum;	
    /* number of imaged inlines */
    if (!sf_getint ("iyn", &ip.yNum))        ip.yNum = rp.is3D ? vp.yNum : 1;	
    /* number of imaged crosslines */
    if (!sf_getfloat ("ito", &ip.zStart))    ip.zStart = dp.zStart;
    /* first imaged time (in ms) */
    if (!sf_getfloat ("ixo", &ip.xStart))    ip.xStart = dp.xStart;
    /* first imaged inline */
    if (!sf_getfloat ("iyo", &ip.yStart))    ip.yStart = dp.yStart;	
    /* first imaged crossline */
    if (!sf_getfloat ("itd", &ip.zStep))     ip.zStep = dp.zStep;
    /* step in imaged times  (in ms) */
    if (!sf_getfloat ("ixd", &ip.xStep))     ip.xStep = dp.xStep;	
    /* step in imaged inlines */
    if (!sf_getfloat ("iyd", &ip.yStep))     ip.yStep = dp.yStep;
    /* step in imaged crosslines */

    checkImageParams ();

    // GATHER PARAMS
    gp.zNum = ip.zNum;
    if (!sf_getint ("dipn" , &gp.dipNum))      gp.dipNum = 1;	
    /* number of dip-angles */
    if (!sf_getint ("sdipn", &gp.sdipNum))     gp.sdipNum = 1;	
    /* number of secondary (azimuth or crossline) angles */
    gp.zStart = ip.zStart;
    if (!sf_getfloat ("dipo",  &gp.dipStart))   gp.dipStart = 0.f;	
    /* first dip-angle */
    if (!sf_getfloat ("sdipo", &gp.sdipStart))  gp.sdipStart = 90.f;
    /* first secondary (azimuth or crossline) angle */
    gp.zStep = ip.zStep;
    if (!sf_getfloat ("dipd",  &gp.dipStep))   gp.dipStep = 1.f;	
    /* step in dip-angle */
    if (!sf_getfloat ("sdipd", &gp.sdipStep))  gp.sdipStep = 1.f;	
    /* step in secondary (azimuth or crossline) angle */

    // Initiate output 

    // image file
    sf_putint (imageFile, "n1", ip.zNum); sf_putint (imageFile, "n2", ip.xNum); sf_putint (imageFile, "n3", ip.yNum); sf_putint (imageFile, "n4", 1);
    sf_putfloat (imageFile, "d1", ip.zStep); sf_putfloat (imageFile, "d2", ip.xStep); sf_putfloat (imageFile, "d3", ip.yStep); 
    sf_putfloat (imageFile, "d4", 1);   
    sf_putfloat (imageFile, "o1", ip.zStart); sf_putfloat (imageFile, "o2", ip.xStart); sf_putfloat (imageFile, "o3", ip.yStart);    
    sf_putfloat (imageFile, "o4", 0);   
    sf_putstring(imageFile, "label1", "Time"); sf_putstring(imageFile, "label2", "Inline"); sf_putstring(imageFile, "label3", "Crossline");
    sf_putstring(imageFile, "unit1", "ms"); sf_putstring(imageFile, "unit2", "m"); sf_putstring(imageFile, "unit3", "m");

    if (rp.isCig) {
	// offset gathers file
    	sf_putint (cigFile, "n1", ip.zNum); sf_putint (cigFile, "n2", rp.hMigNum); sf_putint (cigFile, "n3", ip.xNum); 
	sf_putint (cigFile, "n4", ip.yNum); 
    	sf_putfloat (cigFile, "d1", ip.zStep); sf_putfloat (cigFile, "d2", dp.hStep); sf_putfloat (cigFile, "d3", ip.xStep);
	sf_putfloat (cigFile, "d4", ip.yStep); 
    	sf_putfloat (cigFile, "o1", ip.zStart); sf_putfloat (cigFile, "o2", dp.hStart); sf_putfloat (cigFile, "o3", ip.xStart);
	sf_putfloat (cigFile, "o4", ip.yStart);    
	sf_putstring(cigFile, "label1", "Time"); sf_putstring(cigFile, "label2", "Offset"); 
	sf_putstring(cigFile, "label3", "Inline"); sf_putstring(cigFile, "label4", "Crossline");
	sf_putstring(cigFile, "unit1", "ms"); sf_putstring(cigFile, "unit2", "m"); 
	sf_putstring(cigFile, "unit3", "m"); sf_putstring(cigFile, "unit4", "m");
    }

    if (rp.isDag) {
	// dip-angle gathers file
	sf_putint (dagFile, "n1", ip.zNum); sf_putint (dagFile, "n2", gp.dipNum); sf_putint (dagFile, "n3", gp.sdipNum);
	sf_putint (dagFile, "n4", ip.xNum); sf_putint (dagFile, "n5", ip.yNum);
    	sf_putfloat (dagFile, "d1", ip.zStep); sf_putfloat (dagFile, "d2", gp.dipStep); sf_putfloat (dagFile, "d3", gp.sdipStep);
	sf_putfloat (dagFile, "d4", ip.xStep); sf_putfloat (dagFile, "d5", ip.yStep);    
    	sf_putfloat (dagFile, "o1", ip.zStart); sf_putfloat (dagFile, "o2", gp.dipStart); sf_putfloat (dagFile, "o3", gp.sdipStart);
	sf_putfloat (dagFile, "o4", ip.xStart); sf_putfloat (dagFile, "o5", ip.yStart);    
	sf_putstring(dagFile, "label1", "Time");
	if (rp.isDipAz) {
	    sf_putstring(dagFile, "label2", "Dip angle"); sf_putstring(dagFile, "label3", "Azimuth");
	} else {
	    sf_putstring(dagFile, "label2", "Inline slope"); sf_putstring(dagFile, "label3", "Crossline slope");
	}
	sf_putstring(dagFile, "label4", "Inline"); 	sf_putstring(dagFile, "label5", "Crossline");
	sf_putstring(dagFile, "unit1", "ms"); sf_putstring(dagFile, "unit2", "deg"); sf_putstring(dagFile, "unit3", "deg");
	sf_putstring(dagFile, "unit4", "m"); sf_putstring(dagFile, "unit5", "m");
    }

    if (rp.isSemb) {
	sf_putint (sembFile, "n1", ip.zNum); sf_putint (sembFile, "n2", ip.xNum); sf_putint (sembFile, "n3", ip.yNum); sf_putint (sembFile, "n4", 1);
    	sf_putfloat (sembFile, "d1", ip.zStep); sf_putfloat (sembFile, "d2", ip.xStep); sf_putfloat (sembFile, "d3", ip.yStep); 
	sf_putfloat (sembFile, "d4", 1);   
    	sf_putfloat (sembFile, "o1", ip.zStart); sf_putfloat (sembFile, "o2", ip.xStart); sf_putfloat (sembFile, "o3", ip.yStart);    
	sf_putfloat (sembFile, "o4", 0);   
	sf_putstring(sembFile, "label1", "Time"); sf_putstring (sembFile, "label2", "Inline"); sf_putstring (sembFile, "label3", "Crossline");
	sf_putstring(sembFile, "unit1", "ms"); sf_putstring (sembFile, "unit2", "m"); sf_putstring (sembFile, "unit3", "m");
    }

    // dip-angle gather size
    int dagSize = gp.zNum * gp.dipNum * gp.sdipNum;
    const int sectionSize = dp.zNum * dp.xNum * dp.yNum;

    // velocity trace
    float* velTrace  = sf_floatalloc (vp.zNum);
    // common-offset section
    float* offsetSection = sf_floatalloc (sectionSize);

    // current offset gather
    float* coGather  = sf_floatalloc (dagSize);
    // current offset image
    float* coImage   = sf_floatalloc (ip.zNum);
    // current offset stack of amplitude squares
    float* coImageSq = sf_floatalloc (ip.zNum);
    // full dip-angle gather
    float* mainGather  = sf_floatalloc (dagSize);
    // full image
    float* mainImage   = sf_floatalloc (ip.zNum);
    // full stack of amplitude squares
    float* mainImageSq = sf_floatalloc (ip.zNum);

    // set migrator
    TimeMigratorBase* migrator;
    if (rp.is3D) {
		migrator = new TimeMigrator3D ();
		migrator->initCurveDefiner (true);
    } else {
		migrator = new TimeMigrator2D ();
		migrator->initCurveDefiner (false);
    }
    
	migrator->setImagingParams (&dp, offsetSection, rp.isAA, axis2label, &vp, &ip, &gp);
    migrator->setDataLimits ();

	migrator->getStackTaper (rp.edgeTaper, rp.isDipAz);

    const int fullGatherNum = ip.yNum * ip.xNum;

	int cigind = 1;

    for (int ipy = 0; ipy < ip.yNum; ++ipy) {
		for (int ipx = 0; ipx < ip.xNum; ++ipx, ++cigind) {

	    sf_warning ("gather %d of %d;", cigind, fullGatherNum);	

	    memset (mainGather,   0, dagSize * sizeof (float));
	    memset (mainImage,    0, ip.zNum * sizeof (float));
	    memset (mainImageSq,  0, ip.zNum * sizeof (float));

	    prepareVelocityTrace (ipx, ipy, velTrace);
			
	    for (int ih = 0; ih < rp.hMigNum; ++ih) {

		float curOffset = dp.hStart + ih * dp.hStep;

		memset (coGather,   0, dagSize * sizeof (float));
		memset (coImage,    0, ip.zNum * sizeof (float));
		memset (coImageSq,  0, ip.zNum * sizeof (float));
	
		readOffsetSection (ih, offsetSection);

		Point2D curGather (ipx, ipy);
		migrator->processGather (curGather, curOffset, velTrace, rp.isDipAz, coGather, coImage, coImageSq);

		// add migrated trace to image
		float* ptrToCur  = coImage;	
		float* ptrToMain = mainImage;
		const int size = gp.zNum;// * gp.sdipNum;
		for (int it = 0; it < size; ++it, ++ptrToCur, ++ptrToMain)
		    *ptrToMain += *ptrToCur;

		if (rp.isCig) {
		    // write trace into offset-domain CIG file
		    size_t startPos = (ipx * rp.hMigNum + ih) * ip.zNum * sizeof(float);
		    sf_seek (cigFile, startPos, SEEK_SET);
		    sf_floatwrite (coImage, ip.zNum, cigFile);
		}

		if (rp.isDag) {
		    // add migrated trace to dip-angle gather
		    ptrToCur  = coGather;	
		    ptrToMain = mainGather;
		    for (int it = 0; it < dagSize; ++it, ++ptrToCur, ++ptrToMain)
			*ptrToMain += *ptrToCur;
		}
				
		if (rp.isSemb) {				
		    ptrToCur  = coImageSq;	
		    ptrToMain = mainImageSq;
		    const int size = gp.zNum;// * gp.sdipNum;
		    for (int it = 0; it < size; ++it, ++ptrToCur, ++ptrToMain)
			*ptrToMain += *ptrToCur;
		}
	    }

	    size_t startInd = (ipx + ipy * ip.xNum) * sizeof(float);
		size_t shift = startInd * ip.zNum;
	    sf_seek (imageFile, shift, SEEK_SET);
	    sf_floatwrite (mainImage, ip.zNum, imageFile);

		shift = startInd * dagSize;
	    if (rp.isDag) {
			sf_seek (dagFile, shift, SEEK_SET);
			sf_floatwrite (mainGather, dagSize, dagFile);
	    }			

	    if (rp.isSemb) {
		float* sembTrace = sf_floatalloc (ip.zNum);
		Sembler::getSemblanceForTrace (gp.dipNum * gp.sdipNum, mainImage, mainImageSq, ip.zNum, rp.sembWindow, sembTrace, rp.hMigNum);
		size_t startInd = (ipx + ipy * ip.xNum) * sizeof(float);
		size_t shift = startInd * ip.zNum;
		sf_seek (sembFile, shift, SEEK_SET);
		sf_floatwrite (sembTrace, ip.zNum, sembFile);
		free (sembTrace);
	    }
	}
    }

    free (coGather);
    free (coImage);
    free (coImageSq);

    free (mainGather);
    free (mainImage);

    free (velTrace);
    free (offsetSection);

    sf_fileclose (dagFile);
    sf_fileclose (imageFile);
    sf_fileclose (cigFile);

    sf_fileclose (dataFile);
    sf_fileclose (velFile);

    return 0;
}
