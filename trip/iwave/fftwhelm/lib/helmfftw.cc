#include "helmfftw.hh"
/* Helm solver by FFTW */
/* Yin Huang, Oct 6, 2015 */

namespace TSOpt {

  using RVL::ScalarFieldTraits;
  using RVL::SpaceTest;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::Space;
  using RVL::ProductSpace;
  using RVL::Vector;
  using RVL::Components;
  using RVL::ProtectedDivision;
  using RVL::RnArray;
  using RVL::RVLScale;
  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;
  using RVL::MPISerialFunctionObject;
  using RVL::MPISerialFunctionObjectRedn;
    
#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif

  void HelmFFTWFO::operator()(LocalDataContainer<ireal> & x,
			  LocalDataContainer<ireal> const & y){
    try{
        
      float * indata  = NULL;
      float * work    = NULL;
      float * work1    = NULL;
      float * outdata = NULL;
      fftwf_plan cfg=NULL, icfg=NULL;
        
      IPNT f2c;  // #grid for each axis for FFTW
        
      ContentPackage<ireal, RARR>  & gx =
        dynamic_cast<ContentPackage <ireal, RARR>  &> (x);
      ContentPackage<ireal, RARR> const & gy =
        dynamic_cast<ContentPackage <ireal, RARR> const &> (y);
        
      // precondition - metadata are same dimn
      RARR  & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      //int idxdatum = (int)(datum/d_arr[0]+0.5f);
      int lendom;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);

      if (dimx != dimy) {
          RVLException e;
          e<<"Error: HelmFFTWFO::operator()\n";
          e<<"arguments have different dims:\n";
          e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
          throw e;
      }
        
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
        
      for (int ii=0;ii<dimx;ii++)  {
          s[ii]=max(gsy[ii],gsx[ii]);
          e[ii]=min(gey[ii],gex[ii]);
      }
        //cerr << " idxdatum = " << idxdatum << endl;
        //cerr << " dimx = " << dimx << endl;
        
        lendom = 1;
        for (int ii=0; ii<dimx; ii++) {
            f2c[ii] = n_arr[ii];//(n_arr[ii]+4)%2==1? (n_arr[ii]+5):(n_arr[ii]+4);
            lendom=lendom*f2c[ii];
            //cerr << "f2c["<< ii <<"] = " << f2c[ii] << endl;
            //cerr << "n_arr["<< ii << "] = " << n_arr[ii] << endl;
            //cerr << "length " << ii << " = " << e[ii] - s[ii] +1<< endl;
        }

        float _power=power;
        fftwf_r2r_kind bcf[2][2];
        fftwf_r2r_kind bcb[2][2];

        bcf[0][0] = FFTW_REDFT10;
        bcf[0][1] = FFTW_REDFT11;
        bcf[1][0] = FFTW_RODFT11;
        bcf[1][1] = FFTW_RODFT10;

        bcb[0][0] = FFTW_REDFT01;
        bcb[0][1] = FFTW_REDFT11;
        bcb[1][0] = FFTW_RODFT11;
        bcb[1][1] = FFTW_RODFT01;
        
        //cerr << "lendom = " << lendom<< endl;
        
        IPNT i;
        
        // allocate data arrays
        if (!(indata  = (float *)fftwf_malloc(lendom*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmFFTWOp::apply - failed to allocate " << lendom << " fftwf_complex for input data\n";
            throw e;
        }
        if (!(work    = (float *)fftwf_malloc(lendom*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmFFTWOp::apply - failed to allocate " << lendom << " fftwf_complex for working space\n";
            throw e;
        }
        if (!(work1   = (float *)fftwf_malloc(lendom*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmFFTWOp::apply - failed to allocate " << lendom << " fftwf_complex for working space\n";
            throw e;
        }
        if (!(outdata = (float *)fftwf_malloc(lendom*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmFFTWOp::apply - failed to allocate " << lendom << " fftwf_complex for output data\n";
            throw e;
        }
        
        for (i[1]=0;i[1]<lendom;i[1]++){
            indata[i[1]]=0.0f;
            outdata[i[1]]=0.0f;
        }
//        if (retrieveGlobalRank() == 0) {
//            cerr << "\n before HelmFFTWOp::apply \n";
//            cerr << "\n RARR_MAX_NDIM=" << RARR_MAX_NDIM << endl;
//            cerr << "\n dimx=" << dimx << endl;
//        }        
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
          for (i[0]=s[0];i[0]<=e[0];i[0]++) {
              indata[i[0]-s[0]]=ray._s1[i[0]];
          }
          
          // forward Fourier transform
          if (cfg==NULL) {
              //cfg = fftwf_plan_r2r_1d(f2c[0], indata,  work, FFTW_RODFT01, FFTW_MEASURE);
              cfg = fftwf_plan_r2r_1d(f2c[0], indata,  work, bcf[sbc[0]][ebc[0]], FFTW_MEASURE);
              if (cfg==NULL) fprintf(stderr,"FFTW failure.\n");
          }
          fftwf_execute(cfg);
          
          if (icfg==NULL){
              //icfg = fftwf_plan_r2r_1d( f2c[0], work1, outdata, FFTW_RODFT10, FFTW_MEASURE);
              icfg = fftwf_plan_r2r_1d( f2c[0], work1, outdata, bcb[sbc[0]][ebc[0]], FFTW_MEASURE);
              if (icfg==NULL) fprintf(stderr,"FFTW failure.\n");
          }
          
          float wtz = 2*M_PI/(2*f2c[0]*d_arr[0])*weights[0];
          wtz = wtz * wtz;
          float wz;
          for (i[0]=0; i[0]<f2c[0]; i[0]++) {
              wz = i[0] * i[0];
              work1[i[0]]=work[i[0]]*(pow(1. + wtz * wz, _power));
          }

          fftwf_execute(icfg);
          float wt =  1.0/(2*f2c[0]);

          // copy data back
          for (i[0]=s[0];i[0]<=e[0];i[0]++) {
              rax._s1[i[0]]=wt*outdata[i[0]-s[0]];
          }
      }
#endif
#if RARR_MAX_NDIM > 1
        if (dimx==2) {
            for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                    indata[(i[1]-s[1])*f2c[0] + i[0]-s[0]]=ray._s2[i[1]][i[0]];
                }
            }
            
            // forward Fourier transform
            if (cfg==NULL) {
                cfg = fftwf_plan_r2r_2d(f2c[1], f2c[0], indata, work, bcf[sbc[1]][ebc[1]], bcf[sbc[0]][ebc[0]], FFTW_MEASURE);
                if (cfg==NULL) fprintf(stderr,"FFTW failure.\n");
            }
            fftwf_execute(cfg);
            // inverse Fourier transform
            if (icfg==NULL){
                icfg = fftwf_plan_r2r_2d(f2c[1], f2c[0], work1, outdata, bcb[sbc[1]][ebc[1]], bcb[sbc[0]][ebc[0]], FFTW_MEASURE);
                if (icfg==NULL) fprintf(stderr,"FFTW failure.\n");
            }
            
            float wtx = 2*M_PI/(2*f2c[1]*d_arr[1])*weights[1];
            wtx = wtx * wtx;
            float wtz = 2*M_PI/(2*f2c[0]*d_arr[0])*weights[0];
            wtz = wtz * wtz;
            float wx, wz;
            
            for (i[1]=0; i[1]<f2c[1]; i[1]++) {
                wx = i[1] * i[1];
                for (i[0]=0; i[0]<f2c[0]; i[0]++) {
                    wz = i[0] * i[0];
                    work1[i[1]*f2c[0]+i[0]]=work[i[1]*f2c[0]+i[0]]*(pow(1.+wtx * wx + wtz * wz, _power));
                }
            }
            
            fftwf_execute(icfg);
            float wt =  1.0/(2*(f2c[1])*2*f2c[0]);
            // copy data back
            for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                    rax._s2[i[1]][i[0]]=wt*outdata[(i[1]-s[1])*f2c[0]+i[0]-s[0]];
                }
            }
        }
#endif
#if RARR_MAX_NDIM > 2
        if (dimx==3) {
            for (i[2]=s[2];i[2]<=e[2];i[2]++) {
                for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                        indata[((i[2]-s[2])*f2c[1]+(i[1]-s[1]))*f2c[0]+i[0]-s[0]]=ray._s3[i[2]][i[1]][i[0]];
                    }
                }
            }
            // forward Fourier transform
            if (cfg==NULL) {
                cfg = fftwf_plan_r2r_3d(f2c[2], f2c[1], f2c[0], indata, work, bcf[sbc[2]][ebc[2]], bcf[sbc[1]][ebc[1]], bcf[sbc[0]][ebc[0]],FFTW_MEASURE);
                if (cfg==NULL) fprintf(stderr,"FFTW failure.\n");
            }
            fftwf_execute(cfg);
            // inverse Fourier transform
            if (icfg==NULL){
                icfg = fftwf_plan_r2r_3d(f2c[2], f2c[1], f2c[0], work1, outdata, bcb[sbc[2]][ebc[2]], bcb[sbc[1]][ebc[1]], bcb[sbc[0]][ebc[0]], FFTW_MEASURE);
                if (icfg==NULL) fprintf(stderr,"FFTW failure.\n");
            }
            
            float wty = 2*M_PI/(2*f2c[2]*d_arr[2])*weights[2];
            wty = wty * wty;
            float wtx = 2*M_PI/(2*f2c[1]*d_arr[1])*weights[1];
            wtx = wtx * wtx;
            float wtz = 2*M_PI/(2*f2c[0]*d_arr[0])*weights[0];
            wtz = wtz * wtz;
            float wy, wx, wz;
  
            for (i[2]=0;i[2]<f2c[2];i[2]++) {
                wy = i[2] * i[2];
                for (i[1]=0;i[1]<f2c[1];i[1]++) {
                    wx = i[1] * i[1];
                    for (i[0]=0;i[0]<f2c[0];i[0]++) {
                        wz = i[0] * i[0];
                        work1[(i[2]*f2c[1]+i[1])*f2c[0]+i[0]]=work[(i[2]*f2c[1]+i[1])*f2c[0]+i[0]]*
                        (pow(1. + wty*wy + wtx*wx + wtz*wz, power));
                    }
                }
            }
            
            fftwf_execute(icfg);
            float wt =  1.0/(2*f2c[2]*2*f2c[1]*2*f2c[0]);
            
            // copy data back

            for (i[2]=s[2];i[2]<=e[2];i[2]++) {
                for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                    for (i[0]=s[0]; i[0]<=e[0];i[0]++) {
                        rax._s3[i[2]][i[1]][i[0]]=wt*outdata[((i[2]-s[2])*f2c[1]+(i[1]-s[1]))*f2c[0]+i[0]-s[0]];
                    }
                }
            }
      }
#endif
        fftwf_destroy_plan(icfg);
        fftwf_destroy_plan(cfg);
        
        fftwf_free(indata);
        fftwf_free(outdata);
        fftwf_free(work);
        fftwf_free(work1);
        
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: HelmFFTWFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: HelmFFTWFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from HelmFFTWFO::operator()\n";
      throw e;
    }
        
  }
    

  void GridHelmFFTWOp::apply(const Vector<float> & x,
			 Vector<float> & y) const {
    try {
      // extract components - fine even if only one!
      Components<float> cx(x);
      Components<float> cy(y);

      // detect product structure
      ProductSpace<ireal> const * pdom = NULL;
      pdom = dynamic_cast<ProductSpace<ireal> const *>(&dom);
      int n_fac=1;
      if (pdom) n_fac=pdom->getSize();
      Space<ireal> const * sp = NULL;
   
      // component loop
      for (int j=0; j<n_fac; j++) {
	if (pdom) sp = &((*pdom)[j]);
	else sp = &dom;

	// late tests
	myGridSpace const * gdom = dynamic_cast<myGridSpace const *>(sp);
	if (!gdom) {
	  RVLException e;
	  e<<"Error: GridHelmFFTWOp::apply\n";
	  e<<"  factor "<<j<<" of input space is not a GridSpace\n";
	  e<<"  description:\n";
	  sp->write(e);
	  throw e;	  
	}

	IPNT n_arr;
	RPNT d_arr;
        if (retrieveGlobalRank() == 0) {
	  if (gdom->getGrid().dim != 1 && gdom->getGrid().dim != 2 && gdom->getGrid().dim != 3) {
	    RVLException e;
	    e<<"Error: GridHelmFFTWOp::apply\n";
	    e<<"  current implementation is 2D and 3D only\n";
	    throw e;
	  }
	  get_d(d_arr,gdom->getGrid());
	  get_n(n_arr,gdom->getGrid());
	}
//        if (retrieveGlobalRank() == 0) {
//                cerr << "\n In GridHelmFFTWOp::apply \n";
//            }
	HelmFFTWFO fo(n_arr,d_arr,weights,sbc,ebc,power,datum);
	MPISerialFunctionObject<float> mpifo(fo);
	cy[j].eval(mpifo,cx[j]);    
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled in GridHelmFFTWOp::apply\n";
      throw e;
    }
            
  }
        
  void GridHelmFFTWOp::applyAdj(const Vector<float> & x,
			    Vector<float> & y) const {
    try {
      apply(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled in GridHelmFFTWOp::applyAdj\n";
      throw e;
    }
  }

}
