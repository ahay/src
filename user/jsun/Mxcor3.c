/* OpenMP time- or freq-domain cross-correlation on axes 1,2,3,4 */

/*
  Copyright (C) 2007 Colorado School of Mines
  Copyright (C) 2016 The University of Texas at Austin
  
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

#define rCOR(a,b) (a*b)

#ifdef SF_HAS_COMPLEX_H
#define cCOR(a,b) (crealf(conjf(a)*b))
#else
#define cCOR(a,b) (crealf(sf_cmul(conjf(a),b)))
#endif

int main(int argc, char* argv[])
{
    bool verb,rflg;
    int  axis;

    sf_file Fs,Fr,Fi;    /* I/O files -> source, receiver, image */
    sf_axis a1,a2,a3,a4,aa; /* cube axes */
    int     i1,i2,i3,i4;
    int     n1,n2,n3,n4;

    int nbuf,ibuf;

    float      ****r_us=NULL,****r_ur=NULL,***ii=NULL;
    sf_complex ****c_us=NULL,****c_ur=NULL;

    float scale;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif
    
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint ("axis",&axis)) axis=2;     /* stack axis */
    if(! sf_getint ("nbuf",&nbuf)) nbuf=1;     /* buffer size */

    sf_warning("cross-correlation on axis %d",axis);

    Fs = sf_input ("in" );
    Fr = sf_input ("uu" );
    Fi = sf_output("out");

    rflg = (sf_gettype(Fs) == SF_COMPLEX)?false:true;
    if(rflg) {
	sf_warning("real input");
    } else {
	sf_warning("complex input");
	sf_settype(Fi,SF_FLOAT);
    }

    /* read axes */
    a1=sf_iaxa(Fs,1); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fs,2); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fs,3); if(verb) sf_raxa(a3);
    a4=sf_iaxa(Fs,4); if(verb) sf_raxa(a4);

    aa=sf_maxa(1,0,1); /*length,origin,sampling*/ 
    sf_setlabel(aa,""); 
    sf_setunit (aa,""); 

    n1 = sf_n(a1);
    n2 = sf_n(a2);
    n3 = sf_n(a3);
    n4 = sf_n(a4);

    nbuf = SF_MIN(nbuf,n4);

    switch(axis) {
        case 4:
	    sf_oaxa(Fi,a1,1);
	    sf_oaxa(Fi,a2,2);
	    sf_oaxa(Fi,a3,3);
	    sf_oaxa(Fi,aa,4);
	    ii=sf_floatalloc3(n1,n2,n3); 
	    scale = 1./n4;
	    break;
	case 3:
	    sf_oaxa(Fi,a1,1);
	    sf_oaxa(Fi,a2,2);
	    sf_oaxa(Fi,a4,3);
	    sf_oaxa(Fi,aa,4);
	    ii=sf_floatalloc3(n1,n2,nbuf); 
	    scale = 1./n3;
	    break;
	case 2:
	    sf_oaxa(Fi,a1,1);
	    sf_oaxa(Fi,a3,2);
	    sf_oaxa(Fi,a4,3);
	    sf_oaxa(Fi,aa,4);
	    ii=sf_floatalloc3(n1,n3,nbuf); 
	    scale = 1./n2;
	    break;
	case 1:
	default:
	    sf_oaxa(Fi,a2,1);
	    sf_oaxa(Fi,a3,2);
	    sf_oaxa(Fi,a4,3);
	    sf_oaxa(Fi,aa,4);
	    ii=sf_floatalloc3(n2,n3,nbuf); 
	    scale = 1./n1;
	    break;
    }

    if(axis==4) {
        for        (i3=0; i3<n3; i3++) {
            for    (i2=0; i2<n2; i2++) {
                for(i1=0; i1<n1; i1++) {
                    ii[i3][i2][i1]=0.;
                }
            }
        }
    }

    if(rflg) {
	r_us = sf_floatalloc4  (n1,n2,n3,nbuf);
	r_ur = sf_floatalloc4  (n1,n2,n3,nbuf);
    } else {
	c_us = sf_complexalloc4(n1,n2,n3,nbuf);
	c_ur = sf_complexalloc4(n1,n2,n3,nbuf);
    }

    i4=n4;
    for (; i4 > 0; i4 -= nbuf) {
	if (nbuf > i4) nbuf=i4;
	if(verb) sf_warning("nsiz=%ld nbuf=%ld",i4,nbuf);

	if(rflg) {
	    sf_floatread  (r_us[0][0][0],n1*n2*n3*nbuf,Fs);
	    sf_floatread  (r_ur[0][0][0],n1*n2*n3*nbuf,Fr);
	} else {
	    sf_complexread(c_us[0][0][0],n1*n2*n3*nbuf,Fs);
	    sf_complexread(c_ur[0][0][0],n1*n2*n3*nbuf,Fr);
	}

	switch(axis) {

	    /*------------------------------------------------------------*/
	    case 4:
		
		if(rflg) {
		    
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(i1,i2,i3)				\
    shared( n1,n2,n3,ibuf,ii,r_us,r_ur)
#endif
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[i3][i2][i1] += rCOR(r_us[ibuf][i3][i2][i1],r_ur[ibuf][i3][i2][i1]);
                                }
                            }
			}
		    } /* ibuf */
		    
		} else {
		    
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(i1,i2,i3)				\
    shared( n1,n2,n3,ibuf,ii,c_us,c_ur)
#endif
                        for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[i3][i2][i1] += cCOR(c_us[ibuf][i3][i2][i1],c_ur[ibuf][i3][i2][i1]);
                                }
                            }
			}
		    } /* ibuf */
		    
		}
		break;				
		/*------------------------------------------------------------*/
		    
	    case 3:

		if(rflg) {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2,i3)			   \
    shared( nbuf,n1,n2,n3,ii,r_us,r_ur,scale)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i2=0; i2<n2; i2++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i2][i1]=0;
                            }
                        }
			
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[ibuf][i2][i1] += rCOR(r_us[ibuf][i3][i2][i1],r_ur[ibuf][i3][i2][i1]);
                                }
			    }
			}
			
			for    (i2=0; i2<n2; i2++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i2][i1] *= scale;
                            }
			}
		    } /* ibuf */
		    
		} else {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2,i3)			   \
    shared( nbuf,n1,n2,n3,ii,c_us,c_ur,scale)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i2=0; i2<n2; i2++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i2][i1]=0;
                            }
                        }
			
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[ibuf][i2][i1] += cCOR(c_us[ibuf][i3][i2][i1],c_ur[ibuf][i3][i2][i1]);
                                }
			    }
			}
			
			for    (i2=0; i2<n2; i2++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i2][i1] *= scale;
                            }
			}
		    } /* ibuf */
		    
		}
		
		sf_floatwrite(ii[0][0],n1*n2*nbuf,Fi);
		break;
		/*------------------------------------------------------------*/
		
	    case 2:

		if(rflg) {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2,i3)			   \
    shared( nbuf,n1,n2,n3,ii,r_us,r_ur,scale)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i3=0; i3<n3; i3++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i3][i1]=0;
                            }
                        }
			
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[ibuf][i3][i1] += rCOR(r_us[ibuf][i3][i2][i1],r_ur[ibuf][i3][i2][i1]);
                                }
			    }
			}
			
			for    (i3=0; i3<n3; i3++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i3][i1] *= scale;
                            }
			}
		    } /* ibuf */
		    
		} else {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2,i3)			   \
    shared( nbuf,n1,n2,n3,ii,c_us,c_ur,scale)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i3=0; i3<n3; i3++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i3][i1]=0;
                            }
                        }
			
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[ibuf][i3][i1] += cCOR(c_us[ibuf][i3][i2][i1],c_ur[ibuf][i3][i2][i1]);
                                }
			    }
			}
			
			for    (i3=0; i3<n3; i3++) {
                            for(i1=0; i1<n1; i1++) {
                                ii[ibuf][i3][i1] *= scale;
                            }
			}
		    } /* ibuf */
		    
		}
		
		sf_floatwrite(ii[0][0],n1*n3*nbuf,Fi);
		break;
		/*------------------------------------------------------------*/		
	    case 1:
	    default:

		if(rflg) {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2,i3)			   \
    shared( nbuf,n1,n2,n3,ii,r_us,r_ur,scale)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i3=0; i3<n3; i3++) {
                            for(i2=0; i2<n2; i2++) {
                                ii[ibuf][i3][i2]=0;
                            }
                        }
			
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[ibuf][i3][i2] += rCOR(r_us[ibuf][i3][i2][i1],r_ur[ibuf][i3][i2][i1]);
                                }
			    }
			}
			
			for    (i3=0; i3<n3; i3++) {
                            for(i2=0; i2<n2; i2++) {
                                ii[ibuf][i3][i2] *= scale;
                            }
			}
		    } /* ibuf */
		    
		} else {

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2,i3)			   \
    shared( nbuf,n1,n2,n3,ii,c_us,c_ur,scale)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i3=0; i3<n3; i3++) {
                            for(i2=0; i2<n2; i2++) {
                                ii[ibuf][i3][i2]=0;
                            }
                        }
			
			for        (i3=0; i3<n3; i3++) {
                            for    (i2=0; i2<n2; i2++) {
                                for(i1=0; i1<n1; i1++) {
                                    ii[ibuf][i3][i2] += cCOR(c_us[ibuf][i3][i2][i1],c_ur[ibuf][i3][i2][i1]);
                                }
			    }
			}
			
			for    (i3=0; i3<n3; i3++) {
                            for(i2=0; i2<n2; i2++) {
                                ii[ibuf][i3][i2] *= scale;
                            }
			}
		    } /* ibuf */
		    
		}
		
		sf_floatwrite(ii[0][0],n2*n3*nbuf,Fi);
		break;
		
	} /* n3 */
	
    }
    if(verb) fprintf(stderr,"\n");    
    
    if(axis==4) {
	for        (i3=0; i3<n3; i3++) {
            for    (i2=0; i2<n2; i2++) {
                for(i1=0; i1<n1; i1++) {
                    ii[i3][i2][i1] *=scale;
                }
            }
        }
	sf_floatwrite(ii[0][0],n3*n2*n1,Fi);
    }
    
    /*------------------------------------------------------------*/
    free(**ii); free(*ii); free(ii);
    
    if(rflg) {
	free(***r_us); free(**r_us); free(*r_us); free(r_us);
	free(***r_ur); free(**r_ur); free(*r_ur); free(r_ur);
    } else {
	free(***c_us); free(**c_us); free(*c_us); free(c_us);
	free(***c_ur); free(**c_ur); free(*c_ur); free(c_ur);
    }
    /*------------------------------------------------------------*/
    
    exit (0);
}
