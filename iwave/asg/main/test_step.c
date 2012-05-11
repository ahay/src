/* Test step. */

#include <sgn.h>

// copy of private function from sgn_model.c
int asg_isdyn(int i) {

  // hard-wire!!!
  int m_ndim=3;

  if (m_ndim<4) {
    if (m_ndim>0) {
      if (i==D_P0) return 1;
      if (i==D_V0) return 1;
    }
    if (m_ndim>1) {
      if (i==D_P1) return 1;
      if (i==D_V1) return 1;
    }
    if (m_ndim>2) {
      if (i==D_P2) return 1;
      if (i==D_V2) return 1;
    }
    if (i==D_MV_DYN) return 1;
  }
  return 0;
}

// allocated RDOM - from 2x2x1 dom decomp of demo 11, dom 0
int create_alloc_rdom(RDOM * u) {

  int i;                  // counters
  int ndim=3;             // array dimension
  IPNT gs[RDOM_MAX_NARR]; // start indices
  IPNT n[RDOM_MAX_NARR];  // array lengths

  // initialize to null RDOM
  rd_a_setnull(u);

  // initialize gs, n arrays for zero dim
  for (i=0;i<RDOM_MAX_NARR;i++) {
    gs[i][0]=1; gs[i][1]=1; gs[i][2]=1;
    n[i][0]=0; n[i][1]=0; n[i][2]=0;
  }

  // assign gs and n arrays from example
  gs[0][0] = -4; gs[0][1] = -7; gs[0][2] =  1;
  gs[1][0] =  1; gs[1][1] = -7; gs[1][2] =  1;
  gs[2][0] = -4; gs[2][1] = -7; gs[2][2] =  1;
  gs[3][0] =  0; gs[3][1] = -7; gs[3][2] =  1;
  gs[4][0] =  1; gs[4][1] =-12; gs[4][2] =  1;
  gs[5][0] =  1; gs[5][1] = -8; gs[5][2] =  1;
  gs[6][0] =  1; gs[6][1] = -7; gs[6][2] = -4;
  gs[7][0] =  1; gs[7][1] = -7; gs[7][2] =  0;
  gs[8][0] =  1; gs[8][1] =  0; gs[8][2] =  0;
  gs[9][0] =  0; gs[9][1] =  0; gs[9][2] =  0;
  gs[10][0]=  1; gs[10][1]=-12; gs[10][2]=  1;
  gs[11][0]= -7; gs[11][1]=  0; gs[11][2]=  0;
  gs[12][0]= -8; gs[12][1]=  0; gs[12][2]=  0;
  gs[13][0]=  1; gs[13][1]= -7; gs[13][2]= -4;
  gs[14][0]=  1; gs[14][1]=  0; gs[14][2]=  0;
  gs[15][0]=  0; gs[15][1]=  0; gs[15][2]=  0;
  gs[16][0]= -4; gs[16][1]=-12; gs[16][2]= -4;

  n[0][0] =  58; n[0][1] = 202; n[0][2] = 389;
  n[1][0] =  48; n[1][1] = 202; n[1][2] = 389;
  n[2][0] =  57; n[2][1] = 202; n[2][2] = 389;
  n[3][0] =  49; n[3][1] = 202; n[3][2] = 389;
  n[4][0] =  48; n[4][1] = 211; n[4][2] = 389;
  n[5][0] =  48; n[5][1] = 203; n[5][2] = 389;
  n[6][0] =  48; n[6][1] = 202; n[6][2] = 398;
  n[7][0] =  48; n[7][1] = 202; n[7][2] = 390;
  n[8][0] =  48; n[8][1] =   1; n[8][2] =   1;
  n[9][0] =  49; n[9][1] =   1; n[9][2] =   1;
  n[10][0]=  48; n[10][1]= 212; n[10][2]= 389;
  n[11][0]= 202; n[11][1]=   1; n[11][2]=   1;
  n[12][0]= 203; n[12][1]=   1; n[12][2]=   1;
  n[13][0]=  48; n[13][1]= 202; n[13][2]= 399;
  n[14][0]= 389; n[14][1]=   1; n[14][2]=   1;
  n[15][0]= 390; n[15][1]=   1; n[15][2]=   1;
  n[16][0]=  58; n[16][1]= 212; n[16][2]= 399;

  // create RDOM
  return rd_a_create_s(u,RDOM_MAX_NARR,ndim,gs,n);
}

int create_comp_rdom(RDOM * u) {

  int i;                  // counters
  int err;                // error flag
  IPNT gs[RDOM_MAX_NARR]; // start indices
  IPNT n[RDOM_MAX_NARR];  // array lengths

  // initialize gs, n arrays for zero dim
  for (i=0;i<RDOM_MAX_NARR;i++) {
    gs[i][0]=1; gs[i][1]=1; gs[i][2]=1;
    n[i][0]=0; n[i][1]=0; n[i][2]=0;
  }

  // assign gs and n arrays from example
  gs[0][0] =  1; gs[0][1] = -7; gs[0][2] =  1;
  gs[1][0] =  1; gs[1][1] = -7; gs[1][2] =  1;
  gs[2][0] =  0; gs[2][1] = -7; gs[2][2] =  1;
  gs[3][0] =  0; gs[3][1] = -7; gs[3][2] =  1;
  gs[4][0] =  1; gs[4][1] = -8; gs[4][2] =  1;
  gs[5][0] =  1; gs[5][1] = -8; gs[5][2] =  1;
  gs[6][0] =  1; gs[6][1] = -7; gs[6][2] =  0;
  gs[7][0] =  1; gs[7][1] = -7; gs[7][2] =  0;
  gs[8][0] =  1; gs[8][1] =  0; gs[8][2] =  0;
  gs[9][0] =  0; gs[9][1] =  0; gs[9][2] =  0;
  gs[10][0]=  1; gs[10][1]= -7; gs[10][2]=  1;
  gs[11][0]= -7; gs[11][1]=  0; gs[11][2]=  0;
  gs[12][0]= -8; gs[12][1]=  0; gs[12][2]=  0;
  gs[13][0]=  1; gs[13][1]= -7; gs[13][2]=  1;
  gs[14][0]=  1; gs[14][1]=  0; gs[14][2]=  0;
  gs[15][0]=  0; gs[15][1]=  0; gs[15][2]=  0;
  gs[16][0]=  1; gs[16][1]= -7; gs[16][2]=  1;

  n[0][0] =  48; n[0][1] = 202; n[0][2] = 389;
  n[1][0] =  48; n[1][1] = 202; n[1][2] = 389;
  n[2][0] =  49; n[2][1] = 202; n[2][2] = 389;
  n[3][0] =  49; n[3][1] = 202; n[3][2] = 389;
  n[4][0] =  48; n[4][1] = 203; n[4][2] = 389;
  n[5][0] =  48; n[5][1] = 203; n[5][2] = 389;
  n[6][0] =  48; n[6][1] = 202; n[6][2] = 390;
  n[7][0] =  48; n[7][1] = 202; n[7][2] = 390;
  n[8][0] =  48; n[8][1] =   1; n[8][2] =   1;
  n[9][0] =  49; n[9][1] =   1; n[9][2] =   1;
  n[10][0]=  48; n[10][1]= 202; n[10][2]= 389;
  n[11][0]= 202; n[11][1]=   1; n[11][2]=   1;
  n[12][0]= 203; n[12][1]=   1; n[12][2]=   1;
  n[13][0]=  48; n[13][1]= 202; n[13][2]= 389;
  n[14][0]= 389; n[14][1]=   1; n[14][2]=   1;
  n[15][0]= 390; n[15][1]=   1; n[15][2]=   1;
  n[16][0]=  48; n[16][1]= 202; n[16][2]= 389;

  // reset RDOM
  for (i=0;i<RDOM_MAX_NARR;i++) {
    err=rd_greset_s(u,i,gs[i],n[i]);
    if (err != 0) return err;
  }

  return 0;
}

void assign_rand(RDOM * u) {

  int ndim = 3;
  int j;               // dim counter
  int i;               // grid counter
  int iarr;            // array counter
  int len;             // length of array
  ireal shift;         // shift away from zero, or to zero mean

  srand(19490615);

  // loop over domains
  for (iarr=0;iarr<RDOM_MAX_NARR;iarr++) {
    
    len=1;
    // extract dimensions
    for (j=0;j<ndim;j++) 
      len *= u->_s[iarr]._dims[j].n;

    // assign random numbers (1) between -0.5 and 0.5 with zero mean, for dynamic fields; 
    // (2) between 1.0 and 2.0 for static fields.
    shift=1.0;
    if (asg_isdyn(iarr)) shift = -0.5;
    for (i=0;i<len;i++) 
      u->_s[iarr]._s0[j]=shift + ((ireal)(rand()))/((ireal)(RAND_MAX));
  
  }
}
  
  
int main(int argc, char ** argv) {

  int err=0;
  RDOM rda;
  RDOM rdc;
  SGN_TS_PARS pars;

  // create allocated RDOM
  err = create_alloc_rdom(&rda);

  // dump result
  printf("Allocated ");
  if (!err) rd_a_dump(&rda,stdout);

  // create computational RDOM
  rdc=rda;
  if (!err) err = create_comp_rdom(&rdc);

  // dump result
  printf("Computational ");
  if (!err) rd_a_dump(&rdc,stdout);

  // assign random numbers to allocated rdom
  assign_rand(&rda);

  // dump RDOM - how?

  // assign SGN_TS_PARS data members - use typical values
  pars.dt = 1.0;
  pars.lam[0]=1.0;
  pars.lam[1]=1.0;
  pars.lam[2]=1.0;
  pars.k=5;
  pars.ndim=3;
  pars.psingle=0;
  pars.eflags[0]=1;
  pars.eflags[1]=1;
  pars.eflags[2]=1;

  // loop through updates - p, then v's
  
  sgn_ts3d_210(&rdc, D_P0, &pars);
  
  // dump after P0 update

  sgn_ts3d_210(&rdc, D_V0, &pars);

  // dump after V0 update

  sgn_ts3d_210(&rdc, D_V1, &pars);

  // dump after V1 update

  sgn_ts3d_210(&rdc, D_V2, &pars);

  // dump after V2 update

}
  
  
