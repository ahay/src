/* Convert v(x0,t0) to v(x,z) */
/*
  Copyright (C) 2008 New York University
  
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

#include "fastmarch.h"

#define TOL 1.0e-8
#define NTRYMAX 5

static char solve(int ind,char dir,int b);
static void addtree(int ind); /* adds a node to the binary tree
				  of the "considered" points */ 
static void updtree(int ind); /* updates the binary tree */
static void deltree(void); /* delets the root of the binary tree */
static void checktree(int ind,const char *str);
static char solve(int ind,char dir,int b);
static char solve1pt(int ind,int inda,float ha,char ch);
static char solve2pt(int ind,int inda,int indb,float ha,float hb,char ch);
static float linterp(float f0,float f1,float delta,float x);
static float linterp1(int i,int k,float x);
static float linterp2(float x,float t);
static float fn(float t,float ta,float tb,float xa,float xb,
		 float ha,float hb);
static float lateral(float t,float ua,float ub,float ta,float tb,
		      float ha,float hb);
static float dfn(float t,float ta,float tb,float xa,float xb,
		  float ha,float hb);

static int nx, nz, nt, nx1, nt1;
static float hx, hz, ht;
static int *pup;
static char *ms;
static int count;
static float *x0, *t0, *v, *s;
static int *r; /* r(position in the binary tree)=index of the mesh point (x,y,z) */
static int *p; /* p(index)=position in the binary tree */

void fastmarch_init(int nx_init, int nz_init, int nt_init,
		    float hx_init, float hz_init, float ht_init,
		    float *x0_init, float *t0_init, float *v_init, float *s_init)
/*< initialize >*/
{
    int i,j,ind, nxz;

    nx = nx_init;
    nz = nz_init;
    nt = nt_init;
    nt1 = nt-1;
    nx1 = nx-1;

    hx = hx_init;
    hz = hz_init;
    ht = ht_init;

    x0 = x0_init;
    t0 = t0_init;

    v = v_init;
    s = s_init;

    nxz = nx*nz;

    pup= sf_intalloc(nxz);
    ms= sf_charalloc(nxz);
    r= sf_intalloc(nxz);
    p= sf_intalloc(nxz);

    count=0;

    for( i=0;i<nx;i++ ) {
	*(pup+i)=2;
	*(ms+i)='a';
	ind=i;
	for( j=1;j<nz;j++ ) {
	    ind+=nx;
	    *(pup+ind)=0;
	    *(ms+ind)='u';
	}
    }

}

void fastmarch_close(void)
/*< free allocated storage >*/
{
    free(pup);
    free(ms);
    free(r);
    free(p);
}

void fastmarch(void) 
/*< fast marching time-to-depth conversion >*/
{
  int i,j,ii,jj,ind,m,mm;
  char ch,dir;

  for( i=0;i<nx;i++ ) {
      ind=i+nx;
      ch=solve(ind,'z',1);
      if( ch=='s') {
	  addtree(ind);
	  *(ms+ind)='c';
      }
  }
  
  while( count > 0 ) {
    ind=(*(r+1));
    j=ind/nx;
    i=ind%nx;
    *(ms+ind)='a';
    deltree();
    /*(      printf("accept %i, %i\t",i,j);*/
    /* update the nearest neighbors */
    for( m=0; m<4; m++ ) {
      ch='n';
      dir=(m<2) ? 'x' : 'z';
      mm=2*(m%2)-1;
      if( dir=='x' ) {
	ii=i+mm;
	jj=j;
	if( ii>=0 && ii<nx ) ch='y';
      }
      else if( dir=='z' ) {
	ii=i;
	jj=j+mm;
	if( jj>=0 && jj<nz ) ch='y';
      }
      if( ch=='y' ) {
	ind=ii+jj*nx;
	if( *(ms+ind)=='a' ) ch='n';
	else {
	  ch=solve(ind,dir,mm);
	  if( ch=='s' ) {
	    if( *(ms+ind)=='u' ) {
	      addtree(ind);
	      *(ms+ind)='c';
	    }
	    else updtree(ind);
	    /* printf("(%i,%i) is updated; t0=%.4e\tx0=%.4e\tv=%.4e\t q=%.4e\n",ii,jj,*(t0+ind),*(x0+ind),*(v+ind),*(q+ind));*/
	  }
	}
      }
    }
  }
}



/********************************************************/

static char solve(int ind,char dir,int b) 
{
    int inda=0,indb=0,m0,i,j,m;
    char cha='n',chb='n',ch='n';

    j=ind/nx;
    i=ind%nx;
    if( dir=='x') {inda=ind-b,cha='y';}
    else if( dir=='z' ) {indb=ind-b*nx; chb='y'; }
    for( m0=-1; m0<=1; m0+=2 ) {
	if( dir=='x' ) {
	    m=j+m0;
	    if( m>=0 && m<nz ) {
		indb=ind+m0*nx;
		chb=( *(ms+indb)=='a' ) ? 'y' : 'n';
	    }
	    else chb='n';
	}
	else if( dir=='z' ) {
	    m=i+m0;
	    if( m>=0 && m<nx ) {
		inda=ind+m0;
		cha=( *(ms+inda)=='a' ) ? 'y' : 'n';
	    }
	    else cha='n';
	}
	if( cha=='y' ) {
	    if( chb=='y' ) 
		ch=solve2pt(ind,inda,indb,hx,hz,ch);
	    else if( *(pup+ind)<=1 ) ch=solve1pt(ind,inda,hx,ch);
	}
	else {
	    if( chb=='y' && *(pup+ind)<=1  )
		ch=solve1pt(ind,indb,hz,ch);
	    else ch=(ch!='n') ? ch : 'x';
	}
    }
    /*  printf("Solve: (%i,%i), ch=%c,t=%.4e, x=%.4e, v=%.4e, q=%.4e, pup=%i\n",
	i,j,ch,*(t0+ind),*(x0+ind), *(v+ind), *(q+ind), *(pup+ind));*/
    return ch;
}


/*********************************************************/

static char solve1pt(int ind,int inda,float ha,char ch) 
{
  int i,k,k0;
  float xa,ta,t,st0,st1;
  char ch1='y';

  xa=*(x0+inda);
  ta=*(t0+inda);
  i=floor(xa/hx);
  k=floor(ta/ht);

  if( k<nt1 ) {
    k0=k;
    while( ch1=='y' ) {
      st0=linterp1(i,k0,xa);
      st1=linterp1(i,k0+1,xa);
      t=(ta/ha+(st0*(k0+1)-st1*k0))/(1.0/ha-(st1-st0)/ht);
      if( t>ta && t<=(k0+1)*ht ) {
	ch1='n';
	if( *(pup+ind)<=1 ) {
	  if( t<(*(t0+ind)) ) {
	    *(pup+ind)=1;
	    *(t0+ind)=t;
	    *(x0+ind)=xa;	    
	    *(v+ind)=1.0/linterp(st0,st1,ht,t);;
	    ch='s';
	  }
	  ch=(ch=='s') ? 's' : 'n';
	}
	ch=(ch=='s') ? 's' : 'n';
      }
      else {
	k0++;
	if( k0==nt1 ) {
	  ch1='n';
	  ch=(ch=='s') ? 's' :'f';
	}
      }
    }
  }
  else ch=(ch=='s') ? 's' :'f';
  return ch;
}

/************************************************************/

static char solve2pt(int ind,int inda,int indb,float ha,float hb,char ch) 
{
    int kmin /* ,i,k  ,i1,k1 */;
  float x1,x2;
  float xa,ta,xb,tb,x,t,tmin;
  float d,dp,dp2,df;
  char /* ch1='y', */ ch2='n';
  int ntry=0;

  xa=*(x0+inda);
  ta=*(t0+inda);
  xb=*(x0+indb);
  tb=*(t0+indb);
  tmin=SF_MAX(ta,tb)-TOL;
  kmin=floor(tmin/ht);
  x1=SF_MIN(xa,xb)-TOL;
  x2=SF_MAX(xa,xb)+TOL;
  t=*(t0+ind);
  if( kmin<nt1 ) {
    while( ch2=='n' && ntry<NTRYMAX ) {
      x=lateral(t,xa,xb,ta,tb,ha,hb);
      /*i=SF_MAX(0,SF_MIN(nx1-1,floor(x/hx)));
	k=SF_MAX(0,SF_MIN(nt1-1,floor(t/ht))); */
      d=fn(t,ta,tb,xa,xb,ha,hb);
      dp=2.0*d;
      dp2=dp;
      while( fabs(d)>TOL && fabs(dp2)>fabs(d) ) {
        df=dfn(t,ta,tb,xa,xb,ha,hb);
        t-=d/df;
        dp2=dp;
        dp=d;
        d=fn(t,ta,tb,xa,xb,ha,hb);
      } 
      x=lateral(t,xa,xb,ta,tb,ha,hb);
/*      i1=SF_MAX(0,SF_MIN(nx1-1,floor(x/hx)));
	k1=SF_MAX(0,SF_MIN(nt1-1,floor(t/ht))); */
      if( x>=x1 && x<=x2 && t>=tmin ) ch2='y';
      ntry++;
    }
    if( fabs(d)<=TOL  && t>=tmin  ) {
      if( x>=x1 && x<=x2 ) {
/*	ch1='n'; */
	if( *(pup+ind)<=1 || (*(pup+ind)==2 && t<(*(t0+ind))) ) {
	  ch='s';
	  *(t0+ind)=t;
	  *(x0+ind)=x;
	  *(pup+ind)=2;
	  *(v+ind)=1.0/linterp2(x,t);
	}
	ch=(ch=='s') ? 's' : 'n';
      }
    }
  }
  else
    ch=(ch=='s') ? 's' : 'f';
  return ch;
}

/*********************************************************/

static float linterp1(int i,int k,float x) 
{
  float s0,s1,u=x/hx;
  int ind;

  ind=(k<nt) ? i+nx*k : i+nt1*nx;
  s0=*(s+ind);
  s1=(i<nx1) ? *(s+ind+1) : s0;
  return (s1*(u-i)+s0*(i+1-u));
}

/******************************************************/

static float linterp2(float x,float t) 
{
    float ans,st0,st1;
    int i,k;

    i=SF_MAX(0,SF_MIN(nx1-1,floor(x/hx)));
    k=SF_MAX(0,SF_MIN(nt1-1,floor(t/ht)));
    if( k<nt1 ) {
	st0=linterp1(i,k,x);
	st1=linterp1(i,k+1,x);
	ans=linterp(st0,st1,ht,t);
    }
    else ans=linterp1(i,k,x);
    return ans;
}

/***********************************************************/

static float linterp(float f0,float f1,float delta,float x) 
{
    float u,u0,u1;
    
    u=x/delta;
    u0=(floor(u));
    u1=u0+1.0;
    
    return f0*(u1-u)+f1*(u-u0);
}

/***********************************************************/

static float fn(float t,float ta,float tb,float xa,float xb,
		float ha,float hb) {
    float dta,dtb,x,ss;
    
    dta=(t-ta)/ha;
    dtb=(t-tb)/hb;
    x=lateral(t,xa,xb,ta,tb,ha,hb);
    ss=linterp2(x,t);
    return (dta*dta+dtb*dtb)-ss*ss;
}
/*-------------------------------------------------------*/

static float lateral(float t,float ua,float ub,float ta,float tb,
	      float ha,float hb) 
{
    float dta,dtb;
    
    dta=(t-ta)/(ha*ha);
    dtb=(t-tb)/(hb*hb);
    return (ua*dta+ub*dtb)/(dta+dtb);
}

/*----------------------------------------------------------*/

static float dfn(float t,float ta,float tb,float xa,float xb,
	  float ha,float hb) {
    return (fn(t+TOL,ta,tb,xa,xb,ha,hb)-fn(t,ta,tb,xa,xb,ha,hb))/TOL;
}

/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

static void addtree(int ind) {
  int loc, ptemp;
  int indp, indc;
  char ch;

  count++;
  *(r+count)=ind;
  *(p+ind)=count;
  if( count > 1 ) {
    loc=count;
    indc=*(r+loc);
    indp=*(r+loc/2);
    ch=( (*(t0+indc)) < (*(t0+indp)) ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=*(p+indc);
      *(p+indc)=*(p+indp);
      *(r+loc/2)=indc;
      *(p+indp)=ptemp;
      *(r+loc)=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=*(r+loc);
        indp=*(r+loc/2);
        ch=( (*(t0+indc)) < (*(t0+indp)) ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }
  checktree(ind,"Error is detected in addtree ");
}

/*------------------------------------------------------------------*/

static void updtree(int mind) {
  int loc, lcc, indc, indp,ind, ic=0, ptemp, ic1=0, ic2=0;
  char chu,chd;

  loc=*(p+mind);
  indc=*(r+loc);
  indp=*(r+loc/2);
  chu=( (loc > 1) && ((*(t0+indc)) < (*(t0+indp))) ) ? 'y' : 'n';
  ind=indc;
  lcc=2*loc;

  if( lcc < count )  {
    ic1=*(r+lcc);
    ic2=*(r+lcc+1);

    if( (*(t0+ind)) > (*(t0+ic1)) || (*(t0+ind)) > (*(t0+ic2)) ) {
      if( (*(t0+ic1)) <= (*(t0+ic2)) )  {
        chd='l';
	    ic=ic1;
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc+=1;
      }
    }
    else chd='n';
  }
  else if( lcc == count ) {
    ic=*(r+lcc);

    if( (*(t0+ind)) > (*(t0+ic)) ) { chd='l';}
    else chd='n';
  }
  else chd='n';
  if( chu == 'y' && chd != 'n' ) {
      sf_warning("error detected by updtree, (%i, %i)",mind%nx,mind/nx);
      sf_warning("mind=%i, indc=%i, indp=%i, count=%i, loc=%i, lcc=%i", mind,indc,indp,count,loc,lcc);
      sf_error("t=%.4e, tp=%.4e, tcl=%.4e, tcr=%.4e, chu=%c, chd=%c",
	       *(t0+mind),*(t0+indp),*(t0+ic1),*(t0+ic2),chu,chd);
  }
  while( chu == 'y' ) {
    ptemp=*(p+indc);
    *(p+indc)=*(p+indp);
    *(r+loc/2)=indc;
    *(p+indp)=ptemp;
    *(r+loc)=indp;
    loc=loc/2;
    if( loc > 1 ) {
      indc=*(r+loc);
      indp=*(r+loc/2);
      chu=( (*(t0+indc)) < (*(t0+indp)) ) ? 'y' : 'n';
    }
    else chu='n';
  } /*end while( chu == 'y' ) */
  while( chd != 'n' ) {    
    ptemp=*(p+ind);
    *(p+ind)=*(p+ic);
    *(r+loc)=ic;
    *(p+ic)=ptemp;
    *(r+lcc)=ind;
    loc=lcc;
    lcc=2*loc;
    if( lcc < count )  {
      ic1=*(r+lcc);
      ic2=*(r+lcc+1);
      if( (*(t0+ind)) > (*(t0+ic1)) || (*(t0+ind)) > (*(t0+ic2)) ) {
        if( (*(t0+ic1)) <= (*(t0+ic2)) )  {
          chd='l';
	  ic=ic1;
        }
        else {
          chd='r';
	  ic=ic2;
	  lcc+=1;
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=*(r+lcc);
      if( (*(t0+ind)) > (*(t0+ic)) ) {chd='l';}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
  checktree(mind,"Error is detected in updtree ");
}

/*---------------------------------------------------------------------*/


/* deletes root of the binary tree */
static void deltree(void) 
{
  int loc, ptemp, ind, lcc, ic=0, ic1=0, ic2=0, mind;
  char chd, ch='n';;

  mind=*(r+1);
  /*
  if( mind==53619 ) {
    printf("Deleting root with index %i\n",mind);
    ch='y';
  }
  */
  *(p+(*(r+1)))=0;
  *(r+1)=*(r+count);
  *(p+(*(r+1)))=1;
  count--;
  loc=1;
  ind=*(r+1);
  lcc=2*loc;
  if(ch=='y') printf("parent: loc(%i)=%i, t=%.16e\n",ind,loc,*(t0+ind));
  if( lcc < count )  {
    ic1=*(r+lcc);
    ic2=*(r+lcc+1);
    if(ch=='y') printf("children: loc(%i)=%i, loc(%i)=%i,t1=%.16e, t2=%.16e\n",
		       ic1,lcc,ic2,lcc+1,*(t0+ic1),*(t0+ic2));
    if( (*(t0+ind)) > (*(t0+ic1)) || (*(t0+ind)) > (*(t0+ic2)) ) {
      if( (*(t0+ic1)) <= (*(t0+ic2)) )  {
        chd='l';
	    ic=ic1;
	    if(ch=='y') printf("left\n");
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc++;
	    if(ch=='y') printf("right\n");
      }
    }
    else chd='n';
  }
  else if( lcc == count ) {
    ic=*(r+lcc);
    if( (*(t0+ind)) > (*(t0+ic)) ) {chd='l';}
    else chd='n';
  }
  else chd='n';
  while( chd != 'n' ) {    
    ptemp=*(p+ind);
    *(p+ind)=*(p+ic);
    *(r+loc)=ic;
    *(p+ic)=ptemp;
    *(r+lcc)=ind;
    loc=lcc;
    lcc=2*loc;
    if(ch=='y') printf("parent: loc(%i)=%i, t=%.16e\n",ind,loc,*(t0+ind));
    if( lcc < count )  {
      ic1=*(r+lcc);
      ic2=*(r+lcc+1);
      if(ch=='y') printf("children: loc(%i)=%i, loc(%i)=%i,t1=%.16e, t2=%.16e\n",
		       ic1,lcc,ic2,lcc+1,*(t0+ic1),*(t0+ic2));
      if( (*(t0+ind)) > (*(t0+ic1)) || (*(t0+ind)) > (*(t0+ic2)) ) {
        if( (*(t0+ic1)) <= (*(t0+ic2)) )  {
          chd='l';
	      ic=ic1;
	      if(ch=='y') printf("left\n");
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
	    if(ch=='y') printf("right\n");
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=*(r+lcc);
      if( (*(t0+ind)) > (*(t0+ic)) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
  checktree(mind,"Error is detected in deltree ");
}

/*-------------------------------------------------------*/

static void checktree( int ind, const char * str ) 
{
  int k1, k2, k, i,i1=0,i2=0;

  for( k=1; k<=count; k++ ) {
    k1=k*2;
    k2=k1+1;
	i=*(r+k);
    if( k1 <= count ) {
      i1=*(r+k1);
      if( *(t0+i) > *(t0+i1) ) {
	    i=*(r+k);
	    sf_warning("%s",str);
	    sf_warning("at updating the point (%i,  %i)=%i",ind%nx,ind/nx,ind);
	    sf_warning(" parent: (%i, %i)=%i, %.8e",i%nx,i/nx,i, *(t0+i));
	    sf_error("children: %i, %.8e\t %i, %.8e",i1,*(t0+i1),i2,*(t0+i2));
      }
    }
    if( k2 <= count ) {
      i2=*(r+k2);
      if( *(t0+i) > *(t0+i2) ) {
	    i=*(r+k);
	    sf_warning("%s",str);
	    sf_warning("at updating the point (%i, %i)=%i",ind%nx,ind/nx,ind);
	    sf_warning("parent: (%i, %i)=%i, .%.8e",i%nx,i/nx,i,*(t0+i));
	    sf_error("children: %i, %.8e\t %i, %.8e",i1,*(t0+i1),i2,*(t0+i2));
      }
    }
  } 
}

/********************************************************/		    
