/* memory copy in different cases */
/*
  Copyright (C) 2016 Yangkang Chen
*/

#include<rsf.h>
#include"memcpy.h"

/*real value functions*/
void mcp(float *dst /*destination*/, 
		float *src /*source*/,
		int s1d /*starting index in dst*/,
		int s1s /*starting index in src*/,
		int l1  /*copy length in axis1*/)
/*<memory copy in 1D case>*/
{
	int i1=0;
	for(i1=0;i1<l1;i1++)
	{
		dst[s1d+i1]=src[s1s+i1];
	}
}

void mcp_ad(float *dst /*destination*/, 
		float *src /*source*/,
		int s1d /*starting index in dst*/,
		int s1s /*starting index in src*/,
		int l1  /*copy length in axis1*/)
/*<memory copy and addition in 1D case>*/
{
	int i1=0;
	for(i1=0;i1<l1;i1++)
	{
		dst[s1d+i1]=src[s1s+i1] + dst[s1d+i1] ;
	}
}

void mcp2d(float **dst /*destination*/, 
		float **src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,		
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/)
/*<memory copy in 2D case>*/
{
	int i1=0,i2=0;
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[s2d+i2][s1d+i1]=src[s2s+i2][s1s+i1];
	}
}

void mcp_ad2d(float **dst /*destination*/, 
		float **src /*source (bigger)*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,		
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/)
/*<memory copy and addition in 2D case>*/
{
	int i1,i2;
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{		
		dst[s2d+i2][s1d+i1]=dst[s2d+i2][s1d+i1]+src[s2s+i2][s1s+i1];
	}
}


void mcp3d(float ***dst /*destination*/, 
		float ***src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis2*/)		
/*<memory copy in 3D case>*/
{
	int i1=0,i2=0,i3=0;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[s3d+i3][s2d+i2][s1d+i1]=src[s3s+i3][s2s+i2][s1s+i1];
	}
}


void mcp3d1d(float ***dst /*destination*/, 
		float *src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis3*/,
		int n1 /*source length in axis1*/,
		int n2 /*source length in axis2*/,
		int n3 /*source length in axis3*/)		
/*<memory copy in 3D case>*/
{
	int i1=0,i2=0,i3=0;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[s3d+i3][s2d+i2][s1d+i1]=src[(s3s+i3)*n1*n2+(s2s+i2)*n1+s1s+i1];
	}
}

void mcp3d1d1d(float *dst /*destination*/, 
		float *src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis3*/,
		int n1 /*source length in axis1*/,
		int n2 /*source length in axis2*/,
		int n3 /*source length in axis3*/)		
/*<memory copy in 3D case: 1D to 1D array>*/
{
	int i1=0,i2=0,i3=0;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[(s3d+i3)*l1*l2+(s2d+i2)*l1+s1d+i1]=src[(s3s+i3)*n1*n2+(s2s+i2)*n1+s1s+i1];
	}
}

void mcp3d1d1dint(int *dst /*destination*/, 
		int *src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis3*/,
		int n1 /*source length in axis1*/,
		int n2 /*source length in axis2*/,
		int n3 /*source length in axis3*/)		
/*<memory copy in 3D case: 1D to 1D integer array>*/
{
	int i1=0,i2=0,i3=0;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[(s3d+i3)*l1*l2+(s2d+i2)*l1+s1d+i1]=src[(s3s+i3)*n1*n2+(s2s+i2)*n1+s1s+i1];
	}
}

void mcp_ad1d3d(float *dst /*destination*/, 
		float ***src /*source (bigger)*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,
		int n1 /*source length in axis1*/,
		int n2 /*source length in axis2*/,
		int n3 /*source length in axis3*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis3*/)		
/*<memory copy and addition in 3D case>*/
{
	int i1,i2,i3;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{		
		dst[(s3d+i3)*n1*n2+(s2d+i2)*n1+s1d+i1]=dst[(s3d+i3)*n1*n2+(s2d+i2)*n1+s1d+i1]+src[s3s+i3][s2s+i2][s1s+i1];
	}
}

void mcp_ad3d(float ***dst /*destination*/, 
		float ***src /*source (bigger)*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis3*/)		
/*<memory copy and addition in 3D case>*/
{
	int i1,i2,i3;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{		
		dst[s3d+i3][s2d+i2][s1d+i1]=dst[s3d+i3][s2d+i2][s1d+i1]+src[s3s+i3][s2s+i2][s1s+i1];
	}
}

/*complex value functions*/
void cmcp(kiss_fft_cpx *dst /*destination*/, 
		kiss_fft_cpx *src /*source*/,
		int s1d /*starting index in dst*/,
		int s1s /*starting index in src*/,
		int l1  /*copy length in axis1*/)
/*<memory copy in 1D case>*/
{
	int i1=0;
	for(i1=0;i1<l1;i1++)
	{
		dst[s1d+i1]=src[s1s+i1];
	}
}

void cmcp_ad(kiss_fft_cpx *dst /*destination*/, 
		kiss_fft_cpx *src /*source*/,
		int s1d /*starting index in dst*/,
		int s1s /*starting index in src*/,
		int l1  /*copy length in axis1*/)
/*<memory copy and addition in 1D case>*/
{
	int i1=0;
	for(i1=0;i1<l1;i1++)
	{
		dst[s1d+i1]=sf_cadd(src[s1s+i1], dst[s1d+i1]) ;
	}
}

void cmcp2d(kiss_fft_cpx **dst /*destination*/, 
		kiss_fft_cpx **src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,		
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/)
/*<memory copy in 2D case>*/
{
	int i1=0,i2=0;
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[s2d+i2][s1d+i1]=src[s2s+i2][s1s+i1];
	}
}

void cmcp_ad2d(kiss_fft_cpx **dst /*destination*/, 
		kiss_fft_cpx **src /*source (bigger)*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,		
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/)
/*<memory copy and addition in 2D case>*/
{
	int i1,i2;
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{		
		dst[s2d+i2][s1d+i1]=sf_cadd(dst[s2d+i2][s1d+i1],src[s2s+i2][s1s+i1]);
	}
}

void cmcp3d(kiss_fft_cpx ***dst /*destination*/, 
		kiss_fft_cpx ***src /*source*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis2*/)		
/*<memory copy in 3D case>*/
{
	int i1=0,i2=0,i3=0;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		dst[s3d+i3][s2d+i2][s1d+i1]=src[s3s+i3][s2s+i2][s1s+i1];
	}
}

void cmcp_ad3d(kiss_fft_cpx ***dst /*destination*/, 
		kiss_fft_cpx ***src /*source (bigger)*/,
		int s1d /*starting index 1 in dst*/,
		int s2d /*starting index 2 in dst*/,	
		int s3d /*starting index 3 in dst*/,					
		int s1s /*starting index 1 in src*/,
		int s2s /*starting index 2 in src*/,
		int s3s /*starting index 3 in src*/,		
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*copy length in axis3*/)		
/*<memory copy and addition in 3D case>*/
{
	int i1,i2,i3;
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{		
		dst[s3d+i3][s2d+i2][s1d+i1]=sf_cadd(dst[s3d+i3][s2d+i2][s1d+i1],src[s3s+i3][s2s+i2][s1s+i1]);
	}
}

/*Memory initialization*/
void mi2d(float **din /*input data*/, 
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/)
/*<Memory initialization in 2D case>*/
{
	int i1=0,i2=0;
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		din[i2][i1]=0;
	}
}

void cmi2d(kiss_fft_cpx **din /*input data*/, 
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/)
/*<Memory initialization in 2D case>*/
{
	int i1=0,i2=0;
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		din[i2][i1].r=0;
		din[i2][i1].i=0;
	}
}

void mi3d(float ***din /*input data*/, 
		int l1 /*copy length in axis1*/,
		int l2 /*copy length in axis2*/,
		int l3 /*length in axis3*/)
/*<Memory initialization in 3D case>*/
{
	int i1=0,i2=0,i3=0;
	
	for(i3=0;i3<l3;i3++)
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		din[i3][i2][i1]=0;
	}
}

void cmi3d(kiss_fft_cpx ***din /*input data*/, 
		int l1 /*length in axis1*/,
		int l2 /*length in axis2*/,
		int l3 /*length in axis3*/)
/*<Memory initialization in 3D case>*/
{
	int i1=0,i2=0,i3=0;
	
	for(i3=0;i3<l3;i3++)	
	for(i2=0;i2<l2;i2++)
	for(i1=0;i1<l1;i1++)
	{
		din[i3][i2][i1].r=0;
		din[i3][i2][i1].i=0;
	}
}
