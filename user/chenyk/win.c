/* Local processing window utility function */
/*
  Copyright (C) 2016 Yangkang Chen
*/

#include<rsf.h>
#include"win.h"

void win_weight2d(float **din /*input data*/, 
		int iw1 /*starting window 1 in dst*/,
		int iw2 /*starting window 2 in dst*/,		
		int nw1 /*no of windows 1 in src*/,
		int nw2 /*no of windows 2 in src*/,
		int n1win /*window length 1 in src*/,
		int n2win /*window legnth 2 in src*/,		
		int ov1 /*copy length in axis1*/,
		int ov2 /*copy length in axis2*/)
/*<weights of local window in 2D (checked 100% correct>*/
{
	int i1,i2;

	if(iw2!=0)	
	{	
		for(i1=0;i1<n1win;i1++)
		for(i2=0;i2<ov2;i2++)
			din[i2][i1]=din[i2][i1]*(i2+1)/(ov2+1);
	}
	
	if(iw2!=nw2-1)
	{	
		for(i1=0;i1<n1win;i1++)
		for(i2=0;i2<ov2;i2++)
			din[n2win-ov2+i2][i1]=din[n2win-ov2+i2][i1]*(ov2-i2)/(ov2+1);
	}	
	

	if(iw1!=0)
	{
		/*Upper*/	
		for(i2=0;i2<n2win;i2++)
		for(i1=0;i1<ov1;i1++)
			din[i2][i1]=	din[i2][i1] * (i1+1)/(ov1+1);	
	}
	if(iw1!=nw1-1)
	{
		/*Down*/
		for(i2=0;i2<n2win;i2++)
		for(i1=0;i1<ov1;i1++)
			din[i2][n1win-ov1+i1]=	din[i2][n1win-ov1+i1] * (ov1-i1)/(ov1+1);	
	}

}

void win_weight3d(float ***din /*input data*/, 
		int iw1 /*starting window 1 in dst*/,
		int iw2 /*starting window 2 in dst*/,
		int iw3 /*starting window 3 in dst*/,						
		int nw1 /*no of windows 1 in src*/,
		int nw2 /*no of windows 2 in src*/,
		int nw3 /*no of windows 3 in src*/,		
		int n1win /*window length 1 in src*/,
		int n2win /*window legnth 2 in src*/,		
		int n3win /*window legnth 3 in src*/,				
		int ov1 /*copy length in axis1*/,
		int ov2 /*copy length in axis2*/,
		int ov3 /*copy length in axis3*/)		
/*<weights of local window in 3D >*/
{
	int i1,i2,i3;

	if(iw3!=0)	
	{	
		for(i1=0;i1<n1win;i1++)
		for(i2=0;i2<n2win;i2++)
		for(i3=0;i3<ov3;i3++)
			din[i3][i2][i1]=din[i3][i2][i1]*(i3+1)/(ov3+1);
	}
	
	if(iw3!=nw3-1)
	{	
		for(i1=0;i1<n1win;i1++)	
		for(i2=0;i2<n2win;i2++)
		for(i3=0;i3<ov3;i3++)
			din[n3win-ov3+i3][i2][i1]=din[n3win-ov3+i3][i2][i1]*(ov3-i3)/(ov3+1);
	}
	
	if(iw2!=0)	
	{	
		for(i3=0;i3<n3win;i3++)
		for(i1=0;i1<n1win;i1++)
		for(i2=0;i2<ov2;i2++)
			din[i3][i2][i1]=din[i3][i2][i1]*(i2+1)/(ov2+1);
	}
	
	if(iw2!=nw2-1)
	{	
		for(i3=0;i3<n3win;i3++)	
		for(i1=0;i1<n1win;i1++)
		for(i2=0;i2<ov2;i2++)
			din[i3][n2win-ov2+i2][i1]=din[i3][n2win-ov2+i2][i1]*(ov2-i2)/(ov2+1);
	}	
	

	if(iw1!=0)
	{
		/*Upper*/
		for(i3=0;i3<n3win;i3++)				
		for(i2=0;i2<n2win;i2++)
		for(i1=0;i1<ov1;i1++)
			din[i3][i2][i1]=	din[i3][i2][i1] * (i1+1)/(ov1+1);	
	}
	if(iw1!=nw1-1)
	{
		/*Down*/
		for(i3=0;i3<n3win;i3++)		
		for(i2=0;i2<n2win;i2++)
		for(i1=0;i1<ov1;i1++)
			din[i3][i2][n1win-ov1+i1]=	din[i3][i2][n1win-ov1+i1] * (ov1-i1)/(ov1+1);	
	}

}

