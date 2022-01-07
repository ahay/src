//Diffraction extraction by Fresnel Zone Elimination
//Input DAG file (dip,x,z)
//Input Velocity file (x,z)
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif


int main (int argc, char* argv[]){

  //dimensions relating to dip angle gathers
  int nz, nx, ndip;
  float dz, dx, ddip;
  float z0, x0, dip0;
  //input dip angle gather file
  sf_file gathers_in, dips_in=NULL, fresnel_out;
  
  //parameters relating to the calculated Fresnel Zone
  float width_mask, width_gray;

  //parameters relating to input and output file size
  int gather_size, xz_size;

  //get MADAGASCAR flowing
  sf_init(argc,argv);
  //Declare Inputs
  gathers_in = sf_input("in") ;
  //Dip Angle Gathers are Primary Input, need to be float
    if(SF_FLOAT != sf_gettype (gathers_in)) sf_error("Dip Angle Gathers must be FLOAT");
  /*dip angle gathers*/
  //Read secondary input -- image dips in degrees
  if (NULL != sf_getstring("dip")){/*Image Dips (deg)*/
    dips_in = sf_input("dip");
    if(SF_FLOAT != sf_gettype (dips_in)) sf_error("Dip File must be FLOAT");}
  //Output
  fresnel_out = sf_output("out");

  //Read Input Files to get key Dimensions
  //vertical (3rd axis of dip angle gather)
  if (!sf_histint (gathers_in, "n3", &nz)) sf_error("Need n3=nz in DAG file");
  if (!sf_histfloat (gathers_in, "d3", &dz)) sf_error("Need d3=dz in DAG file");
  if (!sf_histfloat (gathers_in, "o3", &z0)) sf_error("Need o3=oz in DAG file");
  //horizontal (2rd axis of dip angle gather)
  if (!sf_histint (gathers_in, "n2", &nx)) sf_error("Need n2=nx in DAG file");
  if (!sf_histfloat (gathers_in, "d2", &dx)) sf_error("Need d2=dx in DAG file");
  if (!sf_histfloat (gathers_in, "n2", &x0)) sf_error("Need o2=ox in DAG file");
  //dip angle (1st axis of dip angle gather)
  if (!sf_histint (gathers_in, "n1", &ndip)) sf_error("Need n1=ndip in DAG file");
  if (!sf_histfloat (gathers_in, "d1", &ddip)) sf_error("Need d1=ddip in DAG file");
  if (!sf_histfloat (gathers_in, "o1", &dip0)) sf_error("Need 01=odip in DAG file");
  if(!sf_getfloat("mask",&width_mask)) sf_error("Declare mask=");
  if(!sf_getfloat("gray",&width_gray)) sf_error("Declare gray=");

  //read input dips
  xz_size = nx*nz;
  float* dips = sf_floatalloc(xz_size);
  sf_seek(dips_in,0,SEEK_SET);
  sf_floatread(dips,xz_size,dips_in);
  //declare output size
  gather_size = ndip;
  float* fresnel = sf_floatalloc(gather_size);
  
  float local_dip = 0.f;
  float left_dip = 0.f;
  float left_gray = 0.f;
  float right_dip = 0.f;
  float right_gray = 0.f;
  float current_dip = 0.f;
  float yarg = 1 /( width_gray );
  //int omp_thread = 0;
    //looping integers and indices
    int ix,iz,ind,idip;
  ind = 0;
    //loop through vertical axis
  for(iz = 0; iz < nz; iz++){

    //loop through the horizontal
    for (ix = 0; ix < nx; ix++){
      //data is in helix.  which index to call?
      ind = iz * nx + ix;
      //read the local dip data
      local_dip = dips[ind];
       //what is the left extent of the the hard mask
      left_dip = local_dip - width_mask ;
      //right extent?
      right_dip = local_dip + width_mask ;
      //left gray?
      left_gray = left_dip + width_gray;
      //right gray?
      right_gray = right_dip - width_gray;
      //      sf_warning("%f",right_dip);
      //allocate memory for collection of fresnel zones from this horizontal position
      memset (fresnel,0,gather_size * sizeof(float));
      //loop through dip angle for each point's fresnel zone
#ifdef _OPENMP
#pragma omp parallel for schedule(static)	    \
    private(idip,current_dip)		    \
  shared(ndip,ddip,dip0,left_dip,right_dip,left_gray,right_gray,fresnel,yarg)
#endif


      for (idip = 0 ; idip < ndip ; idip++){
/*
#ifdef _OPENMP
	    omp_thread=omp_get_thread_num();
#pragma omp critical
#endif
*/
	current_dip = idip*ddip+dip0;
	if (current_dip > left_dip){//inside of left side
	  if ( current_dip < right_dip ){//inside of right side
	    //in the gray area?
      	    if( current_dip > left_gray){
	      if( current_dip < right_gray){//in the zone
	      fresnel[idip] = 0.0;
	      }else{//right gray area
		fresnel[idip] = 1 - ( right_dip - current_dip ) * yarg;
	      }
	      }else{//left gray area
	      fresnel[idip] = 1 - ( current_dip - left_dip ) * yarg;
	    }
	  }else{
	    fresnel[idip] = 1.0;//wide right
	  }
	}else{//wide left
	  fresnel[idip] = 1.0;
	}
      }//end idip loop
      //write out fresnel zone
      sf_floatwrite(fresnel,gather_size,fresnel_out);
    }//end horizontal axis loop

  }//end vertical axis loop
  //close files
  sf_fileclose ( dips_in );
  sf_fileclose ( gathers_in );
  sf_fileclose ( fresnel_out );
  //free arrays
  free ( fresnel );
  free ( dips ) ;


  return 0;

}//end main
