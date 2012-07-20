// Soft thresholding for 1-D data

#include<rsf.hh>

int main(int argc, char** argv)
{
  // Initialize RSF
  sf_init(argc,argv); 
    
  
  // Get input
  iRSF input;
  int n1;
  input.get("n1",n1);
  
  // Read data
  std::valarray<float> fdata(n1);
  input >> fdata;

  
  // Get parameter
  iRSF par(0);
  float mu;  
  par.get("mu",mu);
  // threshold value


  // Soft thresholding
  for (int i=0; i<n1; i++) {
    if (fdata[i]<=-mu)
      fdata[i]=fdata[i]+mu;
    else if (fdata[i]>=mu)
      fdata[i]=fdata[i]-mu;
    else
      fdata[i]=0;
  }
  

  // Set output
  oRSF output;

  // Write data
  output << fdata;

  exit(0);
}
