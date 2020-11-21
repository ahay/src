/* Generate the header file for all shot and receiver.
 * input: velocity:nz,nx
 * Output: ntrace*15 matrix. The first dimension is 15.data[0]=ishot,data[1]=s,data[2]=r,data[3]=s_z,data[4]=r_z;
 * data[5]=(s+r)/2,data[6]=(r-s)/2;data[14]=tracenow. The other data is used for your personal data.
 * source coordinate s: from min2 to max2 of your velocity axis. The step is d2.
 * The receiver coordinate r: |r-s|\leq maxoffset,r>=min2,r<max2. the step is 2*d2;
 */
/*uploaded by zedong wu*/

#include <rsf.h>
#define LONG_INT long long 
#define LEN_FILE 3000
#define _FILE_OFFSET_BITS 64
int main(int argc, char* argv[])
{
    sf_init(argc,argv);
    sf_file Fin,Fout;
    long long int n1,n2;
    float d1,o1,d2,o2;
    Fin=sf_input ("in" );
    sf_axis a1,a2;
    a1 = sf_iaxa(Fin,1); n1 = sf_n(a1); d1=sf_d(a1);o1=sf_o(a1); 
    a2=sf_iaxa(Fin,2);n2=sf_n(a2);d2=sf_d(a2);o2=sf_o(a2);
    fprintf(stderr,"%lld %lld %f %f %f %f\n",n1,n2,d1,d2,o1,o2);
    int maxoffset;
    if (!sf_getint("maxoffset",&maxoffset)) maxoffset=100000000;
        
    Fout=sf_output("out");
    sf_putint(Fout,"n2",15);
    float temp[15];
    int i,j;
    int ntrace=0;
    for(i=0;i<n2;i++)
    for(j=0;j<n2;j++)
    {
        if((abs(j-i)%2==0)&&(abs(j-i)*d2)<=maxoffset)
        {
            ntrace++;
        }
    }
    sf_putint(Fout,"n1",ntrace);
    ntrace=0;
    for(i=0;i<n2;i++)
    for(j=0;j<n2;j++)
    {
        if((abs(j-i)%2==0)&&(abs(j-i)*d2)<=maxoffset)
        {
            float s=o2+i*d2;
            float r=o2+j*d2;
            float h=(r-s)/2.0;
            float x=(s+r)/2;
            temp[0]=i;
            temp[1]=s;
            temp[2]=r;
            temp[3]=0;
            temp[4]=0;
            temp[5]=x;
            temp[6]=h;
	        temp[14]=ntrace;
	        ntrace++;
            fprintf(stderr,"temp[0]=%f %f %f %f %f ntrace=%f\n",temp[0],s,r,h,x,temp[14]);
            sf_floatwrite(temp,15,Fout);
        }   
    }
    sf_close();
    return 0;
}
