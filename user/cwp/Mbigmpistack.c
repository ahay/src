/* remap and stacks rsf files using mpi 

Assumes that files are commonly named sequentially, e.g.:

File001.rsf
File002.rsf
File003.rsf ...
FileN.rsf

Such that all files can be represented as a prefix, which 
is a printf like statement that will be evaluated for all
files to be included in a range.

For the above example the prefix would be:

prefix="File%03d.rsf" 

The nf, jf, and of parameters specify a range of numbers to evaluate the 
prefix for, giving the program filenames to be used for summing
together.  For example:  

nf=10,of=0,jf=1 --> (0,1,2,3,4,5,6,7,8,9,10)
nf=10,of=5,jf=2 --> (5,7,9,11,13,15,17,19,21,23)

If there are more files than processes, then this program will subdivide
the files onto various processes, and run multiple rounds until
everything is done.

These must be 3D arrays (or 2D ,but with three dimensions), arrays must be
X-Y-Z
a1-a2-a3
*/

#include <rsf.h>
#include <mpi.h>
#include <stdio.h>

void zero_array1(float *array, int n1, int n2, int n3);
void zero_array(float ***array, int n1, int n2, int n3);
int roundupdown(float val);



int roundupdown(float val){
    int lower = floor(val);
    int higher = ceil(val);

    float dhigh = higher - val;

    float dlow  = val - lower;

    if (dlow < dhigh) return lower;
    else return higher;
}
void zero_array1(float *array, int n1, int n2, int n3)
{
    int index;
    int i3, i2, i1;

    for(i3=0; i3 < n3; ++i3){
        for(i2=0; i2 < n2; ++i2){
            for(i1=0; i1 < n1; ++i1){
                index = i1 + i2*n1 + i3*n2*n1;
                array[index] = 0.0f;
            }
        }
    }
}
void zero_array(float ***array, int n1, int n2, int n3)
{
    int i3, i2, i1;

    for(i3 = 0; i3 < n3; ++i3){
        for(i2 = 0; i2 < n2; ++i2){
            for(i1 = 0; i1 < n1; ++i1){
                array[i3][i2][i1] = 0.0;
            } // x
        } // y
    } // w 
}


int main(int argc, char **argv){

    MPI_Init(&argc,&argv);
    
    char *prefix = sf_charalloc(1024);
    char *outName = sf_charalloc(1024);
    char *shotName = sf_charalloc(1024);
    
    int nf,of,jf,nx,ny,nz;
    // index issues
    int inx,iny,inz;
    int idx,idy,idz;
    int iox,ioy,ioz;

    float fdx,fox,fdy,foy,fdz,foz;
   
    float dx,ox,dy,oy,dz,oz;
    int RANK;
    int PROCS;
    
    bool debug,verb,useShots =0;
    int ir, ip, iz, iy, ix;
    
    MPI_Comm_size(MPI_COMM_WORLD,&PROCS);
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
        
    sf_init(argc,argv);
    
    if (! sf_getint("nx",&nx)) sf_error("Must specify nx"); /* origin of files*/
    if (! sf_getint("ny",&ny)) sf_error("Must specify ny"); /* origin of files*/
    if (! sf_getint("nz",&nz)) sf_error("Must specify nz"); /* origin of files*/
    
    if (! sf_getbool("debug",&debug)) debug=false;
    if (! sf_getbool("verb",&verb)) verb=false;
    if (! sf_getfloat("dx",&dx)) sf_error("Must specify dx");
    if (! sf_getfloat("dy",&dy)) sf_error("Must specify dy");
    if (! sf_getfloat("dz",&dz)) sf_error("Must specify dz");
    if (! sf_getfloat("oz",&oz)) sf_error("Must specify oz");
    if (! sf_getfloat("oy",&oy)) sf_error("Must specify oy");
    if (! sf_getfloat("ox",&ox)) sf_error("Must specify ox");
    
    prefix = (char *)sf_getstring("prefix"); /* printf like prefix */
    if (prefix == NULL) sf_error("Must specify prefix");
    
    shotName = (char *)sf_getstring("shots"); /* name of shot file */
    if (shotName != NULL) {
        sf_warning("Found shots file: %s", shotName);
        useShots = 1;
    } else {
        if (! sf_getint("nf",&nf)) sf_error("Must specify how many files to stack"); /* number of files to stack */
        if (! sf_getint("jf",&jf)) sf_error("Must specify jf"); /* delta between files */
        if (! sf_getint("of",&of)) sf_error("Must specify of"); /* origin of files*/
    }
    
    outName    = (char *)sf_getstring("oname"); /* name of output file */
    if (outName== NULL) sf_error("Must specify output name");


    
   
    int **MPI_MAP = NULL;
    int nrounds;
    if (useShots){
        sf_file Fshots = sf_input(shotName);
        sf_axis as = sf_iaxa(Fshots,2);
        int ns = sf_n(as);
        if (verb && RANK == 0) sf_warning("found %d shots",ns); 
        int *shots = sf_intalloc(ns);
        sf_intread(shots,ns,Fshots);
        nrounds = ns/PROCS;
        if (ns % PROCS != 0) nrounds++;
        MPI_MAP = sf_intalloc2(PROCS,nrounds);

        int is = 0;
        for(ir = 0; ir < nrounds; ++ir){
            if (verb&& RANK==0) fprintf(stderr,"ROUND %d .... ", ir);
            for(ip=0; ip < PROCS; ++ip){
                if(is < ns) {
                    MPI_MAP[ir][ip] = shots[is];
                    is += 1;
                } else {
                   MPI_MAP[ir][ip] = -1;
                }
                if (verb&& RANK==0) fprintf(stderr,"%d ", MPI_MAP[ir][ip]);
            }
            if (verb&& RANK==0) fprintf(stderr,"\n");
        }
        if(verb&& RANK == 0) fprintf(stderr, "\n");
        sf_fileclose(Fshots);
        free(shots);
    } else {
        nrounds = nf/PROCS;
        if (nf % PROCS != 0) nrounds++;
        MPI_MAP = sf_intalloc2(PROCS,nrounds);
        int file = of;
        for(ir = 0; ir < nrounds; ++ir){
            if (verb&& RANK==0) fprintf(stderr,"ROUND %d .... ", ir);
           
            for(ip=0; ip < PROCS; ++ip){
                if(file < of+jf*nf) {
                    MPI_MAP[ir][ip] = file;
                    file += jf;
                } else {
                   MPI_MAP[ir][ip] = -1;
                }
                if (verb&& RANK==0) fprintf(stderr,"%d ", MPI_MAP[ir][ip]);
            }
            if (verb&& RANK==0) fprintf(stderr,"\n");
        }
        if(verb&& RANK == 0) fprintf(stderr, "\n");
    }
    
    
    char *filename = sf_charalloc(1024);

    float ***array = sf_floatalloc3(nx,ny,nz);
   // float ***output_array = sf_floatalloc3(nx,ny,nz);

    zero_array(array,nx,ny,nz);
    sf_axis fileX,fileY,fileZ;

    for(ir = 0; ir < nrounds; ++ir){
        
        int fnumber = MPI_MAP[ir][RANK];
        
        if (fnumber < 0) { //I am not reading any files...
            if (verb) sf_warning("%d skipping round %d", RANK,filename,ir);
            continue;
        } else {
            sprintf(filename,prefix,fnumber);
            if (verb) sf_warning("%d reading... %s",RANK,filename);
            sf_file file = sf_input(filename);
            fileX = sf_iaxa(file,1);
            fileY = sf_iaxa(file,2);
            fileZ = sf_iaxa(file,3);
            fox = sf_o(fileX);
            foy = sf_o(fileY);
            foz = sf_o(fileZ);

            fdx = sf_d(fileX);
            fdy = sf_d(fileY);
            fdz = sf_d(fileZ);

            inx = sf_n(fileX);
            iny = sf_n(fileY);
            inz = sf_n(fileZ);

            iox = (int)roundupdown((fox - ox) / dx);
            ioy = (int)roundupdown((foy - oy) / dy);
            ioz = (int)roundupdown((foz - oz) / dz);
            if (iox < 0 || ioy < 0 || ioz < 0) {
                    sf_error("%d got starting indexes less than zero: iox %d ioy %d ioz %d",RANK,iox,ioy,ioz);
                    MPI_Abort(MPI_COMM_WORLD, 15);
            }
            idx = (int)roundupdown(fdx / dx);
            idy = (int)roundupdown(fdy / dy);
            idz = (int)roundupdown(fdz / dz);
            if (idx <= 0 || idy <= 0 || idz <= 0) {
                    sf_warning("got delta indices less than zero: idx %d idy %d",idx,idy);
                    MPI_Abort(MPI_COMM_WORLD,16);
            }
            
            if (debug) sf_warning("%d got indices %d %d %d %d %d %d %d %d %d",RANK,inx,iny,inz,iox,ioy,ioz,idx,idy,idz);
            float ***tarray = sf_floatalloc3(inx,iny,inz);
            zero_array(tarray,inx,iny,inz);

            sf_floatread(tarray[0][0],inx*iny*inz,file);
            if(debug) sf_warning("%d starting updating", RANK);

            int iiz,iiy,iix;
            for(iz = 0; iz < inz; ++iz){
                for(iy = 0; iy < iny; ++iy){
                    for(ix = 0; ix < inx; ++ix){
                        iiz = ioz+iz*idz;
                        iiy = ioy+iy*idy;
                        iix = iox+ix*idx;

                        if ((iiz < 0 || iiz > nz-1) || (iiy < 0 || iiy > ny-1) || (iix < 0 || iix > nx-1)) {
                            fprintf(stderr,"%d %d %d %d %d %d\n",ix,iy,iz,iox+ix*idx,ioy+iy*idy,ioz+idz*iz);
                            fprintf(stderr,"%d %d %d\n",iox,ioy,ioz);
                            fprintf(stderr,"%d %d %d\n",idx,idy,idz);
                            sf_warning("FATAL ERROR OUT OF BOUNDS");
                            MPI_Abort(MPI_COMM_WORLD,233);
                        }
                        array[ioz+iz*idz][ioy+iy*idy][iox+ix*idx] += tarray[iz][iy][ix];
                    }
                }
            }

            if(debug) sf_warning("%d finished updating", RANK);
            
            sf_fileclose(file);
            if(debug) sf_warning("%d close file", RANK);
            free(**tarray); free(*tarray); free(tarray);
            if(debug) sf_warning("%d free junk", RANK);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  
    if(debug) sf_warning("ENTERING REDUCE");
    if(verb && RANK == 0) sf_warning("ENTERING REDUCE");
    float ***hold = sf_floatalloc3(nx,ny,nz);
    zero_array(hold,nx,ny,nz);

    MPI_Reduce(array[0][0],hold[0][0],nx*ny*nz,
               MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);

    if(verb && RANK == 0) sf_warning("Finished reduce");
    if (RANK == 0){
        if (debug) sf_warning("ROOT: Writing");
        sf_axis ax,ay,az;
        ax = sf_maxa(nx,ox,dx);
        ay = sf_maxa(ny,oy,dy);
        az = sf_maxa(nz,oz,dz);

        sf_file Foutput = sf_output(outName);
        sf_oaxa(Foutput,ax,1);
        sf_oaxa(Foutput,ay,2);
        sf_oaxa(Foutput,az,3);
        sf_floatwrite(hold[0][0],nx*ny*nz,Foutput);
        sf_fileclose(Foutput);
        if(debug) sf_warning("ROOT: Finished");
    }
    free(**hold); free(*hold); free(hold);
    free(**array); free(*array); free(array);

    MPI_Finalize();
}
