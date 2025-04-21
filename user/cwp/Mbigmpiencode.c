/* shot encoding with arbitrary phase and amplitude weights using MPI on a distributed cluster 


Use mpiencode if your shots are on the same grid prior to encoding.

Use bigmpiencode if your shots are not on a single grid prior to encoding.
YOUR SHOTS MUST ALL FALL ONTO THE SAME REGULAR GRID.  BigMPIENCODE does not do any
shot interpolation.

Data axes - X, Y, W

*/
#include <rsf.h>
#include <mpi.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159

void zero_array(sf_complex ***array, int n1, int n2, int n3);
void zero_array1(sf_complex *array, int n1, int n2, int n3);
void update_receiver_encoding(sf_complex ***encoding, sf_complex *data, float ow, float dw, float time, float weight, float phase, int nx, int ny, int nw, int ox, int jx, int oy, int jy );
void save_encoding(sf_complex ***encoding, char *prefix, int iencoding, sf_axis ax, sf_axis ay, sf_axis aw);
int roundupdown(float val);



int roundupdown(float val){
    int lower = floor(val);
    int higher = ceil(val);

    float dhigh = higher - val;

    float dlow  = val - lower;

    if (dlow < dhigh) return lower;
    else return higher;
}

void zero_array1(sf_complex *array, int n1, int n2, int n3)
{
    int index;
    int i3, i2, i1;

    for(i3=0; i3 < n3; ++i3){
        for(i2=0; i2 < n2; ++i2){
            for(i1=0; i1 < n1; ++i1){
                index = i1 + i2*n1 + i3*n2*n1;
                array[index] = sf_cmplx(0.0f,0.0f);
            }
        }
    }
}
void zero_array(sf_complex ***array, int n1, int n2, int n3)
{
    int i3, i2, i1;

    for(i3 = 0; i3 < n3; ++i3){
        for(i2 = 0; i2 < n2; ++i2){
            for(i1 = 0; i1 < n1; ++i1){
                array[i3][i2][i1] = sf_cmplx(0.0f,0.0f);
            } // x
        } // y
    } // w 
}

void update_receiver_encoding(sf_complex ***encoding, sf_complex *data, float ow, float dw, float time, float weight, float phase, int nx, int ny, int nw, int ox, int jx, int oy, int jy )
{
    int iw, iy, ix, index;
    float shift;
    sf_complex cs;
    float pshift = -2.0*PI*phase;
    sf_complex cw = sf_cmplx(weight, 0.0f);
    sf_complex cp = sf_cmplx(cosf(pshift),sinf(pshift));


    for(iw = 0; iw < nw; ++iw){
        shift = -2.0*PI*(ow+dw*iw)*time;
        cs = sf_cmplx(cosf(shift),sinf(shift));
        for(iy = 0; iy < ny; ++iy){
            for(ix = 0; ix < nx; ++ix){
                index = ix + iy*nx + iw*nx*ny; // linear index
                encoding[iw][iy*jy+oy][ix*jx+ox] += cw*cs*cp*data[index];
             } // x
        } // y
   } // w 
}

void save_encoding(sf_complex ***encoding, char *prefix, int iencoding, sf_axis ax, sf_axis ay, sf_axis aw)
{

    char encodingName[256];
    sprintf(encodingName,prefix,iencoding);
    sf_file Foutput = sf_output(encodingName);
    sf_settype(Foutput,SF_COMPLEX);
    sf_oaxa(Foutput,ax,1);
    sf_oaxa(Foutput,ay,2);
    sf_oaxa(Foutput,aw,3);
    sf_complexwrite(encoding[0][0],sf_n(ax)*sf_n(ay)*sf_n(aw),Foutput);
    sf_fileclose(Foutput);
}

int main(int argc, char **argv){

    bool verb;
    
    //sf_complex *data = NULL;    //1D array for reading and sending shot records
    
    sf_complex ***receiver_encoding; //output receiver encoding
    
    int   **emap;            // map of encodings to MPI processes
    float **phase;          // array of phase shifts
    float **delays;         // array of time shifts
    float **ampls;          // array of amplitude weights
    //variables that are *xo or *yo are output cube parameters 
    int ns,ne,os,ds,nx,ny,nw,nxi,nyi,iox,ioy,idx,idy;
    float ow,dw,dyi,dxi,oyi,oxi,dy,dx,ox,oy;
    sf_axis ax,ay,aw,as,ae,axi,ayi;
    
    sf_file Fencode, Fshotrecord;

    bool mapShots = false;
 
    char *data_prefix = sf_charalloc(256);     //printf like string
    char *encoding_prefix = sf_charalloc(256); //printf like string for output
    
    int RANK;  // this processes rank
    int PROCS; // total number of processes initialized
    int is,ir,ip;
    
    MPI_Init(&argc,&argv);
    
    MPI_Comm_size(MPI_COMM_WORLD,&PROCS);
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    
    sf_init(argc,argv); //Init RSF

    Fencode = sf_input("encode"); /* encoding file from sfencodemaker */
    
    data_prefix     = (char *)sf_getstring("dprefix"); /* printf like statement that can be evaluated to find the data files corresponding to shot records */
    encoding_prefix = (char *)sf_getstring("eprefix"); /* printf like statement that can be evaluated for the output encodings */
    if (data_prefix == NULL) sf_error("Must enter dprefix!");
    if (encoding_prefix == NULL) sf_error("Must enter eprefix!");

    if (! sf_getint("nx",&nx)) sf_error("Must enter nx"); /* # of output grid x points */
    if (! sf_getint("ny",&ny) )sf_error("Must enter ny"); /* # of output grid y points */
    if (! sf_getfloat("dy",&dy)) sf_error("Must enter dy"); /* dy of output grid points */
    if (! sf_getfloat("dx",&dx) )sf_error("Must enter dx"); /* dx of output grid points */
    if (! sf_getfloat("ox",&ox)) sf_error("Must enter ox"); /* ox of output grid points */
    if (! sf_getfloat("oy",&oy)) sf_error("Must enter oy"); /* ox of output grid points */

   
    char *shotfile = (char *) sf_getstring("shots"); /* shot-file name, dimensions are 1xNS */

    if (shotfile != NULL){
        mapShots = true;
        sf_warning("Will map to shots in shot file %s",shotfile);
    }

    ax = sf_maxa(nx,ox,dx);
    ay = sf_maxa(ny,oy,dy);
    
    if(! sf_getbool("verb",&verb)) verb = false;

    as = sf_iaxa(Fencode,1); 
    ns = sf_n(as); ds = (int)sf_d(as); os = (int)sf_o(as);
    ae = sf_iaxa(Fencode,2);
    ne = sf_n(ae);
    
    sf_file Faxes;
    char tname[256];
    sprintf(tname,data_prefix,os);
    if (verb && RANK==0) sf_warning("will try to read: %s", tname);
    Faxes = sf_input(tname);
    
    aw = sf_iaxa(Faxes,3); nw = sf_n(aw); dw=sf_d(aw); ow=sf_o(aw);
    
    sf_fileclose(Faxes); /* Discard file */
    
    /* Output axes */
    if (verb && RANK==0){
        sf_raxa(as);
        sf_raxa(ae);
        sf_raxa(ax);
        sf_raxa(ay);
        sf_raxa(aw);
    }
    
    int nrounds = ne/PROCS; // How many separate rounds do we need?
    if ( ne % PROCS !=  0) nrounds++; // Catch overflow encodings
    if (verb && RANK == 0) sf_warning("MPI map info: %d encodings; %d procs; %d rounds",ne,PROCS,nrounds);

    int *shots = NULL;

    if (mapShots){

        sf_file Fshots = sf_input(shotfile);
        sf_axis Ashots = sf_iaxa(Fshots,2);
        int ns = sf_n(Ashots);

        shots = sf_intalloc(ns);

        sf_intread(shots,ns,Fshots);

    }

    
    /* Map encodings to processes */
    emap = sf_intalloc2(PROCS,nrounds);
    int te = 0;
    for(ir = 0; ir < nrounds; ++ir){
        if (verb && RANK==0) fprintf(stderr,"ROUND %d .... ", ir);
        for(ip=0; ip < PROCS; ++ip){
            if(te < ne) {
                if (mapShots) {
                    emap[ir][ip] = shots[te];
                    }
                else {
                    emap[ir][ip] = te;
                    }
                ++te;
            } else {
               emap[ir][ip] = -1;
            }

            if (verb && RANK==0) fprintf(stderr,"%d ", emap[ir][ip]);
        }
        if (verb && RANK==0) fprintf(stderr,"\n");
    }
    if(verb && RANK == 0) fprintf(stderr, "\n");
    
    /* Read amplitude and phase encodings */
    phase = sf_floatalloc2(ns,ne);
    ampls = sf_floatalloc2(ns,ne);
    delays = sf_floatalloc2(ns,ne);
    sf_floatread(ampls[0],ns*ne,Fencode);
    sf_floatread(delays[0],ns*ne,Fencode);
    sf_floatread(phase[0],ns*ne,Fencode);

    if (verb) fprintf(stderr,"nx %d ny %d nw %d \n",nx,ny,nw);

    receiver_encoding = sf_complexalloc3(nx,ny,nw);

    int *shot_pars = sf_intalloc(6);
    sf_complex *data = NULL;
    
    if (RANK == 0) {  
        for (ir = 0; ir < nrounds; ++ir){
            
            int encoding = emap[ir][RANK];
            int shot = os;

            zero_array(receiver_encoding,nx,ny,nw);

            for (is = 0; is < ns; ++is){
                if(verb) {
                    fprintf(stderr,"Round %d ---- %d / %d ---- \n", ir, is+1, ns);
                } else {
                fprintf(stderr,"Round %d ---- %d / %d ---- \r", ir, is+1, ns);
                }
                    
                float timeDelay = delays[encoding][is];
                float ampWeight = ampls[encoding][is];
                float phaseShift = phase[encoding][is];

                // Which shot record am I reading?
                char recordName[256];
                sprintf(recordName,data_prefix,shot);
                Fshotrecord = sf_input(recordName); 

                axi = sf_iaxa(Fshotrecord,1); nxi = sf_n(axi); dxi=sf_d(axi); oxi=sf_o(axi);
                ayi = sf_iaxa(Fshotrecord,2); nyi = sf_n(ayi); dyi=sf_d(ayi); oyi=sf_o(ayi);
                
                data = sf_complexalloc(nxi*nyi*nw);
                zero_array1(data,nxi,nyi,nw);

                iox = (int)roundupdown((oxi - ox) / dx);
                ioy = (int)roundupdown((oyi - oy) / dy);
                if (iox < 0 || ioy < 0) {
                    sf_warning("got starting indexes less than zero: iox %d ioy %d",iox,ioy);
                    MPI_Abort(MPI_COMM_WORLD, 15);
                }
                idx = (int)roundupdown(dxi / dx);
                idy = (int)roundupdown(dyi / dy);
                if (idx <= 0 || idy <= 0) {
                    sf_warning("got delta indices less than zero: idx %d idy %d",idx,idy);
                    MPI_Abort(MPI_COMM_WORLD,16);
                }

                shot_pars[0] = nxi; shot_pars[1] = nyi;
                shot_pars[2] = iox; shot_pars[3] = ioy;
                shot_pars[4] = idx; shot_pars[5] = idy;

                if (verb) fprintf(stderr,"nxi %d nyi %d iox %d ioy %d idx %d idy %d\n",nxi,nyi,iox,ioy,idx,idy);

                //Read this shot record
                sf_complexread((sf_complex *)data,nxi*nyi*nw,Fshotrecord); 

                if (verb) fprintf(stderr,"Read data\n");
                
                //Send information about this shot record
                MPI_Bcast((int *)shot_pars,6,MPI_INT,0,MPI_COMM_WORLD);

                if (verb) fprintf(stderr,"Broadcast shot pars\n");

                // Send shot record
                MPI_Bcast((float *)data,nxi*nyi*nw*2,
                    MPI_FLOAT,0,MPI_COMM_WORLD); 

                if (verb) fprintf(stderr,"Broadcast shot record\n");

                // Update my encoding
                update_receiver_encoding(receiver_encoding,
                                    (sf_complex *)data,
                                    ow,dw,timeDelay,ampWeight,phaseShift,
                                    nxi,nyi,nw,
                                    iox,idx,
                                    ioy,idy);

                if (verb) fprintf(stderr,"Added to encoding\n");
                
                sf_fileclose(Fshotrecord);
                if (verb) fprintf(stderr,"Closed record\n");
                shot = shot + ds;
                free(data);
                if (verb) fprintf(stderr,"Freed data\n");
                MPI_Barrier(MPI_COMM_WORLD);
               
            } // shots

            /* Save the encoding */
            save_encoding(receiver_encoding,
                          encoding_prefix,
                          encoding,
                          ax, ay, aw);
            if (verb) fprintf(stderr,"Saved encoding\n");
            
            MPI_Barrier(MPI_COMM_WORLD);
        } // rounds
        if(verb) fprintf(stderr,"\nFinished encoding\n");
    } else {
        for (ir = 0; ir < nrounds; ++ir){
           
            int encoding = emap[ir][RANK]; // which encoding should I write?
            int shot = os;

            if (encoding < 0) { /* We aren't writing to an encoding */
                   for(is = 0; is < ns; ++is){
                        MPI_Bcast((int *)shot_pars,6,MPI_INT,0,MPI_COMM_WORLD);
                        nxi = shot_pars[0]; nyi = shot_pars[1];
                        iox = shot_pars[2]; ioy = shot_pars[3];
                        idx = shot_pars[4]; idy = shot_pars[5];
                        
                        data = sf_complexalloc(nxi*nyi*nw);
                        zero_array1(data,nxi,nyi,nw);

                        MPI_Bcast((float *)data,nxi*nyi*nw*2, 
                            MPI_FLOAT,0,MPI_COMM_WORLD);
                        // Wait
                        free(data);
                        MPI_Barrier(MPI_COMM_WORLD);
                   }
                   MPI_Barrier(MPI_COMM_WORLD);
            }
            else {
                    zero_array(receiver_encoding,nx,ny,nw);

                    for(is = 0; is < ns; ++is){
                        
                        float timeDelay = delays[encoding][is];
                        float ampWeight = ampls[encoding][is];
                        float phaseShift = phase[encoding][is];

                        //Send information about this shot record
                        MPI_Bcast((int *)shot_pars,6,MPI_INT,0,MPI_COMM_WORLD);
                        nxi = shot_pars[0]; nyi = shot_pars[1];
                        iox = shot_pars[2]; ioy = shot_pars[3];
                        idx = shot_pars[4]; idy = shot_pars[5];
                        
                        if(verb) fprintf(stderr,"%d %d %d %d %d %d %d \n",
                                    RANK,nxi,nyi,iox,ioy,idx,idy);

                        data = sf_complexalloc(nxi*nyi*nw);
                        
                        if(verb) fprintf(stderr,"%d allocated\n",RANK);
                        zero_array1(data,nxi,nyi,nw);
                        if(verb) fprintf(stderr,"%d zeroed data\n",RANK);
                        // Get shot record
                        MPI_Bcast((float *)data,nxi*nyi*nw*2,
                                MPI_FLOAT,0,MPI_COMM_WORLD);
                        if(verb) fprintf(stderr,"%d GOT DATA\n",RANK);       
                        // Update my encoding
                        update_receiver_encoding(
                                            receiver_encoding,
                                            (sf_complex *)data,
                                            ow,dw,
                                            timeDelay,ampWeight,phaseShift,
                                            nxi,nyi,nw,
                                            iox,idx,
                                            ioy,idy);
                        if(verb) fprintf(stderr,"%d UPDATED ENCODING\n",RANK);       
                        // Wait 
                        free(data);
                        shot = shot + ds; // get next shot in series

                        MPI_Barrier(MPI_COMM_WORLD);

                    } // shots

                    /* Save my encoding */ 
                    save_encoding(receiver_encoding,
                                  encoding_prefix,
                                  encoding,
                                  ax, ay, aw);
                    
                    // Wait
                    MPI_Barrier(MPI_COMM_WORLD);
            }
        } // rounds
    }

    if(verb) fprintf(stderr,"%d READY TO FREE\n",RANK);
    
    free(**receiver_encoding); free(*receiver_encoding); free(receiver_encoding);
    if(verb) fprintf(stderr,"%d FREE ENCODING\n",RANK);
    free(*phase); free(phase);
    if(verb) fprintf(stderr,"%d FREE PHASE\n",RANK);
    free(*ampls); free(ampls);
    if(verb) fprintf(stderr,"%d FREE AMP\n",RANK);
    free(*emap); free(emap);
    if(verb) fprintf(stderr,"%d FREE MAP\n",RANK);
    free(encoding_prefix);
    free(data_prefix);
    if(verb) fprintf(stderr,"%d FINISHED FREE\n",RANK);
    MPI_Finalize();
    if(verb) fprintf(stderr,"%d FINISHED ALL\n",RANK);
}
