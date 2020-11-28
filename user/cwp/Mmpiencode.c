/* shot encoding with arbitrary phase and amplitude weights using MPI on a distributed cluster 
 * axes are x-y-w

*/

#include <rsf.h>
#include <mpi.h>
#include <stdio.h>

#define PI 3.14159

void zero_array(sf_complex ***array, int n1, int n2, int n3);
void update_receiver_encoding(sf_complex ***encoding, sf_complex *data, float ow, float dw, float time, float* weight, float *phase, int n1, int n2, int n3);
void save_encoding(sf_complex ***encoding, char *prefix, int iencoding, sf_axis ax, sf_axis ay, sf_axis aw);


void zero_array(sf_complex ***array, int n1, int n2, int n3)
{
    int i3, i2, i1;

    for(i3 = 0; i3 < n3; ++i3){
        for(i2 = 0; i2 < n2; ++i2){
            for(i1 = 0; i1 < n1; ++i1){
                array[i3][i2][i1] = sf_cmplx(0.0f,0.0f);
            } // w
        } // y
    } // x 
}

void update_receiver_encoding(sf_complex ***encoding, sf_complex *data, float ow, float dw, float time, float* weight, float *phase, int n1, int n2, int n3)
{
    int i3, i2, i1, index;
    
//    cprint(cp);
    for(i3 = 0; i3 < n3; ++i3){
        float shift = -2.0*PI*(ow+dw*i3)*time;
        float pshift = -2.0*PI*phase[i3];
        sf_complex cp = sf_cmplx(cosf(pshift),sinf(pshift));

        sf_complex cw = sf_cmplx(weight[i3], 0.0f);
        sf_complex cs = sf_cmplx(cosf(shift),sinf(shift));
        sf_complex tot = cw*cp*cs;
        for(i2 = 0; i2 < n2; ++i2){
            for(i1 = 0; i1 < n1; ++i1){
                index = i1 + i2*n1 + i3*n2*n1; // linear index
                encoding[i3][i2][i1] += tot*data[index];
             } // w
        } // y
   } // x 
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
    
    sf_complex *data;    //1D array for reading and sending shot records
    
    sf_complex ***receiver_encoding; //output receiver encoding
    
    int   **emap;            // map of encodings to MPI processes
    float ***phase;          // array of phase shifts
    float **delays;         // array of time delays
    float ***ampls;          // array of amplitude weights
    float *weights;          // subset of ampls that corresponds to a single shot
    float *phase_shifts;          // subset of phase that corresponds to a single shot
    
    int ns,ne,os,ds,nx,ny,nw;
    float dw,ow;
    sf_axis ax,ay,aw,as,ae;
    
    sf_file Fencode, Fshotrecord;
 
    char *data_prefix = sf_charalloc(256);     //printf like string
    char *encoding_prefix = sf_charalloc(256); //printf like string for output
    
    int RANK;  // this processes rank
    int PROCS; // total number of processes initialized

    int ir, is, ip, iw;
    
    MPI_Init(&argc,&argv);
    
    MPI_Comm_size(MPI_COMM_WORLD,&PROCS);
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    
    sf_init(argc,argv); //Init RSF

    Fencode = sf_input("encode"); /* encoding file from sfencodemaker */
    
    data_prefix     = (char *)sf_getstring("dprefix"); /* printf like statement that can be evaluated to find the data files corresponding to shot records */
    encoding_prefix = (char *)sf_getstring("eprefix"); /* printf like statement that can be evaluated for the output encodings */
    if (data_prefix == NULL) sf_error("Must enter dprefix!");
    if (encoding_prefix == NULL) sf_error("Must enter eprefix!");
    
    if(! sf_getbool("verb",&verb)) verb = false;

    as = sf_iaxa(Fencode,1); 
    ns = sf_n(as); ds = (int)sf_d(as); os = (int)sf_o(as);
    ae = sf_iaxa(Fencode,2);
    ne = sf_n(ae);
    
    /* Read the axes from a single shot record for future reference */
    sf_file Faxes;
    char tname[256];
    sprintf(tname,data_prefix,os);
    if (verb && RANK==0) sf_warning("will try to read: %s", tname);
    Faxes = sf_input(tname);
    
    ax = sf_iaxa(Faxes,1); nx = sf_n(ax);
    ay = sf_iaxa(Faxes,2); ny = sf_n(ay);
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
    
    /* Map encodings to processes */
    emap = sf_intalloc2(PROCS,nrounds);
    int te = 0;
    for(ir = 0; ir < nrounds; ++ir){
        if (verb && RANK==0) fprintf(stderr,"ROUND %d .... ", ir);
        for(ip=0; ip < PROCS; ++ip){
            if(te < ne) {
                emap[ir][ip] = te;
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
    delays = sf_floatalloc2(ns,ne);
    phase = sf_floatalloc3(ns,ne,nw);
    ampls = sf_floatalloc3(ns,ne,nw);
    sf_floatread(ampls[0][0],nw*ns*ne,Fencode);
    sf_floatread(phase[0][0],nw*ns*ne,Fencode);
    sf_floatread(delays[0],ns*ne,Fencode);

    weights = sf_floatalloc(nw);
    phase_shifts = sf_floatalloc(nw);

    for(iw = 0; iw < nw; ++iw){
        weights[iw] = 0.0;
        phase_shifts[iw] = 0;
    }

    data              = sf_complexalloc(nx*ny*nw);
    receiver_encoding = sf_complexalloc3(nx,ny,nw);
    
    if (RANK == 0) {  
        for (ir = 0; ir < nrounds; ++ir){
            
            int encoding = emap[ir][RANK];
            int shot = os;

            zero_array(receiver_encoding,nx,ny,nw);

            for (is = 0; is < ns; ++is){
                fprintf(stderr,"Round %d ---- %d / %d ---- \r", ir, is+1, ns);
                    
                float timeDelay = delays[encoding][is];

                for(iw = 0; iw < nw; ++iw){
                    weights[iw] = ampls[iw][encoding][is];
                    phase_shifts[iw] = phase[iw][encoding][is]; 
                }

                // Which shot record am I reading?
                char recordName[256];
                sprintf(recordName,data_prefix,shot);
                Fshotrecord = sf_input(recordName); 

                //Read this shot record
                sf_complexread((sf_complex *)data,nx*ny*nw,Fshotrecord); 

                // Send shot record
                MPI_Bcast((float *)data,nx*ny*nw*2,
                    MPI_FLOAT,0,MPI_COMM_WORLD); 

                // Update my encoding
                update_receiver_encoding(receiver_encoding,
                                    (sf_complex *)data,
                                    ow,dw,timeDelay,weights,phase_shifts,
                                    nx,ny,nw);
                
                MPI_Barrier(MPI_COMM_WORLD);
               
                sf_fileclose(Fshotrecord);
                shot = shot + ds;
            } // shots

            /* Save the encoding */
            save_encoding(receiver_encoding,
                          encoding_prefix,
                          encoding,
                          ax, ay, aw);
            
            MPI_Barrier(MPI_COMM_WORLD);
        } // rounds
    } else {
        for (ir = 0; ir < nrounds; ++ir){
           
            int encoding = emap[ir][RANK]; // which encoding should I write?
            int shot = os;

            if (encoding < 0) { /* We aren't writing to an encoding */
                   for(is = 0; is < ns; ++is){
                        MPI_Bcast((float *)data,nx*ny*nw*2, 
                            MPI_FLOAT,0,MPI_COMM_WORLD);
                        // Wait
                        MPI_Barrier(MPI_COMM_WORLD);
                   }
                   MPI_Barrier(MPI_COMM_WORLD);
            }
            else {
                    zero_array(receiver_encoding,nx,ny,nw);

                    for(is = 0; is < ns; ++is){
                        
                        float timeDelay = delays[encoding][is];

                        for(iw = 0; iw < nw; ++iw){
                            weights[iw] = ampls[iw][encoding][is];
                            phase_shifts[iw] = phase[iw][encoding][is]; 
                        }

                        // Get shot record
                        MPI_Bcast((float *)data,nx*ny*nw*2,
                                MPI_FLOAT,0,MPI_COMM_WORLD);
                        
                        // Update my encoding
                        update_receiver_encoding(
                                            receiver_encoding,
                                            (sf_complex *)data,
                                            ow,dw,
                                            timeDelay,weights,phase_shifts,
                                            nx,ny,nw);

                        // Wait 
                        MPI_Barrier(MPI_COMM_WORLD);

                        shot = shot + ds; // get next shot in series
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
    
    free(**receiver_encoding); free(*receiver_encoding); free(receiver_encoding);
    free(**phase); free(*phase); free(phase);
    free(**ampls); free(*ampls); free(ampls);
    free(*delays); free(delays);
    free(*emap); free(emap);
    free(encoding_prefix);
    free(data_prefix);
    free(data);
    free(weights);
    free(phase_shifts);
    MPI_Finalize();
}
