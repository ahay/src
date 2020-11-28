/* stacks rsf files of the same dimensionality using mpi 

Mode specifies whether to add, multiply, divide or subtract.

mode=0 - add
mode=1 - multiply

If useprefix is set, then:

assume that files are commonly named sequentially, e.g.:

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

This program does not care about dimensionality!  It treats every file
as a 1D array and writes out a 1D array, and then modifies the header
to match the input file size.  

*/

#include <rsf.h>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv){

    MPI_Init(&argc,&argv);
    
    char *prefix = sf_charalloc(1024);
    char *outName = sf_charalloc(1024);
    
    int nf,of=0,jf=1;
    
    int RANK;
    int PROCS;
    
    bool verb, seq;

    int mode;
    char **filenames = NULL;

    int i, j, ir, ip;
    
    MPI_Comm_size(MPI_COMM_WORLD,&PROCS);
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
        
    sf_init(argc,argv);
   
    if (! sf_getint("mode",&mode)) mode=0; /*operation for stack*/
    sf_warning("FOUND %d mode",mode);
    if (! sf_getbool("verb",&verb)) verb=false;

    if (! sf_getbool("seq",&seq)) seq=false; /* not sequentially ordered files*/

    char *shotfile = (char *) sf_getstring("shots");
    if (shotfile != NULL){
        sf_warning("Using shot numbers from shot file: %s", shotfile);
        prefix = (char *)sf_getstring("prefix"); /* printf like prefix */
        if (prefix == NULL) sf_error("Must specify prefix");

        sf_file Fshotnumbers = sf_input(shotfile);
        sf_axis Ashots = sf_iaxa(Fshotnumbers,2);

        nf = sf_n(Ashots);
        int *shotnumbers = sf_intalloc(nf);

        sf_intread(shotnumbers,nf,Fshotnumbers);

        filenames = (char **) sf_alloc((size_t) nf, sizeof(char*));

        for(i = 0; i < nf; ++i){
            char *tfilename = sf_charalloc(1024);
            sprintf(tfilename,prefix,shotnumbers[i]);
            filenames[i] = tfilename;
            if ( verb && RANK == 0) sf_warning("Will use file: %s", tfilename);
        }

        sf_warning("Finished reading shot-file");

    }
    else if (seq){
        sf_warning("Using sequentially ordered input files");
        if (! sf_getint("nf",&nf)) sf_error("Must specify how many files to stack"); /* number of files to stack */
        if (! sf_getint("jf",&jf)) jf=1; /* delta between files */
        if (! sf_getint("of",&of)) of=0; /* origin of files*/

        prefix = (char *)sf_getstring("prefix"); /* printf like prefix */
        if (prefix == NULL) sf_error("Must specify prefix");

        filenames = (char**) sf_alloc ((size_t) nf,sizeof(char*));
        for (i = of; i < nf*jf+of; i = i + jf){
            char *tfilename = sf_charalloc(1024);
            sprintf(tfilename,prefix,i);
            filenames[i] = tfilename;
        }
    } else {
        sf_warning("Using filenames from command line");
        int nfiles = 0;
        for (i = 1; i < argc; ++i){
            if (strchr(argv[i],'=') == NULL) { /* did not find a key=val pair*/
                nfiles++;
            } else {
                continue;
            }
        }
        if (verb && RANK==0) sf_warning("Found %d files", nfiles);
        filenames = (char**) sf_alloc ((size_t) nfiles,sizeof(char*));

        int ifile = 0;
        for(j = 1; j < argc; ++j){
            if (strchr(argv[j],'=') == NULL) { /* did not find a key=val pair*/
                filenames[ifile] = argv[j];
                if (verb && RANK==0) sf_warning("file %d - %s",ifile,argv[j]);
                ifile++;
            }

        }

        nf = nfiles;
        of = 0;
        jf = 1;
    }
    
   
    outName    = (char *)sf_getstring("oname"); /* name of output file */
    if (outName== NULL) sf_error("Must specify output name");
    
    int nrounds = nf/PROCS;
    if (nf % PROCS != 0) nrounds++;
    
    int **MPI_MAP = sf_intalloc2(PROCS,nrounds);
    int file = of;
    for(ir = 0; ir < nrounds; ++ir){
        if (verb && RANK==0) fprintf(stderr,"ROUND %d .... ", ir);
       
        for(ip=0; ip < PROCS; ++ip){
            if(file < of+jf*nf) {
                MPI_MAP[ir][ip] = file;
                file += jf;
            } else {
               MPI_MAP[ir][ip] = -1;
            }
            if (verb && RANK==0) fprintf(stderr,"%d ", MPI_MAP[ir][ip]);
        }
        if (verb && RANK==0) fprintf(stderr,"\n");
    }
    if(verb && RANK == 0) fprintf(stderr, "\n");
    
    
    char *filename = NULL;

    sf_file tfile = sf_input(filenames[0]);
    sf_file output = NULL;
    if (RANK == 0) output = sf_output(outName);
    int fsize = sf_filesize(tfile);
    
    sf_datatype type = sf_gettype(tfile);
    
    float *fsarray = NULL;
    float *frarray = NULL; float *foarray = NULL;
    int   *isarray = NULL;
    int   *irarray = NULL; int   *ioarray = NULL;
    
    switch (type) {
        case SF_FLOAT:
            frarray = sf_floatalloc(fsize);
            fsarray = sf_floatalloc(fsize);
            if (RANK == 0) foarray = sf_floatalloc(fsize);
            for(i = 0; i < fsize; ++i){
                frarray[i] = 0;
                fsarray[i] = 0;
                if (RANK == 0) foarray[i] = 0;
            }
            break;
        case SF_INT:
            irarray = sf_intalloc(fsize);
            isarray = sf_intalloc(fsize);
            if (RANK == 0) ioarray = sf_intalloc(fsize);
            for(i = 0; i < fsize; ++i){
                irarray[i] = 0;
                isarray[i] = 0;
                if (RANK == 0) ioarray[i] = 0;
            }
            break;
        case SF_COMPLEX:
            sf_error("Cannot read complex files");
            break;
            
        case SF_DOUBLE:
        case SF_UCHAR:
        case SF_CHAR:
        case SF_SHORT:
        default:
            sf_error("Unknown data type, must be int or float");
            break;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    for(ir = 0; ir < nrounds; ++ir){
        
        int fnumber = MPI_MAP[ir][RANK];
        
        if (fnumber < 0) {
                  
            MPI_Barrier(MPI_COMM_WORLD);
            
            switch (type){
                case SF_FLOAT:
                    
                    for(i = 0; i < fsize; ++i){
                        fsarray[i] = 0;
                        frarray[i] = 0;
                    }
                    if (mode == 0) {
                        MPI_Reduce(fsarray,frarray, fsize,
                            MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
                    } else {
                        MPI_Reduce(fsarray,frarray, fsize,
                            MPI_FLOAT,MPI_PROD,0,MPI_COMM_WORLD);
                    }
                        
                    if (RANK == 0){
                        for(i = 0; i < fsize; ++i){
                            foarray[i] += frarray[i];
                        }
                    }
                    break;
                case SF_INT:
                    for(i = 0; i < fsize; ++i){
                        isarray[i] = 0;
                    }
                    if (mode == 0){
                        MPI_Reduce(isarray,irarray, fsize,
                            MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
                    } else {
                        MPI_Reduce(isarray,irarray, fsize,
                            MPI_INT,MPI_PROD,0,MPI_COMM_WORLD);
                    }
                        
                    if (RANK == 0){
                        for(i = 0; i < fsize; ++i){
                            ioarray[i] += irarray[i];
                        }
                    }
                    break;
                case SF_DOUBLE:
                case SF_UCHAR:
                case SF_CHAR:
                case SF_SHORT:
                case SF_COMPLEX:
                default:
                    break;
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
        } else {
            filename = filenames[fnumber];
            if (verb) sf_warning("%d reading... %s",RANK,filename);
            sf_file file = sf_input(filename);
            
            if (verb && RANK == 0) 
                {fprintf(stderr,"Round:%d/%d \r", ir,nrounds);}
            
            MPI_Barrier(MPI_COMM_WORLD);
            switch (type){
                case SF_FLOAT:
                    sf_floatread(fsarray,fsize,file);

                    if (mode == 0) {
                        MPI_Reduce(fsarray,frarray, fsize,
                            MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
                    } else {
                        MPI_Reduce(fsarray,frarray, fsize,
                            MPI_FLOAT,MPI_PROD,0,MPI_COMM_WORLD);
                    }
                    if (RANK == 0){
                        for(i = 0; i < fsize; ++i){
                            foarray[i] += frarray[i];
                        }
                    }
                    break;
                case SF_INT:
                    sf_intread(isarray,fsize,file);
                    if (mode == 0) {
                        MPI_Reduce(isarray,irarray, fsize,
                            MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
                    } else {
                        MPI_Reduce(isarray,irarray, fsize,
                            MPI_INT,MPI_PROD,0,MPI_COMM_WORLD);
                    }

                    if (RANK == 0){
                        for(i = 0; i < fsize; ++i){
                            ioarray[i] += irarray[i];
                        }
                    }
                    break;
                    
                case SF_DOUBLE:
                case SF_UCHAR:
                case SF_CHAR:
                case SF_SHORT:
                case SF_COMPLEX:
                default:
                    break;
            }
            
            sf_fileclose(file);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    if (RANK == 0) {
        sf_warning("Finished reducing, writing out to %s",outName);
        //sprintf(tfilename,prefix,of);
        sf_fileflush(output,tfile);
   
        switch(type){
            case SF_FLOAT:
                sf_warning("Writing %d floats",fsize);
                sf_floatwrite(&foarray[0],fsize,output);
                free(foarray);
                break;
            case SF_INT:
                sf_warning("Writing %d ints",fsize);
                sf_intwrite(&ioarray[0],fsize,output);
                free(ioarray);
                break;
                
            case SF_DOUBLE:
            case SF_UCHAR:
            case SF_CHAR:
            case SF_SHORT:
            case SF_COMPLEX:
            default:
                break;
        }
        
        sf_fileclose(output);

        sf_warning("Done");       
    } 

    MPI_Barrier(MPI_COMM_WORLD);

    switch(type){
        case SF_FLOAT:
            free(fsarray);
            free(frarray);
            break;
        case SF_INT:
            free(isarray);
            free(irarray);
            break;
        default:
            break;
    }

    sf_fileclose(tfile);

    if (seq) free(prefix);

    /*for(i = 0; i < nf; ++i){
        free(filenames[i]);
    }*/
    free(filenames);
    free(outName);   
    MPI_Finalize();
}
