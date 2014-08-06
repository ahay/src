#include "iwinfo.hh"
#include "grid.h"
#include "acdpml_defn.hh"
#include "gtest/gtest.h"
#include "acdpml_defn.hh"
#include "traceio.h"
#include "istate.hh"

using namespace std;

using TSOpt::IWaveEnvironment;


//#define GTEST_VERBOSE

IOKEY IWaveInfo::iwave_iokeys[]
= {
    {"csq",    0, true,  true },
    {"data",   1, false, true },
    {"movie",  1, false, false},
    {"init",   1, true,  false},
    {"",       0, false, false}
};

namespace {
    
    using RVL::parse;
    using RVL::RVLException;
    using TSOpt::IWaveEnvironment;
    using TSOpt::IWaveTree;
    using TSOpt::IWaveSampler;
    using TSOpt::IWaveSim;
    using TSOpt::TASK_RELN;
    using TSOpt::IOTask;
    
    void create_fixed_2D_data(string fnm,
                              int nt,
                              float dt,
                              float ot,
                              int nrec,
                              int ntr_in_rec,
                              float sx0,
                              float dsx,
                              float sz,
                              float rx0,
                              float drx,
                              float rz,
                              int scalel,
                              int scalco) {
        
        string CWPROOT = getenv("CWPROOT");
        stringstream cmd;
        cmd << CWPROOT + "/bin/sunull nt=" <<nt<<" ntr="<<ntr_in_rec*nrec<<" dt="<<dt;
        cmd << "| " + CWPROOT + "/bin/sushw key=delrt a="<<ot;
        cmd << "| " + CWPROOT + "/bin/sushw key=sx a="<<sx0<<" c="<<dsx<<" j="<<ntr_in_rec;
        cmd << "| " + CWPROOT + "/bin/sushw key=selev a=" << (-1.0f)*sz;
        cmd << "| " + CWPROOT + "/bin/sushw key=gx a="<<rx0<<" b="<<drx<<" j="<<ntr_in_rec;
        cmd <<" | " + CWPROOT + "/bin/sushw key=gelev a="<<(-1.0f)*rz;
        cmd << "| " + CWPROOT + "/bin/sushw key=scalel a="<<scalel;
        cmd << "| " + CWPROOT + "/bin/sushw key=scalco a="<<scalco;
        cmd << "| " + CWPROOT + "/bin/suchw key1=offset key2=gx key3=sx b=1 c=-1 > "<<fnm<<"\n";
        
        if (system(cmd.str().c_str())) {
            RVLException e;
            e<<"Error: create_fixed_2D_data\n";
            e<<"  failed system call on command \n";
            e<<cmd.str();
            throw e;
        }
        
    }
    
    void create_hfile(string hfile, string dfile, grid g, float val, bool var=false) {
        
        if (hfile.size()>120) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"filename "<<hfile<<" longer than 120 chars\n";
            throw e;
        }
        
        char * fname=(char *)malloc(128*sizeof(char));
        FILE * fp = NULL;
        
        strcpy(fname,hfile.c_str());
        
        fp = iwave_fopen(&fname,"w",NULL,stderr);
        if (!fp) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"file testgrid.rsf not opened\n";
            throw e;
        }
        
        fprint_grid(fp,g);
        
        fprintf(fp,"data_format=native_float\n");
        
        fprintf(fp,"data_type = csq\n");
        fprintf(fp,"in=%s\n",dfile.c_str());
        
        fflush(fp);
        
        iwave_fclose(fp);
        
        strcpy(fname,dfile.c_str());
        float * buf = (float *)malloc(get_datasize_grid(g)*sizeof(float));
        if (var) {
            for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val/((float)(i+1));
        }
        else {
            for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val;
        }
        //      FILE *
        fp = iwave_fopen(&fname,"w",NULL,stderr);
        if (!fp) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"file "<<fname<<" not opened\n";
            throw e;
        }
        for (int i=0;i<get_panelnum_grid(g);i++)
            fwrite(buf,sizeof(float),get_datasize_grid(g),fp);
        
        free(buf);
        fflush(fp);
        
        iwave_fclose(fp);
        free(fname);
    }
    
    void create_gauss(grid g, float lam) {
        
        char * fname=(char *)malloc(128*sizeof(char));
        FILE * fp = NULL;
        
        strcpy(fname,"gauss.rsf");
        
        fp = iwave_fopen(&fname,"w",NULL,stderr);
        if (!fp) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"file gauss.rsf not opened\n";
            throw e;
        }
        
        fprint_grid(fp,g);
        
        fprintf(fp,"data_format=native_float\n");
        
        fprintf(fp,"data_type = csq\n");
        fprintf(fp,"in=%s\n","gauss.rsf@");
        
        fflush(fp);
        
        iwave_fclose(fp);
        
        strcpy(fname,"gauss.rsf@");
        
        // compute center
        RPNT ctr;
        float radsq=0.0f;
        
        for (int i=0;i<g.dim;i++) {
            ctr[i] = g.axes[i].o + 0.5 * g.axes[i].d * g.axes[i].n;
            radsq += (0.5 * g.axes[i].d * g.axes[i].n)*(0.5 * g.axes[i].d * g.axes[i].n);
        }
        float x;
        float y;
        float z;
        float rsq;
        float * buf = (float *)malloc(get_datasize_grid(g)*sizeof(float));
        if (g.dim==2) {
            //	cerr<<"g.axes[0].n="<<g.axes[0].n<<" g.axes[1].n="<<g.axes[1].n<<endl;
            //	cerr<<"g.axes[0].o="<<g.axes[0].o<<" g.axes[1].o="<<g.axes[1].o<<endl;
            //	cerr<<"g.axes[0].d="<<g.axes[0].d<<" g.axes[1].d="<<g.axes[1].d<<endl;
            //	cerr<<"radsq="<<radsq<<endl;
            for (int j=0;j<g.axes[1].n;j++) {
                for (int i=0;i<g.axes[0].n;i++) {
                    z=g.axes[0].o + i*g.axes[0].d;
                    x=g.axes[1].o + j*g.axes[1].d;
                    rsq=(z-ctr[0])*(z-ctr[0])+(x-ctr[1])*(x-ctr[1]);
                    //	    cerr<<"rsq/radsq="<<rsq/radsq<<endl;
                    buf[i+j*g.axes[0].n]=exp(-lam*lam*rsq/radsq);
                }
            }
        }
        else if (g.dim==3) {
            for (int k=0;k<g.axes[2].n;k++) {
                for (int j=0;j<g.axes[1].n;j++) {
                    for (int i=0;i<g.axes[0].n;i++) {
                        z=g.axes[0].o + i*g.axes[0].d;
                        x=g.axes[1].o + j*g.axes[1].d;
                        y=g.axes[2].o + k*g.axes[2].d;
                        buf[i+j*g.axes[0].n]=exp(-lam*lam*((z-ctr[0])*(z-ctr[0])+(x-ctr[1])*(x-ctr[1])+(y-ctr[2])*(y-ctr[2])));
                    }
                }
            }  
        }
        else {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"  wake up and smell the roses\n";
            throw e;
        }
        //      FILE * 
        fp = iwave_fopen(&fname,"w",NULL,stderr);
        if (!fp) {
            RVLException e;
            e<<"Error: create_hfile\n";
            e<<"file "<<fname<<" not opened\n";
            throw e;
        }
        for (int i=0;i<get_panelnum_grid(g);i++) 
            fwrite(buf,sizeof(float),get_datasize_grid(g),fp);
        
        free(buf);
        fflush(fp);
        
        iwave_fclose(fp);
        free(fname);
    }

    class ACDPMLTest : public ::testing::Test {
    public:
        
        string hfile;
        string dfile;
        string hfilem;
        string dfilem;
        string hfilei;
        string dfilei;
        
        IWaveInfo ic;
        
        ACDPMLTest(): ic() {
            
            int nt;
            float dt;
            float ot;
            int nrec;
            int ntr_in_rec;
            float sx0;
            float dsx;
            float sz;
            float rx0;
            float drx;
            float rz;
            int scalel;
            int scalco;
            
            // first segy - single gather
            nt = 101;
            dt = 0.004;
            ot = 0.0;
            nrec = 3;
            ntr_in_rec = 101;
            sx0=3300;
            dsx=100;
            sz=40;
            rx0=100;
            drx=20;
            rz=50;
            scalel=0;
            scalco=0;
            
            string dn = "data.su";
            
            create_fixed_2D_data(dn,nt,dt,ot,
                                 nrec,ntr_in_rec,
                                 sx0,dsx,sz,
                                 rx0,drx,rz,
                                 scalel,scalco);
            
            // segcond segy - source pulses
            nt = 51;
            dt = 0.004;
            ot = -100.0;
            nrec = 3;
            ntr_in_rec = 1;
            sx0=3300;
            dsx=100;
            sz=40;
            rx0=0.0;
            drx=0.0;
            rz=20;
            scalel=0;
            scalco=0;
            
            string wn = "wavelet.su";
            
            create_fixed_2D_data(wn,nt,dt,ot,
                                 nrec,ntr_in_rec,
                                 sx0,dsx,sz,
                                 rx0,drx,rz,
                                 scalel,scalco);
            
            // hfile etc. = simple 2D
            // hfile1 etc. = x, z reversed
            // hfile2      = gdim=3, n3=3 etc.
            // hfile3      = gdim=3, h=1, z=2, x=3
            hfile = "csq_4layer.rsf";
            dfile = "csq_4layer.rsf@";
            hfilem= "movie.rsf";
            dfilem= "movie.rsf@";
            hfilei= "init.rsf";
            dfilei= "init.rsf@";

       
            
            ofstream qs("initfile");
            qs<<"INPUT DATA FOR IWAVE\n";
            qs<<"------------------------------------------------------------------------\n";
            qs<<"FD:\n";
            qs<<"\n";
            qs<<"         order = 1           scheme half-order\n";
            qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
            qs<<"          cmin = 1.0         min velocity - checked\n";
            qs<<"          cmax = 2.0         max velocity - checked\n";
            qs<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
            qs<<"         fpeak = 0.010       nominal central frequency \n";
            qs<<"            dt = 0.01";
            qs<<"\n";
            qs<<"------------------------------------------------------------------------\n";
            qs<<"Model info:\n";
            qs<<"\n";
            qs<<"           csq = csq_4layer.rsf\n";
            qs<<"          init = gauss.rsf\n";
            qs<<"         movie = movie.rsf\n";
            qs<<"\n";
            qs<<"------------------------------------------------------------------------\n";
            qs<<"PML info:\n";
            qs<<"\n";
            qs<<"           LZ = 40\n";
            qs<<"           LX = 40\n";
            qs<<"         pmlampl = 40.0   set to 0.0 turn off pml\n";
            qs<<"\n";
            qs<<"------------------------------------------------------------------------\n";
            qs<<"MPI info:\n";
            qs<<"\n";
            qs<<"       mpi_np1 = 1      n_doms along axis 1\n";
            qs<<"       mpi_np2 = 1      n_doms along axis 2\n";
            qs<<"       mpi_np3 = 1      n_doms along axis 3\n";
            qs<<"       partask = 1      task parallelization\n";
            qs<<"\n";
            qs<<"------------------------------------------------------------------------\n";
            qs<<"Output info:\n";
            qs<<"\n";
            qs<<"     printact = 1           per-time-step verbosity level\n";
            qs<<"                            0 - none\n";
            qs<<"                            1 - time step index\n";
            qs<<"                            2 - internal time step info\n";
            qs<<"                            > 5: dump everything\n";
            qs<<"      dump_pi = 1           dump parallel/dom. decomp info\n";
            qs<<"     dump_lda = 1           dump grid data for allocated arrays\n";
            qs<<"     dump_ldc = 1           dump grid data for computational arrays\n";
            qs<<"     dump_lds = 1           dump grid data for send arrays\n";
            qs<<"     dump_ldr = 1           dump grid data for receive arrays\n";
            qs<<"    dump_term = 0           dump terminator data\n";
            qs<<"    dump_pars = 0           print parameter table in IWaveOp\n";
            qs<<"   dump_steqs = 0           print major steqs in IWaveOp\n";
            qs.flush();
            qs.close();
            
            grid g;
            g.dim=2;
            g.gdim=2;
            g.axes[0].n=416;
            g.axes[1].n=400;//800;
            //g.axes[2].n=1;
            g.axes[0].d=0.025;
            g.axes[1].d=0.025;
            //g.axes[2].d=100.0;
            g.axes[0].o=0.0;
            g.axes[1].o=0.0;
            //g.axes[2].o=-100.0;
            g.axes[0].id = 0;
            g.axes[1].id = 1;
            //g.axes[2].id = 3;
            
            float val = 1.5;
            create_hfile(hfile, dfile, g, val);
            
            g.dim=2;
            g.gdim=2;
            g.axes[0].n=416;
            g.axes[1].n=400;//800;
            g.axes[0].d=0.025;
            g.axes[1].d=0.025;
            g.axes[0].o=0.0;
            g.axes[1].o=0.0;
            g.axes[0].id = 0;
            g.axes[1].id = 1;
            
            val = 1.0;
            create_hfile(hfilei, dfilei, g, val);

            
            g.dim=2;
            g.gdim=3;
            g.axes[0].n=416;
            g.axes[1].n=400;//800;
            g.axes[2].n=20;
            g.axes[0].d=0.025;
            g.axes[1].d=0.025;
            g.axes[2].d=0.5;
            g.axes[0].o=0.0;
            g.axes[1].o=0.0;
            g.axes[2].o=0.01;
            g.axes[0].id = 0;
            g.axes[1].id = 1;
            g.axes[2].id = 2;
            
            val = 0.0;
            create_hfile(hfilem, dfilem, g, val);
            
            g.dim=2;
            g.gdim=3;
            g.axes[0].n=416;
            g.axes[1].n=400;//800;
            g.axes[2].n=1;
            g.axes[0].d=0.025;
            g.axes[1].d=0.025;
            g.axes[2].d=2.0;
            g.axes[0].o=0.0;
            g.axes[1].o=0.0;
            g.axes[2].o=0.0;
            g.axes[0].id = 0;
            g.axes[1].id = 1;
            g.axes[2].id = 2;
            
            float lam = 10;
            create_gauss(g,lam);
            
        }
    };
    TEST_F(ACDPMLTest, setup_environment) {
        try {
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=initfile");
            PARARRAY * par = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &par, &stream);
            delete [] argv;
            delete [] argvv;
#ifdef GTEST_VERBOSE
            ps_printall(*par,stderr);
#endif

            // check one of the main par entries
            string csq="";
            parse(*par,"csq",csq);
            EXPECT_EQ("csq_4layer.rsf",csq);
            
            IWAVE state0;
            int err;
            IWaveInfo ic;
            // construct IWAVE state using acd_modelinit
            err=iwave_construct(&state0,par,stream,ic);
            
            if (err) {
                fprintf(stream,"ERROR: main from iwave_construct. ABORT\n");
                iwave_destroy(&state0,ic.get_mdest());
            }
            
            iwave_destroy(&state0,ic.get_mdest());
            
            ps_delete(&par);
            fclose(stream);
            
#ifndef GTEST_VERBOSE
            //      unlink("parfile");
#endif
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }

TEST_F(ACDPMLTest, dryrun_sim_fwd_ord0) {
     try {
     
     // fake command line environment
     int argc = 2;
     char * argvv = new char[128];
     char ** argv = new char*[2];
     argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
     strcpy(argv[1],"par=initfile");
     PARARRAY * par = NULL;
     FILE * stream = NULL;
     IWaveEnvironment(argc, argv, 0, &par, &stream);
     delete [] argv;
     delete [] argvv;
     
     // build order zero IWaveTree, check
     int order=0;
     bool fwd = true;
     int printact=0; int snaps=0;
     
     bool dryrun=true;
     ofstream drystr("dryrun_sim_fwd_ord0");
     IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
     sim->run();
     drystr.close();
     
     delete sim;
     ps_delete(&par);
     fclose(stream);
     }
     catch (RVLException & e) {
     e.write(cerr);
     exit(1);
     }
}
    
TEST_F(ACDPMLTest, run_sim_fwd_ord0) {
    try {
        
        // fake command line environment
        int argc = 2;
        char * argvv = new char[128];
        char ** argv = new char*[2];
        argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
        strcpy(argv[1],"par=initfile");
        PARARRAY * par = NULL;
        FILE * stream = NULL;
        IWaveEnvironment(argc, argv, 0, &par, &stream);
        delete [] argv;
        delete [] argvv;
        
        // build order zero IWaveTree, check
        int order=0;
        bool fwd = true;
        int printact=1; int snaps=0;
        
        bool dryrun=false;
        ofstream drystr("run_sim_fwd_ord0");
        IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
        sim->run();
        drystr.close();
        
        delete sim;
        ps_delete(&par);
        fclose(stream);
    }
    catch (RVLException & e) {
        e.write(cerr);
        exit(1);
    }
}

}

int xargc;
char **xargv;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int err = RUN_ALL_TESTS();
    iwave_fdestroy();
    return err;
}
