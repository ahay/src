#include <except.hh>
#include <utility.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

const char * sdoc[] = { 
  "Usage: planewave.x CWPROOT= src= pwhdr= pwsrc= ",
  "       nt= ot= nx= dx= ox= zs= zr=",
  "       nts= ots= nxs= dxs= oxs=",
  "       np= dp= op=",
  " ",
  "Purpose: create a pair of (a) SU data file with zero samples and ",
  "correct trace headers for plane wave modeling output, and (b) SU data",
  "file with source traces for plane wave modeling input",
  "and corresponding headers, for use in 2D synthetic plane wave studies",
  "",
  "Notes:",
  "",
  "1.This utility delegates to SU commands, hence depends on proper",
  "installation of SU. In particular, CWPROOT must be defined in the",
  "execution environment and point to the SU root directory. CWPROOT",
  "passed as parameter in order to work in the brain-dead environment",
  "created by Flow",
  "",
  "2. suconv creates a minimal time axis of the source traces, based on",
  "min and max slowness, pivot position for impulse plane waves,",
  "horizontal extent of the source array, and the length of the",
  "source wavelet.",
  "",
  "3. Time step (dt) for output files is that of input file (source pulse)",
  "",
  "4. TO DO: add option for taking receiver positions from file, to",
  "permit nonuniform trace spacing",
  "",
  "5. TO DO; 3D-ize",
  "",
  "Required parameters:",
  "  CWPROOT        = path to SU root directory",
  "  src            = filename for source wavelet or pulse (single trace)",
  "  pwhdr          = filename for multiple plane wave data file (zero samples) [OUTPUT]",
  "  pwsrc          = filename for multiple plane wave source file = impulse pw traces ",
  "                   convolved with src [OUTPUT]",
  "  nt             = number of time samples in pwhdr traces",
  "  ot             = time of first sample in pwhdr traces (ms)",
  "  nx             = number of traces in each pwhdr gather",
  "  dx             = x increment between traces in pwhdr gather (m)",
  "  ox             = x coord of first trace in pwhdr gather (m)",
  "  zr             = depth of receiver points (pwhdr)",
  "  zs             = depth of source points (pwsrc)",
  "  nxs            = number of traces in each pwsrc gather",
  "  dxs            = x increment between traces in pwsrc gather (m) (must be > 0)",
  "  oxs            = x coord of first trace in pwsrc gather (m)",
  "  np             = number of plane waves = number of pwhdr, pwsrc gathers",
  "  dp             = slowness increment between plane waves (ms/m)",
  "  op             = slowness of first gather (ms/m)",
  NULL};  

using RVL::RVLException;
using RVL::valparse;
using RVL::ProtectedDivision;
int xargc_;
char **xargv_;

int main(int argc, char ** argv) {

  try {

    xargc_=argc; xargv_=argv;
    requestdoc(1);

    PARARRAY * par = ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }
         
    ps_printall(*par,stderr);

    // test for presence of SU
    std::string cwp = valparse<std::string>(*par,"CWPROOT");
    int nt     = valparse<int>(*par,"nt");
    int nx     = valparse<int>(*par,"nx");
    int nxs    = valparse<int>(*par,"nxs");
    int np     = valparse<int>(*par,"np");
    float ot   = valparse<float>(*par,"ot");
    float ox   = valparse<float>(*par,"ox");
    float oxs  = valparse<float>(*par,"oxs");
    float op   = valparse<float>(*par,"op");
    float zs   = valparse<float>(*par,"zs");
    float zr   = valparse<float>(*par,"zr");
    float dx   = valparse<float>(*par,"dx");
    float dxs  = valparse<float>(*par,"dxs");
    if (!(dxs>0)) {
      RVLException e;
      e<<"Error: planewave\n";
      e<<"  source trace spacing (dxs) must be > 0\n";
      throw e;
    }
    float dp   = valparse<float>(*par,"dp");

    // extract time step from source wavelet, convert unit from mus to s
    std::string src = valparse<std::string>(*par,"src");
    FILE * fp = fopen(src.c_str(),"r");
    segy tr;
    fgettr(fp,&tr);
    // dt in ms
    float dt = tr.dt*1.e-3;

    // note: dp = ms/m, so ms/trace needed by suplane (dip1) is dp*dxs
    // compute time axis for impulse response
    float pmin = iwave_min(op,op+(np-1)*dp);
    float pmax = iwave_max(op,op+(np-1)*dp);
    // min time of plane wave with slope pmin
    float ptmin = pmin*dxs*nxs/2;
    float ptmax = pmax*dxs*nxs/2;
    float tmin = iwave_min(ptmin,iwave_min(-ptmin,iwave_min(ptmax,-ptmax)));
    float tmax = iwave_max(ptmin,iwave_max(-ptmin,iwave_max(ptmax,-ptmax)));
    cerr<<"dt="<<dt<<" pmin="<<pmin<<" pmax="<<pmax<<" ptmin="<<ptmin<<" ptmax="<<ptmax<<" tmin="<<tmin<<" tmax="<<tmax<<endl;
    float ots  = dt*(int(tmin/dt)-1);
    int nts    = int((tmax-tmin)/dt) + 1;
    // center pivot index, for suplane

    // set up paths to SU commands
    string sunull=cwp;
    string suplane=cwp;
    string sushw=cwp;
    string suconv=cwp;
    sunull=sunull + "/bin/sunull";
    suplane=suplane + "/bin/suplane";
    sushw=sushw + "/bin/sushw";
    suconv=suconv + "/bin/suconv";
      
    // identify outputs
    std::string pwsrc = valparse<std::string>(*par,"pwsrc");
    std::string pwhdr = valparse<std::string>(*par,"pwhdr");

    // preclean
    string cln = "/bin/rm -f " + pwsrc + " " + pwhdr + ";";
    if (system(cln.c_str()) == -1)
      abort();
    // loop over plane waves
    // note need to convert from ms to s in sunull, suplane commands
    for (int ip=0; ip<np; ip++) {
      // create header file
      std::stringstream hcmd;
      hcmd << sunull <<" nt="<<nt<<" ntr="<<nx<<" dt="<<dt*.001;
      hcmd <<"| " << sushw <<" key=selev a="<<-zs*1000<<" | " << sushw << " key=gelev a="<<-zr*1000;
      hcmd <<"| " << sushw <<" key=sx a="<<1000.0*(op+ip*dp)<<" | " << sushw <<" key=gx a="<<1000.0*ox<<" b="<<1000.0*dx;
      hcmd <<"| " << sushw <<" key=tracl a="<<ip*nx<<" b=1 | " << sushw <<" key=tracr a=0 b=1";
      hcmd <<"| " << sushw <<" key=scalco a=-1000 | "<< sushw <<" key=scalel a=-1000 | " << sushw << " key=delrt a="<<ot;
      hcmd <<" >> " << pwhdr <<";";
      cerr << "ip="<<ip<<"\ncmd="<<hcmd.str()<<"\n";
      if (system(hcmd.str().c_str()) == -1)
        abort();
      std::stringstream scmd;
      scmd << suplane << " nt="<<nts<<" ntr="<<nxs<<" dt="<< dt*.001 <<" npl=1 ct1="<<nts/2+1<<" len1="<<nxs<<" dip1="<<(op+ip*dp)*dxs << " cx1="<<nxs/2+1;
      scmd << "| " << sushw << " key=gelev a="<<-zs*1000 << " | " << sushw <<" key=scalel a=-1000 ";
      scmd << "| " << sushw <<" key=sx a="<<1000.0*(op+ip*dp)<<" | " << sushw <<" key=gx a="<<1000.0*oxs<<" b="<<1000.0*dxs;
      scmd << "| " << sushw <<" key=tracl a="<<ip*nxs<<" b=1 | " << sushw <<" key=tracr a=0 b=1";
      scmd << "| " << sushw <<" key=scalco a=-1000 | " << sushw <<" key=offset a=0 | " << sushw <<" key=delrt a=" << ots;
      scmd << "| " << suconv << " sufile=" << src << " >> " << pwsrc <<";";
      cerr << "ip="<<ip<<"\ncmd="<<scmd.str()<<"\n";
      if (system(scmd.str().c_str()) == -1)
        abort();
    }

    fclose(fp);

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
