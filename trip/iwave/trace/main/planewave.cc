#include <except.hh>
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
  "Purpose: create a pair of (a) data file with zero samples and ",
  "correct trace headers, and (b) data file with plane wave traces",
  "and corresponding headers, for use in 2D synthetic plane wave studies",
  "",
  "Required parameters:",
  "  CWPROOT        = path to SU root directory",
  "  src            = source wavelet or pulse (single trace)",
  "  pwhdr          = multiple plane wave data file (zero samples) [OUTPUT]",
  "  pwsrc          = multiple plane wave source file, impulse pw traces convolved with src [OUTPUT]",
  "  nt             = number of time samples in pwhdr traces",
  "  ot             = time of first sample in pwhdr traces (ms)",
  "  nx             = number of traces in each pwhdr gather",
  "  dx             = x increment between traces in pwhdr gather (m)",
  "  ox             = x coord of first trace in pwhdr gather (m)",
  "  zs             = depth of source points (pwsrc)",
  "  zr             = depth of receiver points (pwhdr)",
  "  nts            = number of time samples in pwsrc traces",
  "  ots            = time of first sample in pwsrc traces (ms)",
  "  nxs            = number of traces in each pwsrc gather",
  "  dxs            = x increment between traces in pwsrc gather (m)",
  "  oxs            = x coord of first trace in pwsrc gather (m)",
  "  np             = number of plane waves = number of pwhdr, pwsrc gathers",
  "  dp             = slowness increment between plane waves (ms/m)",
  "  op             = slowness of first gather (ms/m)",
  NULL};  

using RVL::RVLException;
using RVL::valparse;
int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

    xargc=argc; xargv=argv;
    requestdoc(1);

    PARARRAY * par = ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }
         
    // extract dt
    std::string src = valparse<std::string>(*par,"src");
    std::string pwsrc = valparse<std::string>(*par,"pwsrc");
    std::string pwhdr = valparse<std::string>(*par,"pwhdr");

    FILE * fp = fopen(src.c_str(),"r");
    segy tr;
    fgettr(fp,&tr);
  
    //
    string cwp = valparse<string>(*par,"CWPROOT");
    int nt     = valparse<int>(*par,"nt");
    int nts    = valparse<int>(*par,"nts");
    int nx     = valparse<int>(*par,"nx");
    int nxs    = valparse<int>(*par,"nxs");
    int np     = valparse<int>(*par,"np");
    float ot   = valparse<float>(*par,"ot");
    float ots  = valparse<float>(*par,"ots");
    float ox   = valparse<float>(*par,"ox");
    float oxs  = valparse<float>(*par,"oxs");
    float op   = valparse<float>(*par,"op");
    float zs   = valparse<float>(*par,"zs");
    float zr   = valparse<float>(*par,"zr");
    float dx   = valparse<float>(*par,"dx");
    float dxs  = valparse<float>(*par,"dxs");
    float dp   = valparse<float>(*par,"dp");
    // note: dp = ms/m, so ms/trace needed by suplane (dip1) is dp*dxs

    // set up path
    string sunull=cwp;
    string suplane=cwp;
    sunull=sunull + "/bin/sunull";
    suplane=suplane + "/bin/suplane";

    // preclean
    string cln = "/bin/rm -f " + pwsrc + " " + pwhdr + ";";
    system(cln.c_str());
    // loop over plane waves
    for (int ip=0; ip<np; ip++) {
      // create header file
      std::stringstream hcmd;
      hcmd << "sunull nt="<<nt<<" ntr="<<nx<<" dt="<<tr.dt*1.e-6;
      hcmd <<"| sushw key=selev a="<<-zs<<" | sushw key=gelev a="<<-zr;
      hcmd <<"| sushw key=sx a="<<1000.0*(op+ip*dp)<<" | sushw key=gx a="<<1000.0*ox<<" b="<<1000.0*dx;
      hcmd <<"| sushw key=tracl a="<<ip*nx<<" b=1 | sushw key=tracr a=0 b=1";
      hcmd <<"| sushw key=scalco a=-1000 | sushw key=delrt a="<<ot;
      hcmd <<" >> " << pwhdr <<";";
      cerr << "ip="<<ip<<"\ncmd="<<hcmd.str()<<"\n";
      system(hcmd.str().c_str());
      std::stringstream scmd;
      scmd << "suplane nt="<<nts<<" ntr="<<nxs<<" dt="<<tr.dt*1.e-6<<" npl=1 len1="<<nxs<<" dip1="<<(op+ip*dp)*dxs;
      scmd << "| sushw key=sx a="<<1000.0*(op+ip*dp)<<" | sushw key=gx a="<<1000.0*oxs<<" b="<<1000.0*dxs;
      scmd << "| sushw key=tracl a="<<ip*nxs<<" b=1 | sushw key=tracr a=0 b=1";
      scmd << "| sushw key=scalco a=-1000 | sushw key=offset a=0 | sushw key=delrt a=" << ots;
      scmd << "| suconv sufile=" << src << " >> " << pwsrc <<";";
      cerr << "ip="<<ip<<"\ncmd="<<scmd.str()<<"\n";
      system(scmd.str().c_str());
    }
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
