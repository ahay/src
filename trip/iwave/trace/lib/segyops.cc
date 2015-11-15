#include "segyops.hh"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

float mutefun(float t) {
    t = MIN(1.0f,MAX(t,0.0f));
    return 3*t*t-2*t*t*t;
}
float cosfun(float t){
    return 0.5 - cos(M_PI * t)/2.0f;
}
float cosfun2(float t){
    return 0.5 + cos(M_PI * t)/2.0f;
}
namespace TSOpt {
    
    using RVL::LocalDataContainer;
    using RVL::RVLException;
    
    void SEGYLinMute::operator()(LocalDataContainer<float> & x,
                                 LocalDataContainer<float> const & y) {
        
        try {
            
            segytrace & cpout = dynamic_cast<segytrace &>(x);
            segytrace const & cpin = dynamic_cast<segytrace const &>(y);
            
            segy & trout = cpout.getMetadata();
            segy const & trin = cpin.getMetadata();
            
            Value val;
            float gx;
            float x;
            float dt;
            float t0;
            float wm;
            int nt;
            
            string name="ns";
            gethdval(&trin,(char *)(name.c_str()),&val);
            nt=vtoi(hdtype((char *)(name.c_str())),val);
            
            name="dt";
            gethdval(&trin,(char *)(name.c_str()),&val);
            dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
            name="delrt";
            gethdval(&trin,(char *)(name.c_str()),&val);
            t0=vtof(hdtype((char *)(name.c_str())),val);
            
            name="gx";
            gethdval(&trin,(char *)(name.c_str()),&val);
            gx = vtof(hdtype((char *)(name.c_str())),val);
            
            name="sx";
            gethdval(&trin,(char *)(name.c_str()),&val);
            x = gx - vtof(hdtype((char *)(name.c_str())),val);
            
            // adjust mute so that w=0 works
            if (mute_type) {
                wm = MAX(dt/s,w);
                float wtmp = mutefun((gx - s)/wm) * mutefun((tm - gx)/wm);
                for (int i=0;i<nt;i++)
                    trout.data[i] = trin.data[i]*wtmp;
            }
            else {
                wm=MAX(dt,w);
                //      cerr<<"mute at offset="<<x<<" slope="<<s<<" zotime="<<tm<<" width="<<wm<<endl;
                for (int i=0;i<nt;i++)
                    //	trout.data[i] = trin.data[i]*mutefun((t0+i*dt-s*fabs(x)-tm)/wm);
                    // this version: tm = time AFTER first sample
                    trout.data[i] = trin.data[i]*mutefun((i*dt+t0-s*fabs(x)-tm)/wm);
            }
        }
        catch (bad_cast) {
            RVLException e;
            e<<"Error: SEGYLinMute::operator()\n";
            e<<"input LDC not segytrace\n";
            throw e;
        }
    }
    
    
    void SEGYTaperMute::operator()(LocalDataContainer<float> & x,
                                 LocalDataContainer<float> const & y) {
        
        try {
            
            segytrace & cpout = dynamic_cast<segytrace &>(x);
            segytrace const & cpin = dynamic_cast<segytrace const &>(y);
            
            segy & trout = cpout.getMetadata();
            segy const & trin = cpin.getMetadata();
            
            Value val;
            float gx;
            float sx;
            float x;
            float dt;
            float t0;
            float wm;
            float wt;
            int nt;
            
            string name="ns";
            gethdval(&trin,(char *)(name.c_str()),&val);
            nt=vtoi(hdtype((char *)(name.c_str())),val);
            
            name="dt";
            gethdval(&trin,(char *)(name.c_str()),&val);
            dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
            name="delrt";
            gethdval(&trin,(char *)(name.c_str()),&val);
            t0=vtof(hdtype((char *)(name.c_str())),val);
            
            name="gx";
            gethdval(&trin,(char *)(name.c_str()),&val);
            gx = vtof(hdtype((char *)(name.c_str())),val);
            
            name="sx";
            gethdval(&trin,(char *)(name.c_str()),&val);
            sx = vtof(hdtype((char *)(name.c_str())),val);
            x = gx - sx; 
            
            int itw = (int) (tw/dt + 0.1);
            itw = MAX(1,itw);

            float sxt_wt = 1.0f;
            wt = MAX(1.0,sx_width);
            if (sx < sx_min) sxt_wt = cosfun2((sx_min-sx)/wt);
            if (sx > sx_max) sxt_wt = cosfun2((sx-sx_max)/wt);
//cerr << " sxt_wt = " << sxt_wt << endl;
            // adjust mute so that w=0 works
            if (mute_type) {
                wm = MAX(dt/s,w);
                float wtmp = mutefun((gx - s)/wm) * mutefun((tm - gx)/wm);
                wtmp = wtmp * sxt_wt;
                wt=MAX(1.0,width);
                float wttmp=1.0;
                if (taper_type){ // offset
                    if (x < taper_min) wttmp = cosfun2((taper_min-x)/wt);
                    else if (x > taper_max) wttmp = cosfun2((x - taper_max)/wt);
#pragma ivdep                    
                    for (int i=0;i<nt-itw;i++)
                        trout.data[i] = trin.data[i]*wtmp*wttmp;
#pragma ivdep
                    for (int i=nt-itw;i<nt;i++)
                        trout.data[i] = trin.data[i]*wtmp*wttmp*cosfun(float(nt-i)/(itw+0.0f));

                }
                else { // taper geophone position
                    if (gx < taper_min) wttmp = cosfun2((taper_min-gx)/wt);
                    else if (gx > taper_max) wttmp = cosfun2((gx - taper_max)/wt);
#pragma ivdep                    
                    for (int i=0;i<nt-itw;i++)
                        trout.data[i] = trin.data[i]*wtmp*wttmp;
#pragma ivdep
                    for (int i=nt-itw;i<nt;i++)
                        trout.data[i] = trin.data[i]*wtmp*wttmp*cosfun(float(nt-i)/(itw+0.0f));
                }
            }
            else {
                wm=MAX(dt,w);
                float wttmp=1.0;
                wt=MAX(1.0,width);
                if (taper_type){ // offset
                    if (x < taper_min) wttmp = cosfun2((taper_min-x)/wt);
                    else if (x > taper_max) wttmp = cosfun2((x - taper_max)/wt);
		    wttmp = wttmp * sxt_wt;
#pragma ivdep                    
                    for (int i=0;i<nt-itw;i++)
                        trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm);
#pragma ivdep
                    for (int i=nt-itw;i<nt;i++)
                        trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm)*cosfun(float(nt-i)/(itw+0.0f));
                }
                else { // taper geophone position
                    if (gx < taper_min) wttmp = cosfun2((taper_min-gx)/wt);
                    else if (gx > taper_max) wttmp = cosfun2((gx - taper_max)/wt);
                    wttmp = wttmp * sxt_wt;
#pragma ivdep
                    for (int i=0;i<nt-itw;i++)
                        trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm);
#pragma ivdep
                    for (int i=nt-itw;i<nt;i++)
                        trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm)*cosfun(float(nt-i)/(itw+0.0f));
                }
            }
        }
        catch (bad_cast) {
            RVLException e;
            e<<"Error: SEGYTaperMute::operator()\n";
            e<<"input LDC not segytrace\n";
            throw e;
        }
    }
    
    void SEGYFwdInt::operator()(LocalDataContainer<float> & x,
                                LocalDataContainer<float> const & y) {
        
        try {
            
            segytrace & cpout = dynamic_cast<segytrace &>(x);
            segytrace const & cpin = dynamic_cast<segytrace const &>(y);
            
            segy & trout = cpout.getMetadata();
            segy const & trin = cpin.getMetadata();
            
            Value val;
            float dt;
            int nt;
            
            string name="ns";
            gethdval(&trin,(char *)(name.c_str()),&val);
            nt=vtoi(hdtype((char *)(name.c_str())),val);
            gethdval(&trout,(char *)(name.c_str()),&val);
            if (nt != vtoi(hdtype((char *)(name.c_str())),val)) {
                RVLException e;
                e<<"Error: SEGYFwdInt::operator()\n";
                e<<"input, output incompatible\n";
                throw e;
            }
            name="dt";
            gethdval(&trin,(char *)(name.c_str()),&val);
            dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
            memcpy(trout.data,trin.data,nt*sizeof(float));
            trout.data[0]=0.0f;
            for (int j=0;j<nint;j++) {
#pragma ivdep
                for (int i=1;i<nt;i++) trout.data[i]
                    = trout.data[i-1]+trin.data[i-1]*dt;
            }
        }
        catch (bad_cast) {
            RVLException e;
            e<<"Error: SEGYFwdInt::operator()\n";
            e<<"input LDC not segytrace\n";
            throw e;
        }
    }
    
    void SEGYAdjInt::operator()(LocalDataContainer<float> & x,
                                LocalDataContainer<float> const & y) {
        
        try {
            
            segytrace & cpout = dynamic_cast<segytrace &>(x);
            segytrace const & cpin = dynamic_cast<segytrace const &>(y);
            
            segy & trout = cpout.getMetadata();
            segy const & trin = cpin.getMetadata();
            
            Value val;
            float dt;
            int nt;
            
            string name="ns";
            gethdval(&trin,(char *)(name.c_str()),&val);
            nt=vtoi(hdtype((char *)(name.c_str())),val);
            gethdval(&trout,(char *)(name.c_str()),&val);
            if (nt != vtoi(hdtype((char *)(name.c_str())),val)) {
                RVLException e;
                e<<"Error: SEGYAdjInt::operator()\n";
                e<<"input, output incompatible\n";
                throw e;
            }
            name="dt";
            gethdval(&trin,(char *)(name.c_str()),&val);
            dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
            memcpy(trout.data,trin.data,nt*sizeof(float));
            trout.data[nt-1]=0.0f;
            for (int j=0;j<nint;j++) {
#pragma ivdep
                for (int i=nt-2;i>-1;i--) trout.data[i] 
                    = trout.data[i+1]+trin.data[i+1]*dt;
            }
        }
        catch (bad_cast) {
            RVLException e;
            e<<"Error: SEGYAdjInt::operator()\n";
            e<<"input LDC not segytrace\n";
            throw e;
        }
    }
    
    
}
