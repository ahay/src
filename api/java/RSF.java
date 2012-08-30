package rsf;

public class RSF{

    public RSF(String[] args){      
        String[] targs = null;
        if (args.length == 0){
            targs = new String[]{"java"};
        }  else {
            targs = new String[args.length+1];
            targs[0] = "java";
            for(int i = 0; i < args.length; ++i){
                targs[i+1] = args[i];
                //System.err.printf("targs: %d %s\n",i+1,args[i]);
            }
        }
        m8r.sf_init(targs.length,targs);
    }

    public boolean getBool(String key, boolean fallback){
        boolean found = m8r.sf_getbool(key);
        return found;
    }

    public int getInt(String key, int fallback){
        int[] value = new int[]{0};
        boolean found = m8r.sf_getint(key,value);
        if (found) return value[0];
        else return fallback;
    }

    public float getFloat(String key, float fallback){
        float[] value = new float[]{0.0f};
        boolean found = m8r.sf_getfloat(key,value);
        if (found) return value[0];
        else return fallback;
    }

    public String getString(String key, String fallback){
        String value = m8r.sf_getstring(key);
        if (value == null) return fallback;
        else return value;
    }
}
